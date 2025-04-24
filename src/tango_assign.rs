use std::borrow::Cow;
use std::collections::HashMap;
use std::io;
use std::path::Path;
use std::sync::RwLock;

use anyhow::{bail, Context};
use clap::Args;
use csv::{StringRecord, WriterBuilder};
use extended_rational::{Rational, URational};
use meteor_macros::PolarsSchema;
use rayon::iter::{ParallelBridge, ParallelIterator};
use serde::{Deserialize, Serialize};
use std::sync::mpsc::{Receiver, Sender};

use crate::csv::filter::FilterableRecord;
use crate::preprocessed_taxonomy::{
    with_some_ncbi_or_newick_taxonomy, with_some_taxonomy, PreprocessedTaxonomy,
    PreprocessedTaxonomyArgs,
};
use crate::taxonomy::formats::ncbi::{NamesAssoc, NcbiTaxonomy};
use crate::taxonomy::tree::node_id::NodeIdSet;
use crate::taxonomy::tree::postorder_ann::{PostorderAnnPool, PostorderAnnSliceMut};
use crate::taxonomy::tree::walk::{LcaCache, LcaMapError};
use crate::taxonomy::{LabeledTaxonomy, NodeId, RootedTree, Taxonomy};
use crate::tool::blast::blastout::BlastOutFullRecord;
use crate::util::clap::ExistingFilePath;
use crate::util::io::MaybeGzDecoder;
use crate::util::{self, writing_new_file_or_stdout};


#[derive(Clone, Debug, Default)]
pub struct TaxonMatchCount {
    is_leaf: bool,
    descendant_leaves: u32,
    matches: u32,
}

impl TaxonMatchCount {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn clear_matches(&mut self) {
        self.matches = 0;
    }
}

type TaxonomyMatchCountPool = PostorderAnnPool<TaxonMatchCount>;
type TaxonomyMatchCountSliceMut<'p> = PostorderAnnSliceMut<'p, TaxonMatchCount>;

fn initialize<Tax: Taxonomy>(tax: &Tax, ann_pool: &mut TaxonomyMatchCountPool) {
    let mut slice = ann_pool.slice_root_mut();

    for i in slice.indices() {
        let &mut (node, ref mut ann) = &mut slice[i];

        let is_leaf = tax.is_leaf(node);
        ann.is_leaf = is_leaf;

        let descendant_leaves = if is_leaf {
            ann.descendant_leaves = 1;
            1
        } else {
            ann.descendant_leaves
        };

        if node != tax.get_root() {
            let parent_index = slice.parent_index(i).unwrap();
            let (_, parent_ann) = &mut slice[parent_index];

            parent_ann.descendant_leaves += descendant_leaves;
        }
    }
}

fn annotate_match_count<'p, Tax: LabeledTaxonomy>(
    taxonomy: &Tax,
    lca_cache: &LcaCache,
    reads_node_ids: &NodeIdSet,
    ann_pool: &'p mut TaxonomyMatchCountPool,
) -> anyhow::Result<TaxonomyMatchCountSliceMut<'p>> {
    let reads_lca_id = match lca_cache.lca_reduce(reads_node_ids) {
        Some(Ok(reads_lca_id)) => reads_lca_id,
        Some(Err(LcaMapError::NodeIdToRange(node_id))) => {
            bail!(
                "Could not map node {node_id} labeled as {} to a LCA range",
                taxonomy.some_label_of(node_id).unwrap_or("")
            );
        }
        Some(Err(LcaMapError::IndexToNodeId(idx))) => {
            bail!("Could not map LCA index {idx:?} back to a NodeId");
        }
        None => {
            bail!("Empty list of NodeIds provided to compute LCA")
        }
    };

    let Some(mut ann_slice) = ann_pool.slice_clade_mut(taxonomy, reads_lca_id) else {
        bail!("Could not locate first and last leaf from the annotation pool for LCA = {reads_lca_id}");
    };

    for node_index in ann_slice.indices() {
        let &mut (node, ref mut ann) = &mut ann_slice[node_index];

        let matches = if ann.is_leaf {
            let matches = reads_node_ids.contains(node) as u32;
            ann.matches = matches;
            matches
        } else {
            ann.matches
        };

        if node != reads_lca_id {
            let parent_index = ann_slice.parent_index(node_index).unwrap();
            let (_, parent_ann) = &mut ann_slice[parent_index];

            parent_ann.matches += matches;
        }
    }

    Ok(ann_slice)
}

#[derive(Clone, Debug, Default)]
struct TaxonPrecisionRecall {
    tps: u32,
    fps: u32,
    //tns: u32,
    fns: u32,
}

impl TaxonPrecisionRecall {
    pub fn from_match_count(lca_ann: &TaxonMatchCount, ann: &TaxonMatchCount) -> Self {
        TaxonPrecisionRecall {
            tps: ann.matches,
            fps: ann.descendant_leaves - ann.matches, // = ann.nonmatches
            //tns: lca_ann.nonmatches - ann.nonmatches;
            fns: lca_ann.matches - ann.matches,
        }
    }

    pub fn penalty(&self, q: URational) -> Rational {
        let tps = self.tps as i32;
        let fps = self.fps as i32;
        let fns = self.fns as i32;

        let qn = q.numerator() as i32;
        let qd = q.denominator() as i32;

        let num = qn * fns + qd * fps - qn * fps;

        let denom = if tps != 0 {
            qd * tps
        } else {
            qd
        };

        Rational::new(num as i64, denom as i64)

        /*
        // With rationals (2x slower)

        // q = qn/qd
        let one = Rational::one();
        if tps != 0 {
            let fns_ratio = Rational::new(fns, tps);
            let fps_ratio = Rational::new(fps, tps);
            q * fns_ratio + (one - q) * fps_ratio
            // qn/qd fn/tp + (1 - qn/qd) fp/tp
            // (qn fn) / (qd tp) + (qd fp - qn fp) / (qd tp)
            // (qn fn + qd fp - qn fp) / (qd tp)
        } else {
            let fns = Rational::from(fns);
            let fps = Rational::from(fps);
            q * fns + (one - q) * fps
            // qn/qd fn + (1 - qn/qd) fp
            // (qn fn)/qd + (qd - qn)/qd fp
            // (qn fn)/qd + (qd fp - qn fp)/qd
            // (qn fn + qd fp - qn fp) / qd
        }
        */
    }
}

fn assign_reads_and_clear(
    anns: &mut PostorderAnnSliceMut<TaxonMatchCount>,
    q: URational,
) -> (Vec<NodeId>, Rational) {
    let mut min_nodes = Vec::new();
    let mut min_penalty = Rational::new(-1, 1);

    let lca_ann = anns.lca().clone();

    for (node_id, ann) in anns.iter_mut() {
        let penalty = TaxonPrecisionRecall::from_match_count(&lca_ann, ann).penalty(q);
        ann.clear_matches();

        debug_assert!(!penalty.is_infinity());

        if min_penalty.is_negative() || penalty < min_penalty {
            min_penalty = penalty;
            min_nodes.clear();
            min_nodes.push(node_id);
        } else if penalty == min_penalty {
            min_nodes.push(node_id);
        }
    }

    (min_nodes, min_penalty)
}

/// Produce a metagenomic assignment using the TANGO algorithm.
#[derive(Args)]
pub struct TangoAssignArgs {
    #[clap(flatten)]
    taxonomy: PreprocessedTaxonomyArgs,

    /// Path to the preprocessed reads file (see `meteor preprocess-blastout`). It should be a
    /// file containing the parsed output of a mapping program, in the (tab-separated) format
    ///
    /// [read_id]    [species_id_1];...;[species_id_n]
    preprocessed_reads: ExistingFilePath,

    /// Output path (use '-' for STDOUT)
    #[clap(short, long, default_value = "-")]
    output: String,

    /// Parameter that allows balancing the taxonomic assignment between precision (q = 0) and
    /// recall (q = 1), with q = 0.5 (default) corresponding to the F-measure (harmonic mean of
    /// precision and recall)
    #[clap(short, default_value_t = 0.5)]
    q: f32,
}

struct ReadsCsvHeader {
    //query_id_col: String,
    subjects_id_col: String,
    weights_col: Option<String>,
}

fn peek_reads_csv_header<R>(csv_reader: &mut csv::Reader<R>) -> anyhow::Result<ReadsCsvHeader>
where
    R: io::Read,
{
    let mut record = StringRecord::new();
    csv_reader.read_record(&mut record)?;
    assert!(record.len() == 2);
    let _query_id_col = record[0].to_owned();
    let subjects_col = &record[1];

    if let Some((subjects_id_col, weights_id)) = subjects_col.split_once('/') {
        Ok(ReadsCsvHeader {
            //query_id_col,
            subjects_id_col: subjects_id_col.to_owned(),
            weights_col: Some(weights_id.to_owned()),
        })
    } else {
        Ok(ReadsCsvHeader {
            //query_id_col,
            subjects_id_col: subjects_col.to_owned(),
            weights_col: None,
        })
    }
}

struct QueryAssignments {
    query_id: String,
    assigned_nodes: Vec<NodeId>,
    penalty: f32,
}

fn produce_assignments<R, Tax, F>(
    mut csv_reader: csv::Reader<R>,
    header: &ReadsCsvHeader,
    tax: &Tax,
    read_to_taxid: F,
    q: URational,
    sender: Sender<QueryAssignments>,
) -> anyhow::Result<()>
where
    R: io::Read + Send,
    Tax: LabeledTaxonomy + Sync,
    F: Fn(&str) -> Option<NodeId> + Send + Sync,
{
    let lca_cache = LcaCache::compute(tax);

    let mut ann_pool = TaxonomyMatchCountPool::new(tax, tax.get_root());
    initialize(tax, &mut ann_pool);

    let skipped_queries = RwLock::new(HashMap::new());
    let skipped_subjects = RwLock::new(HashMap::new());

    csv_reader.records().par_bridge().try_for_each_with(
        (sender, ann_pool),
        |(sender_copy, ann_pool), record| {
            let record = record?;

            assert!(record.len() == 2);

            let query_id = &record[0];
            let subject_ids = &record[1];

            if subject_ids.is_empty() {
                return Ok(());
            }

            let mut num_subject_ids = 0;

            let read_to_taxid_warn = |name| {
                let taxid = read_to_taxid(name);

                if read_to_taxid(name).is_none() {
                    *skipped_subjects
                        .write()
                        .unwrap()
                        .entry(name.to_owned())
                        .or_insert(0) += 1;
                }

                taxid
            };

            let subject_ids_iter = subject_ids.split(';').inspect(|_| {
                num_subject_ids += 1;
            });

            let subject_node_ids: NodeIdSet = if header.weights_col.is_some() {
                subject_ids_iter
                    .flat_map(|item| read_to_taxid_warn(item.split_once('/')?.0))
                    .collect()
            } else {
                subject_ids_iter.flat_map(read_to_taxid_warn).collect()
            };

            if subject_node_ids.is_empty() {
                skipped_queries
                    .write()
                    .unwrap()
                    .insert(query_id.to_owned(), num_subject_ids);
            } else {
                let mut anns = annotate_match_count(tax, &lca_cache, &subject_node_ids, ann_pool)?;

                let (assigned_nodes, penalty) = assign_reads_and_clear(&mut anns, q);

                sender_copy.send(QueryAssignments {
                    query_id: query_id.to_owned(),
                    assigned_nodes,
                    penalty: f32::from(penalty),
                })?;
            }

            anyhow::Ok(())
        },
    )?;

    eprintln!();

    let skipped_subjects = skipped_subjects.into_inner().unwrap();
    for (name, count) in &skipped_subjects {
        eprintln!(
            "Warning: No leaf taxid matching {name:?} found in the taxonomy, skipped {count} times"
        );
    }

    let skipped_queries = skipped_queries.into_inner().unwrap();
    for (name, num_subjects) in &skipped_queries {
        eprintln!(
            "Warning: Skipped query {name} as none out of {num_subjects} subjects were recognized",
        );
    }

    eprintln!(
        "Skipped {} queries and {} subjects",
        skipped_queries.len(),
        skipped_subjects.values().sum::<i32>(),
    );

    Ok(())
}

fn ncbi_taxonomy_lookup_taxid_leaf<Names: 'static + NamesAssoc + Send>(
    tax: &NcbiTaxonomy<Names>,
) -> impl Fn(&str) -> Option<NodeId> + '_ {
    |name: &str| {
        let node_id = tax.fixup_node(name.parse().ok()?)?;

        if tax.is_leaf(node_id) {
            Some(node_id)
        } else {
            None
        }
    }
}

fn labeled_taxonomy_lookup_taxid_leaf<Tax>(tax: &Tax) -> impl Fn(&str) -> Option<NodeId> + '_
where
    Tax: LabeledTaxonomy,
{
    |name: &str| {
        let mut iter = tax.nodes_with_label(name).filter(|&node| tax.is_leaf(node));

        let first = iter.next();

        if first.is_some() && iter.next().is_some() {
            eprintln!("\nWarning: Multiple taxids matched {name:?} in the taxonomy, returning the first one. Consider using staxid column in the BLAST output if possible.");
        }

        first
    }
}

fn load_reads_and_produce_assignments(
    reads_path: &Path,
    taxonomy: &PreprocessedTaxonomy,
    q: URational,
    sender: Sender<QueryAssignments>,
) -> anyhow::Result<()> {
    let reader = MaybeGzDecoder::from_path(reads_path)?;

    let mut csv_reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(reader);

    let header = peek_reads_csv_header(&mut csv_reader)?;

    macro_rules! produce_assignments_with {
        ($tax:expr, $lookup_taxid:expr) => {
            produce_assignments(csv_reader, &header, $tax, $lookup_taxid, q, sender)
        };
    }

    with_some_ncbi_or_newick_taxonomy!(&taxonomy.tree,
        // TODO: does it really need to be a leaf?
        ncbi: tax => {
            if header.subjects_id_col == BlastOutFullRecord::STAXID.0 {
                produce_assignments_with!(tax, ncbi_taxonomy_lookup_taxid_leaf(tax))?;
            } else if header.subjects_id_col == BlastOutFullRecord::SSCINAME.0 {
                produce_assignments_with!(tax, labeled_taxonomy_lookup_taxid_leaf(tax))?;
            } else {
                anyhow::bail!("Unable to map read subject ID's ({}) to NCBI Taxonomy ID's", &header.subjects_id_col)
            }
        },
        newick: tax => {
            eprintln!("Warning: Unable to ensure that Newick taxonomy labels match the search subject ID key ({})", &header.subjects_id_col);

            produce_assignments_with!(tax, labeled_taxonomy_lookup_taxid_leaf(tax))?;
        }
    );

    Ok(())
}

#[derive(FilterableRecord, Serialize, Deserialize, PolarsSchema)]
#[filter(polars_literal)]
pub struct AssignmentRecord<'a> {
    #[filter(as_string)]
    #[serde(borrow)]
    pub query_id: Cow<'a, str>,
    pub assigned_taxid: u32,
    #[filter(as_string)]
    #[serde(borrow)]
    pub assigned_name: Cow<'a, str>,
    #[filter(as_string)]
    #[serde(borrow)]
    pub assigned_rank: Cow<'a, str>,
    pub penalty: f32,
}

fn write_assignments(
    output: &str,
    taxonomy: &PreprocessedTaxonomy,
    assignments: Receiver<QueryAssignments>,
) -> anyhow::Result<()> {
    with_some_taxonomy!(&taxonomy.tree, tax => {
        writing_new_file_or_stdout!(output, writer => {
            let writer = writer.context("Error creating assignments file")?;

            let mut csv_writer = WriterBuilder::new()
                .delimiter(b'\t')
                .has_headers(true)
                .from_writer(writer);

            let mut count = 0;

            for record in assignments {
                for node in record.assigned_nodes {
                    let assigned_name = tax.labels_of(node).next().unwrap_or_else(|| {
                        eprintln!("\nWarning: No label found for taxid {node}");
                        ""
                    });

                    let rank = tax
                        .find_rank(node)
                        .and_then(|rank_sym| tax.rank_sym_str(rank_sym))
                        .unwrap_or("");

                    let penalty = record.penalty;

                    csv_writer.serialize(AssignmentRecord {
                        query_id: Cow::Borrowed(&record.query_id),
                        assigned_taxid: node.into(),
                        assigned_name: Cow::Borrowed(assigned_name),
                        assigned_rank: Cow::Borrowed(rank),
                        penalty,
                    })?;
                }

                count += 1;

                if count % 100 == 0 {
                    eprint!("\rProcessed {count} query sequences");
                }
            }
            eprintln!("\rProcessed {count} query sequences");

            csv_writer.flush()?;
        });
    });

    Ok(())
}

pub fn tango_assign(args: TangoAssignArgs) -> anyhow::Result<()> {
    let taxonomy = args.taxonomy.deserialize()?;

    eprintln!("Loaded taxonomy");

    let taxonomy = &taxonomy;
    let q = Rational::from(args.q);

    let q = if q.is_negative() {
        anyhow::bail!("-q {q}: q must be positive or 0");
    } else {
        q.unsigned()
    };

    let (sender, receiver) = std::sync::mpsc::channel(); // unbounded channel

    //let start_time = std::time::Instant::now();

    let (rl, rr) = rayon::join(
        move || load_reads_and_produce_assignments(args.preprocessed_reads.as_path(), taxonomy, q, sender),
        move || write_assignments(&args.output, taxonomy, receiver),
    );

    //eprintln!("\nTime: {:?}", start_time.elapsed());

    // Times for ES_AV.LargeContigs:
    // base               66.34s
    // +NodeIdSet        121.91s
    // +pool              15.45s
    // -rational           5.45s
    // +postorder_index    5.20s

    util::ignore_broken_pipe_anyhow(rr)?;
    rl?;

    Ok(())
}
