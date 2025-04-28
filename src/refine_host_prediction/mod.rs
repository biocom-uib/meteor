use std::sync::Arc;

use clap::{Args, ValueEnum};
use lending_iterator::{
    higher_kinded_types::{HKTRef, HKT},
    LendingIterator,
};
use metagenomic_matches::{MergeCrisprMatchArgs, MetagenomicMatches};
use string_interner::DefaultSymbol;

use crate::{
    csv::stream::{self as csv_stream, CsvReaderIterExt},
    preprocessed_taxonomy::{PreprocessedTaxonomyArgs, with_some_taxonomy},
    tango_assign::AssignmentRecord,
    taxonomy::{
        LabeledTaxonomy, RootedTree, Taxonomy,
        tree::{NodeId, postorder_ann::PostorderAnnPool, walk::{LcaCache, RootedTreeWalk}},
    },
    util::{
        clap::{ExistingFilePath, OutputFileFlag},
        interners::{Interner, Resolver},
    },
};


mod metagenomic_matches;

#[cfg(feature = "refine_host_prediction_polars")]
mod polars_impl;
#[cfg(feature = "refine_host_prediction_polars")]
use polars_impl::refine_host_prediction_with_tax_impl;

#[cfg(not(feature = "refine_host_prediction_polars"))]
mod csv_impl;
#[cfg(not(feature = "refine_host_prediction_polars"))]
use csv_impl::refine_host_prediction_with_tax_impl;


#[derive(Copy, Clone, PartialEq, Eq, Default, ValueEnum)]
pub enum TreeRelation {
    /// Keep a prediction P if an assignment which is more general than P exists.
    Ascendants = 0b10,
    /// Keep a prediction P if an assignment which is more specific than P exists.
    Descendants = 0b01,
    /// Keep a prediction P if an assignment which is either more general than P exists.
    #[default]
    Related = 0b11,
    /// Keep a prediction P only if an assignment equal to P exists (not recommended).
    Exact = 0b00,
}

impl TreeRelation {
    pub fn include_ascendants(self) -> bool {
        (self as u32 & Self::Ascendants as u32) != 0
    }

    pub fn include_descendants(self) -> bool {
        (self as u32 & Self::Descendants as u32) != 0
    }
}

#[derive(Copy, Clone, PartialEq, Eq, Default, ValueEnum)]
pub enum Aggregation {
    /// Semicolon-separated list of matching contig names
    List,
    /// Number of matching contigs
    #[default]
    Count,
}

#[derive(Args, Clone)]
pub struct RefinementOptions {
    /// Keep prediction class names
    #[clap(long)]
    keep_names: bool,

    /// Move up to a specific rank, and drop less specific predictions
    #[clap(long)]
    rank: Option<String>,

    /// How to filter predictions based on metagenomic assignments.
    #[clap(long, value_enum, default_value_t)]
    keep_if: TreeRelation,

    /// How to report matching metagenomic contigs.
    #[clap(long, value_enum, default_value_t)]
    agg_contigs: Aggregation,

    /// How to report matching CRISPRs.
    #[clap(long, value_enum, default_value_t, requires = "crispr_match")]
    agg_crispr: Aggregation,
}

/// Refine a host prediction using metagenomic data.
#[derive(Args)]
pub struct RefineHostPredictionArgs {
    /// The metagenomic assignment (see `meteor tango-assign`)
    #[clap(long)]
    assignment: ExistingFilePath,

    /// The merged host prediction (see `meteor merge-host-prediction`)
    #[clap(long)]
    host_prediction: ExistingFilePath,

    #[clap(flatten)]
    refinement_options: RefinementOptions,

    #[clap(flatten)]
    output: OutputFileFlag,

    #[clap(flatten)]
    taxonomy: PreprocessedTaxonomyArgs,

    #[clap(flatten)]
    crispr_args: Option<MergeCrisprMatchArgs>,

    //#[clap(flatten)]
    //trna_args: Option<MergeTrnaMatchArgs>,

    //#[clap(flatten)]
    //metilase_args: Option<MergeMetilaseMatchArgs>,
}

#[repr(transparent)]
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct VirusSymbol(DefaultSymbol);

#[repr(transparent)]
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct HostSymbol(DefaultSymbol);

#[derive(Default)]
struct TaxonEvidence {
    // populated only if refinement_options.include_ascendants()
    ascendant_assignments: u32,
    ascendant_crispr_matches: u32,

    // populated only if refinement_options.include_descendants()
    descendant_assignments: u32,
    descendant_crispr_matches: u32,

    // always populated
    assigned_contigs: Vec<HostSymbol>,
    crispr_matches: Vec<(HostSymbol, VirusSymbol)>,
}

struct MetagenomicEvidence<VI: Resolver, HI: Resolver> {
    ann_pool: Option<PostorderAnnPool<TaxonEvidence>>,
    virus_interner: VI,
    host_interner: HI,
    has_crispr: bool,
}

impl<VI: Resolver, HI: Resolver> MetagenomicEvidence<VI, HI> {
    fn new(virus_interner: VI, host_interner: HI) -> Self {
        Self {
            ann_pool: None, //PostorderAnnPool::new(tax, tax.get_root()),
            virus_interner,
            host_interner,
            has_crispr: false,
        }
    }
}

impl<VI, HI> MetagenomicEvidence<VI, HI>
where
    VI: Resolver<Symbol = DefaultSymbol>,
    HI: Interner<Symbol = DefaultSymbol, ValueHKT = HKTRef<str>>,
{
    fn load_assignments(
        &mut self,
        tax: &impl RootedTree,
        lca_cache: &LcaCache,
        args: &RefineHostPredictionArgs,
        matches: &MetagenomicMatches,
    ) -> anyhow::Result<()> {
        self.has_crispr = matches.crispr.is_some();

        csv_stream::tsv_reader_builder(true)
            .from_path(args.assignment.as_path())?
            .into_lending_iter()
            .into_deserialize::<HKT!(AssignmentRecord<'_>)>(None)
            .try_for_each(|record| {
                let record = &record?;

                let _: Option<()> =
                    try {
                        let assigned_node = tax.fixup_node(record.assigned_taxid)?;

                        let mut slice = PostorderAnnPool::new_or_expand_for(
                            &mut self.ann_pool,
                            tax,
                            lca_cache,
                            assigned_node,
                        )
                        .ok()?;

                        let ann = slice.lca_mut();

                        let contig_sym = self.host_interner.get_or_intern(&record.query_id);
                        ann.assigned_contigs.push(HostSymbol(contig_sym));

                        if let Some(crispr_matches) = &matches.crispr {
                            if let Some(crisprs) = crispr_matches.0.get(contig_sym) {
                                ann.crispr_matches.extend(crisprs.map(|virus_sym| {
                                    (HostSymbol(contig_sym), VirusSymbol(virus_sym))
                                }));
                            }
                        }
                    };

                anyhow::Ok(())
            })?;

        let Some(mut root_slice) = self.ann_pool.as_mut().map(|pool| pool.slice_root_mut()) else {
            return Ok(());
        };

        for i in root_slice.indices() {
            let (_, ann) = &mut root_slice[i];

            let assignments = ann.assigned_contigs.len() as u32 + ann.descendant_assignments;
            let crisprs = ann.crispr_matches.len() as u32 + ann.descendant_crispr_matches;

            if let Some(parent_index) = root_slice.parent_index(i) {
                let (_, parent_ann) = &mut root_slice[parent_index];
                parent_ann.descendant_assignments += assignments;
                parent_ann.descendant_crispr_matches += crisprs;
            }
        }

        for i in root_slice.indices().rev() {
            if let Some(parent_index) = root_slice.parent_index(i) {
                let (_, parent_ann) = &mut root_slice[parent_index];

                let assignments =
                    parent_ann.assigned_contigs.len() as u32 + parent_ann.ascendant_assignments;
                let crisprs =
                    parent_ann.crispr_matches.len() as u32 + parent_ann.ascendant_crispr_matches;

                let (_, ann) = &mut root_slice[i];
                ann.ascendant_assignments += assignments;
                ann.ascendant_crispr_matches += crisprs;
            }
        }

        Ok(())
    }
}

struct PredictionEvidence {
    assigned_contigs: Vec<HostSymbol>,
    crispr_matches: Vec<HostSymbol>,
}

impl<VI: Resolver, HI: Resolver> MetagenomicEvidence<VI, HI> {
    pub fn filter_crispr_matches(
        virus_sym: Option<VirusSymbol>,
        crispr_matches: &[(HostSymbol, VirusSymbol)],
    ) -> impl Iterator<Item = HostSymbol> + '_ {

        // If virus_sym is none, we assume there is no CRISPR info.
        // Otherwise it would have been already interned.

        virus_sym
            .into_iter()
            .flat_map(|virus_sym| {
                crispr_matches
                    .iter()
                    .filter(move |(_host, virus)| *virus == virus_sym)
                    .map(|(host, _virus)| *host)
            })
    }

    pub fn refine_and_enrich<Tax: Taxonomy>(
        &self,
        tax: &Tax,
        virus_sym: Option<VirusSymbol>,
        node: NodeId,
    ) -> Option<PredictionEvidence> {

        let ann_slice = self.ann_pool.as_ref()?.slice_clade(tax, node)?;
        let node_index = ann_slice.global_lca_index();

        let mut assigned_contigs = Vec::new();
        let mut crispr_matches = Vec::new();

        let node_ann = ann_slice.lca();
        assigned_contigs.extend(&node_ann.assigned_contigs);
        crispr_matches.extend(Self::filter_crispr_matches(virus_sym, &node_ann.crispr_matches));

        if node_ann.descendant_assignments > 0 {
            // skip(1) because we already have ours
            for (_desc, desc_ann) in ann_slice.iter().rev().skip(1) {
                assigned_contigs.extend(&desc_ann.assigned_contigs);
                crispr_matches.extend(Self::filter_crispr_matches(virus_sym, &desc_ann.crispr_matches));
            }
        }

        if node_ann.ascendant_assignments > 0 {
            if let Some(pool) = &self.ann_pool {
                let root_slice = pool.slice_root();

                let mut i = node_index;

                // intentionally skipping node_index, as we already have them
                while let Some(parent_index) = root_slice.parent_index(i) {
                    i = parent_index;

                    let (_asc, asc_ann) = &root_slice[i];
                    assigned_contigs.extend(&asc_ann.assigned_contigs);
                    crispr_matches.extend(Self::filter_crispr_matches(virus_sym, &asc_ann.crispr_matches));
                }
            }
        }

        assigned_contigs.sort();
        assigned_contigs.dedup();

        crispr_matches.sort();
        crispr_matches.dedup();

        if !assigned_contigs.is_empty() || !crispr_matches.is_empty() {
            Some(PredictionEvidence { assigned_contigs, crispr_matches })
        } else {
            None
        }
    }
}

fn adjust_prediction_to_rank<Tax>(tax: &Tax, rank_sym: Tax::RankSym, node: NodeId) -> Option<NodeId>
where
    Tax: Taxonomy,
{
    if tax.find_rank(node) == Some(rank_sym) {
        Some(node)
    } else {
        tax
            .strict_ancestors(node)
            .find(move |ancestor| tax.find_rank(*ancestor) == Some(rank_sym))
    }
}

fn adjust_prediction_to_rank_with_names<'a, Tax>(
    tax: &'a Tax,
    rank_sym: Tax::RankSym,
    node: NodeId,
    class_name: Option<&'a str>,
) -> Option<(NodeId, Option<&'a str>)>
where
    Tax: LabeledTaxonomy
{
    let ancestor_at_rank = adjust_prediction_to_rank(tax, rank_sym, node)?;

    let ancestor_class_name = if node == ancestor_at_rank {
        class_name
    } else {
        tax.some_label_of(ancestor_at_rank)
    };

    Some((ancestor_at_rank, ancestor_class_name))
}

pub fn refine_host_prediction_with_tax<Tax>(
    args: RefineHostPredictionArgs,
    tax: Tax,
) -> anyhow::Result<()>
where
    Tax: LabeledTaxonomy + Send + Sync + 'static,
{
    use string_interner::DefaultStringInterner;

    let tax = Arc::new(tax);
    let lca_cache = LcaCache::compute(&*tax);

    let rank_sym = if let Some(rank) = &args.refinement_options.rank {
        if let Some(rank_sym) = tax.lookup_rank_sym(rank) {
            Some(rank_sym)
        } else {
            anyhow::bail!("Rank not found in the taxonomy: {}", rank)
        }
    } else {
        None
    };

    let mut virus_interner = DefaultStringInterner::new();
    let mut host_interner = DefaultStringInterner::new();

    let matches = {
        let mut matches = MetagenomicMatches::new();

        if let Some(crispr_args) = &args.crispr_args {
            eprintln!("Loading CRISPR matches");
            matches.load_crispr_matches_using(
                crispr_args,
                &mut host_interner,
                &mut virus_interner,
            )?;
        }

        matches
    };

    let evidence = {
        let mut evidence = MetagenomicEvidence::new(virus_interner, host_interner);

        eprintln!("Loading metagenomic assignment");
        evidence.load_assignments(&*tax, &lca_cache, &args, &matches)?;
        evidence
    };

    eprintln!("Refining...");

    refine_host_prediction_with_tax_impl(&args, tax, rank_sym, evidence)
}

pub fn refine_host_prediction(args: RefineHostPredictionArgs) -> anyhow::Result<()> {
    eprintln!("Loading taxonomy");
    let taxonomy = args.taxonomy.deserialize()?;

    with_some_taxonomy!(taxonomy.tree, tax => {
        refine_host_prediction_with_tax(args, tax)?;
    });

    Ok(())
}
