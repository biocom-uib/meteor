use std::{fs::File, io, marker::PhantomData, path::Path};

use clap::Args;
use itertools::Itertools;
use lending_iterator::{LendingIterator, HKT};
use serde::{Deserialize, Serialize};
use string_interner::{DefaultStringInterner, DefaultSymbol};

use crate::{
    csv::stream::CsvReaderIterExt,
    fetch::{self, string_viruses::ProteinLinksRecord},
    preprocessed_taxonomy::{with_some_ncbi_or_newick_taxonomy, PreprocessedTaxonomyArgs},
    tango_assign::AssignmentRecord,
    taxonomy::{
        formats::ncbi::NcbiTaxonomy,
        tree::{
            node_id::NodeIdSet,
            walk::{LcaCache, RootedTreeWalk},
        },
        NodeId, RootedTree, Taxonomy,
    },
    util::{
        self,
        interned_mapping::InternedMapping,
        interners::{PairInterner, Resolver},
        io::maybe_gzdecoder,
        writing_new_file_or_stdout,
    },
};

use super::ProteinVirusTaxidMapping;

struct StringProteinAliases<Entry> {
    mapping: InternedMapping<DefaultStringInterner, PairInterner>,
    phantom: PhantomData<Entry>,
}

#[derive(Copy, Clone)]
struct ProteinStringId<S = DefaultSymbol> {
    string_id: S,
    source: S,
}

#[derive(Copy, Clone)]
struct ProteinAlias<S = DefaultSymbol> {
    alias: S,
    source: S,
}

impl<Entry> StringProteinAliases<Entry> {
    pub fn get_source_sym(&self, source: &str) -> Option<DefaultSymbol> {
        self.mapping.value_interner().1.get(source)
    }
}

#[allow(dead_code)]
impl StringProteinAliases<ProteinStringId> {
    pub fn read_tsv<R: io::Read>(reader: R) -> anyhow::Result<Self> {
        let mapping = InternedMapping::read_grouped_tsv_with::<SourceRepeaterHKT, _, _>(
            reader,
            false,
            |record| {
                // (alias, (string_id, source))
                Ok((
                    &record[1],
                    SourceRepeater {
                        value: &record[0],
                        source: record[2].split_whitespace(),
                    },
                ))
            },
        )?;

        Ok(StringProteinAliases {
            mapping,
            phantom: PhantomData,
        })
    }

    pub fn find_string_ids(&self, alias: &str) -> impl Iterator<Item = ProteinStringId> + '_ {
        self.mapping
            .lookup_syms(alias)
            .map(|(string_id, source)| ProteinStringId { string_id, source })
    }

    pub fn resolve_string_id(&self, string_id: ProteinStringId) -> Option<ProteinStringId<&str>> {
        self
            .mapping
            .value_interner()
            .resolve((string_id.string_id, string_id.source))
            .map(|(string_id, source)| ProteinStringId { string_id, source })
    }

    pub fn resolve_alias_name(&self, alias: DefaultSymbol) -> Option<&str> {
        self.mapping.key_interner().resolve(alias)
    }

    pub fn resolve_alias(&self, ProteinAlias { alias, source }: ProteinAlias) -> Option<ProteinAlias<&str>> {
        let alias = self.resolve_alias_name(alias)?;

        let source = self
            .mapping
            .value_interner()
            .1
            .resolve(source)?;

        Some(ProteinAlias { alias, source })
    }
}

struct SourceRepeater<'a> {
    value: &'a str,
    source: std::str::SplitWhitespace<'a>,
}

type SourceRepeaterHKT = HKT!(SourceRepeater<'_>);

impl<'a> Iterator for SourceRepeater<'a> {
    type Item = (&'a str, &'a str);

    fn next(&mut self) -> Option<Self::Item> {
        let source = self.source.next()?;
        Some((self.value, source))
    }
}

#[allow(dead_code)]
impl StringProteinAliases<ProteinAlias> {
    pub fn read_tsv<R: io::Read>(reader: R) -> anyhow::Result<Self> {
        let mapping = InternedMapping::read_grouped_tsv_with::<SourceRepeaterHKT, _, _>(
            reader,
            false,
            |record| {
                // (string_id, (alias, source))
                Ok((
                    &record[0],
                    SourceRepeater {
                        value: &record[1],
                        source: record[2].split_whitespace(),
                    },
                ))
            },
        )?;

        Ok(StringProteinAliases {
            mapping,
            phantom: PhantomData,
        })
    }

    pub fn find_aliases(&self, string_id: &str) -> impl Iterator<Item = ProteinAlias> + '_ {
        self.mapping
            .lookup_syms(string_id)
            .map(|(alias, source)| ProteinAlias { alias, source })
    }

    pub fn resolve_alias_name(&self, alias: DefaultSymbol) -> Option<&str> {
        self.mapping.value_interner().0.resolve(alias)
    }

    pub fn resolve_alias(&self, alias: ProteinAlias) -> Option<ProteinAlias<&str>> {
        self
            .mapping
            .value_interner()
            .resolve((alias.alias, alias.source))
            .map(|(alias, source)| ProteinAlias { alias, source })
    }

    pub fn resolve_string_id(&self, ProteinStringId { string_id, source }: ProteinStringId) -> Option<ProteinStringId<&str>> {
        let string_id = self.mapping.key_interner().resolve(string_id)?;

        let source = self
            .mapping
            .value_interner()
            .1
            .resolve(source)?;

        Some(ProteinStringId { string_id, source })
    }
}

#[derive(Args)]
pub struct StringVirusesInteractionsArgs {
    #[clap(flatten)]
    taxonomy: PreprocessedTaxonomyArgs,

    /// TANGO metagenomic assignment output (see `meteor tango-assign --help`)
    metagenomic_assignment: String,

    /// Does the presence of a higher-ranked node in the metagenomic assignment imply the presence
    /// of any of its descendants?
    #[clap(long)]
    include_descendants: bool,

    /// Path to the UniProt - Viral contig name mapping. See `meteor ppi match-viral-proteins`
    /// present.
    #[clap(long)]
    uniprot_virus_mapping: String,

    /// Directory for STRING-Viruses downloads (see `meteor fetch string-viruses --help`)
    #[clap(long)]
    stringdb_dir: String,

    /// Output path. The result is a tab-separated values with columns
    ///   protein_accession_a, protein_accession_b, combined_score, virus_name, virus_taxid, host_taxid, edge_kind
    /// where edge_kind is always one of virus-virus, virus-host, host-virus or host-host
    ///
    /// Use '-' to print to standard output.
    #[clap(short, long, default_value = "-")]
    output: String,
}


struct VirusAliasHitInfo<'a> {
    accession: &'a str,
    hit_virus_name: &'a str,
}

enum NodeKind<'a, F>
where
    F: FnOnce() -> Option<VirusAliasHitInfo<'a>>
{
    Host(&'a str),
    MaybeVirus(F),
}

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum EdgeKind {
    VirusVirus,
    VirusHost,
    HostVirus,
    HostHost,
}

#[derive(Serialize)]
pub struct PredictedInteraction<S = String> {
    pub protein_accession_a: S,
    pub protein_accession_b: S,
    pub combined_score: i32,
    pub virus_name: Option<S>,
    pub virus_taxid: Option<u32>,
    pub host_taxid: Option<u32>,
    pub edge_kind: EdgeKind,
}

struct PredictionContext<Tax> {
    taxonomy: Tax,
    lca_cache: LcaCache,
    assigned_taxids: NodeIdSet,
    uniprot_virus_mapping: ProteinVirusTaxidMapping,
    string_uniprot_mapping: StringProteinAliases<ProteinAlias>,
    preferred_source_sym: DefaultSymbol,
}

impl<Names: 'static> PredictionContext<NcbiTaxonomy<Names>> {
    fn inspect_node_from_string_id<'b, 'a: 'b>(
        &'a self,
        string_id: &'b str,
    ) -> Option<(
        NodeId,
        NodeKind<'a, impl FnOnce() -> Option<VirusAliasHitInfo<'a>> + 'b>,
    )> {
        let Some((taxid, _)) = fetch::string_viruses::parse_string_id(string_id) else {
            eprintln!("Warning: Could not parse StringDB identifier: {string_id}");
            return None;
        };

        let taxid = self.taxonomy.fixup_node(taxid)?;

        let mut aliases = self.string_uniprot_mapping.find_aliases(string_id);

        let result = if self.assigned_taxids.contains(taxid) {
            NodeKind::Host(
                self.string_uniprot_mapping.resolve_alias_name(
                    aliases
                        .find(|alias| alias.source == self.preferred_source_sym)?
                        .alias,
                )?,
            )
        } else {
            NodeKind::MaybeVirus(move || {
                let mut candidate = None;

                for alias in aliases {
                    let Some(accession) = self.string_uniprot_mapping.resolve_alias_name(alias.alias) else {
                        continue;
                    };

                    let is_preferred_source = alias.source == self.preferred_source_sym;

                    for (hit_virus_name, hit_taxid) in self.uniprot_virus_mapping.0.lookup(accession) {
                        let Some(hit_taxid) = self.taxonomy.fixup_node(hit_taxid) else {
                            eprintln!("UniProt BLAST search result taxid not found in the taxonomy: {hit_taxid}");
                            continue;
                        };

                        if self.taxonomy.are_independent(taxid, hit_taxid) {
                            match self.lca_cache.lca_node(taxid, hit_taxid).ok().and_then(|lca| self.taxonomy.find_rank_str(lca)) {
                                Some("genus") => {
                                },
                                r => {
                                    eprintln!("Warning: Virus TaxID mismatch between StringDB and UniProt BLAST search: {string_id} vs {hit_taxid} (LCA rank: {r:?}");
                                    eprintln!("STRING: {}", self.taxonomy.strict_ancestors(taxid).map(|n| n.to_string()).join(" - "));
                                    eprintln!("UniProt: {}", self.taxonomy.strict_ancestors(hit_taxid).map(|n| n.to_string()).join(" - "));
                                    continue;
                                }
                            }
                        }

                        if is_preferred_source {
                            candidate = Some(VirusAliasHitInfo {
                                accession,
                                hit_virus_name,
                            });
                            break;
                        } else if candidate.is_none() {
                            candidate = Some(VirusAliasHitInfo {
                                accession,
                                hit_virus_name,
                            });
                        }
                    }
                }

                candidate
            })
        };

        Some((taxid, result))
    }

    fn process_protein_links_record<'a>(
        &'a self,
        record: ProteinLinksRecord<&str>,
    ) -> Option<PredictedInteraction<&'a str>> {
        let (taxid_a, node_kind_a) = self.inspect_node_from_string_id(record.protein_id_a)?;
        let (taxid_b, node_kind_b) = self.inspect_node_from_string_id(record.protein_id_b)?;

        use EdgeKind::*;
        use NodeKind::*;

        match (node_kind_a, node_kind_b) {
            (Host(acc_a), Host(acc_b)) => {
                if taxid_a != taxid_b {
                    eprintln!(
                        "Warning: Host-host interaction TaxID mismatch: {} - {}",


                        record.protein_id_a, record.protein_id_b
                    );
                    None
                } else {
                    Some(PredictedInteraction {
                        protein_accession_a: acc_a,
                        protein_accession_b: acc_b,
                        combined_score: record.combined_score,
                        virus_name: None,
                        virus_taxid: None,
                        host_taxid: Some(taxid_a.into()),
                        edge_kind: HostHost,
                    })
                }
            }
            (Host(acc_a), MaybeVirus(best_hit_b)) => {
                let best_hit_b = best_hit_b()?;

                Some(PredictedInteraction {
                    protein_accession_a: acc_a,
                    protein_accession_b: best_hit_b.accession,
                    combined_score: record.combined_score,
                    virus_name: Some(best_hit_b.hit_virus_name),
                    virus_taxid: Some(taxid_b.0),
                    host_taxid: Some(taxid_a.0),
                    edge_kind: HostVirus,
                })
            }
            (MaybeVirus(best_hit_a), Host(acc_b)) => {
                let best_hit_a = best_hit_a()?;

                Some(PredictedInteraction {
                    protein_accession_a: best_hit_a.accession,
                    protein_accession_b: acc_b,
                    combined_score: record.combined_score,
                    virus_name: Some(best_hit_a.hit_virus_name),
                    virus_taxid: Some(taxid_a.0),
                    host_taxid: Some(taxid_b.0),
                    edge_kind: VirusHost,
                })
            }
            (MaybeVirus(best_hit_a), MaybeVirus(best_hit_b)) => {
                let best_hit_a = best_hit_a()?;
                let best_hit_b = best_hit_b()?;

                if taxid_a != taxid_b {
                    eprintln!(
                        "Warning: Virus-virus interaction TaxID mismatch: {} - {}",
                        record.protein_id_a, record.protein_id_b
                    );
                    None
                } else {
                    Some(PredictedInteraction {
                        protein_accession_a: best_hit_a.accession,
                        protein_accession_b: best_hit_b.accession,
                        combined_score: record.combined_score,
                        virus_name: Some(best_hit_a.hit_virus_name),
                        virus_taxid: Some(taxid_a.0),
                        host_taxid: None,
                        edge_kind: VirusVirus,
                    })
                }
            }
        }
    }
}


fn read_assignment_taxids<Tax: Taxonomy>(
    assignment: &Path,
    taxonomy: &Tax,
    include_descendants: bool,
) -> anyhow::Result<NodeIdSet> {
    let mut taxids = NodeIdSet::new();

    let mut csv_reader_iter = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(assignment)?
        .into_lending_iter()
        .into_deserialize::<HKT!(AssignmentRecord<'_>)>(None);

    while let Some(record) = csv_reader_iter.next() {
        let taxid = NodeId(record?.assigned_taxid);

        taxids.extend(taxonomy.strict_ancestors(taxid));

        if include_descendants {
            taxids.extend(taxonomy.postorder_descendants(taxid));
        } else {
            taxids.insert(taxid);
        }
    }

    Ok(taxids)
}

fn read_string_uniprot_mapping(stringdb_dir: &Path) -> anyhow::Result<StringProteinAliases<ProteinAlias>> {
    use fetch::string_viruses::FileName;

    maybe_gzdecoder!(FileName::ProteinAliases.open_read(stringdb_dir)?, rdr => {
        StringProteinAliases::<ProteinAlias>::read_tsv(rdr)
    })
}

fn read_uniprot_virus_mapping(path: &Path) -> anyhow::Result<ProteinVirusTaxidMapping> {
    let file = File::open(path)?;

    ProteinVirusTaxidMapping::read_tsv(file)
}

pub fn string_viruses_interactions(args: StringVirusesInteractionsArgs) -> anyhow::Result<()> {
    let stringdb_dir = args.stringdb_dir.as_ref();

    let context = {
        eprintln!("Reading the reference taxonomy");
        let taxonomy = with_some_ncbi_or_newick_taxonomy!(args.taxonomy.deserialize()?.tree,
            ncbi: taxonomy => taxonomy.without_names(),
            newick: _ => {
                anyhow::bail!("Newick taxonomies are not supported in this utility");
            },
        );

        let lca_cache = LcaCache::compute(&taxonomy);

        eprintln!("Reading the metagenomic assignment from {}", &args.metagenomic_assignment);
        let assigned_taxids = read_assignment_taxids(
            args.metagenomic_assignment.as_ref(),
            &taxonomy,
            args.include_descendants
        )?;

        eprintln!(
            "Reading the UniProt - viral contig name mapping from {}",
            &args.uniprot_virus_mapping
        );
        let uniprot_virus_mapping = read_uniprot_virus_mapping(args.uniprot_virus_mapping.as_ref())?;


        eprintln!(
            "Reading the StringDB - UniProt mapping in {}",
            &args.stringdb_dir
        );

        let string_uniprot_mapping = read_string_uniprot_mapping(stringdb_dir)?;

        let preferred_source_sym = string_uniprot_mapping
            .get_source_sym("UniProtKB-EI")
            .expect("Could not find StringDB alias source: UniProtKB-EI");

        PredictionContext {
            taxonomy,
            lca_cache,
            assigned_taxids,
            uniprot_virus_mapping,
            string_uniprot_mapping,
            preferred_source_sym,
        }
    };

    let context = &context;

    let (sender, receiver) = std::sync::mpsc::channel(); // unbounded

    eprintln!("Processing Viruses-StringDB network links...");

    let (r1, r2) = rayon::join(
        move || {
            use fetch::string_viruses::FileName::ProteinLinks;

            maybe_gzdecoder!(ProteinLinks.open_read(stringdb_dir)?, protein_links_reader => {
                let mut records = csv::ReaderBuilder::new()
                    .has_headers(true)
                    .delimiter(b' ')
                    .from_reader(protein_links_reader)
                    .into_lending_iter()
                    .into_deserialize::<HKT!(fetch::string_viruses::ProteinLinksRecord<&str>)>(None);

                while let Some(record) = records.next() {
                    if let Some(interaction) = context.process_protein_links_record(record?) {
                        if sender.send(interaction).is_err() {
                            eprintln!("Could not send interaction -- channel closed?");
                            return anyhow::Ok(());
                        };
                    }
                }
            });

            anyhow::Ok(())
        },
        move || {
            writing_new_file_or_stdout!(&args.output, writer => {
                let mut csv_writer = csv::WriterBuilder::new()
                    .delimiter(b'\t')
                    .has_headers(true)
                    .from_writer(writer?);

                for interaction in receiver {
                    if let Err(e) = csv_writer.serialize(interaction) {
                        if util::io::is_broken_pipe(&e) {
                            break;
                        } else {
                            return Err(e.into());
                        }
                    }
                }
            });

            anyhow::Ok(())
        },
    );

    r2?;
    r1?;

    Ok(())
}
