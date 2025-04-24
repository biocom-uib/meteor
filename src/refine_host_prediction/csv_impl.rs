use std::sync::{mpsc::Receiver, Arc};
use std::io::Write;

use csv::StringRecord;
use lending_iterator::higher_kinded_types::HKTRef;
use itertools::{Either, Itertools};
use rayon::iter::{ParallelBridge, ParallelIterator};
use string_interner::DefaultSymbol;

use crate::taxonomy::LabeledTaxonomy;
use crate::{util, writing_new_file_or_stdout};
use crate::{taxonomy::{NodeId, tree::walk::RootedTreeWalk}, util::interners::Interned};

use super::{Aggregation, HostSymbol, MetagenomicEvidence, RefineHostPredictionArgs, RefinementOptions, VirusSymbol};

struct HostPredictionRecord {
    record: StringRecord,
    virus_name_idx: usize,
    merged_taxids_idx: usize,
    merged_class_names_idx: Option<usize>,
}

impl HostPredictionRecord {
    fn virus_name(&self) -> &str {
        &self.record[self.virus_name_idx]
    }

    fn predicted_taxids(&self) -> impl Iterator<Item = Option<u32>> + '_ {
        self.record[self.merged_taxids_idx].split(';').map(|taxid| taxid.parse().ok())
    }

    fn predicted_class_names(&self) -> Option<impl Iterator<Item = &str> + '_> {
        Some(self.record[self.merged_class_names_idx?].split(';'))
    }

    fn predictions(&self) -> impl Iterator<Item = (Option<u32>, Option<&str>)> + '_ {
        let taxids = self.predicted_taxids();

        let class_names = if let Some(class_names) = self.predicted_class_names() {
            itertools::Either::Left(class_names.map(Some))
        } else {
            itertools::Either::Right(std::iter::repeat(None))
        };

        taxids.zip(class_names)
    }
}

struct RefinedRecord {
    virus_name: String,
    predicted_taxid: NodeId,
    predicted_class_name: Option<String>,
    assigned_contigs: String,
    crispr_matches: Option<String>
}

impl RefinedRecord {
    fn write_header<W: Write>(
        has_class_names: bool,
        has_crispr: bool,
        writer: &mut csv::Writer<W>,
    ) -> csv::Result<()> {

        writer.write_field("virus_name")?;
        writer.write_field("predicted_taxid")?;

        if has_class_names {
            writer.write_field("predicted_class_name")?;
        }

        writer.write_field("assigned_contigs")?;

        if has_crispr {
            writer.write_field("crispr_matches")?;
        }

        writer.write_record(None::<&[u8]>)?;

        Ok(())
    }

    fn write<W: Write>(&self, writer: &mut csv::Writer<W>) -> csv::Result<()> {
        writer.write_field(&self.virus_name)?;
        writer.write_field(self.predicted_taxid.to_string())?;

        if let Some(predicted_class_name) = &self.predicted_class_name {
            writer.write_field(predicted_class_name)?;
        }

        writer.write_field(&self.assigned_contigs)?;

        if let Some(crispr_matches) = &self.crispr_matches {
            writer.write_field(crispr_matches)?;
        }

        writer.write_record(None::<&[u8]>)?;

        Ok(())
    }
}

fn load_host_prediction(
    args: &RefineHostPredictionArgs,
) -> anyhow::Result<(
    bool, /* has_class_names */
    impl Iterator<Item = csv::Result<HostPredictionRecord>>,
)> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(args.host_prediction.as_path())?;

    let headers = reader.headers()?;

    use crate::merge_outputs::host_prediction::columns;

    let Some(virus_name_idx) = headers.iter().position(|name| name == columns::VIRUS_NAME) else {
        anyhow::bail!(
            "Could not find column '{}' in the host prediction",
            columns::VIRUS_NAME
        )
    };

    let Some(merged_taxids_idx) = headers
        .iter()
        .position(|name| name == columns::MERGED_TAXIDS)
    else {
        anyhow::bail!(
            "Could not find column '{}' in the host prediction",
            columns::MERGED_TAXIDS
        )
    };

    let mut merged_class_names_idx = headers
        .iter()
        .position(|name| name == columns::MERGED_CLASS_NAMES);

    if !args.refinement_options.keep_names {
        merged_class_names_idx = None;
    }

    let records = reader.into_records().map(move |result| {
        result.map(move |record| HostPredictionRecord {
            record,
            virus_name_idx,
            merged_taxids_idx,
            merged_class_names_idx,
        })
    });

    Ok((merged_class_names_idx.is_some(), records))
}

fn agg_node_ids<'a, IN>(agg: Aggregation, nodes: IN) -> String
where
    IN: ExactSizeIterator<Item = &'a str>,
{
    match agg {
        Aggregation::List => nodes.into_iter().join(";"),
        Aggregation::Count => nodes.into_iter().len().to_string(),
    }
}

fn record_predictions<'a, Tax>(
    tax: &'a Tax,
    rank_sym: Option<Tax::RankSym>,
    record: &'a HostPredictionRecord,
) -> impl Iterator<Item = (NodeId, Option<&'a str>)> + 'a
where
    Tax: LabeledTaxonomy
{
    let predictions = record.predictions().filter_map(move |(taxid, class_name)| {
        let taxid = taxid?;

        let Some(node) = tax.fixup_node(taxid) else {
            eprintln!("Could not locate node {taxid} in the taxonomy");
            return None;
        };

        Some((node, class_name))
    });

    let Some(rank_sym) = rank_sym else {
        return Either::Left(predictions);
    };

    let mut predictions_vec = predictions.filter_map(move |(node, class_name)| {
        if tax.find_rank(node) == Some(rank_sym) {
            Some((node, class_name))
        } else {
            let ancestor_at_rank = tax
                .strict_ancestors(node)
                .find(move |ancestor| tax.find_rank(*ancestor) == Some(rank_sym))?;

            let ancestor_class_name = if class_name.is_some() {
                tax.some_label_of(ancestor_at_rank)
            } else {
                None
            };

            Some((ancestor_at_rank, ancestor_class_name))
        }
    })
    .collect_vec();

    predictions_vec.sort_by(|x, y| x.0.cmp(&y.0));
    predictions_vec.dedup_by(|x, y| x.0 == y.0);

    Either::Right(predictions_vec.into_iter())
}


fn refine_record<'a, Tax, VI, HI>(
    evidence: &'a MetagenomicEvidence<VI, HI>,
    opts: &'a RefinementOptions,
    tax: &'a Tax,
    rank_sym: Option<Tax::RankSym>,
    record: &'a HostPredictionRecord,
) -> impl Iterator<Item = RefinedRecord> + 'a
where
    Tax: LabeledTaxonomy,
    VI: Interned<Symbol = DefaultSymbol, ValueHKT = HKTRef<str>> + Send + Sync,
    HI: Interned<Symbol = DefaultSymbol, ValueHKT = HKTRef<str>> + Send + Sync,
{
    let virus_name = record.virus_name();
    let virus_sym = evidence.virus_interner.get(virus_name).map(VirusSymbol);

    record_predictions(tax, rank_sym, record).filter_map(move |(node, class_name)| {
        let pred_evidence = evidence.refine_and_enrich(tax, virus_sym, node)?;

        let contigs = pred_evidence
            .assigned_contigs
            .iter()
            .map(|&HostSymbol(sym)| evidence.host_interner.resolve(sym).unwrap());

        let crisprs = pred_evidence
            .crispr_matches
            .iter()
            .map(|&HostSymbol(sym)| evidence.host_interner.resolve(sym).unwrap());

        Some(RefinedRecord {
            virus_name: virus_name.to_owned(),
            predicted_taxid: node,
            predicted_class_name: class_name.map(str::to_owned),
            assigned_contigs: agg_node_ids(opts.agg_contigs, contigs),
            crispr_matches: evidence
                .has_crispr
                .then(|| agg_node_ids(opts.agg_crispr, crisprs)),
        })
    })
}

fn write_records(
    has_class_names: bool,
    has_crispr: bool,
    receiver: Receiver<Vec<RefinedRecord>>,
    writer: impl Write,
) -> anyhow::Result<()> {
    let mut writer = csv::WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_writer(writer);

    RefinedRecord::write_header(has_class_names, has_crispr, &mut writer)?;

    let mut count = 0;

    for records in receiver.iter() {
        for record in records {
            record.write(&mut writer)?;
        }

        count += 1;

        if count % 100 == 0 {
            eprint!("\rProcessed {count} viral contigs");
        }
    }

    eprintln!("\rProcessed {count} viral contigs");

    writer.flush()?;

    Ok(())
}

pub fn refine_host_prediction_with_tax_impl<Tax, VI, HI>(
    args: &RefineHostPredictionArgs,
    tax: Arc<Tax>,
    evidence: MetagenomicEvidence<VI, HI>,
) -> anyhow::Result<()>
where
    Tax: LabeledTaxonomy + Send + Sync,
    VI: Interned<Symbol = DefaultSymbol, ValueHKT = HKTRef<str>> + Send + Sync,
    HI: Interned<Symbol = DefaultSymbol, ValueHKT = HKTRef<str>> + Send + Sync,
{
    let (sender, receiver) = std::sync::mpsc::channel(); // unbounded channel

    let evidence = &evidence;
    let has_crispr = evidence.has_crispr;

    let (has_class_names, host_prediction) = load_host_prediction(args)?;

    let rank_sym =
        if let Some(rank) = &args.refinement_options.rank {
            Some(tax.lookup_rank_sym(rank).ok_or_else(|| {
                anyhow::anyhow!("Rank not found in the taxonomy: {}", rank)
            })?)
        } else {
            None
        };

    let (rp, rc) = rayon::join(
        move || {
            host_prediction.par_bridge().try_for_each_with(
                sender,
                move |sender, record| {
                    let record = record?;
                    let records = refine_record(evidence, &args.refinement_options, &*tax, rank_sym, &record)
                        .collect();
                    sender.send(records)?;
                    anyhow::Ok(())
                },
            )
        },
        move || {
            writing_new_file_or_stdout!(args.output.as_ref(), writer => {
                write_records(has_class_names, has_crispr, receiver, writer?)
            })
        },
    );

    util::ignore_broken_pipe_anyhow(rc)?;
    rp?;

    Ok(())
}
