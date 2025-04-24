use std::sync::Arc;

use anyhow::Context;
use clap::Args;
use polars::{
    datatypes::DataType,
    io::SerWriter,
    prelude::{CsvWriter, LazyFrame},
    series::Series,
};

use crate::{
    csv::{
        filter::Filter,
        polars::{lazyframe_from_file, ExprExt},
    },
    merge_outputs::NamedLazyFrame,
    preprocessed_taxonomy::{with_some_taxonomy, PreprocessedTaxonomyArgs},
    taxonomy::LabeledTaxonomy,
    tool::{
        phabox::{self, CherryPredictionRecord, CherryPredictionRecordFilter},
        vpf_class::{VpfClassHostRank, VpfClassRecordFilter},
    },
    util::{
        self,
        clap::{ExistingDirPath, OutputFileFlag},
        writing_new_file_or_stdout,
    },
};

use super::{MergeCommonArgs, MergeOutputContext};

pub use super::columns;

#[derive(Args)]
struct MergeVpfClassArgs {
    /// Path to the VPF-Class output directory
    #[clap(long, help_heading = "VPF-Class", required = false)]
    vpf_class: ExistingDirPath,

    /// Filters to apply to the VPF-Class output
    #[clap(long, help_heading = "VPF-Class")]
    vpf_class_filter: Vec<VpfClassRecordFilter>,

    /// Ranks to use from VPF-Class (family and/or genus)
    #[clap(long, value_enum, default_values_t = [VpfClassHostRank::Genus], value_delimiter = ',', help_heading = "VPF-Class")]
    vpf_class_ranks: Vec<VpfClassHostRank>,
}

#[derive(Args)]
pub struct MergeCherryArgs {
    /// Path to the PhaBOX output directory
    #[clap(long, help_heading = "PhaBOX-CHERRY", required = false)]
    cherry: ExistingDirPath,

    /// Filters to apply to the CHERRY output
    #[clap(long, help_heading = "PhaBOX-CHERRY")]
    cherry_filter: Vec<CherryPredictionRecordFilter>,
}

/// Combine host predictions from multiple tools
#[derive(Args)]
pub struct MergeHostPredictionArgs {
    #[clap(flatten)]
    taxonomy: PreprocessedTaxonomyArgs,

    #[clap(flatten)]
    output: OutputFileFlag,

    #[clap(flatten)]
    merge_args: MergeCommonArgs,

    #[clap(flatten)]
    vpf_class_args: Option<MergeVpfClassArgs>,

    #[clap(flatten)]
    cherry_args: Option<MergeCherryArgs>,
}

fn load_vpf_class<Tax>(
    args: &MergeVpfClassArgs,
    ctx: &MergeOutputContext<Tax>,
) -> anyhow::Result<LazyFrame>
where
    Tax: LabeledTaxonomy + Sync + Send + 'static,
{
    let output_dir = args.vpf_class.as_ref();

    let ranks = &args.vpf_class_ranks;

    let filters = &args.vpf_class_filter;

    super::load_vpf_class_for_merge(ctx, output_dir, ranks, filters, true)
}

fn load_cherry<Tax>(
    args: &MergeCherryArgs,
    ctx: &MergeOutputContext<Tax>,
) -> anyhow::Result<LazyFrame>
where
    Tax: LabeledTaxonomy + Sync + Send + 'static,
{
    use polars::lazy::dsl::{col, lit};

    let output_dir = args.cherry.as_ref();

    let output_file = phabox::locate_cherry_prediction(output_dir)?;

    let tax = ctx.taxonomy();

    let mut df = lazyframe_from_file::<CherryPredictionRecord>(&output_file)?;

    df = Filter::apply_polars(df, &args.cherry_filter)
        .filter(col("Pred").neq(lit("filtered")))
        .select([
            col("Accession").alias(columns::VIRUS_NAME),
            col("Pred").alias(columns::CLASS_NAMES),
        ])
        .with_column(
            col(columns::CLASS_NAMES)
                .downcast_map_apply(Series::str, move |class_name| {
                    phabox::pred_to_taxid(&*tax, class_name)
                })
                .cast(DataType::List(DataType::UInt32.boxed()))
                .alias(columns::TAXIDS),
        );

    if ctx.args().keep_names {
        df = df.with_column(col(columns::CLASS_NAMES).cast(DataType::List(DataType::String.boxed())));
    } else {
        df = df.drop([col(columns::CLASS_NAMES)]);
    }

    Ok(df)
}

fn merge_host_prediction_with_tax<Tax>(
    args: MergeHostPredictionArgs,
    tax: Tax,
) -> anyhow::Result<()>
where
    Tax: LabeledTaxonomy + Send + Sync + 'static,
{
    let tax = Arc::new(tax);

    let ctx = MergeOutputContext::new(&args.merge_args, tax);

    let mut predictions = Vec::new();

    if let Some(vpf_class_args) = &args.vpf_class_args {
        eprintln!("Loading VPF-Class output");
        let vpf_class = load_vpf_class(vpf_class_args, &ctx)?;
        predictions.push(NamedLazyFrame::new("vpf_class", vpf_class));
    }

    if let Some(cherry_args) = &args.cherry_args {
        eprintln!("Loading PhaBOX-CHERRY output");
        let cherry = load_cherry(cherry_args, &ctx)?;
        predictions.push(NamedLazyFrame::new("cherry", cherry));
    }

    let merged = super::merge_virus_predictions(&ctx, predictions);

    let merged = super::prepare_merged_virus_predictions_for_writing(merged);

    let mut merged = merged.collect()?;

    eprintln!("Loaded & merged data");

    writing_new_file_or_stdout!(args.output.as_ref(), writer => {
        let writer = writer.context("Error creating merged taxonomy prediction")?;

        eprintln!("Writing results...");

        util::ignore_broken_pipe(
            CsvWriter::new(writer)
                .include_header(true)
                .with_separator(b'\t')
                .finish(&mut merged)
        )?;
    });

    Ok(())
}

pub fn merge_host_prediction(args: MergeHostPredictionArgs) -> anyhow::Result<()> {
    eprintln!("Loading taxonomy");

    let tax = args.taxonomy.deserialize()?.tree;

    with_some_taxonomy!(tax, tax => {
        merge_host_prediction_with_tax(args, tax)?;
    });

    Ok(())
}
