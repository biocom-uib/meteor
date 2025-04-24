use std::sync::Arc;

use anyhow::Context;
use clap::Args;
use polars::io::SerWriter;
use polars::prelude::{CsvWriter, LazyFrame};
use polars::series::Series;

use crate::csv::filter::Filter;
use crate::csv::polars::{lazyframe_from_file, ExprExt};
use crate::merge_outputs::{MergeOutputContext, NamedLazyFrame};
use crate::preprocessed_taxonomy::{with_some_taxonomy, PreprocessedTaxonomyArgs};
use crate::taxonomy::LabeledTaxonomy;
use crate::tool::genomad::{self, GenomadVirusRecordFilter};
use crate::tool::phagenus::{self, PhaGenusRecord, PhaGenusRecordFilter};
use crate::tool::vpf_class::{VpfClassRecordFilter, VpfClassVirusRank};
use crate::util::clap::{ExistingDirPath, ExistingFilePath, OutputFileFlag};
use crate::util::{self, writing_new_file_or_stdout};

use super::MergeCommonArgs;

pub use super::columns;

#[derive(Args)]
struct MergeVpfClassArgs {
    /// Path to the VPF-Class output directory
    #[clap(long, help_heading = "VPF-Class", required = false)]
    vpf_class: ExistingDirPath,

    /// Filters to apply to the VPF-Class output
    #[clap(long, help_heading = "VPF-Class", requires = "vpf_class")]
    vpf_class_filter: Vec<VpfClassRecordFilter>,

    /// Ranks to use from VPF-Class (family and/or genus)
    #[clap(long, help_heading = "VPF-Class", requires = "vpf_class", value_enum, default_values_t = [VpfClassVirusRank::Genus], value_delimiter = ',')]
    vpf_class_ranks: Vec<VpfClassVirusRank>,
}

#[derive(Args)]
struct MergeGenomadArgs {
    /// Path to the geNomad output directory
    #[clap(long, help_heading = "geNomad", required = false)]
    genomad: ExistingDirPath,

    /// Filters to apply to the geNomad output
    #[clap(long, help_heading = "geNomad", requires = "genomad")]
    genomad_filter: Option<GenomadVirusRecordFilter>,
}

#[derive(Args)]
#[group(required = false)]
pub struct MergePhaGenusArgs {
    /// Path to the PhaGenus output csv file
    #[clap(long, help_heading = "PhaGenus", required = false)]
    phagenus: ExistingFilePath,

    /// Filters to apply to the PhaGenus output
    #[clap(long, help_heading = "PhaGenus", requires = "phagenus")]
    phagenus_filter: Vec<PhaGenusRecordFilter>,
}

/// Combine virus taxonomy annotations from multiple tools
#[derive(Args)]
pub struct MergeVirusTaxonomyArgs {
    #[clap(flatten)]
    taxonomy: PreprocessedTaxonomyArgs,

    #[clap(flatten)]
    output: OutputFileFlag,

    #[clap(flatten)]
    merge_args: MergeCommonArgs,

    #[clap(flatten)]
    vpf_class_args: Option<MergeVpfClassArgs>,

    #[clap(flatten)]
    genomad_args: Option<MergeGenomadArgs>,

    #[clap(flatten)]
    phagenus_args: Option<MergePhaGenusArgs>,
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

fn load_genomad<Tax>(
    args: &MergeGenomadArgs,
    ctx: &MergeOutputContext<Tax>,
) -> anyhow::Result<LazyFrame>
where
    Tax: LabeledTaxonomy + Sync + Send + 'static,
{
    use polars::lazy::dsl::{col, lit};

    let output_dir = args.genomad.as_ref();

    let mut df = genomad::virus_summary_lazyframe_from_dir(output_dir)?;

    let tax = ctx.taxonomy();

    df = Filter::apply_polars(df, &args.genomad_filter)
        .filter(col("topology").neq(lit("Provirus")))
        .select([col("seq_name").alias(columns::VIRUS_NAME), col("taxonomy")])
        .with_column(
            col("taxonomy")
                .downcast_map_apply(Series::str, move |lineage| {
                    genomad::lineage_to_taxid(&*tax, lineage)
                })
                .alias(columns::TAXIDS),
        );

    df = if ctx.args.keep_names {
        df.with_column(
            col("taxonomy")
                .downcast_map_apply_str(genomad::lineage_most_specific)
        )
        .rename(["taxonomy"], [columns::CLASS_NAMES], true)
    } else {
        df.drop(["taxonomy"])
    };

    df = df.group_by([columns::VIRUS_NAME]).agg([col("*")]);

    Ok(df)
}

fn load_phagenus<Tax>(
    args: &MergePhaGenusArgs,
    ctx: &MergeOutputContext<Tax>,
) -> anyhow::Result<LazyFrame>
where
    Tax: LabeledTaxonomy + Sync + Send + 'static,
{
    use polars::lazy::dsl::col;

    let output_file = args.phagenus.as_ref();

    let rank_sym = ctx
        .tax
        .lookup_rank_sym("genus")
        .ok_or_else(|| anyhow::anyhow!("Could not find rank \"genus\" in the taxonomy"))?;

    let mut df = lazyframe_from_file::<PhaGenusRecord>(output_file)?;

    df = Filter::apply_polars(df, &args.phagenus_filter).select([
        col("contigs").alias(columns::VIRUS_NAME),
        col("predict_label").alias(columns::CLASS_NAMES),
    ]);

    // {CLASS_NAMES}: str (name1/name2...)

    df = phagenus::explode_predict_label(df, columns::CLASS_NAMES).with_column(
        super::class_names_to_taxid_list(col(columns::CLASS_NAMES), ctx.taxonomy(), Some(rank_sym))
            .alias(columns::TAXIDS),
    );

    // {CLASS_NAMES}: str (name_i), {TAXIDS}: list[u32]

    if ctx.args().keep_names {
        df = df
            .with_column(super::repeat_class_names_for_taxids(
                col(columns::CLASS_NAMES),
                col(columns::TAXIDS),
            ))
            .group_by([columns::VIRUS_NAME])
            .agg([col("*").explode()]);
    } else {
        df = df.drop([columns::CLASS_NAMES]).explode([columns::TAXIDS]);
    }

    Ok(df)
}

fn merge_virus_taxonomy_with_tax<Tax>(args: MergeVirusTaxonomyArgs, tax: Tax) -> anyhow::Result<()>
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
    };

    if let Some(genomad_args) = &args.genomad_args {
        eprintln!("Loading geNomad output");
        let genomad = load_genomad(genomad_args, &ctx)?;
        predictions.push(NamedLazyFrame::new("genomad", genomad));
    };

    if let Some(phagenus_args) = &args.phagenus_args {
        eprintln!("Loading PhaGenus output");
        let phagenus = load_phagenus(phagenus_args, &ctx)?;
        predictions.push(NamedLazyFrame::new("phagenus", phagenus));
    };

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

pub fn merge_virus_taxonomy(args: MergeVirusTaxonomyArgs) -> anyhow::Result<()> {
    eprintln!("Loading taxonomy");

    let tax = args.taxonomy.deserialize()?.tree;

    with_some_taxonomy!(tax, tax => {
        merge_virus_taxonomy_with_tax(args, tax)?;
    });

    Ok(())
}
