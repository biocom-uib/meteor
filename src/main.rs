#![feature(
    associated_type_defaults,
    cmp_minmax,
    extend_one,
    hash_raw_entry,
    hash_set_entry,
    impl_trait_in_assoc_type,
    io_error_more,
    iter_next_chunk,
    iterator_try_reduce,
    once_cell_try,
    never_type,
    result_flattening,
    try_blocks,
    type_alias_impl_trait,
)]

use clap::{builder::styling, Parser, Subcommand};

extern crate self as meteor;

pub mod taxonomy;
pub(crate) mod preprocessed_taxonomy;

pub(crate) mod util;
pub(crate) mod csv;

#[cfg(feature = "cache")]
mod cache;
mod crispr_match;
mod fetch;
#[cfg(feature = "ppin")]
mod ppin;
mod get_lineage;
mod merge_outputs;
mod preprocess_blastout;
mod preprocess_taxonomy;
mod refine_host_prediction;
mod refine_vpf_class;
mod tango_assign;
mod tool;

#[derive(Subcommand)]
enum Commands {
    Fetch(fetch::FetchArgs),
    PreprocessTaxonomy(preprocess_taxonomy::PreprocessTaxonomyArgs),
    PreprocessBlastout(preprocess_blastout::PreprocessBlastOutArgs),
    TangoAssign(tango_assign::TangoAssignArgs),
    MergeHostPrediction(merge_outputs::host_prediction::MergeHostPredictionArgs),
    MergeVirusTaxonomy(merge_outputs::virus_taxonomy::MergeVirusTaxonomyArgs),
    RefineHostPrediction(refine_host_prediction::RefineHostPredictionArgs),
    RefineVpfClass(refine_vpf_class::RefineVpfClassArgs),
    CrisprMatch(crispr_match::CrisprMatchArgs),
    GetLineage(get_lineage::GetLineageArgs),
    #[cfg(feature = "ppin")]
    Ppin(ppin::PpinArgs),
}

const CLI_STYLES: styling::Styles = styling::Styles::styled()
    .header(styling::AnsiColor::BrightWhite.on_default().bold().underline())
    .usage(styling::AnsiColor::BrightWhite.on_default().bold().underline())
    .literal(styling::AnsiColor::BrightYellow.on_default().bold())
    .placeholder(styling::AnsiColor::Red.on_default())
    .error(styling::AnsiColor::BrightRed.on_default().bold())
    .valid(styling::AnsiColor::BrightBlue.on_default().bold())
    .invalid(styling::AnsiColor::BrightBlue.on_default().bold());

/// METEOR: Metagenome and Metavirome Joint Analysis
#[derive(Parser)]
#[command(author, version, about, styles = CLI_STYLES)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();

    match args.command {
        Commands::Fetch(args) => {
            fetch::fetch(args)?;
        }
        Commands::PreprocessTaxonomy(args) => {
            preprocess_taxonomy::preprocess_taxonomy(args)?;
        }
        Commands::GetLineage(args) => {
            get_lineage::get_lineage(args)?;
        }
        Commands::PreprocessBlastout(args) => {
            preprocess_blastout::preprocess_blastout(args)?;
        }
        Commands::TangoAssign(args) => {
            tango_assign::tango_assign(args)?;
        }
        Commands::MergeHostPrediction(args) => {
            merge_outputs::host_prediction::merge_host_prediction(args)?;
        }
        Commands::MergeVirusTaxonomy(args) => {
            merge_outputs::virus_taxonomy::merge_virus_taxonomy(args)?;
        }
        Commands::RefineHostPrediction(args) => {
            refine_host_prediction::refine_host_prediction(args)?;
        }
        Commands::RefineVpfClass(args) => {
            refine_vpf_class::refine_vpf_class(args)?;
        }
        #[cfg(feature = "ppin")]
        Commands::Ppin(args) => {
            ppin::ppi(args)?;
        }
        Commands::CrisprMatch(args) => {
            crispr_match::crispr_match(args)?;
        }
    }

    Ok(())
}
