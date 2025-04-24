use clap::{Args, Subcommand};

pub mod ncbi_taxonomy;
#[cfg(feature = "ppin")]
pub mod string_viruses;
#[cfg(feature = "ppin")]
pub mod uniprot;

#[derive(Subcommand)]
enum Commands {
    NcbiTaxonomy(ncbi_taxonomy::NcbiTaxonomyFetchArgs),
    #[cfg(feature = "ppin")]
    StringViruses(string_viruses::StringVirusesFetchArgs),
    #[cfg(feature = "ppin")]
    Uniprot(uniprot::UniProtFetchArgs),
}

/// Download external files to be used by Meteor.
#[derive(Args)]
pub struct FetchArgs {
    #[clap(subcommand)]
    command: Commands,
}

pub fn fetch(args: FetchArgs) -> anyhow::Result<()> {
    match args.command {
        Commands::NcbiTaxonomy(args) => {
            ncbi_taxonomy::fetch(&args)?;
        }
        #[cfg(feature = "ppin")]
        Commands::StringViruses(args) => {
            string_viruses::fetch(&args)?;
        }
        #[cfg(feature = "ppin")]
        Commands::Uniprot(args) => {
            uniprot::fetch(&args)?;
        }
    }

    Ok(())
}

