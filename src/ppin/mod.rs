use std::io;

use clap::{Args, Subcommand};
use serde::{Serialize, Deserialize};
use string_interner::DefaultStringInterner;

use crate::util::{interned_mapping::InternedMapping, interners::{PairInterner, TrivialInterner}};

mod uniprot_xml_parser;
mod prodigal_uniprot_pipeline;
mod prepare_uniprot_blastdb;
mod string_viruses;


#[derive(Serialize, Deserialize)]
pub struct VirusProteinTaxidRecord<S = String> {
    uniprot_accession: S,
    virus_name: S,
    virus_taxid: u32,
}

pub struct ProteinVirusTaxidMapping(
    InternedMapping<
        DefaultStringInterner,
        PairInterner<DefaultStringInterner, TrivialInterner<u32>>,
    >,
);

impl ProteinVirusTaxidMapping {
    pub fn read_tsv<R: io::Read>(reader: R) -> anyhow::Result<Self> {
        let mapping = InternedMapping::read_tsv_with(reader, true, |record| {
            Ok((&record[0], (&record[1], record[2].parse()?)))
        })?;

        Ok(Self(mapping))
    }

    pub fn write_tsv<W: io::Write>(&self, writer: W) -> csv::Result<()> {
        self.0.write_tsv_with(
            writer,
            &["uniprot_accession", "virus_name", "virus_taxid"],
            |csv_writer, uniprot_accession, (virus_name, virus_taxid)| {
                csv_writer.serialize(VirusProteinTaxidRecord {
                    uniprot_accession,
                    virus_name,
                    virus_taxid,
                })
            },
        )
    }
}

#[derive(Subcommand)]
enum Commands {
    PrepareUniprotBlastdb(prepare_uniprot_blastdb::PrepareUniProtBlastDBArgs),
    MatchViralProteins(prodigal_uniprot_pipeline::MatchViralProteinsArgs),
    StringVirusesInteractions(string_viruses::StringVirusesInteractionsArgs),
}

/// Commands that can be used to obtain protein-protein interaction networks.
/// The subcommands to use to generate a virus-host PPI network are the following (in order):
/// - meteor ppi prepare-uniprot-blastdb ...
/// - meteor ppi match-viral-proteins ...
/// - meteor ppi string-viruses-interactions ...
#[derive(Args)]
pub struct PpinArgs {
    #[clap(subcommand)]
    command: Commands,
}

pub fn ppi(args: PpinArgs) -> anyhow::Result<()> {
    match args.command {
        Commands::PrepareUniprotBlastdb(args) => {
            prepare_uniprot_blastdb::prepare_uniprot_blastdb(args)?;
        }
        Commands::MatchViralProteins(args) => {
            prodigal_uniprot_pipeline::match_viral_proteins(args)?;
        }
        Commands::StringVirusesInteractions(args) => {
            string_viruses::string_viruses_interactions(args)?;
        }
    }

    Ok(())
}
