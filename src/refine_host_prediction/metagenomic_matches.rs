use std::fs::File;

use crate::{
    crispr_match::HostVirusMapping,
    util::{
        clap::ExistingFilePath,
        interners::{BorrowedInterner, Interner, PlaceholderInterner, PlaceholderInternerFor},
    },
};
use clap::Args;
use lending_iterator::higher_kinded_types::HKTRef;
use string_interner::{DefaultStringInterner, DefaultSymbol};


#[derive(Args)]
pub struct MergeCrisprMatchArgs {
    /// CRISPR match data (see `meteor crispr-match`)
    #[clap(long, help_heading = "CRISPR", required = false)]
    pub crispr_match: ExistingFilePath,
}

pub struct MetagenomicMatches {
    pub crispr: Option<
        HostVirusMapping<
            PlaceholderInternerFor<DefaultStringInterner>,
            PlaceholderInternerFor<DefaultStringInterner>,
        >,
    >,
}

impl MetagenomicMatches {
    pub fn new() -> Self {
        Self { crispr: None }
    }

    pub fn load_crispr_matches_using<HI, VI>(
        &mut self,
        args: &MergeCrisprMatchArgs,
        host_interner: &mut HI,
        virus_interner: &mut VI,
    ) -> anyhow::Result<()>
    where
        HI: Interner<Symbol = DefaultSymbol, ValueHKT = HKTRef<str>>,
        VI: Interner<Symbol = DefaultSymbol, ValueHKT = HKTRef<str>>,
    {
        let crispr = HostVirusMapping::read_tsv_using(
            BorrowedInterner::new(host_interner),
            BorrowedInterner::new(virus_interner),
            File::open(args.crispr_match.as_path())?,
        )?
        .replace_interners(|_, _| (PlaceholderInterner::new(), PlaceholderInterner::new()));

        self.crispr = Some(crispr);
        Ok(())
    }
}
