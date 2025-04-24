use std::{borrow::Cow, io, path::Path};

use anyhow::Context;
use clap::Args;
use lending_iterator::higher_kinded_types::HKTRef;
use serde::{Deserialize, Serialize};
use string_interner::DefaultStringInterner;
use tempdir::TempDir;

use crate::{
    crispr_match::minced_spacers::MincedSpacersPipeline,
    util::{self, cli_tools, interned_mapping::InternedMapping, interners::{Interner, Resolver}, writing_new_file_or_stdout},
};


const DEFAULT_PERC_IDENTITY: i32 = 95;

#[derive(Default)]
pub struct VirusHostMapping<
    VI: Resolver = DefaultStringInterner,
    HI: Resolver = VI,
>(pub InternedMapping<VI, HI>);

impl<VI, HI> VirusHostMapping<VI, HI>
where
    VI: Interner<ValueHKT = HKTRef<str>>,
    HI: Interner<ValueHKT = HKTRef<str>>,
{
    pub const HEADER: [&'static str; 2] = ["virus_name", "host_name"];

    pub fn read_tsv_using<R: io::Read>(virus_interner: VI, host_interner: HI, reader: R) -> anyhow::Result<Self> {
        let mut mapping = InternedMapping::new(virus_interner, host_interner);
        mapping.extend_from_tsv(reader, true)?;
        Ok(Self(mapping))
    }

    pub fn write_tsv<W: io::Write>(&self, writer: W) -> csv::Result<()> {
        self.0.write_tsv(writer, &Self::HEADER)
    }
}

impl<VI, HI> VirusHostMapping<VI, HI>
where
    VI: Default + Interner<ValueHKT = HKTRef<str>>,
    HI: Default + Interner<ValueHKT = HKTRef<str>>,
{
    pub fn read_tsv<R: io::Read>(reader: R) -> anyhow::Result<Self> {
        Self::read_tsv_using(VI::default(), HI::default(), reader)
    }
}

impl<VI: Resolver, HI: Resolver> VirusHostMapping<VI, HI> {
    pub fn replace_interners<F, VI2, HI2>(self, f: F) -> VirusHostMapping<VI2, HI2>
    where
        VI2: Resolver<Symbol = VI::Symbol>,
        HI2: Resolver<Symbol = HI::Symbol>,
        F: FnOnce(VI, HI) -> (VI2, HI2)
    {
        VirusHostMapping(self.0.replace_interners(f))
    }
}

#[derive(Default)]
pub struct HostVirusMapping<HI: Resolver, VI: Resolver = HI>(pub InternedMapping<HI, VI>);

impl<HI, VI> HostVirusMapping<HI, VI>
where
    HI: Interner<ValueHKT = HKTRef<str>>,
    VI: Interner<ValueHKT = HKTRef<str>>,
{
    pub const HEADER: [&'static str; 2] = ["host_name", "virus_name"];

    pub fn read_tsv_using<R: io::Read>(host_interner: HI, virus_interner: VI, reader: R) -> anyhow::Result<Self> {
        let mut mapping = InternedMapping::new(host_interner, virus_interner);

        mapping.extend_from_tsv_with(reader, true, |record| {
            Ok((&record[1], &record[0]))
        })?;

        Ok(Self(mapping))
    }

    pub fn write_tsv<W: io::Write>(&self, writer: W) -> csv::Result<()> {
        self.0.write_tsv_with(
            writer,
            &Self::HEADER,
            |csv_writer, host_name, virus_name| csv_writer.write_record([virus_name, host_name]),
        )
    }
}

impl<HI, VI> HostVirusMapping<HI, VI>
where
    HI: Default + Interner<ValueHKT = HKTRef<str>>,
    VI: Default + Interner<ValueHKT = HKTRef<str>>,
{
    pub fn read_tsv<R: io::Read>(reader: R) -> anyhow::Result<Self> {
        Self::read_tsv_using(HI::default(), VI::default(), reader)
    }
}


impl<HI: Resolver, VI: Resolver> HostVirusMapping<HI, VI> {
    pub fn replace_interners<F, HI2, VI2>(self, f: F) -> HostVirusMapping<HI2, VI2>
    where
        HI2: Resolver<Symbol = HI::Symbol>,
        VI2: Resolver<Symbol = VI::Symbol>,
        F: FnOnce(HI, VI) -> (HI2, VI2)
    {
        HostVirusMapping(self.0.replace_interners(f))
    }
}

#[derive(Serialize, Deserialize)]
pub struct CrisprMatchRecord<'a> {
    #[serde(borrow)]
    pub virus_name: Cow<'a, str>,
    #[serde(borrow)]
    pub host_name: Cow<'a, str>,
}

mod minced_spacers;

/// Find matching CRISPR spacers between a metaviromic and a metagenomic sample.
#[derive(Args)]
pub struct CrisprMatchArgs {
    /// Directory to place intermediate files in, instead of a temporary directory. Existing files
    /// will not be re-created.
    #[clap(long, short = 'd')]
    work_dir: Option<String>,

    /// Path to the MinCED .jar file.
    #[clap(long = "minced", value_name = "MINCED_JAR", env = "METEOR_MINCED_JAR")]
    minced_jar: String,

    /// Installation prefix of the NCBI BLAST+ toolkit to use.
    #[clap(long, help_heading = Some("BLAST Options"), env = "METEOR_BLAST_PREFIX")]
    blast_prefix: Option<String>,

    /// Number of threads (blastp -num_threads option).
    #[clap(long, help_heading = Some("BLAST Options"), env = "METEOR_BLAST_NUM_THREADS")]
    num_threads: Option<i32>,

    /// Minimum hit identity % for the blastn search (0-100).
    #[clap(long, help_heading = Some("BLAST Options"), env = "METEOR_BLAST_PERC_IDENTITY", default_value_t = DEFAULT_PERC_IDENTITY)]
    perc_identity: i32,

    /// Do not actually run any external tools, print the command instead.
    #[clap(long)]
    dry_run: bool,

    /// FASTA file of the metaviromic sample.
    viral_seqs: String,

    /// FASTA file of the metagenomic sample.
    metagenomic_seqs: String,

    /// Output path. The result is a tab-separated values file containing two columns: virus_name,
    /// host_name. Use - to print to standard output.
    #[clap(long, short, default_value = "-")]
    output: String,
}

pub fn crispr_match(args: CrisprMatchArgs) -> anyhow::Result<()> {
    let Some(java) = cli_tools::Java::resolve()? else {
        anyhow::bail!("Could not locate a Java executable in PATH");
    };

    let Some(minced) = minced_spacers::MincedTool::resolve(java, args.minced_jar.as_ref()) else {
        anyhow::bail!("Could not locate the provided MinCED .jar file");
    };

    let blast_prefix = args.blast_prefix.as_ref().map(|s| s.as_ref());

    let Some(makeblastdb) = cli_tools::BlastTool::resolve(blast_prefix, "makeblastdb")? else {
        anyhow::bail!("Could not find a makeblastdb executable");
    };

    let Some(blastn) = cli_tools::BlastTool::resolve(blast_prefix, "blastn")? else {
        anyhow::bail!("Could not find a blastn executable");
    };

    let viral_seqs = args.viral_seqs.as_ref();
    let metagenomic_seqs = args.metagenomic_seqs.as_ref();

    let mut temp_dir = None;

    let work_dir = if let Some(work_dir) = &args.work_dir {
        Path::new(work_dir)
    } else {
        temp_dir.insert(TempDir::new("crispr-match")?).path()
    };

    let run_pipeline = || {
        let mut pipeline = MincedSpacersPipeline::new(minced, makeblastdb, blastn, work_dir);
        pipeline.set_blastn_num_threads(args.num_threads);
        pipeline.set_dry_run(args.dry_run);
        pipeline.match_spacers(viral_seqs, metagenomic_seqs, args.perc_identity)?;
        pipeline.collect_virus_host_mapping()
            .context("Collecting the virus-host mapping from the BLAST output")
    };

    let mapping = match run_pipeline() {
        Ok(mapping) => mapping,
        Err(e) => {
            if let Some(temp_dir) = temp_dir.take() {
                let path = temp_dir.into_path();
                eprintln!("Intentionally not deleting temporary work directory for inspection due to a previous error: {path:?}");
            }
            return Err(e);
        }
    };

    writing_new_file_or_stdout!(&args.output, output_file => {
        util::ignore_broken_pipe(mapping.write_tsv(output_file?))?;
    });

    Ok(())
}
