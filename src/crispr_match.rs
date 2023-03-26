use std::{io, path::Path};

use clap::Args;
use tempdir::TempDir;

use crate::{
    crispr_match::minced_spacers::MincedSpacersPipeline,
    util::{cli_tools, interned_mapping::InternedMultiMapping, writing_new_file_or_stdout},
};


const DEFAULT_PERC_IDENTITY: i32 = 95;

#[derive(Default)]
pub struct VirusHostMapping(pub InternedMultiMapping);

impl VirusHostMapping {
    pub fn read_tsv<R: io::Read>(reader: R) -> anyhow::Result<Self> {
        Ok(Self(InternedMultiMapping::read_tsv(reader, true)?))
    }

    pub fn write_tsv<W: io::Write>(&self, writer: W) -> csv::Result<()> {
        self.0.write_tsv(writer, &["virus_name", "host_name"])
    }
}

#[derive(Default)]
pub struct HostVirusMapping(pub InternedMultiMapping);

impl HostVirusMapping {
    pub fn read_tsv<R: io::Read>(reader: R) -> anyhow::Result<Self> {
        let mapping = InternedMultiMapping::<_>::read_tsv_with(reader, true, |record| {
            Ok((&record[1], &record[0]))
        })?;

        Ok(Self(mapping))
    }

    pub fn write_tsv<W: io::Write>(&self, writer: W) -> csv::Result<()> {
        self.0.write_tsv_with(
            writer,
            &["virus_name", "host_name"],
            |csv_writer, host_name, virus_name| csv_writer.write_record([virus_name, host_name]),
        )
    }
}

mod minced_spacers;

#[derive(Args)]
pub struct CrisprMatchArgs {
    /// Directory to place intermediate files in, instead of a temporary directory. Existing files
    /// will not be re-created.
    #[clap(long, short = 'd')]
    work_dir: Option<String>,

    /// Path to the MinCED .jar file.
    #[clap(long = "minced", value_name = "MINCED_JAR")]
    minced_jar: String,

    /// Installation prefix of the NCBI BLAST+ toolkit to use.
    #[clap(long, help_heading = Some("BLAST Options"))]
    blast_prefix: Option<String>,

    /// Number of threads (blastp -num_threads option).
    #[clap(long, help_heading = Some("BLAST Options"))]
    num_threads: Option<i32>,

    /// Minimum hit identity % for the blastn search (0-100).
    #[clap(long, help_heading = Some("BLAST Options"), default_value_t = DEFAULT_PERC_IDENTITY)]
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
        let mut pipeline = MincedSpacersPipeline::new(minced, blastn, work_dir);
        pipeline.set_blastn_num_threads(args.num_threads);
        pipeline.set_dry_run(args.dry_run);
        pipeline.match_spacers(viral_seqs, metagenomic_seqs, args.perc_identity)?;
        pipeline.collect_virus_host_mapping()
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
        mapping.write_tsv(output_file?)?;
    });

    Ok(())
}
