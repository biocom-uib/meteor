use std::{path::{PathBuf, Path, self}, process::{Command, Stdio}, io, fs::{File, self}};

use anyhow::Context;
use clap::Args;
use regex::Regex;
use serde::Deserialize;
use tempdir::TempDir;

use crate::util::{cli_tools::{self, CliTool}, interned_mapping::InternedMapping, writing_new_file_or_stdout, self};

use super::ProteinVirusTaxidMapping;

const FAA_FILE_NAME: &str = "predicted_viral_proteins.faa";

const PROTEIN_VIRUS_NAME_REGEX: &str = r#"_\d+$"#;
const DEFAULT_PRODIGAL_PROCEDURE: &str = "meta";

const DEFAULT_EVALUE_THRESHOLD: f32 = 1e-3;

const BLASTOUT_FMT: &str = "6 qseqid saccver staxid";
const BLASTOUT_FILE_NAME: &str = "uniprot_viral_search.blastout";


#[derive(Deserialize)]
struct ProteinProteinTaxidRecord<S = String> {
    predicted_name: S,
    uniprot_accession: S,
    taxid: u32,
}

pub struct ProdigalBlastpPipeline {
    prodigal: cli_tools::Prodigal,
    prodigal_procedure: String,
    blastp: cli_tools::BlastTool,
    work_dir: PathBuf,
    blastp_num_threads: Option<i32>,
    dry_run: bool,
}

impl ProdigalBlastpPipeline {
    pub fn new(prodigal: cli_tools::Prodigal, blastp: cli_tools::BlastTool, work_dir: PathBuf) -> Self {
        Self {
            prodigal,
            prodigal_procedure: DEFAULT_PRODIGAL_PROCEDURE.to_owned(),
            blastp,
            work_dir,
            blastp_num_threads: None,
            dry_run: false,
        }
    }

    pub fn set_prodigal_procedure(&mut self, procedure: &str) {
        self.prodigal_procedure = procedure.to_owned();
    }

    pub fn set_blastp_num_threads(&mut self, num_threads: Option<i32>) {
        self.blastp_num_threads = num_threads;
    }

    pub fn set_dry_run(&mut self, dry_run: bool) {
        self.dry_run = dry_run;
    }

    fn prodigal_cmd(&self, viral_seqs: &Path) -> io::Result<Command> {
        let mut cmd = self.prodigal.new_command();

        cmd.current_dir(&self.work_dir)
            .arg("-a")
            .arg(FAA_FILE_NAME)
            .arg("-i")
            .arg(path::absolute(viral_seqs)?)
            .args(["-p", &self.prodigal_procedure])
            .stdout(Stdio::null()); // silence prodigal flatfile output

        Ok(cmd)
    }

    fn blastp_cmd(&self, uniprot_blastdb: &Path, evalue: f32) -> io::Result<Command> {
        let mut cmd = self.blastp.new_command();

        if let Some(num_threads) = self.blastp_num_threads {
            cmd.args(["-num_threads", &num_threads.to_string()]);
        }

        cmd.current_dir(&self.work_dir)
            .arg("-db")
            .arg(path::absolute(uniprot_blastdb)?)
            .arg("-evalue")
            .arg(evalue.to_string())
            .arg("-query")
            .arg(FAA_FILE_NAME)
            .args(["-outfmt", BLASTOUT_FMT])
            .arg("-out")
            .arg(BLASTOUT_FILE_NAME);

        Ok(cmd)
    }

    pub fn find_and_match_proteins(
        &self,
        viral_seqs: &Path,
        uniprot_blastdb: &Path,
        evalue: f32,
    ) -> anyhow::Result<()> {

        if self.work_dir.join(FAA_FILE_NAME).exists() {
            eprintln!("{FAA_FILE_NAME} already exists in the working directory, skipping MinCED. \
                Delete the file to re-generate it");

        } else {
            let mut prodigal_cmd = self.prodigal_cmd(viral_seqs)?;

            eprintln!("Running {prodigal_cmd:?}");

            if !self.dry_run {
                let prodigal_status = prodigal_cmd.status().context("Running Prodigal")?;

                if !prodigal_status.success() {
                    anyhow::bail!("Prodigal finished with non-zero exit code: {}", prodigal_status.code().unwrap_or(0));
                }
            }
        }

        if self.work_dir.join(BLASTOUT_FILE_NAME).exists() {
            eprintln!("{BLASTOUT_FILE_NAME} already exists in the working directory, skipping BLAST \
                search. Delete the file to re-generate it");

        } else {
            let mut blastp_cmd = self.blastp_cmd(uniprot_blastdb, evalue)?;

            eprintln!("Running {blastp_cmd:?}");

            if !self.dry_run {
                let blastp_status = blastp_cmd.status().context("Running blastp")?;

                if !blastp_status.success() {
                    anyhow::bail!("blastp finished with non-zero exit code: {}", blastp_status.code().unwrap_or(0));
                }
            }
        }

        Ok(())
    }
}

impl ProdigalBlastpPipeline {
    pub fn collect_viral_protein_mapping(&self) -> anyhow::Result<ProteinVirusTaxidMapping> {
        let protein_virus_name_regex =
            Regex::new(PROTEIN_VIRUS_NAME_REGEX).expect("Error compiling PROTEIN_VIRUS_NAME_REGEX");

        let file_reader = File::open(self.work_dir.join(BLASTOUT_FILE_NAME))?;

        let mapping = InternedMapping::read_tsv_with(file_reader, false, |record| {
            let record = record.deserialize::<ProteinProteinTaxidRecord<&str>>(None)?;

            if let Some(suffix_match) = protein_virus_name_regex.find(record.predicted_name) {
                let virus_name = &record.predicted_name[..suffix_match.start()];
                Ok((record.uniprot_accession, (virus_name, record.taxid)))
            } else {
                anyhow::bail!(
                    "Could not parse virus name from blast output: {:?}",
                    record.predicted_name
                );
            }
        })
        .context("Reading BLAST+ blastp output")?;

        Ok(ProteinVirusTaxidMapping(mapping))
    }
}


/// Extract viral proteins and search them in a UniProt BLAST+ database.
#[derive(Args)]
pub struct MatchViralProteinsArgs {
    /// Directory to use to prepare the inputs for makeblastdb. A temporary directory is used by
    /// default.
    #[clap(long)]
    work_dir: Option<String>,

    /// Do not actually run any external tools, print the command instead.
    #[clap(long)]
    dry_run: bool,

    /// Path to the prodigal executable if not in PATH.
    #[clap(long, help_heading = Some("Prodigal Options"), env = "METEOR_PRODIGAL_BIN")]
    prodigal_bin: Option<String>,

    /// Prodigal procedure (Prodigal 2.x -p flag) to use.
    #[clap(long, help_heading = Some("Prodigal Options"), default_value = "meta")]
    prodigal_procedure: String,

    /// Installation prefix of the NCBI BLAST+ toolkit to use.
    #[clap(long, help_heading = Some("BLAST Options"), env = "METEOR_BLAST_PREFIX")]
    blast_prefix: Option<String>,

    /// Location of the UniProt BLAST+ database. See `meteor ppi prepare-uniprot-blastdb`.
    #[clap(long, help_heading = Some("BLAST Options"))]
    uniprot_blastdb: String,

    /// E-value threshold for the blastp search.
    #[clap(long, help_heading = Some("BLAST Options"), env = "METEOR_BLAST_EVALUE_THRESHOLD", default_value_t = DEFAULT_EVALUE_THRESHOLD)]
    evalue: f32,

    /// Number of threads (blastp -num_threads option).
    #[clap(long, help_heading = Some("BLAST Options"), env = "METEOR_BLAST_NUM_THREADS")]
    num_threads: Option<i32>,

    /// Input FASTA of viral contigs.
    viral_seqs: String,

    /// Output path. The result is a tab-separated values file containing three columns:
    /// virus_name, uniprot_accession and taxid. Use '-' to print to standard output.
    #[clap(short, long, default_value = "-")]
    output: String,
}

pub fn match_viral_proteins(args: MatchViralProteinsArgs) -> anyhow::Result<()> {
    let prodigal_bin = args.prodigal_bin.as_ref().map(|s| s.as_ref());

    let Some(prodigal) = cli_tools::Prodigal::resolve(prodigal_bin)? else {
        anyhow::bail!("Could not locate a Prodigal executable");
    };

    let blast_prefix = args.blast_prefix.as_ref().map(|s| s.as_ref());

    let Some(blastp) = cli_tools::BlastTool::resolve(blast_prefix, "blastp")? else {
        anyhow::bail!("Could not find a blastp executable");
    };

    let mut temp_dir = None;

    let work_dir = if let Some(work_dir) = &args.work_dir {
        Path::new(work_dir)
    } else {
        temp_dir.insert(TempDir::new("match-viral-proteins")?).path()
    };

    fs::create_dir_all(work_dir)?;

    let run_pipeline = || {
        let mut pipeline = ProdigalBlastpPipeline::new(prodigal, blastp, work_dir.to_path_buf());
        pipeline.set_dry_run(args.dry_run);
        pipeline.set_prodigal_procedure(&args.prodigal_procedure);
        pipeline.set_blastp_num_threads(args.num_threads);
        pipeline.find_and_match_proteins(args.viral_seqs.as_ref(), args.uniprot_blastdb.as_ref(), args.evalue)?;
        pipeline.collect_viral_protein_mapping()
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
