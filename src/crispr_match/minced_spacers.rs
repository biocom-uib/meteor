use std::{path::{Path, PathBuf, self}, process::Command, io, fs::File};

use anyhow::Context;
use regex::Regex;
use serde::{Deserialize, Serialize};

use crate::util::{cli_tools::{self, CliTool, BlastTool}, interned_mapping::InternedMultiMapping};

use super::VirusHostMapping;

const CRISPRS_FILE_NAME: &str = "crisprs.txt";
const SPACERS_FILE_NAME: &str = "crisprs_spacers.fa"; // hardcoded in minced since v0.1.5

const BLASTOUT_FILE_NAME: &str = "spacer_search.blastout";

const BLASTOUT_FMT: &str = "6 qaccver saccver";

const DEFAULT_PERC_IDENTITY: i32 = 95;
const SPACER_SUFFIX_REGEX: &str = r#"_CRISPR_\d+_spacer_\d+$"#;

#[derive(Deserialize)]
struct SpacersBlastRecord<'a> {
    spacer_name: &'a str,
    host_name: &'a str,
}

#[derive(Serialize, Deserialize)]
pub struct CrisprHostPrediction<'a> {
    pub virus_name: &'a str,
    pub host_name: &'a str,
}


pub struct MincedTool {
    java: cli_tools::Java,
    minced_jar: PathBuf,
}

impl MincedTool {
    pub fn resolve(java: cli_tools::Java, minced_jar: &Path) -> Option<Self> {
        if minced_jar.exists() {
            Some(MincedTool { java, minced_jar: minced_jar.to_path_buf() })
        } else {
            None
        }
    }
}

impl CliTool for MincedTool {
    fn new_command(&self) -> Command {
        let mut cmd = self.java.new_command();
        cmd.args(["-jar", &self.minced_jar.to_string_lossy()]);
        cmd
    }
}

pub struct MincedSpacersPipeline {
    minced: MincedTool,
    blastn: BlastTool,
    work_dir: PathBuf,
    blastn_num_threads: Option<i32>,
    dry_run: bool,
}

impl MincedSpacersPipeline {
    pub fn new(minced: MincedTool, blastn: BlastTool, work_dir: &Path) -> Self {
        MincedSpacersPipeline {
            minced,
            blastn,
            work_dir: work_dir.to_path_buf(),
            blastn_num_threads: None,
            dry_run: false,
        }
    }

    pub fn set_blastn_num_threads(&mut self, num_threads: Option<i32>) {
        self.blastn_num_threads = num_threads;
    }

    pub fn set_dry_run(&mut self, dry_run: bool) {
        self.dry_run = dry_run;
    }

    fn minced_cmd(&self, metagenomic_seqs: &Path) -> io::Result<Command> {
        let mut cmd = self.minced.new_command();

        cmd.current_dir(&self.work_dir)
            .arg("-spacers")
            .arg(path::absolute(metagenomic_seqs)?)
            .arg(CRISPRS_FILE_NAME);

        Ok(cmd)
    }

    fn blastn_cmd(&self, perc_identity: i32, subject: &Path) -> io::Result<Command> {
        let mut cmd = self.blastn.new_command();

        let query_path = self.work_dir.join(SPACERS_FILE_NAME);

        if let Some(num_threads) = self.blastn_num_threads {
            cmd.args(["-num_threads", &num_threads.to_string()]);
        }

        cmd.current_dir(&self.work_dir)
            .args(["-task", "blastn-short"])
            .arg("-perc_identity")
            .arg(perc_identity.to_string())
            .arg("-query")
            .arg(query_path)
            .arg("-subject")
            .arg(path::absolute(subject)?)
            .args(["-outfmt", BLASTOUT_FMT])
            .arg("-out")
            .arg(BLASTOUT_FILE_NAME);

        Ok(cmd)
    }

    pub fn match_spacers(
        &self,
        viral_seqs: &Path,
        metagenomic_seqs: &Path,
        perc_identity: Option<i32>,
    ) -> anyhow::Result<()> {
        if self.work_dir.join(SPACERS_FILE_NAME).exists() {
            eprintln!("{SPACERS_FILE_NAME} already exists in the working directory, skipping
                MinCED. Delete the file to re-generate it");

        } else {
            let mut minced_cmd = self.minced_cmd(metagenomic_seqs)?;

            eprintln!("Running {minced_cmd:?}");

            if !self.dry_run {
                let minced_status = minced_cmd.status().context("Running MinCED")?;

                if !minced_status.success() {
                    anyhow::bail!("MinCED finished with non-zero exit code: {}", minced_status.code().unwrap_or(0));
                }
            }
        }

        if self.work_dir.join(BLASTOUT_FILE_NAME).exists() {
            eprintln!("{BLASTOUT_FILE_NAME} already exists in the working directory, skipping BLAST
                search. Delete the file to re-generate it");

        } else {
            let perc_identity = perc_identity.unwrap_or(DEFAULT_PERC_IDENTITY);
            let mut blastn_cmd = self.blastn_cmd(perc_identity, viral_seqs)?;

            eprintln!("Running {blastn_cmd:?}");

            if !self.dry_run {
                let blastn_status = blastn_cmd.status().context("Running blastn")?;

                if !blastn_status.success() {
                    anyhow::bail!("blastn finished with non-zero exit code: {}", blastn_status.code().unwrap_or(0));
                }
            }
        }

        Ok(())
    }

    pub fn collect_virus_host_mapping(&self) -> anyhow::Result<VirusHostMapping> {
        let spacer_suffix_regex =
            Regex::new(SPACER_SUFFIX_REGEX).expect("Error compiling SPACER_SUFFIX_REGEX");

        let file_reader = File::open(self.work_dir.join(BLASTOUT_FILE_NAME))?;

        let mapping = InternedMultiMapping::read_tsv_with(file_reader, false, |record| {
            let record = record.deserialize::<SpacersBlastRecord>(None)?;

            if let Some(suffix_match) = spacer_suffix_regex.find(record.spacer_name) {
                let virus_name = &record.spacer_name[..suffix_match.start()];

                Ok((virus_name, record.host_name))
            } else {
                anyhow::bail!(
                    "Could not parse spacer name from blast output: {:?}",
                    record.spacer_name
                );
            }
        })
        .context("Reading BLAST+ blastn output")?;

        Ok(VirusHostMapping(mapping))
    }
}