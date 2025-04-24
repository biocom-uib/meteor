use std::{
    fs::{self, File},
    io,
    path::{Path, PathBuf},
    time::Duration,
};

use clap::{Args, ValueEnum};
use flate2::read::GzDecoder;
use serde::{Serialize, Deserialize};

use crate::util::{progress_monitor, io::MaybeGzDecoder};

// viruses.string-db.org does not support https (23/03/2023)
const STRING_VIRUSES_BASE_URL: &str = "http://viruses.string-db.org/download";
const STRING_VIRUSES_VERSION: &str = "v10.5";


#[derive(Serialize, Deserialize)]
pub struct ProteinLinksRecord<S = String> {
    pub protein_id_a: S,
    pub protein_id_b: S,
    pub combined_score: i32,
}

#[derive(ValueEnum, Debug, Copy, Clone, PartialEq, Eq)]
pub enum FileName {
    #[value(name = "protein.links")]
    ProteinLinks,
    #[value(name = "protein.links.detailed")]
    ProteinLinksDetailed,
    #[value(name = "protein.links.full")]
    ProteinLinksFull,
    #[value(name = "protein.sequences")]
    ProteinSequences,
    #[value(name = "species")]
    Species,
    #[value(name = "protein.aliases")]
    ProteinAliases,
}

impl FileName {
    pub fn into_full_file_name(self) -> String {
        use FileName::*;

        match self {
            ProteinLinks => format!("protein.links.{STRING_VIRUSES_VERSION}.txt.gz"),
            ProteinLinksDetailed => {
                format!("protein.links.detailed.{STRING_VIRUSES_VERSION}.txt.gz")
            }
            ProteinLinksFull => format!("protein.links.full.{STRING_VIRUSES_VERSION}.txt.gz"),
            ProteinSequences => format!("protein.sequences.{STRING_VIRUSES_VERSION}.fa.gz"),
            Species => format!("species.{STRING_VIRUSES_VERSION}.txt"),
            ProteinAliases => format!("protein.aliases.{STRING_VIRUSES_VERSION}.txt.gz"),
        }
    }

    pub fn into_local_path(self, download_dir: &Path) -> PathBuf {
        download_dir.join(self.into_full_file_name())
    }

    pub fn is_gz_compressed(self) -> bool {
        self != Self::Species
    }

    pub fn open_read(self, download_dir: &Path) -> io::Result<MaybeGzDecoder<File>> {
        let file = File::open(self.into_local_path(download_dir))?;

        if self.is_gz_compressed() {
            Ok(MaybeGzDecoder::GzDecoder(GzDecoder::new(file)))
        } else {
            Ok(MaybeGzDecoder::Reader(file))
        }
    }
}

pub fn parse_string_id(string_id: &str) -> Option<(u32, &str)> {
    let (taxid, name) = string_id.split_once('.')?;

    Some((taxid.parse().ok()?, name))
}

/// Fetch STRING Viruses exported data
#[derive(Args)]
pub struct StringVirusesFetchArgs {
    /// File names to download (not including version or extension).
    #[clap(num_args = 1.., value_enum)]
    pub files: Vec<FileName>,

    /// Always re-download existing files
    #[clap(long, short)]
    pub force_download: bool,

    /// Directory to download Viruses StringDB files into
    #[clap(long, short = 'd')]
    pub download_dir: String,
}

pub fn fetch(args: &StringVirusesFetchArgs) -> anyhow::Result<()> {
    let work_dir = Path::new(&args.download_dir);

    fs::create_dir_all(work_dir)?;

    let client = reqwest::blocking::Client::new();

    for &file in &args.files {
        let file_name = file.into_full_file_name();
        let local_partial_path = work_dir.join(format!("{file_name}.part"));
        let local_path = work_dir.join(&file_name);

        if args.force_download || !local_path.exists() {
            let url = format!("{STRING_VIRUSES_BASE_URL}/{file_name}");

            eprintln!("Downloading {url:?} to {local_path:?}");

            let response = client.get(&url).send()?.error_for_status()?;

            let local_file = File::create(&local_partial_path)?;

            let r = if let Some(file_size) = response.content_length() {
                progress_monitor::monitor_copy(
                    response,
                    local_file,
                    Duration::from_secs(1),
                    progress_monitor::default_progress_callback(file_size),
                    progress_monitor::default_finish_callback(),
                )
            } else {
                progress_monitor::monitor_copy(
                    response,
                    local_file,
                    Duration::from_secs(1),
                    progress_monitor::unknown_progress_callback(),
                    progress_monitor::default_finish_callback(),
                )
            };

            eprintln!();

            if r.is_err() {
                fs::remove_file(&local_partial_path)?;
            }

            r?;

            fs::rename(&local_partial_path, &local_path)?;
        } else {
            println!("File {file_name:?} is already present as {local_path:?}, skipping (delete it or use --force-download to re-download)");
        }
    }

    Ok(())
}
