use std::{fs::File, io::{BufReader, BufWriter}, path::Path};

use anyhow::Context;
use clap::{Args, ValueEnum};
use serde::{de::Error, Deserialize, Serialize, Deserializer};

use crate::taxonomy::formats::{newick::NewickTaxonomy, ncbi::{self, NcbiTaxonomy}};

#[derive(Deserialize, Serialize)]
pub enum SomeTaxonomy {
    NcbiTaxonomyWithSingleClassNames(NcbiTaxonomy<ncbi::SingleClassNames>),
    NcbiTaxonomyWithManyNames(NcbiTaxonomy<ncbi::AllNames>),
    NewickTaxonomy(NewickTaxonomy),
}

#[derive(Deserialize, Serialize)]
pub struct PreprocessedTaxonomy {
    #[serde(deserialize_with = "deserialize_check_version")]
    version: u64,
    pub tree: SomeTaxonomy,
    pub ordered_ranks: Option<Vec<String>>,
    pub contraction_ranks: Option<Vec<String>>,
}

macro_rules! with_some_taxonomy {
    ($some_taxonomy:expr, $tax:pat => $body:expr $(,)?) => {
        match $some_taxonomy {
            crate::preprocessed_taxonomy::SomeTaxonomy::NcbiTaxonomyWithSingleClassNames($tax) => $body,
            crate::preprocessed_taxonomy::SomeTaxonomy::NcbiTaxonomyWithManyNames($tax) => $body,
            crate::preprocessed_taxonomy::SomeTaxonomy::NewickTaxonomy($tax) => $body,
        }
    };
}

macro_rules! with_some_ncbi_or_newick_taxonomy {
    (
        $some_taxonomy:expr,
        ncbi: $ncbi_tax:pat => $ncbi_body:expr,
        newick: $newick_tax:pat => $newick_body:expr
        $(,)?
    ) => {
        match $some_taxonomy {
            crate::preprocessed_taxonomy::SomeTaxonomy::NcbiTaxonomyWithSingleClassNames($ncbi_tax) => $ncbi_body,
            crate::preprocessed_taxonomy::SomeTaxonomy::NcbiTaxonomyWithManyNames($ncbi_tax) => $ncbi_body,
            crate::preprocessed_taxonomy::SomeTaxonomy::NewickTaxonomy($newick_tax) => $newick_body,
        }
    };
}

pub(crate) use with_some_ncbi_or_newick_taxonomy;
pub(crate) use with_some_taxonomy;

#[derive(ValueEnum, Debug, Default, Copy, Clone)]
pub enum PreprocessedTaxonomyFormat {
    #[default]
    Bincode,
    #[cfg(feature = "taxonomy-serialize-cbor")]
    Cbor,
    #[cfg(feature = "taxonomy-serialize-json")]
    Json,
}

fn deserialize_check_version<'de, D: Deserializer<'de>>(d: D) -> Result<u64, D::Error> {
    let version = u64::deserialize(d)?;

    if version != PreprocessedTaxonomy::FORMAT_VERSION {
        let msg = format!("Preprocessed taxonomy version mismatch (found: {}, current: {}). Re-generate it using preprocess-taxonomy.", version, PreprocessedTaxonomy::FORMAT_VERSION);
        return Err(D::Error::custom(msg));
    }

    Ok(version)
}

impl PreprocessedTaxonomy {
    pub const FORMAT_VERSION: u64 = 2
        + NcbiTaxonomy::<ncbi::SingleClassNames>::FORMAT_VERSION as u64
        + NcbiTaxonomy::<ncbi::AllNames>::FORMAT_VERSION as u64
        + NewickTaxonomy::FORMAT_VERSION as u64;

    pub fn new(
        tree: SomeTaxonomy,
        ordered_ranks: Option<Vec<String>>,
        contraction_ranks: Option<Vec<String>>,
    ) -> Self {
        Self {
            version: Self::FORMAT_VERSION,
            tree,
            ordered_ranks,
            contraction_ranks,
        }
    }

    pub fn deserialize_with_format<P: AsRef<Path>>(
        path: P,
        format: PreprocessedTaxonomyFormat,
    ) -> anyhow::Result<Self> {
        let reader = BufReader::new(File::open(path)?);

        match format {
            PreprocessedTaxonomyFormat::Bincode => Ok(bincode::deserialize_from(reader)?),

            #[cfg(feature = "taxonomy-serialize-cbor")]
            PreprocessedTaxonomyFormat::Cbor => Ok(serde_cbor::from_reader(reader)?),

            #[cfg(feature = "taxonomy-serialize-json")]
            PreprocessedTaxonomyFormat::Json => Ok(serde_json::from_reader(reader)?),
        }
    }

    pub fn serialize_with_format<P: AsRef<Path>>(
        &self,
        path: P,
        format: PreprocessedTaxonomyFormat,
    ) -> anyhow::Result<()> {
        let writer = BufWriter::new(File::create(path)?);

        match format {
            PreprocessedTaxonomyFormat::Bincode => {
                Ok(bincode::serialize_into(writer, self)?)
            }
            #[cfg(feature ="taxonomy-serialize-cbor")]
            PreprocessedTaxonomyFormat::Cbor => {
                Ok(serde_cbor::to_writer(writer, self)?)
            }
            #[cfg(feature ="taxonomy-serialize-json")]
            PreprocessedTaxonomyFormat::Json => {
                Ok(serde_json::to_writer(writer, self)?)
            }
        }
    }
}

#[derive(Args)]
pub struct PreprocessedTaxonomyArgs {
    /// Path to the preprocessed taxonomy (from preprocess-taxonomy)
    #[clap(help_heading = "Taxonomy", long = "taxonomy", env = "METEOR_TAXONOMY")]
    preprocessed_taxonomy: String,

    /// Format in which the preprocessed taxonomy was serialized.
    #[clap(help_heading = "Taxonomy", long, value_enum, default_value_t, env = "METEOR_TAXONOMY_FORMAT")]
    taxonomy_format: PreprocessedTaxonomyFormat,
}

impl PreprocessedTaxonomyArgs {
    pub fn deserialize(&self) -> anyhow::Result<PreprocessedTaxonomy> {
        PreprocessedTaxonomy::deserialize_with_format(
            &self.preprocessed_taxonomy,
            self.taxonomy_format
        )
        .context("Could not deserialize the preprocessed taxonomy")
    }
}
