use std::{borrow::Cow, path::{Path, PathBuf}};

use polars::prelude::{LazyCsvReader, NullValues};
use serde::{Deserialize, Serialize};

use crate::{csv::{filter::FilterableRecord, polars::PolarsSchema}, taxonomy::LabeledTaxonomy};


#[derive(Debug, Clone, Serialize, Deserialize, FilterableRecord, PolarsSchema)]
#[filter(alias_field, alias_filter, polars_literal)]
#[polars(set_csv_options = csv_options)]
pub struct CherryPredictionRecord<'a> {
    #[filter(as_string, rename = "Accession")]
    #[polars(rename = "Accession")]
    #[serde(borrow)]
    accession: Cow<'a, str>,
    #[filter(rename = "Length")]
    #[polars(rename = "Length")]
    length: u32,
    #[filter(as_string, rename = "Pred")]
    #[polars(rename = "Pred")]
    #[serde(borrow)]
    pred: Cow<'a, str>,
    #[filter(rename = "Score")]
    #[polars(rename = "Score")]
    score: f32,
    #[filter(as_string, rename = "Type")]
    #[polars(rename = "Type")]
    #[serde(borrow)]
    r#type: Cow<'a, str>,
}

fn csv_options(reader: LazyCsvReader) -> LazyCsvReader {
    // yes, the typo is intentional
    let null_values = NullValues::Named(vec![
        ("Pred".into(), "unknown".into()),
        ("Type".into(), "-".into()),
    ]);

    crate::csv::polars::default_tsv_options(reader)
        .with_null_values(Some(null_values))
        .with_separator(b',')
}

pub fn locate_cherry_prediction(output_dir: &Path) -> anyhow::Result<PathBuf> {
    if !output_dir.is_dir() {
        anyhow::bail!("{output_dir:?} is not a directory");
    }

    let mut result = output_dir.to_path_buf();

    result.push("cherry_prediction.csv");

    if !result.is_file() {
        anyhow::bail!("{result:?} does not exist");
    }

    Ok(result)
}

pub fn pred_to_taxid<Tax>(tax: &Tax, pred: &str) -> Option<u32>
where
    Tax: LabeledTaxonomy,
{
    let mut node_ids = tax.nodes_with_label(pred);

    let node_id = node_ids.next()?;

    if node_ids.next().is_some() {
        eprintln!("Found multiple nodes with name {pred} in the taxonomy");
    }

    Some(node_id.into())
}
