use std::{
    borrow::Cow, fs, path::{Path, PathBuf}
};

use anyhow::Context;
use polars::prelude::{LazyCsvReader, LazyFrame, NullValues};
use serde::{Deserialize, Serialize};

use crate::{
    csv::{filter::FilterableRecord, polars::{lazyframe_from_file, PolarsSchema}},
    taxonomy::LabeledTaxonomy,
};

#[derive(Debug, Clone, Serialize, Deserialize, FilterableRecord, PolarsSchema)]
#[filter(alias_field, alias_filter, polars_literal)]
#[polars(set_csv_options = csv_options)]
pub struct GenomadVirusRecord<'a> {
    #[filter(as_string)]
    #[serde(borrow)]
    pub seq_name: Cow<'a, str>,
    pub length: u32,
    #[filter(skip)]
    #[serde(borrow)]
    pub topology: Cow<'a, str>,
    #[filter(skip)]
    #[serde(borrow)]
    pub coordinates: Cow<'a, str>,
    pub n_genes: u32,
    #[filter(skip)]
    #[serde(borrow)]
    pub genetic_code: Cow<'a, str>,
    pub virus_score: f32,
    #[filter(skip)]
    #[serde(borrow)]
    pub fdr: Cow<'a, str>,
    pub n_hallmarks: u32,
    pub marker_enrichment: f32,
    #[filter(skip)]
    #[serde(borrow)]
    pub taxonomy: Cow<'a, str>,
}

fn csv_options(reader: LazyCsvReader) -> LazyCsvReader {
    // yes, the typo is intentional
    let null_values = NullValues::Named(vec![
        ("taxonomy".into(), "Unclassified".into()),
    ]);

    crate::csv::polars::default_tsv_options(reader)
        .with_null_values(Some(null_values))
}

pub fn locate_summary_dir(output_dir: &Path) -> anyhow::Result<PathBuf> {
    if !output_dir.is_dir() {
        anyhow::bail!("{output_dir:?} is not a directory");
    }

    let mut summary_dirs = Vec::new();

    for entry in fs::read_dir(output_dir)? {
        let entry = entry?;
        let entry_path = entry.path();

        let Ok(dir_name) = entry.file_name().into_string() else {
            anyhow::bail!("could not convert {entry_path:?} base name to std::String");
        };

        if entry_path.is_dir() && dir_name.ends_with("_summary") {
            summary_dirs.push(entry_path);
        }
    }

    match &mut *summary_dirs {
        [] => anyhow::bail!("no *_summary directories found in {output_dir:?}"),
        [summary_dir] => Ok(std::mem::take(summary_dir)),
        _ => anyhow::bail!(
            "found multiple *_summary directories in {output_dir:?}: {summary_dirs:?}"
        ),
    }
}

pub fn locate_virus_summary_file(summary_dir: &Path) -> anyhow::Result<PathBuf> {
    let mut summary_files = Vec::new();

    for entry in fs::read_dir(summary_dir)? {
        let entry = entry?;
        let entry_path = entry.path();

        let Ok(dir_name) = entry.file_name().into_string() else {
            anyhow::bail!("could not convert {summary_dir:?} base name to std::String");
        };

        if entry_path.is_file() && dir_name.ends_with("_virus_summary.tsv") {
            summary_files.push(entry_path);
        }
    }

    match &mut *summary_files {
        [] => anyhow::bail!("no *_virus_summary.tsv files found in {summary_dir:?}"),
        [summary_file] => Ok(std::mem::take(summary_file)),
        _ => anyhow::bail!(
            "found multiple *_summary directories in {summary_dir:?}: {summary_files:?}"
        ),
    }
}

pub fn virus_summary_lazyframe_from_dir(output_dir: &Path) -> anyhow::Result<LazyFrame> {
    let virus_summary = locate_summary_dir(output_dir)
        .and_then(|summary_dir| locate_virus_summary_file(&summary_dir))
        .with_context(|| "Could not locate {output_dir:?}/*_summary/*_virus_summary.tsv file")?;

    Ok(lazyframe_from_file::<GenomadVirusRecord>(&virus_summary)?)
}

pub fn lineage_most_specific(lineage: &str) -> Option<&str> {
    lineage.rsplit(';').find(|name| !name.is_empty())
}

pub fn lineage_to_taxid<Tax>(tax: &Tax, lineage: &str) -> Option<u32>
where
    Tax: LabeledTaxonomy,
{
    for name in lineage.rsplit(';') {
        if name.is_empty() {
            continue;
        }

        let mut node_ids = tax.nodes_with_label(name);

        if let Some(node_id) = node_ids.next() {
            if node_ids.next().is_some() {
                eprintln!("Found multiple nodes with name {name} in the taxonomy");
            }

            return Some(node_id.into());
        }
    }

    None
}


#[cfg(test)]
pub(crate) mod tests {
    use std::sync::Arc;

    use polars::series::Series;

    use crate::{
        csv::polars::{lazyframe_from_string_strip, ExprExt},
        taxonomy::formats::newick,
    };

    use super::GenomadVirusRecord;

    pub const VIRUS_SUMMARY: &str = "\
      seq_name  \t length \t topology            \t coordinates \t n_genes \t genetic_code \t virus_score \t fdr \t n_hallmarks \t marker_enrichment \t taxonomy
        v1      \t 17430  \t No terminal repeats \t NA          \t 20      \t 11           \t 0.9835      \t NA  \t 10          \t 33.3694           \t S1;p2;c3;o4;f5;;
        v2      \t 17700  \t No terminal repeats \t NA          \t 25      \t 11           \t 0.9835      \t NA  \t 1           \t 36.6676           \t S1;p2;c3;o4;f5;;
        v3      \t 11041  \t No terminal repeats \t NA          \t 13      \t 15           \t 0.9835      \t NA  \t 4           \t 20.6106           \t S1;p2;c3;o4;f5;g6;s9
        v4      \t 20278  \t No terminal repeats \t NA          \t 23      \t 11           \t 0.9835      \t NA  \t 7           \t 38.6188           \t S1;p2;c3;o4;f5;;
        v5      \t 18846  \t No terminal repeats \t NA          \t 29      \t 11           \t 0.9835      \t NA  \t 7           \t 47.6528           \t S1;p2;c3;o4;f5;;
        v6      \t 12467  \t No terminal repeats \t NA          \t 11      \t 11           \t 0.9835      \t NA  \t 0           \t 17.1491           \t S1;p2;c3;o4;f5;;
        v7      \t 17476  \t No terminal repeats \t NA          \t 21      \t 11           \t 0.9835      \t NA  \t 7           \t 34.3338           \t S1;p2;c3;o4;f5;;s19
        v8      \t 17462  \t No terminal repeats \t NA          \t 21      \t 11           \t 0.9834      \t NA  \t 10          \t 33.9309           \t S1;p2;c3;o4;f5;;
        v9      \t 18768  \t No terminal repeats \t NA          \t 20      \t 11           \t 0.9834      \t NA  \t 6           \t 33.7222           \t ;;;;;;
    ";

    #[test]
    fn test_taxonomy_to_taxid() -> anyhow::Result<()> {
        use polars::lazy::dsl::col;

        let tax = Arc::new(newick::tests::sample_taxonomy());

        let df = lazyframe_from_string_strip::<GenomadVirusRecord>(VIRUS_SUMMARY)?;

        let taxid = df
            .with_column(
                col("taxonomy").downcast_map_apply(Series::str, move |lineage| {
                    super::lineage_to_taxid(&*tax, lineage)
                }),
            )
            .collect()?
            .column("taxonomy")?
            .u32()?
            .to_vec();

        let expected: Vec<Option<u32>> = vec![
            Some(5),
            Some(5),
            Some(9),
            Some(5),
            Some(5),
            Some(5),
            Some(19),
            Some(5),
            None,
        ];

        assert_eq!(&taxid, &expected);

        Ok(())
    }
}
