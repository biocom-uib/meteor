use std::borrow::Cow;

use polars::prelude::{LazyCsvReader, LazyFrame, NullValues};
use serde::{Deserialize, Serialize};

use crate::csv::{filter::FilterableRecord, polars::PolarsSchema};

#[derive(Debug, Clone, Serialize, Deserialize, FilterableRecord, PolarsSchema)]
#[filter(alias_field, alias_filter, polars_literal)]
#[polars(set_csv_options = csv_options)]
pub struct PhaGenusRecord<'a> {
    #[filter(as_string)]
    contigs: Cow<'a, str>,
    #[filter(as_string)]
    predict_label: Cow<'a, str>,
    uncertainty: f32,
}

fn csv_options(reader: LazyCsvReader) -> LazyCsvReader {
    // yes, the typo is intentional
    let null_values = NullValues::Named(vec![
        ("predict_label".into(), "no predication".into()),
        ("uncertainty".into(),   "NA".into()),
    ]);

    crate::csv::polars::default_tsv_options(reader)
        .with_null_values(Some(null_values))
}

//pub fn locate_virus_output(output_dir: &Path) -> anyhow::Result<PathBuf> {
//    if !output_dir.is_dir() {
//        anyhow::bail!("{output_dir:?} is not a directory");
//    }

//    let mut result = fs::canonicalize(output_dir)?;

//    // The output is often .../<name>/<name>.csv
//    let base_name = result
//        .file_name()
//        .ok_or_else(|| anyhow::anyhow!("{result:?} has no base name"))?
//        .to_owned();

//    result.push(base_name);
//    result.set_extension("csv");

//    if !result.is_file() {
//        anyhow::bail!("{result:?} does not exist");
//    }

//    Ok(result)
//}

/// Inptut columns: `..., {predict_label_col}: str (name1/name2/...), ...`
///
/// Output columns: `..., {predict_label_col}: str (name_i), ...`
pub fn explode_predict_label(df: LazyFrame, predict_label_col: &str) -> LazyFrame {
    use polars::lazy::dsl::{col, lit};

    df.with_column(col(predict_label_col).str().split(lit("/")))
        .explode([col(predict_label_col)])
}

/**
 * Inptut columns:
 *   - ..., {predict_label_col}: str, ...
 *
 * Output columns:
 *   - ..., {predict_label_col}: [str], {taxids_new_col}: [u32], ...
 *
 * where
 *     {predict_label_col}.list().len() == {taxids_new_col}.list().len()
 *
*/
//pub fn repeat_predict_label_with_taxids<TaxPtr>(
//    df: LazyFrame,
//    tax: TaxPtr,
//    rank_sym: <TaxPtr::Target as Taxonomy>::RankSym,
//    predict_label_col: &str,
//    taxids_new_col: &str,
//) -> LazyFrame
//where
//    TaxPtr: Deref<Target: LabeledTaxonomy> + Send + Sync + 'static,
//{
//    use polars::lazy::dsl::col;

//    df.with_column(
//        col(predict_label_col)
//            .downcast_map_apply_to_chunked(Series::str, move |class_name| {
//                let result = tax.nodes_with_label_and_rank(rank_sym, class_name);

//                if let Some(result) = result {
//                    Some(UInt32Chunked::new_vec(PlSmallStr::EMPTY, tree::nodeid_to_u32_vec(result.into())))
//                } else {
//                    Some(UInt32Chunked::from_iter([None]))
//                }
//            })
//            .alias(taxids_new_col),
//    )
//    .with_column(col(predict_label_col).repeat_by(col(predict_label_col).list().len()))
//}

#[cfg(test)]
pub(crate) mod tests {
    use std::sync::Arc;

    use itertools::Itertools;
    use polars::series::Series;

    use crate::{
        csv::polars::{lazyframe_from_string_strip, ExprExt},
        taxonomy::{
            formats::newick::{self},
            tree::nodeid_to_u32_vec,
            LabeledTaxonomy, Taxonomy,
        },
    };

    use super::PhaGenusRecord;

    pub const OUTPUT: &str = "\
        contigs \t predict_label    \t uncertainty
        v1      \t g6/s11           \t 0.024633251
        v2      \t no predication   \t NA
        v3      \t no predication   \t NA
        v4      \t no predication   \t NA
        v5      \t no predication   \t NA
        v6      \t g17              \t 0.0139116105
        v7      \t xxx              \t 0.0139116105
        v8      \t no predication   \t NA
        v9      \t f22/g31          \t 1.0442458e-08
    ";

    #[test]
    fn test_predict_label_to_taxids() -> anyhow::Result<()> {
        use polars::lazy::dsl::col;

        let tax = Arc::new(newick::tests::sample_taxonomy());
        let rank_sym = tax.lookup_rank_sym("genus").unwrap();

        let output = lazyframe_from_string_strip::<PhaGenusRecord>(OUTPUT)?.drop(["uncertainty"]);

        let grouped = super::explode_predict_label(output, "predict_label")
            .with_column(
                col("predict_label")
                    .downcast_map_apply_to_vec(Series::str, move |label| {
                        let taxids = tax.nodes_with_label_and_rank(rank_sym, label)?;
                        Some(nodeid_to_u32_vec(taxids.into()))
                    })
                    .alias("taxids"),
            )
            .explode(["taxids"])
            .group_by_stable(["contigs"])
            .agg([col("*")])
            .collect()?;

        let taxids = grouped
            .column("taxids")?
            .list()?
            .amortized_iter()
            .map(|taxid| {
                let taxid = taxid?;
                let taxid = taxid
                    .as_ref()
                    .u32()
                    .expect("Each taxid list should be an array of u32's");
                Some(taxid.to_vec())
            })
            .collect_vec();

        assert_eq!(
            taxids,
            vec![
                Some(vec![Some(6), Some(11)]),
                Some(vec![None]),
                Some(vec![None]),
                Some(vec![None]),
                Some(vec![None]),
                Some(vec![Some(17)]),
                Some(vec![None]),
                Some(vec![None]),
                Some(vec![Some(22), Some(31)]),
            ]
        );

        Ok(())
    }
}
