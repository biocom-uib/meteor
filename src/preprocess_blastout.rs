use anyhow::Context;
use clap::{Args, ValueEnum};

use polars::datatypes::DataType;
use polars::lazy::prelude::Expr;
use polars::prelude::{
    col, lit, CsvWriter, LazyFrame, SerWriter, UniqueKeepStrategy
};

use crate::csv::filter::Filter;
use crate::tool::blast::blastout::{BlastOutFilter, BlastOutFmt};
use crate::util::{self, writing_new_file_or_stdout};


#[derive(ValueEnum, Debug, Default, Copy, Clone)]
pub enum WeightColAgg {
    Max,
    Mean,
    Median,
    Min,
    #[default]
    Count,
    Product,
    Sum,
}

impl WeightColAgg {
    fn into_expr_fn(self) -> fn(Expr) -> Expr {
        use WeightColAgg::*;

        match self {
            Max => Expr::max,
            Mean => Expr::mean,
            Median => Expr::median,
            Min => Expr::min,
            Count => Expr::count,
            Product => Expr::product,
            Sum => Expr::sum,
        }
    }
}

/// Load the output of NCBI's BLAST+ (blastn) to produce a suitable input file for the assign
/// subcommand.
#[derive(Args)]
pub struct PreprocessBlastOutArgs {
    /// Path to the BLAST+ output file. STDIN is not supported.
    blastout_path: String,

    /// Preprocessed output file. Use '-' to write to stdout.
    #[clap(short, long, default_value = "-")]
    output: String,

    /// BLAST+ output format. Only 6 and 7 are supported at the moment. Columns can be specified
    /// just like in BLAST+, like '7 qaccver saccver pident bitscore'. Irrelevant columns can be
    /// ignored using _ instead, e.g. '7 qaccver saccver _ bitscore'.
    #[clap(long)]
    blast_outfmt: String,

    /// Column name (as specified in outfmt) to use as query identifier.
    #[clap(long)]
    query_id_col: String,

    /// Column name (as specified in outfmt) to use as subject identifier.
    #[clap(long)]
    subject_id_col: String,

    /// Column name (as specified in outfmt) to use as weight.
    #[clap(long)]
    weight_col: Option<String>,

    /// Aggregation function to apply to the selected weight column (if needed).
    #[clap(long, value_enum, default_value_t)]
    weight_col_agg: WeightColAgg,

    /// Filters to apply before processing. Example: --filter 'evalue<=1e-3'
    #[clap(long)]
    filter: Vec<BlastOutFilter>,
}

fn aggregate_subjects(subject_id: Expr) -> Expr {
    subject_id.str().join(";", true)
}

fn group_blast_hits(
    hits: LazyFrame,
    query_id_col: &str,
    subject_id_col: &str,
    stable: bool,
) -> LazyFrame {
    let hits = hits
        .drop_nulls(Some(vec![col(query_id_col), col(subject_id_col)]))
        .select([
            col(query_id_col),
            col(subject_id_col).cast(DataType::String),
        ]);

    let hits = if stable {
        hits.unique_stable(None, UniqueKeepStrategy::First)
            .group_by_stable([col(query_id_col)])
    } else {
        hits.unique(None, UniqueKeepStrategy::First)
            .group_by([col(query_id_col)])
    };

    hits.agg([aggregate_subjects(col(subject_id_col))])
}

fn group_blast_hits_with_weights(
    hits: LazyFrame,
    query_id_col: &str,
    subject_id_col: &str,
    weight_col: &str,
    weight_col_agg: impl FnOnce(Expr) -> Expr,
    stable: bool,
) -> LazyFrame {
    let zipped_subject_col = format!("{subject_id_col}/{weight_col}");

    let hits = hits
        .drop_nulls(Some(vec![
            col(query_id_col),
            col(subject_id_col),
            col(weight_col),
        ]))
        .select([
            col(query_id_col),
            col(subject_id_col).cast(DataType::String),
            col(weight_col),
        ]);

    let hits = if stable {
        hits.group_by_stable([col(query_id_col), col(subject_id_col)])
    } else {
        hits.group_by([col(query_id_col), col(subject_id_col)])
    };

    hits.agg([weight_col_agg(col(weight_col))])
        .select([
            col(query_id_col),
            (col(subject_id_col) + lit("/") + col(weight_col).cast(DataType::String))
                .alias(&zipped_subject_col),
        ])
        .group_by([col(query_id_col)])
        .agg([aggregate_subjects(col(&zipped_subject_col))])
}

pub fn preprocess_blastout(args: PreprocessBlastOutArgs) -> anyhow::Result<()> {
    let df = {
        let format = args
            .blast_outfmt
            .parse::<BlastOutFmt>()
            .context("Error parsing blast outfmt specifier")?;

        let mut df = format.load_lazyframe_from_path(args.blastout_path.as_ref())
            .context("Error loading blast output")?;

        eprintln!("Loaded BLAST+ output with schema {:?}", df.collect_schema()?);

        Filter::apply_polars(df, &args.filter)
    };

    let grouped = if let Some(weight_col) = &args.weight_col {
        let weight_col_agg = args.weight_col_agg.into_expr_fn();

        group_blast_hits_with_weights(
            df,
            &args.query_id_col,
            &args.subject_id_col,
            weight_col,
            weight_col_agg,
            false,
        )
    } else {
        group_blast_hits(df, &args.query_id_col, &args.subject_id_col, false)
    };

    let mut grouped = grouped.collect()?;

    //for col_name in grouped.get_column_names_owned() {
    //    if let DataType::List(inner_dtype) = grouped.column(&col_name)?.dtype() {
    //        let inner_dtype_is_string = inner_dtype.is_string();
    //        grouped.try_apply(&col_name, |series| {
    //            let mut series = Cow::Borrowed(series);
    //            if !inner_dtype_is_string {
    //                series = Cow::Owned(series.cast(&DataType::List(DataType::String.boxed()))?);
    //            }
    //            Ok(series.list()?.join_literal(";", true)?.into_series())
    //        })?;
    //    }
    //}

    writing_new_file_or_stdout!(&args.output, writer => {
        let writer = writer.context("Error creating grouped hits file")?;

        util::ignore_broken_pipe(
            CsvWriter::new(writer)
                .with_separator(b'\t')
                .finish(&mut grouped)
        )?;
    });

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use itertools::Itertools;

    use crate::tool::blast::blastout::{self, BlastOutFmt};

    #[test]
    fn test_group_blast_hits() {
        let fmt = BlastOutFmt::from_str(blastout::tests::BLAST_OUT_FMT).unwrap();

        let df = fmt
            .load_lazyframe_from_static(blastout::tests::BLAST_OUT)
            .unwrap();

        let result = super::group_blast_hits(df, "qaccver", "saccver", true)
            .collect()
            .unwrap();

        let columns = result.get_columns();

        assert_eq!(columns.len(), 2);
        assert_eq!(columns[0].name().as_str(), "qaccver");
        assert_eq!(columns[1].name().as_str(), "saccver");

        let qaccver = Vec::from_iter(columns[0].str().unwrap());
        assert_eq!(qaccver, vec![Some("c1"), Some("c2")]);

        let saccver = columns[1]
            .str()
            .unwrap()
            .iter()
            .map(|group| {
                Some(group?.split(';').map(str::to_owned).sorted().collect_vec())
            })
            .collect_vec();

        assert_eq!(
            saccver,
            vec![
                Some(vec!["CP034340.1".to_owned(), "CP034345.1".to_owned()]),
                Some(vec!["CP034343.1".to_owned(), "CP034345.1".to_owned()]),
            ],
        );
    }
}
