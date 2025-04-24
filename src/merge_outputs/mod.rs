pub mod taxid_reduction;

pub mod host_prediction;
pub mod virus_taxonomy;

use std::{
    ops::Deref,
    path::Path,
    sync::{Arc, LazyLock},
};

use clap::Args;
use itertools::Itertools;
use polars::{
    datatypes::DataType,
    frame::DataFrame,
    lazy::dsl::Expr,
    prelude::{Column, IntoLazy, LazyFrame, ListNameSpaceExtension},
    series::Series,
};
use polars_arrow::array::PrimitiveArray;
use taxid_reduction::TaxIdReduction;

use crate::{
    csv::{filter::Filter, polars::ExprExt},
    taxonomy::{
        tree::{
            self,
            walk::{LazyLcaCache, LcaCache},
        },
        LabeledTaxonomy, NodeId, Taxonomy,
    },
    tool::vpf_class::{self, VpfClassRank, VpfClassRecordFilter},
};


pub mod columns {
    pub const VIRUS_NAME: &str = "virus_name";
    pub const TAXIDS: &str = "taxids";
    pub const CLASS_NAMES: &str = "class_names";

    pub const MERGED_TAXIDS: &str = "merged_taxids";
    pub const MERGED_CLASS_NAMES: &str = "merged_class_names";
}

use columns::*;

#[derive(Args)]
pub struct MergeCommonArgs {
    /// Reduce the set of possible taxids with a chosen strategy.
    #[clap(long, value_enum)]
    pub reduce_taxids: Option<TaxIdReduction>,

    /// Keep names in addition to taxids.
    #[clap(long)]
    pub keep_names: bool,
}

pub struct MergeOutputContext<'a, Tax: Taxonomy + Send + Sync + 'static> {
    args: &'a MergeCommonArgs,
    tax: Arc<Tax>,
    lca_cache: Arc<LazyLcaCache<Tax>>,
}

impl<'a, Tax> Clone for MergeOutputContext<'a, Tax>
where
    Tax: Taxonomy + Send + Sync + 'static,
{
    fn clone(&self) -> Self {
        Self {
            args: self.args,
            tax: Arc::clone(&self.tax),
            lca_cache: Arc::clone(&self.lca_cache),
        }
    }
}

impl<'a, Tax> MergeOutputContext<'a, Tax>
where
    Tax: Taxonomy + Send + Sync + 'static,
{
    /// Precomputes the LCA cache concurrently if needed
    pub fn new(args: &'a MergeCommonArgs, tax: Arc<Tax>) -> MergeOutputContext<'a, Tax> {
        let lca_cache = precompute_lca_cache_concurrently(args, Arc::clone(&tax));

        MergeOutputContext { args, tax, lca_cache }
    }

    pub fn args(&self) -> &'a MergeCommonArgs {
        self.args
    }

    pub fn taxonomy(&self) -> Arc<Tax> {
        Arc::clone(&self.tax)
    }

    pub fn lca_cache(&self) -> Arc<LazyLcaCache<Tax>> {
        Arc::clone(&self.lca_cache)
    }
}


/// Maps `class_names: str` to:
///   - If a mapping exists, `{class_names}: list[u32] ([taxid1, ...])`
///   - If no mappings exist, `{class_names}: list[u32] ([null])`
pub fn class_names_to_taxid_list<TaxPtr, Tax>(
    class_names: Expr,
    tax: TaxPtr,
    rank_sym: Option<Tax::RankSym>,
) -> Expr
where
    TaxPtr: Deref<Target = Tax> + Send + Sync + 'static,
    Tax: LabeledTaxonomy,
{
    if let Some(rank_sym) = rank_sym {
        class_names.downcast_map_apply_to_chunk(Series::str, move |class_name| {
            let result = tax.nodes_with_label_and_rank(rank_sym, class_name);

            if let Some(result) = result {
                Some(PrimitiveArray::from_vec(tree::nodeid_to_u32_vec(result.into())))
            } else {
                Some(PrimitiveArray::from_iter([None]))
            }
        })
    } else {
        class_names.downcast_map_apply_to_chunk(Series::str, move |class_name| {
            let mut iter = tax.nodes_with_label(class_name).map(u32::from).peekable();

            if iter.peek().is_some() {
                Some(PrimitiveArray::from_values(iter))
            } else {
                Some(PrimitiveArray::from_iter([None]))
            }
        })
    }
}

/// Input columns: `..., {class_names}: str, {taxids}: [u32], ...`
/// Output columns: `..., {class_names}: [str], {taxids}: [u32], ...`
/// where
///     `{class_names}.list().len() == {class_names}.list().len()`
pub fn repeat_class_names_for_taxids(class_names: Expr, taxids: Expr) -> Expr {
    class_names.repeat_by(taxids.list().len())
}

/// To be used after combining multiple VPF-Class output files
/// Always uses `TaxIdReduction::KeepMostSpecific`.
pub fn reduce_vpf_class_taxids<Tax>(ctx: &MergeOutputContext<Tax>, df: LazyFrame) -> LazyFrame
where
    Tax: Taxonomy + Send + Sync + 'static,
{
    use polars::lazy::dsl::col;

    if ctx.args().keep_names {
        let indices = TaxIdReduction::most_specific_indices(col(TAXIDS), ctx.lca_cache());

        df.with_column(indices.alias("indices")).select([
            col(VIRUS_NAME),
            col(CLASS_NAMES).list().gather(col("indices"), false),
            col(TAXIDS).list().gather(col("indices"), false),
        ])
    } else {
        df.with_column(
            TaxIdReduction::KeepMostSpecific.apply(col(TAXIDS), ctx.lca_cache()),
        )
    }
}

/// Columns:
///   - If `ctx.args().keep_names`: `virus_name: str, class_names: list[str], taxids: list[str]`
///   - Otherwise: `virus_name: str, taxids: list[str]`
///
/// Where
///     `{class_names}.list().len() == {class_names}.list().len()`
pub fn load_vpf_class_for_merge<Tax>(
    ctx: &MergeOutputContext<Tax>,
    output_dir: &Path,
    ranks: &[impl Into<VpfClassRank> + Copy],
    filters: &[VpfClassRecordFilter],
    reduce_taxids: bool,
) -> anyhow::Result<LazyFrame>
where
    Tax: LabeledTaxonomy + Send + Sync + 'static,
{
    let load_rank = |rank: VpfClassRank| -> anyhow::Result<LazyFrame> {
        use polars::lazy::dsl::col;

        let tax = ctx.taxonomy();

        let rank_sym = rank.lookup_sym(&*tax)?;

        let mut df = vpf_class::lazyframe_from_rank(output_dir, rank)?;

        df = Filter::apply_polars(df, filters)
            .select([col(VIRUS_NAME), col("class_name").alias(CLASS_NAMES)])
            .with_column(
                class_names_to_taxid_list(col(CLASS_NAMES), tax, Some(rank_sym)).alias(TAXIDS),
            );

        if ctx.args().keep_names {
            df = df
                .select([col(VIRUS_NAME), col(CLASS_NAMES), col(TAXIDS)])
                .with_column(repeat_class_names_for_taxids(col(CLASS_NAMES), col(TAXIDS)));
        } else {
            df = df.drop([CLASS_NAMES]);
        }

        df = df
            .group_by([VIRUS_NAME])
            .agg([col("*").explode()]);

        Ok(df)
    };

    let mut result = ranks
        .iter()
        .map(|&rank| load_rank(rank.into()))
        .process_results(|iter| {
            let cols: &[&str] = if ctx.args().keep_names {
                &[CLASS_NAMES, TAXIDS]
            } else {
                &[TAXIDS]
            };
            iter.reduce(|df1, df2| vpf_class::join_lists_by_virus_name(df1, df2, cols))
        })?
        .ok_or_else(|| anyhow::anyhow!("No ranks specified for VPF-Class predictions"))?;

    if reduce_taxids {
        result = reduce_vpf_class_taxids(ctx, result)
    }

    Ok(result)
}


pub fn precompute_lca_cache_concurrently<Tax>(
    args: &MergeCommonArgs,
    tax: Arc<Tax>,
) -> Arc<LazyLcaCache<Tax>>
where
    Tax: Taxonomy + Send + Sync + 'static,
{
    let lca_cache = {
        let tax = Arc::clone(&tax);
        Arc::new(LcaCache::compute_lazy(tax))
    };

    if args.reduce_taxids.is_some() {
        let lca_cache = Arc::clone(&lca_cache);

        rayon::spawn(move || {
            LazyLock::force(&*lca_cache);
        });
    }

    lca_cache
}

#[derive(Clone)]
pub struct NamedLazyFrame {
    name: String,
    df: LazyFrame,
}

impl NamedLazyFrame {
    pub fn new(name: impl Into<String>, df: LazyFrame) -> Self {
        Self { name: name.into(), df }
    }
}

pub fn merge_virus_predictions<Tax>(
    ctx: &MergeOutputContext<Tax>,
    dfs: Vec<NamedLazyFrame>,
) -> LazyFrame
where
    Tax: LabeledTaxonomy + Send + Sync + 'static,
{
    use polars::lazy::dsl::{col, cols};

    let init = DataFrame::new(vec![
        Column::new_empty(VIRUS_NAME.into(), &DataType::String),
        Column::new_empty(
            MERGED_TAXIDS.into(),
            &DataType::List(DataType::UInt32.boxed()),
        ),
    ])
    .unwrap();

    let mut result = dfs.into_iter().fold(init.lazy(), |acc, named_df| {
        let NamedLazyFrame { name: prefix, df } = named_df;

        let taxids_col = format!("{prefix}_{TAXIDS}");
        let class_names_col = format!("{prefix}_{CLASS_NAMES}");

        let df = if ctx.args().keep_names {
            df.select([
                col(VIRUS_NAME),
                col(TAXIDS).alias(&taxids_col),
                col(CLASS_NAMES).alias(&class_names_col),
            ])
        } else {
            df.select([col(VIRUS_NAME), col(TAXIDS).alias(&taxids_col)])
        };

        crate::csv::polars::full_join_builder(acc, df)
            .on([col(VIRUS_NAME)])
            .finish()
            .with_column(
                crate::csv::polars::col_list_union(MERGED_TAXIDS, &taxids_col)
                    .alias(MERGED_TAXIDS),
            )
    });

    if let Some(reduction) = ctx.args().reduce_taxids {
        result = result.with_column(reduction.apply(col(MERGED_TAXIDS), ctx.lca_cache()));
    }

    if ctx.args().keep_names {
        let tax = ctx.taxonomy();

        let taxid_to_label = col("").downcast_map_apply(Series::u32, move |taxid| {
            tax.some_label_of(NodeId(taxid)).map(str::to_owned)
        });

        result = result
            .with_column(
                col(MERGED_TAXIDS)
                    .list()
                    .eval(taxid_to_label, true)
                    .alias(MERGED_CLASS_NAMES),
            )
            .select([
                cols([VIRUS_NAME, MERGED_TAXIDS, MERGED_CLASS_NAMES]),
                col("*").exclude([VIRUS_NAME, MERGED_TAXIDS, MERGED_CLASS_NAMES]),
            ]);
    }

    result
}

pub fn prepare_merged_virus_predictions_for_writing(df: LazyFrame) -> LazyFrame {
    use polars::lazy::dsl::{col, dtype_col, lit};

    df
        .with_column(
            dtype_col(&DataType::List(DataType::UInt32.boxed()))
                .list()
                .eval(col("").cast(DataType::String), false),
        )
        .with_column(
            dtype_col(&DataType::List(DataType::String.boxed()))
                .list()
                .eval(col("").fill_null(lit("")), false)
                .list()
                .join(lit(";"), false),
        )
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use polars::{error::PolarsResult, prelude::IntoLazy};

    use crate::{
        csv::polars::{series_opt_vec, series_vec},
        merge_outputs::{MergeOutputContext, NamedLazyFrame},
        taxonomy::{formats::newick, Taxonomy},
    };

    use super::MergeCommonArgs;

    #[test]
    fn test_class_names_to_taxid_list() -> PolarsResult<()> {
        use polars::lazy::dsl::col;

        let tax = Arc::new(newick::tests::ambiguous_taxonomy());
        let genus_sym = tax.lookup_rank_sym("genus").unwrap();

        let df = polars::df!(
            "virus_name" => ["v1", "v2", "v3"],
            "class_name" => ["x", "y", "7"],
        )?;

        let result_none = df
            .clone()
            .lazy()
            .with_column(
                super::class_names_to_taxid_list(col("class_name"), Arc::clone(&tax), None)
                    .list()
                    .sort(Default::default())
                    .alias("taxids"),
            )
            .collect()?;

        let expected_none = polars::df!(
            "virus_name" => ["v1", "v2", "v3"],
            "class_name" => ["x", "y", "7"],
            "taxids" => series_vec([
                vec![Some(5), Some(6u32)],
                vec![None],
                vec![Some(7)]
            ]),
        )?;

        assert_eq!(result_none, expected_none);

        let result_some_genus = df
            .lazy()
            .with_column(
                super::class_names_to_taxid_list(col("class_name"), tax, Some(genus_sym))
                    .list()
                    .sort(Default::default())
                    .alias("taxids"),
            )
            .collect()?;

        let expected_some_genus = polars::df!(
            "virus_name" => ["v1", "v2", "v3"],
            "class_name" => ["x", "y", "7"],
            "taxids" => series_vec([
                vec![Some(6u32)],
                vec![None],
                vec![Some(7)]
            ]),
        )?;

        assert_eq!(result_some_genus, expected_some_genus);

        Ok(())
    }

    #[test]
    fn test_repeat_class_names_for_taxids() -> PolarsResult<()> {
        use polars::lazy::dsl::col;

        let taxids = series_vec([
            vec![None, Some(1u32)],
            vec![],
            vec![Some(2)],
        ]);

        let df = polars::df!(
            "class_names" => ["a", "b", "c"],
            "taxids" => &taxids,
        )?;

        let expected = polars::df!(
            "class_names" => series_vec([
                vec!["a"; 2],
                vec![],
                vec!["c"],
            ]),
            "taxids" => &taxids,
        )?;

        let repeated = df
            .lazy()
            .with_column(super::repeat_class_names_for_taxids(
                col("class_names"),
                col("taxids"),
            ))
            .collect()?;

        assert_eq!(repeated, expected);

        Ok(())
    }

    #[test]
    fn test_merge_virus_predictions() -> PolarsResult<()> {
        let tax = Arc::new(newick::tests::ambiguous_taxonomy());

        let df1 = polars::df!(
            "virus_name" => ["v1", "v2", "v3"],
            "class_names" => series_vec([vec!["x", "x"], vec!["12"], vec!["30"]]),
            "taxids" => series_vec([vec![5u32, 6], vec![12], vec![30]]),
        )?;

        let df2 = polars::df!(
            "virus_name" => ["v1", "v2", "v4"],
            "class_names" => series_vec([vec!["z"], vec!["11", "y"], vec!["31"]]),
            "taxids" => series_vec([vec![None], vec![Some(11u32), None], vec![Some(31)]]),
        )?;

        let dfs = vec![
            NamedLazyFrame::new("df1", df1.clone().lazy()),
            NamedLazyFrame::new("df2", df2.clone().lazy()),
        ];

        let args = MergeCommonArgs {
            reduce_taxids: Some(super::TaxIdReduction::KeepMostSpecific),
            keep_names: false,
        };

        let ctx = MergeOutputContext::new(&args, tax);

        let result = super::merge_virus_predictions(&ctx, dfs.clone())
            .sort(["virus_name"], Default::default())
            .collect()?;

        let expected = polars::df!(
            "virus_names" => ["v1", "v2", "v3", "v4"],
            "merged_taxids" => series_vec([
                vec![Some(6)],
                vec![Some(12), Some(11)],
                vec![Some(30)],
                vec![Some(31)],
            ]),
            "df1_taxids" => series_opt_vec([
                Some(vec![Some(5u32), Some(6)]),
                Some(vec![Some(12)]),
                Some(vec![Some(30)]),
                None,
            ]),
            "df2_taxids" => series_opt_vec([
                Some(vec![None]),
                Some(vec![Some(11), None]),
                None,
                Some(vec![Some(31)]),
            ])
        )?;

        assert_eq!(result, expected);

        let args = MergeCommonArgs {
            keep_names: true,
            ..args
        };

        let ctx = MergeOutputContext { args: &args, ..ctx };

        let result_with_names = super::merge_virus_predictions(&ctx, dfs)
            .sort(["virus_name"], Default::default())
            .collect()?;

        let expected_with_names = polars::df!(
            "virus_names" => ["v1", "v2", "v3", "v4"],
            "merged_taxids" => series_vec([
                vec![Some(6)],
                vec![Some(12), Some(11)],
                vec![Some(30)],
                vec![Some(31)],
            ]),
            "merged_class_names" => series_vec([
                vec!["x"],
                vec!["12", "11"],
                vec!["30"],
                vec!["31"],
            ]),
            "df1_taxids" => series_opt_vec([
                Some(vec![Some(5u32), Some(6)]),
                Some(vec![Some(12)]),
                Some(vec![Some(30)]),
                None,
            ]),
            "df1_class_names" => series_opt_vec([
                Some(vec!["x", "x"]),
                Some(vec!["12"]),
                Some(vec!["30"]),
                None,
            ]),
            "df2_taxids" => series_opt_vec([
                Some(vec![None]),
                Some(vec![Some(11), None]),
                None,
                Some(vec![Some(31)]),
            ]),
            "df2_class_names" => series_opt_vec([
                Some(vec!["z"]),
                Some(vec!["11", "y"]),
                None,
                Some(vec!["31"]),
            ]),
        )?;

        assert_eq!(result_with_names, expected_with_names);

        Ok(())
    }
}
