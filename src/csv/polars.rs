use std::{borrow::Cow, path::Path, sync::Arc};

use polars::{
    chunked_array::ChunkedArray,
    datatypes::{
        self, ArrayFromIter, DataType, FalseT, ListChunked, NumericNative, PlSmallStr,
        PolarsDataType, PolarsNumericType, StringType,
    },
    error::PolarsResult,
    lazy::dsl::Expr,
    prelude::{
        arity::unary_elementwise, ChunkedCollectIterExt, IntoColumn, JoinBuilder, JoinCoalesce, JoinType, JoinValidation, LazyCsvReader, LazyFileListReader, LazyFrame, Schema
    },
    series::{IntoSeries, Series},
};

use polars_arrow::{
    array::{Array, PrimitiveArray},
    types::NativeType,
};
use polars_plan::{dsl::GetOutput, plans::ScanSources};

pub use meteor_macros::PolarsSchema;

pub trait PolarsCompatible {
    type DataType: PolarsDataType;
    const DTYPE: DataType;
}

impl PolarsCompatible for String {
    type DataType = datatypes::StringType;
    const DTYPE: DataType = DataType::String;
}

impl PolarsCompatible for &str {
    type DataType = datatypes::StringType;
    const DTYPE: DataType = DataType::String;
}

impl PolarsCompatible for Cow<'_, str> {
    type DataType = datatypes::StringType;
    const DTYPE: DataType = DataType::String;
}

impl PolarsCompatible for u32 {
    type DataType = datatypes::UInt32Type;
    const DTYPE: DataType = DataType::UInt32;
}

impl PolarsCompatible for u64 {
    type DataType = datatypes::UInt64Type;
    const DTYPE: DataType = DataType::UInt64;
}

impl PolarsCompatible for i32 {
    type DataType = datatypes::Int32Type;
    const DTYPE: DataType = DataType::Int32;
}

impl PolarsCompatible for i64 {
    type DataType = datatypes::Int64Type;
    const DTYPE: DataType = DataType::Int64;
}

impl PolarsCompatible for f32 {
    type DataType = datatypes::Float32Type;
    const DTYPE: DataType = DataType::Float32;
}

impl PolarsCompatible for f64 {
    type DataType = datatypes::Float64Type;
    const DTYPE: DataType = DataType::Float64;
}

pub trait PolarsSchema {
    const ALL_FIELDS: &[(&str, &DataType)];

    fn polars_schema() -> Schema;

    fn set_csv_options(reader: LazyCsvReader) -> LazyCsvReader {
        default_tsv_options(reader)
    }
}

pub fn default_tsv_options(reader: LazyCsvReader) -> LazyCsvReader {
    reader.with_has_header(true).with_separator(b'\t')
}

pub fn schema_subset<'s>(
    full_schema: &Schema,
    cols: &[Option<&'s str>],
) -> Result<Schema, &'s str> {
    let mut schema = Vec::new();

    for &col in cols {
        if let Some(col) = col {
            if let Some(dtype) = full_schema.get(col) {
                schema.push((col.into(), dtype.clone()));
            } else {
                return Err(col);
            }
        } else {
            let ignored = format!("__ignored{}", schema.len());
            schema.push((ignored.into(), DataType::String));
        }
    }

    Ok(schema.into_iter().collect())
}

pub fn lazyframe_from<Record: PolarsSchema>(sources: ScanSources) -> PolarsResult<LazyFrame> {
    let reader = LazyCsvReader::new_with_sources(sources)
        .with_schema(Some(Arc::new(Record::polars_schema())));

    Record::set_csv_options(reader).finish()
}

pub fn lazyframe_from_file<Record: PolarsSchema>(path: &Path) -> PolarsResult<LazyFrame> {
    lazyframe_from::<Record>(ScanSources::Paths(Arc::new([path.to_owned()])))
}

/// removes spaces
#[cfg(test)]
pub fn lazyframe_from_string_strip<Record: PolarsSchema>(string: &str) -> PolarsResult<LazyFrame> {
    use regex::{Regex, RegexBuilder};

    let leading = RegexBuilder::new("^ *").multi_line(true).build().unwrap();
    let middle = Regex::new(r" *\t *").unwrap();
    let trailing = RegexBuilder::new(" *$").multi_line(true).build().unwrap();

    let string = leading.replace_all(string, "");
    let string = middle.replace_all(&string, "\t");
    let string = trailing.replace_all(&string, "");
    let bytes_buf = string.into_owned().into_bytes().into();

    lazyframe_from::<Record>(ScanSources::Buffers(Arc::new([bytes_buf])))
}

pub fn full_join_builder(df1: LazyFrame, df2: LazyFrame) -> JoinBuilder {
    df1.join_builder()
        .with(df2)
        .how(JoinType::Full)
        .validate(JoinValidation::ManyToMany)
        .coalesce(JoinCoalesce::CoalesceColumns)
}

// concat two list columns that can be null
pub fn col_list_concat(col1: &str, col2: &str) -> Expr {
    use polars::lazy::dsl::{col, when};

    when(col(col1).is_null()).then(col(col2)).otherwise(
        when(col(col2).is_null())
            .then(col(col1))
            .otherwise(polars::prelude::concat_list([col(col1), col(col2)]).unwrap()),
    )
}

// union two list columns that can be null
pub fn col_list_union(col1: &str, col2: &str) -> Expr {
    use polars::lazy::dsl::{col, when};

    when(col(col1).is_null()).then(col(col2)).otherwise(
        when(col(col2).is_null())
            .then(col(col1))
            .otherwise(col(col1).list().union(col(col2))),
    )
}

#[cfg(test)]
pub fn series_vec<T, Phantom>(lists: impl IntoIterator<Item = T>) -> Vec<Series>
where
    Phantom: ?Sized,
    Series: polars::prelude::NamedFrom<T, Phantom>,
{
    use itertools::Itertools;
    use polars::prelude::NamedFrom;

    lists
        .into_iter()
        .map(|list| Series::new(PlSmallStr::EMPTY, list))
        .collect_vec()
}

#[cfg(test)]
pub fn series_opt_vec<T, Phantom>(lists: impl IntoIterator<Item = Option<T>>) -> Vec<Option<Series>>
where
    Phantom: ?Sized,
    Series: polars::prelude::NamedFrom<T, Phantom>,
{
    use itertools::Itertools;
    use polars::prelude::NamedFrom;

    lists
        .into_iter()
        .map(|opt_list| Some(Series::new(PlSmallStr::EMPTY, opt_list?)))
        .collect_vec()
}

pub trait ExprExt: Sized {
    fn downcast_map_to<T, U, C, F>(self, cast: C, output_dtype: DataType, f: F) -> Self
    where
        C: Fn(&Series) -> PolarsResult<&ChunkedArray<T>> + Send + Sync + 'static,
        F: Fn(&ChunkedArray<T>) -> ChunkedArray<U> + Send + Sync + 'static,
        T: PolarsDataType,
        U: PolarsDataType,
        ChunkedArray<U>: IntoSeries;

    fn downcast_map<T, U, C, F>(self, cast: C, f: F) -> Self
    where
        C: Fn(&Series) -> PolarsResult<&ChunkedArray<T>> + Send + Sync + 'static,
        F: Fn(&ChunkedArray<T>) -> ChunkedArray<U> + Send + Sync + 'static,
        T: PolarsDataType<IsNested = FalseT>,
        U: PolarsDataType<IsNested = FalseT>,
        ChunkedArray<U>: IntoSeries,
    {
        self.downcast_map_to(cast, U::get_dtype(), f)
    }

    fn downcast_map_to_list<T, C, F>(self, cast: C, output_inner: DataType, f: F) -> Self
    where
        C: Fn(&Series) -> PolarsResult<&ChunkedArray<T>> + Send + Sync + 'static,
        F: Fn(&ChunkedArray<T>) -> ListChunked + Send + Sync + 'static,
        T: PolarsDataType,
    {
        self.downcast_map_to(cast, DataType::List(output_inner.boxed()), f)
    }

    fn downcast_map_lists<T, U, IC, IF>(self, inner_cast: IC, inner_f: IF) -> Self
    where
        IC: Fn(&Series) -> PolarsResult<&ChunkedArray<T>> + Send + Sync + 'static,
        IF: Fn(&ChunkedArray<T>) -> ChunkedArray<U> + Send + Sync + 'static,
        T: PolarsDataType,
        U: PolarsDataType,
        ChunkedArray<U>: IntoSeries;

    fn downcast_map_apply<T, B, C, F>(self, cast: C, f: F) -> Self
    where
        C: Fn(&Series) -> PolarsResult<&ChunkedArray<T>> + Send + Sync + 'static,
        F: Fn(T::Physical<'_>) -> Option<B> + Send + Sync + 'static,
        T: PolarsDataType<IsNested = FalseT>,
        B: PolarsCompatible<DataType: PolarsDataType<IsNested = FalseT>>,
        <B::DataType as PolarsDataType>::Array: ArrayFromIter<Option<B>>,
        ChunkedArray<B::DataType>: IntoSeries,
    {
        self.downcast_map(cast, move |chunked| {
            unary_elementwise::<T, B::DataType, _>(chunked, |t| f(t?))
        })
    }

    fn downcast_map_apply_str<F>(self, f: F) -> Self
    where
        F: Fn(&str) -> Option<&str> + Send + Sync + 'static,
    {
        self.downcast_map(Series::str, move |chunked| {
            unary_elementwise::<StringType, StringType, _>(chunked, |t| f(t?))
        })
    }

    fn downcast_map_apply_to_vec<T, B, C, F>(self, cast: C, f: F) -> Self
    where
        C: Fn(&Series) -> PolarsResult<&ChunkedArray<T>> + Send + Sync + 'static,
        F: Fn(T::Physical<'_>) -> Option<Vec<B>> + Send + Sync + 'static,
        T: PolarsDataType<IsNested = FalseT>,
        B: PolarsCompatible + NativeType + NumericNative,
        B::DataType: PolarsNumericType<Native = B, Array = PrimitiveArray<B>>,
        ChunkedArray<B::DataType>: IntoSeries,
    {
        self.downcast_map_apply_to_chunk(cast, move |t| Some(PrimitiveArray::<B>::from_vec(f(t)?)))
    }

    fn downcast_map_apply_to_chunk<T, B, C, F>(self, cast: C, f: F) -> Self
    where
        C: Fn(&Series) -> PolarsResult<&ChunkedArray<T>> + Send + Sync + 'static,
        F: Fn(T::Physical<'_>) -> Option<PrimitiveArray<B>> + Send + Sync + 'static,
        T: PolarsDataType<IsNested = FalseT>,
        B: NumericNative + PolarsCompatible,
    {
        self.downcast_map_to_list(cast, B::DataType::get_dtype(), move |chunked| {
            chunked
                .iter()
                .map(|t| Some::<Box<dyn Array>>(Box::new(f(t?)?)))
                .collect_ca_trusted_with_dtype(
                    PlSmallStr::EMPTY,
                    DataType::List(B::DataType::get_dtype().boxed()),
                )
        })
    }
}

pub trait ExprTryExt: ExprBinExt {
    fn downcast_try_map_to<T, U, C, F>(self, cast: C, output_dtype: DataType, f: F) -> Self
    where
        C: Fn(&Series) -> PolarsResult<&ChunkedArray<T>> + Send + Sync + 'static,
        F: Fn(&ChunkedArray<T>) -> PolarsResult<ChunkedArray<U>> + Send + Sync + 'static,
        T: PolarsDataType,
        U: PolarsDataType,
        ChunkedArray<U>: IntoSeries;

    fn downcast_try_map2_to<T1, T2, U, C1, C2, F>(
        self,
        cast1: C1,
        col2_name: &str,
        cast2: C2,
        output_dtype: DataType,
        f: F,
    ) -> Self
    where
        C1: Fn(&Series) -> PolarsResult<&ChunkedArray<T1>> + Send + Sync + 'static,
        C2: Fn(&Series) -> PolarsResult<&ChunkedArray<T2>> + Send + Sync + 'static,
        F: Fn(&ChunkedArray<T1>, &ChunkedArray<T2>) -> PolarsResult<ChunkedArray<U>>
            + Send
            + Sync
            + 'static,
        T1: PolarsDataType,
        T2: PolarsDataType,
        U: PolarsDataType,
        ChunkedArray<U>: IntoSeries;
}

pub trait ExprBinExt: ExprExt {
    fn downcast_map2_to<T1, T2, U, C1, C2, F>(
        self,
        cast1: C1,
        col2_name: &str,
        cast2: C2,
        output_dtype: DataType,
        f: F,
    ) -> Self
    where
        C1: Fn(&Series) -> PolarsResult<&ChunkedArray<T1>> + Send + Sync + 'static,
        C2: Fn(&Series) -> PolarsResult<&ChunkedArray<T2>> + Send + Sync + 'static,
        F: Fn(&ChunkedArray<T1>, &ChunkedArray<T2>) -> ChunkedArray<U> + Send + Sync + 'static,
        T1: PolarsDataType,
        T2: PolarsDataType,
        U: PolarsDataType,
        ChunkedArray<U>: IntoSeries;

    fn downcast_map2<T1, T2, U, C1, C2, F>(
        self,
        cast1: C1,
        col2_name: &str,
        cast2: C2,
        f: F,
    ) -> Self
    where
        C1: Fn(&Series) -> PolarsResult<&ChunkedArray<T1>> + Send + Sync + 'static,
        C2: Fn(&Series) -> PolarsResult<&ChunkedArray<T2>> + Send + Sync + 'static,
        F: Fn(&ChunkedArray<T1>, &ChunkedArray<T2>) -> ChunkedArray<U> + Send + Sync + 'static,
        T1: PolarsDataType,
        T2: PolarsDataType,
        U: PolarsDataType<IsNested = FalseT>,
        ChunkedArray<U>: IntoSeries,
    {
        self.downcast_map2_to(cast1, col2_name, cast2, U::get_dtype(), f)
    }

    fn downcast_map2_to_list<T1, T2, C1, C2, F>(
        self,
        cast1: C1,
        col2_name: &str,
        cast2: C2,
        output_inner: DataType,
        f: F,
    ) -> Self
    where
        C1: Fn(&Series) -> PolarsResult<&ChunkedArray<T1>> + Send + Sync + 'static,
        C2: Fn(&Series) -> PolarsResult<&ChunkedArray<T2>> + Send + Sync + 'static,
        F: Fn(&ChunkedArray<T1>, &ChunkedArray<T2>) -> ListChunked + Send + Sync + 'static,
        T1: PolarsDataType,
        T2: PolarsDataType,
    {
        self.downcast_map2_to(
            cast1,
            col2_name,
            cast2,
            DataType::List(output_inner.boxed()),
            f,
        )
    }

    fn downcast_map2_apply_to_vec<T1, T2, B, C1, C2, F>(
        self,
        cast1: C1,
        col2_name: &str,
        cast2: C2,
        f: F,
    ) -> Self
    where
        C1: Fn(&Series) -> PolarsResult<&ChunkedArray<T1>> + Send + Sync + 'static,
        C2: Fn(&Series) -> PolarsResult<&ChunkedArray<T2>> + Send + Sync + 'static,
        F: Fn(T1::Physical<'_>, T2::Physical<'_>) -> Option<Vec<B>> + Send + Sync + 'static,
        T1: PolarsDataType<IsNested = FalseT>,
        T2: PolarsDataType<IsNested = FalseT>,
        B: PolarsCompatible + NativeType + NumericNative,
        B::DataType: PolarsNumericType<Native = B, Array = PrimitiveArray<B>>,
        ChunkedArray<B::DataType>: IntoSeries,
    {
        self.downcast_map2_apply_to_chunk(cast1, col2_name, cast2, move |t1, t2| {
            Some(PrimitiveArray::<B>::from_vec(f(t1, t2)?))
        })
    }

    fn downcast_map2_apply_to_chunk<T1, T2, B, C1, C2, F>(
        self,
        cast1: C1,
        col2_name: &str,
        cast2: C2,
        f: F,
    ) -> Self
    where
        C1: Fn(&Series) -> PolarsResult<&ChunkedArray<T1>> + Send + Sync + 'static,
        C2: Fn(&Series) -> PolarsResult<&ChunkedArray<T2>> + Send + Sync + 'static,
        F: Fn(T1::Physical<'_>, T2::Physical<'_>) -> Option<PrimitiveArray<B>>
            + Send
            + Sync
            + 'static,
        T1: PolarsDataType<IsNested = FalseT>,
        T2: PolarsDataType<IsNested = FalseT>,
        B: NumericNative + PolarsCompatible,
    {
        self.downcast_map2_to_list(
            cast1,
            col2_name,
            cast2,
            B::DataType::get_dtype(),
            move |chunked1, chunked2| {
                chunked1
                    .iter()
                    .zip(chunked2.iter())
                    .map(|(t1, t2)| Some::<Box<dyn Array>>(Box::new(f(t1?, t2?)?)))
                    .collect_ca_trusted_with_dtype(
                        PlSmallStr::EMPTY,
                        DataType::List(B::DataType::get_dtype().boxed()),
                    )
            },
        )
    }
}

impl ExprExt for Expr {
    fn downcast_map_to<T, U, C, F>(self, cast: C, output_dtype: DataType, f: F) -> Self
    where
        C: Fn(&Series) -> PolarsResult<&ChunkedArray<T>> + Send + Sync + 'static,
        F: Fn(&ChunkedArray<T>) -> ChunkedArray<U> + Send + Sync + 'static,
        T: PolarsDataType,
        U: PolarsDataType,
        ChunkedArray<U>: IntoSeries,
    {
        self.map(
            move |column| {
                Ok(Some(
                    f(cast(column.as_materialized_series())?).into_column(),
                ))
            },
            GetOutput::from_type(output_dtype),
        )
    }

    fn downcast_map_lists<T, U, IC, IF>(self, inner_cast: IC, inner_f: IF) -> Self
    where
        IC: Fn(&Series) -> PolarsResult<&ChunkedArray<T>> + Send + Sync + 'static,
        IF: Fn(&ChunkedArray<T>) -> ChunkedArray<U> + Send + Sync + 'static,
        T: PolarsDataType,
        U: PolarsDataType,
        ChunkedArray<U>: IntoSeries,
    {
        self.map_list(
            move |column| {
                let result = column
                    .as_materialized_series()
                    .list()?
                    .try_apply_amortized(|group| {
                        let inner = inner_cast(group.as_ref())?;
                        Ok(inner_f(inner).into_series())
                    })?;

                Ok(Some(result.into_column()))
            },
            GetOutput::from_type(DataType::List(U::get_dtype().boxed())),
        )
    }
}

impl ExprBinExt for Expr {
    fn downcast_map2_to<T1, T2, U, C1, C2, F>(
        self,
        cast1: C1,
        col2_name: &str,
        cast2: C2,
        output_dtype: DataType,
        f: F,
    ) -> Self
    where
        C1: Fn(&Series) -> PolarsResult<&ChunkedArray<T1>> + Send + Sync + 'static,
        C2: Fn(&Series) -> PolarsResult<&ChunkedArray<T2>> + Send + Sync + 'static,
        F: Fn(&ChunkedArray<T1>, &ChunkedArray<T2>) -> ChunkedArray<U> + Send + Sync + 'static,
        T1: PolarsDataType,
        T2: PolarsDataType,
        U: PolarsDataType,
        ChunkedArray<U>: IntoSeries,
    {
        use polars::prelude::col;

        self.map_many(
            move |columns| {
                let col1 = cast1(columns[0].as_materialized_series())?;
                let col2 = cast2(columns[1].as_materialized_series())?;

                Ok(Some(f(col1, col2).into_column()))
            },
            &[col(col2_name)],
            GetOutput::from_type(output_dtype),
        )
    }
}

impl ExprTryExt for Expr {
    fn downcast_try_map_to<T, U, C, F>(self, cast: C, output_dtype: DataType, f: F) -> Self
    where
        C: Fn(&Series) -> PolarsResult<&ChunkedArray<T>> + Send + Sync + 'static,
        F: Fn(&ChunkedArray<T>) -> PolarsResult<ChunkedArray<U>> + Send + Sync + 'static,
        T: PolarsDataType,
        U: PolarsDataType,
        ChunkedArray<U>: IntoSeries,
    {
        self.map(
            move |column| {
                Ok(Some(
                    f(cast(column.as_materialized_series())?)?.into_column(),
                ))
            },
            GetOutput::from_type(output_dtype),
        )
    }

    fn downcast_try_map2_to<T1, T2, U, C1, C2, F>(
        self,
        cast1: C1,
        col2_name: &str,
        cast2: C2,
        output_dtype: DataType,
        f: F,
    ) -> Self
    where
        C1: Fn(&Series) -> PolarsResult<&ChunkedArray<T1>> + Send + Sync + 'static,
        C2: Fn(&Series) -> PolarsResult<&ChunkedArray<T2>> + Send + Sync + 'static,
        F: Fn(&ChunkedArray<T1>, &ChunkedArray<T2>) -> PolarsResult<ChunkedArray<U>>
            + Send
            + Sync
            + 'static,
        T1: PolarsDataType,
        T2: PolarsDataType,
        U: PolarsDataType,
        ChunkedArray<U>: IntoSeries,
    {
        use polars::prelude::col;

        self.map_many(
            move |columns| {
                let col1 = cast1(columns[0].as_materialized_series())?;
                let col2 = cast2(columns[1].as_materialized_series())?;

                Ok(Some(f(col1, col2)?.into_column()))
            },
            &[col(col2_name)],
            GetOutput::from_type(output_dtype),
        )
    }
}
