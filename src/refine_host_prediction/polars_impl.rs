use std::sync::Arc;

use itertools::Itertools;
use lending_iterator::higher_kinded_types::HKTRef;
use polars::{
    chunked_array::StructChunked,
    datatypes::{
        ArrowDataType, BooleanChunked, DataType, Field, PlSmallStr, PolarsDataType, StringChunked,
        StringType, UInt32Type, UInt32Chunked,
    },
    error::PolarsResult,
    prelude::{
        ChunkedBuilder, CsvWriterOptions, LazyCsvReader, LazyFileListReader, LazyFrame,
        PrimitiveChunkedBuilder, SerializeOptions, StringChunkedBuilder,
    },
    series::{IntoSeries, Series},
};
use polars_arrow::{array::BooleanArray, bitmap::MutableBitmap};
use polars_plan::dsl::GetOutput;
use string_interner::DefaultSymbol;

use crate::{
    taxonomy::Taxonomy,
    util::{self, interners::{Interned, Resolver}},
};

use super::{Aggregation, HostSymbol, MetagenomicEvidence, RefineHostPredictionArgs, VirusSymbol};

trait SeriesBuilder {
    type DataType: PolarsDataType;

    fn with_capacity(name: PlSmallStr, capacity: usize) -> Self;
    fn append_list(&mut self, count: usize, f: impl FnOnce() -> String);
    fn append_null(&mut self);
    fn finish(self) -> Series;
}

struct CountingAggregator(PrimitiveChunkedBuilder<UInt32Type>);

impl SeriesBuilder for CountingAggregator {
    type DataType = UInt32Type;

    fn with_capacity(name: PlSmallStr, capacity: usize) -> Self {
        Self(PrimitiveChunkedBuilder::new(name, capacity))
    }

    fn append_list(&mut self, count: usize, _f: impl FnOnce() -> String) {
        self.0.append_value(count as u32);
    }

    fn append_null(&mut self) {
        self.0.append_null();
    }

    fn finish(self) -> Series {
        self.0.finish().into_series()
    }
}

struct ListAggregator(StringChunkedBuilder);

impl SeriesBuilder for ListAggregator {
    type DataType = StringType;

    fn with_capacity(name: PlSmallStr, capacity: usize) -> Self {
        Self(StringChunkedBuilder::new(name, capacity))
    }

    fn append_list(&mut self, _count: usize, f: impl FnOnce() -> String) {
        self.0.append_value(f());
    }

    fn append_null(&mut self) {
        self.0.append_null();
    }

    fn finish(self) -> Series {
        self.0.finish().into_series()
    }
}

const KEEP: Field = Field {
    name: PlSmallStr::from_static("keep"),
    dtype: DataType::Boolean,
};
const PRESENT_CONTIGS: Field = Field {
    name: PlSmallStr::from_static("present_contigs"),
    dtype: DataType::String,
};
const CRISPR_MATCHES: Field = Field {
    name: PlSmallStr::from_static("crispr_matches"),
    dtype: DataType::String,
};

fn enrichment_fields<VI, HI>(evidence: &MetagenomicEvidence<VI, HI>) -> Vec<Field>
where
    VI: Resolver,
    HI: Resolver,
{
    let mut fields = Vec::with_capacity(3);

    fields.extend([KEEP, PRESENT_CONTIGS]);

    if evidence.has_crispr {
        fields.push(CRISPR_MATCHES);
    }

    fields
}

fn enrichment_dtype<VI, HI>(evidence: &MetagenomicEvidence<VI, HI>) -> DataType
where
    VI: Resolver,
    HI: Resolver,
{
    DataType::Struct(enrichment_fields(evidence))
}

fn refine_and_enrich_chunked<Tax, VI, HI, ContigAgg, CrisprAgg>(
    evidence: &MetagenomicEvidence<VI, HI>,
    tax: &Tax,
    virus_names: &StringChunked,
    taxids: &UInt32Chunked,
) -> PolarsResult<StructChunked>
where
    Tax: Taxonomy,
    ContigAgg: SeriesBuilder,
    CrisprAgg: SeriesBuilder,
    VI: Interned<Symbol = DefaultSymbol, ValueHKT = HKTRef<str>>,
    HI: Interned<Symbol = DefaultSymbol, ValueHKT = HKTRef<str>>,
{
    let len = taxids.len();

    let mut keep_validity = MutableBitmap::with_capacity(len);
    let mut keep = MutableBitmap::with_capacity(len);

    let mut contigs = ContigAgg::with_capacity(PRESENT_CONTIGS.name, len);
    let mut crisprs = CrisprAgg::with_capacity(CRISPR_MATCHES.name, len);

    for (virus_name, taxid) in virus_names.iter().zip(taxids) {
        let Some(node) = taxid.and_then(|taxid| tax.fixup_node(taxid)) else {
            keep_validity.push(false);
            keep.push(false);
            contigs.append_null();
            crisprs.append_null();
            continue;
        };

        let virus_sym = virus_name
            .and_then(|virus_name| evidence.virus_interner.get(virus_name))
            .map(VirusSymbol);

        if let Some(pred_evidence) = evidence.refine_and_enrich(tax, virus_sym, node) {
            keep_validity.push(true);
            keep.push(true);

            let contig_list = || {
                pred_evidence
                    .assigned_contigs
                    .iter()
                    .flat_map(|HostSymbol(contig_sym)| evidence.host_interner.resolve(*contig_sym))
                    .join(";")
            };

            contigs.append_list(pred_evidence.assigned_contigs.len(), contig_list);

            let crispr_list = || {
                pred_evidence
                    .crispr_matches
                    .iter()
                    .flat_map(|HostSymbol(host_sym)| evidence.virus_interner.resolve(*host_sym))
                    .join(";")
            };

            crisprs.append_list(pred_evidence.crispr_matches.len(), crispr_list);
        } else {
            keep_validity.push(true);
            keep.push(false);
            contigs.append_null();
            crisprs.append_null();
        }
    }

    let mut series = if evidence.has_crispr {
        Vec::with_capacity(3)
    } else {
        Vec::with_capacity(2)
    };

    let keep = BooleanArray::new(
        ArrowDataType::Boolean,
        keep.freeze(),
        Some(keep_validity.freeze()),
    );

    series.push(BooleanChunked::with_chunk(KEEP.name, keep).into_series());

    series.push(contigs.finish());

    if evidence.has_crispr {
        series.push(crisprs.finish());
    }

    StructChunked::from_series(PlSmallStr::EMPTY, &series)
}


fn load_host_prediction(args: &RefineHostPredictionArgs) -> anyhow::Result<(LazyFrame, bool)> {
    use polars::lazy::dsl::{col, lit};

    let mut prediction = {
        let pathbuf = args.host_prediction.as_path().to_path_buf();

        let reader = LazyCsvReader::new_paths(Arc::new([pathbuf]));

        crate::csv::polars::default_tsv_options(reader)
            .with_chunk_size(5000)
            .finish()?
    };

    use crate::merge_outputs::host_prediction::columns;

    let prediction_has_class_names = args.refinement_options.keep_names
        && prediction
            .collect_schema()?
            .contains(columns::MERGED_CLASS_NAMES);

    if prediction_has_class_names {
        prediction = prediction
            .select([
                col(columns::VIRUS_NAME),
                col(columns::MERGED_TAXIDS).alias("predicted_taxid"),
                col(columns::MERGED_CLASS_NAMES).alias("predicted_class_name"),
            ])
            .with_column(col("predicted_class_names").str().split(lit(";")))
    } else {
        prediction = prediction.select([
            col(columns::VIRUS_NAME),
            col(columns::MERGED_TAXIDS).alias("predicted_taxid"),
        ]);
    }

    prediction = prediction.with_column(
        col("predicted_taxid")
            .str()
            .split(lit(";"))
            .cast(DataType::List(DataType::UInt32.boxed())),
    );

    Ok((prediction, prediction_has_class_names))
}

fn adjust_predictions_to_rank<Tax>(
    tax: Arc<Tax>,
    rank_sym: Tax::RankSym,
    prediction: LazyFrame,
    has_class_names: bool,
) -> LazyFrame
where
    Tax: LabeledTaxonomy + Send + Sync + 'static,
{
    use polars::lazy::dsl::col;

    if has_class_names {
        unimplemented!()
    } else {
        prediction.with_column(
            col("predicted_taxids").map(|series| {

            })
        )
    }
}


pub(super) fn refine_host_prediction_with_tax_impl<Tax, VI, HI>(
    args: &RefineHostPredictionArgs,
    tax: Arc<Tax>,
    evidence: MetagenomicEvidence<VI, HI>,
) -> anyhow::Result<()>
where
    Tax: LabeledTaxonomy + Send + Sync + 'static,
    VI: Interned<Symbol = DefaultSymbol, ValueHKT = HKTRef<str>> + Send + Sync + 'static,
    HI: Interned<Symbol = DefaultSymbol, ValueHKT = HKTRef<str>> + Send + Sync + 'static,
{
    use polars::lazy::dsl::col;

    let output_type = GetOutput::from_type(enrichment_dtype(&evidence));

    let opts = args.refinement_options.clone();

    // TODO: Run in parallel. Requires https://github.com/pola-rs/polars/issues/7243
    // Until then, disabling feature = refine_host_prediction_polars falls back to csv + rayon.

    let (prediction, has_class_names) = load_host_prediction(args)?;

    let prediction = if let Some(rank) = &args.refinement_options.rank {
        let rank_sym = tax.lookup_rank_sym(rank).ok_or_else(|| {
            anyhow::anyhow!("Rank not found in the taxonomy: {}", rank)
        })?;

        adjust_predictions_to_rank(tax, rank_sym, prediciton, has_class_names);
    } else {
        if has_class_names {
            prediction.explode(["predicted_taxid", "predicted_class_name"])
        } else {
            prediction.explode(["predicted_taxid"])
        };
    };

    prediction
        .with_column(
            col("predicted_taxid")
                .map_many(
                    move |series| {
                        let taxids = series[0].u32()?;
                        let virus_name = series[1].str()?;

                        macro_rules! refine_and_enrich_chunked_with {
                            ($contig_agg:ty, $crispr_agg:ty) => {
                                refine_and_enrich_chunked::<_, _, _, $contig_agg, $crispr_agg>(&evidence, &*tax, virus_name, taxids)
                            }
                        }

                        use Aggregation::*;

                        let result = match (opts.agg_contigs, opts.agg_crispr) {
                            (List,  List) => refine_and_enrich_chunked_with!(ListAggregator, ListAggregator)?,
                            (List,  Count) => refine_and_enrich_chunked_with!(ListAggregator, CountingAggregator)?,
                            (Count, List) => refine_and_enrich_chunked_with!(CountingAggregator, ListAggregator)?,
                            (Count, Count) => refine_and_enrich_chunked_with!(CountingAggregator, CountingAggregator)?,
                        };

                        Ok(Some(result.into_series()))
                    },
                    &[col("virus_name")],
                    output_type
                )
                .alias("evidence"),
        )
        .unnest([col("evidence")])
        .filter(col("keep"))
        .drop(["keep"]);

    util::ignore_broken_pipe({
        prediction.with_streaming(true).sink_csv(
            &args.output,
            CsvWriterOptions {
                include_header: true,
                serialize_options: SerializeOptions {
                    separator: b'\t',
                    ..Default::default()
                },
                ..Default::default()
            },
        )
    })?;

    Ok(())
}
