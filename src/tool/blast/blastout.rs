use std::ops::Not;
use std::path::Path;
use std::str::FromStr;
use std::sync::Arc;

use csv::StringRecord;
use itertools::Itertools;
use ::polars::lazy::prelude::LazyCsvReader;
use ::polars::prelude::{
    col, LazyFileListReader, LazyFrame, NullValues, Schema,
};
use polars_plan::plans::ScanSources;

use crate::csv::filter::{FilterableRecord, FromStrField};
use crate::csv::polars::{self, PolarsSchema};


#[allow(dead_code)]
#[derive(FilterableRecord, PolarsSchema)]
#[filter(alias_field, alias_filter, polars_literal, as_string = String)]
pub struct BlastOutFullRecord {
    pub qseqid: String,
    pub qgi: String,
    pub qacc: String,
    pub qaccver: String,
    pub qlen: i32,
    pub sseqid: String,
    pub sallseqid: String,
    pub sgi: String,
    pub sallgi: String,
    pub sacc: String,
    pub saccver: String,
    pub sallacc: String,
    pub slen: i32,
    pub qstart: i32,
    pub qend: i32,
    pub sstart: i32,
    pub send: i32,
    pub qseq: String,
    pub sseq: String,
    pub evalue: f32,
    pub bitscore: f32,
    pub score: f32,
    pub length: i32,
    pub pident: f32,
    pub nident: i32,
    pub mismatch: i32,
    pub positive: i32,
    pub gapopen: i32,
    pub gaps: i32,
    pub ppos: f32,
    pub frames: String,
    pub qframe: String,
    pub sframe: String,
    pub btop: String,
    pub staxid: u32,
    pub ssciname: String,
    pub scomname: String,
    pub sblastname: String,
    pub sskingdom: String,
    pub staxids: String,
    pub sscinames: String,
    pub scomnames: String,
    pub sblastnames: String,
    pub sskingdoms: String,
    pub stitle: String,
    pub salltitles: String,
    pub sstrand: String,
    pub qcovs: f32,
    pub qcovhsp: f32,
    pub qcovus: f32,
}

pub type BlastOutField = BlastOutFullRecordField;
pub type BlastOutFilter = BlastOutFullRecordFilter;

impl BlastOutFullRecord {
    pub const DEFAULT_BLASTN_COLUMNS: &'static [&'static str] = &[
        Self::QACCVER.0,
        Self::SACCVER.0,
        Self::PIDENT.0,
        Self::LENGTH.0,
        Self::MISMATCH.0,
        Self::GAPOPEN.0,
        Self::QSTART.0,
        Self::QEND.0,
        Self::SSTART.0,
        Self::SEND.0,
    ];
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BlastOutFmt {
    // TSV without header
    Six {
        delimiter: char,
        present_columns: Vec<String>,
    },

    // like Six, but with comments
    Seven {
        delimiter: char,
        present_columns: Vec<String>,
    },
}

impl BlastOutFmt {
    pub fn comment_char(&self) -> Option<char> {
        match &self {
            Self::Six { .. } => None,
            Self::Seven { .. } => Some('#'),
        }
    }

    pub fn comment_prefix(&self) -> Option<String> {
        Some(self.comment_char()?.to_string())
    }

    pub fn delimiter(&self) -> char {
        match &self {
            Self::Six { delimiter, .. } => *delimiter,
            Self::Seven { delimiter, .. } => *delimiter,
        }
    }

    pub fn present_columns(&self) -> &[String] {
        match &self {
            Self::Six { present_columns, .. } => present_columns,
            Self::Seven { present_columns, .. } => present_columns,
        }
    }

    pub fn column_index(&self, column: &str) -> Option<usize> {
        self.present_columns().iter().position(|c| c == column)
    }

    pub fn csv_reader_builder(&self) -> csv::ReaderBuilder {
        let mut builder = csv::ReaderBuilder::new();

        builder
            .delimiter(self.delimiter() as u8)
            .has_headers(false)
            .comment(self.comment_char().map(|c| c as u8));

        builder
    }

    pub fn to_schema(&self) -> anyhow::Result<Schema> {
        let names = self
            .present_columns()
            .iter()
            .map(|col| col.starts_with('_').not().then_some(col.as_str()))
            .collect_vec();

        polars::schema_subset(&BlastOutFullRecord::polars_schema(), &names)
            .map_err(|unknown| anyhow::anyhow!("Unrecognized BLAST+ column name {unknown}"))
    }

    pub fn load_lazyframe_from_sources(&self, sources: ScanSources) -> anyhow::Result<LazyFrame> {
        let delim = self.delimiter();
        let comment_prefix = self.comment_prefix();

        let schema = self.to_schema()?;

        let result = LazyCsvReader::new_with_sources(sources)
            .with_has_header(false)
            .with_separator(delim as u8)
            .with_comment_prefix(comment_prefix.map(Into::into))
            .with_schema(Some(Arc::new(schema)))
            .with_null_values(Some(NullValues::AllColumnsSingle("N/A".into())))
            .finish()?
            .select([col("*").exclude(["_*"])]);

        Ok(result)
    }

    pub fn load_lazyframe_from_path(&self, path: &Path) -> anyhow::Result<LazyFrame> {
        self.load_lazyframe_from_sources(ScanSources::Paths(Arc::new([path.to_path_buf()])))
    }

    pub fn load_lazyframe_from_static(&self, buf: &'static str) -> anyhow::Result<LazyFrame> {
        self.load_lazyframe_from_sources(ScanSources::Buffers(Arc::new([buf.into()])))
    }

    pub fn load_lazyframe_from_string(&self, buf: String) -> anyhow::Result<LazyFrame> {
        self.load_lazyframe_from_sources(ScanSources::Buffers(Arc::new([buf.into()])))
    }
}

impl FromStr for BlastOutFmt {
    type Err = anyhow::Error;

    fn from_str(fmt: &str) -> Result<Self, Self::Err> {
        let split_fmt = fmt.split_whitespace().collect_vec();

        fn parse_blast_outfmt_6_or_7(mut fmt_args: &[&str]) -> anyhow::Result<(char, Vec<String>)> {
            let delim = if let Some(delim_str) = fmt_args.first().and_then(|s| s.strip_prefix("delim=")) {
                if let Ok(delim_chr) = delim_str.chars().exactly_one() {
                    fmt_args = &fmt_args[1..];
                    delim_chr
                } else {
                    anyhow::bail!("Error parsing delimiter specifier")
                }
            } else {
                '\t'
            };

            let present_columns = if fmt_args.is_empty() {
                BlastOutFullRecord::DEFAULT_BLASTN_COLUMNS
                    .iter()
                    .map(|&s| s.to_owned())
                    .collect()
            } else {
                fmt_args.iter().map(|&s| s.to_owned()).collect()
            };

            Ok((delim, present_columns))
        }

        match split_fmt.first().and_then(|s| s.parse::<u32>().ok()) {
            Some(6) => {
                let (delimiter, present_columns) = parse_blast_outfmt_6_or_7(&split_fmt[1..])?;
                Ok(BlastOutFmt::Six {
                    delimiter,
                    present_columns,
                })
            }
            Some(7) => {
                let (delimiter, present_columns) = parse_blast_outfmt_6_or_7(&split_fmt[1..])?;
                Ok(BlastOutFmt::Seven {
                    delimiter,
                    present_columns,
                })
            }
            Some(n) => anyhow::bail!("Unsupported outfmt type: {n} (only 6 and 7 are supported)"),
            None => anyhow::bail!("Unable to parse outfmt type"),
        }
    }
}

#[expect(dead_code)]
pub fn filters_into_string_record_predicate(
    filters: Vec<BlastOutFilter>,
    outfmt: &BlastOutFmt,
) -> Option<impl Fn(&StringRecord) -> Result<bool, <BlastOutField as FromStrField>::Err>> {
    fn and<P, E>(preds: Vec<P>) -> impl Fn(&StringRecord) -> Result<bool, E>
    where
        P: Fn(&StringRecord) -> Result<bool, E>,
    {
        move |sr| {
            for pred in preds.iter() {
                if !pred(sr)? {
                    return Ok(false);
                }
            }
            Ok(true)
        }
    }

    let preds = filters
        .into_iter()
        .map(|f| f.into_string_record_predicate(outfmt.present_columns()))
        .collect::<Option<_>>()?;

    Some(and(preds))
}

//pub fn load_blastout(
//    path: impl AsRef<Path>,
//    blast_outfmt: &BlastOutFmt,
//    wanted_columns: Option<Vec<&str>>,
//) -> anyhow::Result<LazyFrame> {
//    let delim = blast_outfmt.delimiter();
//    let comment_prefix = blast_outfmt.comment_prefix();
//    let present_columns = blast_outfmt.present_columns();

//    let schema = blast_outfmt.to_schema()?;

//    let result = LazyCsvReader::new(path)
//        .with_has_header(false)
//        .with_separator(delim as u8)
//        .with_comment_prefix(comment_prefix.map(Into::into))
//        .with_schema(Some(Arc::new(schema)))
//        .with_null_values(Some(NullValues::AllColumnsSingle("N/A".into())))
//        .finish()?;

//    let wanted_columns: Vec<&str> = wanted_columns.unwrap_or_else(|| {
//        present_columns
//            .iter()
//            .filter(|p| !p.starts_with('_'))
//            .map(|s| s.as_ref())
//            .collect()
//    });

//    let result = result.select(wanted_columns.iter().map(|name| col(*name)).collect_vec());

//    Ok(result)
//}

#[cfg(test)]
pub(crate) mod tests {
    use std::str::FromStr;

    use itertools::Itertools;

    use crate::tool::blast::blastout::BlastOutFmt;

    pub const BLAST_OUT_FMT: &str = "7 _ qaccver _ saccver staxid pident evalue bitscore";

    pub const BLAST_OUT: &str = "\
        # BLASTN 2.13.0+
        # Query: ES_AV_contig-100_0 length_81420 read_count_843037
        # RID: 492V5FC2013
        # Database: nt
        # Fields: query id, query acc.ver, subject id, subject acc.ver, subject tax id, % identity, evalue, bit score
        # 3051 hits found
        c1	c1	gi|1785199060|gb|CP034345.1|	CP034345.1	1816183	79.053	0.0	3487
        c1	c1	gi|1785199060|gb|CP034345.1|	CP034345.1	1816183	82.548	0.0	2294
        c1	c1	gi|1785199060|gb|CP034345.1|	CP034345.1	1816183	84.000	0.0	1794
        c1	c1	gi|1785199060|gb|CP034340.1|	CP034340.1	1816183	82.927	0.0	1753
        c1	c1	gi|1785199060|gb|CP034345.1|	CP034345.1	1816183	82.151	0.0	1487
        c1	c1	gi|1785199060|gb|CP034345.1|	CP034345.1	1816183	79.718	0.0	1040
        c2	c2	gi|1785199060|gb|CP034345.1|	CP034345.1	1816183	83.040	0.0	953
        c2	c2	gi|1785199060|gb|CP034345.1|	CP034345.1	1816183	75.866	0.0	704
        c2	c2	gi|1785199060|gb|CP034345.1|	CP034345.1	1816183	85.584	0.0	682
        c2	c2	gi|1785199060|gb|CP034345.1|	CP034345.1	1816183	75.912	9.80e-143	527
        c2	c2	gi|1785199060|gb|CP034345.1|	CP034345.1	1816183	78.596	6.29e-95	368
        c2	c2	gi|1785199060|gb|CP034345.1|	CP034345.1	1816183	79.604	1.06e-87	344
        c2	c2	gi|1785199060|gb|CP034343.1|	CP034343.1	1816183	73.453	2.43e-44	200
        c2	c2	gi|1785199060|gb|CP034345.1|	CP034345.1	1816183	97.500	6.94e-25	135
    ";

    #[test]
    fn test_blastout_from_str() {
        let fmt = BlastOutFmt::from_str("7 _ qaccver _ saccver staxid pident evalue bitscore").unwrap();

        assert_eq!(
            fmt,
            BlastOutFmt::Seven {
                delimiter: '\t',
                present_columns: vec![
                    "_".to_owned(),
                    "qaccver".to_owned(),
                    "_".to_owned(),
                    "saccver".to_owned(),
                    "staxid".to_owned(),
                    "pident".to_owned(),
                    "evalue".to_owned(),
                    "bitscore".to_owned()
                ]
            }
        );
    }

    #[test]
    fn test_load_lazyframe() {
        let df = BlastOutFmt::from_str(BLAST_OUT_FMT)
            .unwrap()
            .load_lazyframe_from_static(BLAST_OUT)
            .unwrap();

        let schema = df.clone().collect_schema().unwrap();

        let col_names = schema
            .iter_names()
            .map(|s| s.as_str())
            .filter(|s| !s.starts_with('_'))
            .collect_vec();

        assert_eq!(
            col_names,
            vec!["qaccver", "saccver", "staxid", "pident", "evalue", "bitscore"]
        );
    }
}
