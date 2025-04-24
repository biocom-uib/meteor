use clap::builder::{PossibleValue, TypedValueParser, ValueParserFactory};
use csv::StringRecord;
use polars::lazy::dsl::Expr;
use polars::prelude::{LazyFrame, Literal};
use regex::Regex;
use std::convert::Infallible;
use std::error::Error;
use std::ffi::OsStr;
use std::fmt::Debug;
use std::marker::PhantomData;
use std::num::{ParseFloatError, ParseIntError};
use std::str::FromStr;
use std::sync::LazyLock;
use thiserror::Error;

pub use meteor_macros::FilterableRecord;

#[derive(PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Debug)]
pub enum Op {
    Eq,
    Neq,
    Lt,
    Leq,
    Gt,
    Geq,
}

impl Op {
    pub fn apply<T: ?Sized + PartialOrd>(self, a: &T, b: &T) -> bool {
        match self {
            Op::Eq => a == b,
            Op::Neq => a != b,
            Op::Lt => a < b,
            Op::Leq => a <= b,
            Op::Gt => a > b,
            Op::Geq => a >= b,
        }
    }

    pub fn apply_polars(self, a: Expr, b: Expr) -> Expr {
        match self {
            Op::Eq => a.eq(b),
            Op::Neq => a.neq(b),
            Op::Lt => a.lt(b),
            Op::Leq => a.lt_eq(b),
            Op::Gt => a.gt(b),
            Op::Geq => a.gt_eq(b),
        }
    }
}

impl FromStr for Op {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use Op::*;

        match s {
            "<" => Ok(Lt),
            "<=" => Ok(Leq),
            ">" => Ok(Gt),
            ">=" => Ok(Geq),
            "=" | "==" => Ok(Eq),
            "!=" => Ok(Neq),
            _ => Err(()),
        }
    }
}

#[derive(Debug, Error)]
pub enum FieldParseError {
    #[error("Error parsing floating point number")]
    FloatParseError(#[from] ParseFloatError),

    #[error("Error parsing point number")]
    IntParseError(#[from] ParseIntError),

    #[error("Unknown field {:?}. Valid names are: {:?}", .0, .1)]
    UnknownField(String, &'static [&'static str]),
}

impl From<Infallible> for FieldParseError {
    fn from(value: Infallible) -> Self {
        match value {}
    }
}

pub trait FromStrField: Sized {
    fn field_names() -> &'static [&'static str];

    fn field_name(&self) -> &str;

    fn parse_apply_op(other: &str, op: Op, this: &Self) -> Result<bool, Self::Err>;

    type Err: Error + 'static = FieldParseError;
    fn try_from_parts(key: &str, value: &str) -> Result<Self, Self::Err>;
}


#[derive(Debug, Error)]
pub enum FilterParseError<FieldErr: Error> {
    #[error("Could not parse filter {:?}. The format should be 'field <comparison> value', e.g. 'evalue>=1e-3.'", .0)]
    Split(String),

    #[error("Could not parse filter operator {:?} in {:?}. Supported operators are: <, <=, >, >=, ==, !=.", .1, .0)]
    Op(String, String),

    #[error("Error building filter from parts of {:?}: {}", .0, .1)]
    Field(String, #[source] FieldErr),
}

#[derive(Debug, Clone)]
pub struct Filter<Field> {
    pub op: Op,
    pub field_value: Field,
}

impl<Field> Filter<Field> {
    pub fn new(op: Op, field_value: Field) -> Self {
        Filter { op, field_value }
    }
}

impl<Field: FromStrField> Filter<Field> {
    pub fn parse_many<I, S>(
        filter_strs: I,
    ) -> Result<Vec<Filter<Field>>, FilterParseError<Field::Err>>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        filter_strs
            .into_iter()
            .map(|f| f.as_ref().parse())
            .collect()
    }

    pub fn into_string_record_predicate(
        self,
        columns: &[impl AsRef<str>],
    ) -> Option<impl for<'a> Fn(&'a StringRecord) -> Result<bool, Field::Err>> {
        let i = columns
            .iter()
            .position(|col| col.as_ref() == self.field_value.field_name())?;

        Some(move |sr: &StringRecord| {
            FromStrField::parse_apply_op(&sr[i], self.op, &self.field_value)
        })
    }
}

impl<Field: FromStrField + Literal + Clone> Filter<Field> {
    pub fn to_polars_expr(&self) -> Expr {
        use polars::lazy::dsl::col;

        self.op.apply_polars(
            col(self.field_value.field_name()),
            self.field_value.clone().lit(),
        )
    }

}

impl<Field: FromStrField + Literal + Clone + 'static> Filter<Field> {
    pub fn apply_polars<'a, I>(df: LazyFrame, filters: I) -> LazyFrame
    where
        I: IntoIterator<Item = &'a Self>,
    {
        let combined = filters
            .into_iter()
            .map(|filter| filter.to_polars_expr())
            .reduce(Expr::and);

        if let Some(filter) = combined {
            df.filter(filter)
        } else {
            df
        }
    }
}

impl<Field: FromStrField> FromStr for Filter<Field> {
    type Err = FilterParseError<Field::Err>;

    fn from_str(filter_str: &str) -> Result<Self, FilterParseError<Field::Err>> {
        static RE: LazyLock<Regex> =
            LazyLock::new(|| Regex::new(r"^(\w+)\s*([!><=]{1,2})\s*(\S+)$").unwrap());

        let cap = RE
            .captures(filter_str)
            .ok_or_else(|| FilterParseError::Split(filter_str.to_owned()))?;

        let key = &cap[1];

        let op = cap[2]
            .parse()
            .map_err(|()| FilterParseError::Op(filter_str.to_owned(), cap[2].to_owned()))?;

        let value_str = &cap[3];

        let field_value = Field::try_from_parts(key, value_str)
            .map_err(|e| FilterParseError::Field(filter_str.to_owned(), e))?;

        Ok(Filter { field_value, op })
    }
}

#[derive(Debug, Copy, Clone)]
pub struct FilterValueParser<Field>(PhantomData<fn() -> Field>);

impl<Field> TypedValueParser for FilterValueParser<Field>
where
    Field: FromStrField<Err: Send + Sync> + Clone + Send + Sync + 'static,
{
    type Value = Filter<Field>;

    fn parse_ref(
        &self,
        cmd: &clap::Command,
        arg: Option<&clap::Arg>,
        value: &OsStr,
    ) -> Result<Self::Value, clap::Error> {
        clap::builder::StringValueParser::new()
            .try_map(|filter_str| filter_str.parse())
            .parse_ref(cmd, arg, value)
    }

    fn possible_values(&self) -> Option<Box<dyn Iterator<Item = PossibleValue> + '_>> {
        Some(Box::new(
            Field::field_names()
                .iter()
                .map(|&field| PossibleValue::new(field)),
        ))
    }
}

impl<Field> ValueParserFactory for Filter<Field> {
    type Parser = FilterValueParser<Field>;

    fn value_parser() -> Self::Parser {
        FilterValueParser(PhantomData)
    }
}

pub trait FilterableRecord {
    type Field: 'static;

    fn apply_filter(&self, filter: &Filter<Self::Field>) -> bool;
}

impl<Field> Filter<Field> {
    pub fn apply_all<R>(r: &R, filters: &[Self]) -> bool
    where
        R: FilterableRecord<Field = Field>
    {
        filters.iter().all(|f| r.apply_filter(f))
    }
}
