use thiserror::Error;

use crate::taxonomy::unlabeled_intmap::TaxonomyBuildError;
use std::{io, num::ParseIntError};

#[derive(Debug, Error)]
pub enum DmpError {
    #[error("Tree consistency error")]
    TreeBuildError(#[from] TaxonomyBuildError),

    #[error("Error parsing taxid: {}", .0)]
    TaxIdParseError(#[from] ParseIntError),

    #[error("Missing .dmp field: {}", .0)]
    MissingField(&'static str),

    #[error("error reading .dmp file")]
    ReadError(#[from] io::Error),
}

macro_rules! parse_fields {
    ($line:expr => { $(let $fields:pat = next as $field_names:expr;)+ }; $($body:tt)+) => {{
        let mut line_iter = $line.strip_suffix("\t|").unwrap_or($line).split("\t|\t");

        $(
            let $fields = line_iter.next().ok_or(DmpError::MissingField($field_names))?;
        )+

        $($body)+
    }};
}

pub(crate) use parse_fields;
