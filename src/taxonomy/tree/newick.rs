use std::{error::Error, io::{self, Read, Write}, string::FromUtf8Error};

use newick_rs::{newick, SimpleTree};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum NewickLoadError {
    #[error("IO error loading newick")]
    IoError(#[from] io::Error),

    #[error("Encoding error of newick data")]
    EncodingError(#[from] FromUtf8Error),

    #[error("Newick parse error: {}", .0)]
    ParseError(#[source] Box<dyn Error + Send + Sync>),
}


#[cfg(test)]
macro_rules! newick_literal {
    ([ $($trees:tt),+ ] $name:expr) => {
        newick_rs::SimpleTree {
            name: $name.to_string(),
            children: vec![ $( newick_literal!($trees) ),+ ],
            length: None
        }
    };

    (($sub:expr)) => {
        $sub
    };

    ($name:expr) => {
        newick_rs::SimpleTree {
            name: $name.to_string(),
            children: vec![],
            length: None
        }
    };
}

#[cfg(test)]
pub(crate) use newick_literal;

pub fn read_newick_simple_tree<R: Read>(mut reader: R) -> Result<SimpleTree, NewickLoadError> {
    let contents = io::read_to_string(&mut reader)?;

    let tree: SimpleTree = newick::from_newick(&contents)
        .map_err(|err| NewickLoadError::ParseError(Box::new(err.to_owned())))?;

    Ok(tree)
}

pub fn write_newick_simple_tree<W: Write>(
    tree: &SimpleTree,
    mut writer: W,
) -> Result<(), io::Error> {
    writer.write_all(newick::to_newick(tree).as_bytes())
}
