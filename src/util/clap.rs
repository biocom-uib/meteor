use std::{fmt::{Debug, Display}, io::{self, ErrorKind}, path::{Path, PathBuf}};

use clap::{builder::{TypedValueParser, ValueParserFactory}, Args};

#[derive(Clone)]
pub struct ExistingFilePath(PathBuf);

impl ExistingFilePath {
    pub fn as_path(&self) -> &Path {
        self.as_ref()
    }
}

impl<T: ?Sized> AsRef<T> for ExistingFilePath
where
    PathBuf: AsRef<T>
{
    fn as_ref(&self) -> &T {
        self.0.as_ref()
    }
}

impl Debug for ExistingFilePath {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        <PathBuf as Debug>::fmt(&self.0, f)
    }
}

impl Display for ExistingFilePath {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0.to_string_lossy())
    }
}

impl ValueParserFactory for ExistingFilePath {
    type Parser = impl TypedValueParser<Value = ExistingFilePath>;

    fn value_parser() -> Self::Parser {
        clap::builder::StringValueParser::new().try_map(|path| {
            let path = PathBuf::from(path);

            let meta = std::fs::metadata(&path)?;

            if meta.is_file() {
                Ok(ExistingFilePath(path))
            } else if path.is_dir() {
                let msg = format!("{} is a directory", path.to_string_lossy());
                Err(io::Error::new(ErrorKind::IsADirectory, msg))
            } else {
                let msg = format!("{} does not exist or is not a file", path.to_string_lossy());
                Err(io::Error::new(ErrorKind::NotFound, msg))
            }
        })
    }
}

#[derive(Clone)]
pub struct ExistingDirPath(PathBuf);

impl ExistingDirPath {
    pub fn as_path(&self) -> &Path {
        self.as_ref()
    }
}

impl<T: ?Sized> AsRef<T> for ExistingDirPath
where
    PathBuf: AsRef<T>
{
    fn as_ref(&self) -> &T {
        self.0.as_ref()
    }
}

impl Debug for ExistingDirPath {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        <PathBuf as Debug>::fmt(&self.0, f)
    }
}

impl Display for ExistingDirPath {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0.to_string_lossy())
    }
}

impl ValueParserFactory for ExistingDirPath {
    type Parser = impl TypedValueParser<Value = ExistingDirPath>;

    fn value_parser() -> Self::Parser {
        clap::builder::StringValueParser::new().try_map(|path| {
            let path = PathBuf::from(path);

            let meta = std::fs::metadata(&path)?;

            if meta.is_dir() {
                Ok(ExistingDirPath(path))
            } else {
                let msg = format!("{} is not a directory", path.to_string_lossy());
                Err(io::Error::new(ErrorKind::NotADirectory, msg))
            }
        })
    }
}

#[derive(Args)]
pub struct OutputFileArg {
    /// Output file. Use '-' to write to the standard output instead.
    output: String,
}

impl<T: ?Sized> AsRef<T> for OutputFileArg
where
    String: AsRef<T>,
{
    fn as_ref(&self) -> &T {
        self.output.as_ref()
    }
}

#[derive(Args)]
pub struct OutputFileFlag {
    /// Output file. Use '-' to write to the standard output instead.
    #[clap(short, long, default_value = "-")]
    output: String,
}

impl<T: ?Sized> AsRef<T> for OutputFileFlag
where
    String: AsRef<T>,
{
    fn as_ref(&self) -> &T {
        self.output.as_ref()
    }
}

#[macro_export]
macro_rules! writing_new_file_or_stdout {
    ($path:expr, $writer:pat => $body:expr $(,)?) => {{
        let path: &str = $path;

        if path == "-" {
            let $writer = Ok::<_, std::io::Error>(std::io::stdout());
            $body
        } else {
            let $writer = std::fs::File::create(path);
            $body
        }
    }};
}

pub use writing_new_file_or_stdout;
