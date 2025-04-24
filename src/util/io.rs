use std::{io::Read, path::Path};

use flate2::read::GzDecoder;

#[cfg(target_family = "unix")]
fn is_broken_pipe_impl(mut err: &(dyn std::error::Error + 'static)) -> bool {
    use polars::error::PolarsError;

    loop {
        if let Some(io_err) = err.downcast_ref::<std::io::Error>() {
            if io_err.kind() == std::io::ErrorKind::BrokenPipe {
                return true;
            }

        } else if let Some(cause) = err.source() {
            err = cause;

        // currently (0.43) polars does not report IO { error } as .cause()
        } else if let Some(PolarsError::IO { error: io_err, msg: _ }) = err.downcast_ref::<PolarsError>() {
            err = io_err.as_ref();

        } else {
            break;
        }
    }

    false
}

#[cfg(not(target_family = "unix"))]
fn is_broken_pipe_impl(result: &(dyn std::error::Error + 'static)) -> bool {
    false
}

pub fn is_broken_pipe<E>(err: &E) -> bool
where
    E: std::error::Error + 'static
{
    is_broken_pipe_impl(err)
}

pub fn ignore_broken_pipe<E>(result: Result<(), E>) -> Result<(), E>
where
    E: std::error::Error + 'static
{
    if let Err(err) = &result {
        if is_broken_pipe(err) {
            return Ok(())
        }
    }

    result
}

pub fn is_broken_pipe_anyhow<E>(err: &E) -> bool
where
    E: AsRef<dyn std::error::Error + 'static>,
{
    is_broken_pipe_impl(err.as_ref())
}

pub fn ignore_broken_pipe_anyhow<E>(result: Result<(), E>) -> Result<(), E>
where
    E: AsRef<dyn std::error::Error + 'static>,
{
    if let Err(err) = &result {
        if is_broken_pipe_anyhow(err) {
            return Ok(())
        }
    }

    result
}

#[allow(clippy::large_enum_variant, dead_code)]
pub enum MaybeGzDecoder<R> {
    GzDecoder(GzDecoder<R>),
    Reader(R),
}

impl MaybeGzDecoder<std::fs::File> {
    pub fn from_path(path: &Path) -> std::io::Result<Self> {
        let file = std::fs::File::open(path)?;

        if let Some(ext) = path.extension().and_then(|ext| ext.to_str()) {
            if ext.to_lowercase() == "gz" {
                return Ok(Self::GzDecoder(GzDecoder::new(file)))
            }
        }

        Ok(Self::Reader(file))
    }
}

impl<R: Read> Read for MaybeGzDecoder<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            MaybeGzDecoder::GzDecoder(r) => r.read(buf),
            MaybeGzDecoder::Reader(r) => r.read(buf),
        }
    }
}

#[allow(unused_macros)]
macro_rules! maybe_gzdecoder {
    ($reader_expr:expr, $reader_pat:pat => $body:expr $(,)?) => {{
        match $reader_expr {
            crate::util::io::MaybeGzDecoder::GzDecoder($reader_pat) => $body,
            crate::util::io::MaybeGzDecoder::Reader($reader_pat) => $body,
        }
    }}
}

#[allow(unused_imports)]
pub(crate) use maybe_gzdecoder;
