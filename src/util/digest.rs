use std::{
    fmt::Debug,
    fs::File,
    io::{self, BufReader},
    ops::Deref,
    path::Path,
};

use base64ct::{Base64, Encoding};
use sha2::Digest;

pub type Sha256Digest = impl Deref<Target = [u8]> + Clone + Debug;

#[define_opaque(Sha256Digest)]
pub fn sha256_file(path: &Path) -> io::Result<Box<Sha256Digest>> {
    let mut hasher = sha2::Sha256::new();

    let reader = File::open(path)?;
    let mut reader = BufReader::new(reader);

    io::copy(&mut reader, &mut hasher)?;
    let mut result = Box::default();

    hasher.finalize_into(&mut *result);

    Ok(result)
}

pub fn digest_to_base64(hash: &[u8]) -> String {
    Base64::encode_string(hash)
}
