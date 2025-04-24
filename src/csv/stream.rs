use std::{io, marker::PhantomData};

use csv::StringRecord;
use lending_iterator::prelude::*;
use serde::Deserialize;

pub fn tsv_reader_builder(has_header: bool) -> csv::ReaderBuilder {
    let mut builder = csv::ReaderBuilder::new();

    builder
        .has_headers(has_header)
        .delimiter(b'\t');

    builder
}

pub struct CsvReaderIter<R> {
    buf_record: StringRecord,
    reader: csv::Reader<R>,
}

pub struct DeserializedCsvReaderIter<T: HKT, R> {
    inner: CsvReaderIter<R>,
    headers: Option<StringRecord>,
    #[expect(clippy::type_complexity)]
    phantom: PhantomData<for<'a> fn(&'a StringRecord) -> Option<Feed<'a, T>>>,
}

impl<R: io::Read> CsvReaderIter<R> {
    pub fn new(reader: csv::Reader<R>) -> Self {
        Self {
            buf_record: StringRecord::new(),
            reader,
        }
    }

    pub fn into_deserialize<T: HKT>(
        self,
        headers: Option<StringRecord>,
    ) -> DeserializedCsvReaderIter<T, R> {
        DeserializedCsvReaderIter {
            inner: self,
            headers,
            phantom: PhantomData,
        }
    }

    pub fn into_deserialize_ref<T>(
        self,
        headers: Option<StringRecord>,
    ) -> DeserializedCsvReaderIter<HKTRef<T>, R> {
        DeserializedCsvReaderIter {
            inner: self,
            headers,
            phantom: PhantomData,
        }
    }
}

pub trait CsvReaderIterExt {
    type Iter: LendingIterator;

    fn into_lending_iter(self) -> Self::Iter;
}

impl<R: io::Read> CsvReaderIterExt for csv::Reader<R> {
    type Iter = CsvReaderIter<R>;

    fn into_lending_iter(self) -> Self::Iter {
        CsvReaderIter::new(self)
    }
}

#[gat]
impl<R: io::Read> LendingIterator for CsvReaderIter<R> {
    type Item<'next> = csv::Result<&'next StringRecord>;

    fn next(&mut self) -> Option<Self::Item<'_>> {
        match self.reader.read_record(&mut self.buf_record) {
            Err(e) => Some(Err(e)),
            Ok(true) => Some(Ok(&self.buf_record)),
            Ok(false) => None,
        }
    }
}

#[gat]
impl<T, R> LendingIterator for DeserializedCsvReaderIter<T, R>
where
    T: HKT,
    for<'de> Feed<'de, T>: Deserialize<'de>,
    R: io::Read,
{
    type Item<'next> = csv::Result<Feed<'next, T>>;

    fn next(&mut self) -> Option<Self::Item<'_>> {
        let headers = self.headers.as_ref();

        let next = self.inner.next()?;

        Some(next.and_then(|record| record.deserialize(headers)))
    }
}
