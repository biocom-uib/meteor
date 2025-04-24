use std::{collections::{HashMap, HashSet}, io};

use csv::StringRecord;
use lending_iterator::prelude::*;

use crate::csv::stream::CsvReaderIter;

use super::interners::{Interned, Interner, Resolver, InternedValue};

pub struct InternedMapping<KeyInterner, ValueInterner>
where
    KeyInterner: Resolver,
    ValueInterner: Resolver,
{
    key_interner: KeyInterner,
    value_interner: ValueInterner,

    mapping: HashMap<KeyInterner::Symbol, HashSet<ValueInterner::Symbol>>,
}

impl<KeyInterner, ValueInterner> Default for InternedMapping<KeyInterner, ValueInterner>
where
    KeyInterner: Default + Resolver,
    ValueInterner: Default + Resolver,
{
    fn default() -> Self {
        Self {
            key_interner: Default::default(),
            value_interner: Default::default(),
            mapping: Default::default(),
        }
    }
}

pub trait IntoValues<'a, I: Resolver + 'a>: IntoIterator<Item = InternedValue<'a, I>> {}

impl<'a, I, Iter> IntoValues<'a, I> for Iter
where
    I: Resolver + 'a,
    Iter: IntoIterator<Item = InternedValue<'a, I>>,
{
}

impl<KeyInterner, ValueInterner> InternedMapping<KeyInterner, ValueInterner>
where
    KeyInterner: Resolver,
    ValueInterner: Resolver,
{
    pub fn new(key_interner: KeyInterner, value_interner: ValueInterner) -> Self {
        Self {
            key_interner,
            value_interner,
            mapping: Default::default(),
        }
    }

    pub fn get(
        &self,
        key_sym: KeyInterner::Symbol,
    ) -> Option<impl Iterator<Item = ValueInterner::Symbol> + '_> {
        self.mapping
            .get(&key_sym)
            .map(|set| set.iter().copied())
    }

    pub fn key_interner(&self) -> &KeyInterner {
        &self.key_interner
    }

    pub fn value_interner(&self) -> &ValueInterner {
        &self.value_interner
    }

    pub fn replace_interners<KeyInterner2, ValueInterner2, F>(
        self,
        f: F,
    ) -> InternedMapping<KeyInterner2, ValueInterner2>
    where
        KeyInterner2: Resolver<Symbol = KeyInterner::Symbol>,
        ValueInterner2: Resolver<Symbol = ValueInterner::Symbol>,
        F: FnOnce(KeyInterner, ValueInterner) -> (KeyInterner2, ValueInterner2),
    {
        let (key_interner, value_interner) = f(self.key_interner, self.value_interner);

        InternedMapping {
            key_interner,
            value_interner,
            mapping: self.mapping,
        }
    }
}

impl<KeyInterner, ValueInterner> InternedMapping<KeyInterner, ValueInterner>
where
    KeyInterner: Interner<ValueHKT = HKTRef<str>>,
    ValueInterner: Interner,
{
    pub fn add<'a, 'b>(&mut self, key: InternedValue<'a, KeyInterner>, value: InternedValue<'b, ValueInterner>)
    where
        KeyInterner: 'a,
        ValueInterner: 'b,
    {
        let key_sym = self.key_interner.get_or_intern(key);
        let value_sym = self.value_interner.get_or_intern(value);

        self.mapping.entry(key_sym).or_default().insert(value_sym);
    }
}

impl<KeyInterner, ValueInterner> InternedMapping<KeyInterner, ValueInterner>
where
    KeyInterner: Interned,
    ValueInterner: Resolver,
{
    pub fn lookup_syms<'b, 'a: 'b>(
        &'a self,
        key: InternedValue<'b, KeyInterner>,
    ) -> impl Iterator<Item = ValueInterner::Symbol> + 'a
    where
        KeyInterner: 'b,
    {
        self
            .key_interner
            .get(key)
            .and_then(|key_sym| self.mapping.get(&key_sym))
            .into_iter()
            .flatten()
            .copied()
    }

    pub fn lookup<'b, 'a: 'b>(
        &'a self,
        key: InternedValue<'b, KeyInterner>,
    ) -> impl Iterator<Item = InternedValue<'a, ValueInterner>> + 'a {
        self.lookup_syms(key)
            .filter_map(|value_sym| self.value_interner.resolve(value_sym))
    }
}

impl<KeyInterner, ValueInterner> InternedMapping<KeyInterner, ValueInterner>
where
    KeyInterner: Interner,
    ValueInterner: Interner,
{
    #[apply(Gat!)]
    pub fn try_add_many<'a, Iter, E>(&'a mut self, mut iter: Iter) -> Result<(), E>
    where
        Iter: for<'n> LendingIterator<
                Item<'n> = Result<(InternedValue<'n, KeyInterner>, InternedValue<'n, ValueInterner>), E>,
            > + 'a,
    {
        let mut entry_cache_key_sym = None;
        let mut entry_cache_values: Option<&mut HashSet<ValueInterner::Symbol>> = None;

        while let Some(result) = iter.next() {
            let (key, value) = result?;
            let key_sym = self.key_interner.get_or_intern(key);
            let value_sym = self.value_interner.get_or_intern(value);

            entry_cache_values = match entry_cache_values {
                Some(values) if entry_cache_key_sym == Some(key_sym) => {
                    values.insert(value_sym);
                    Some(values)
                },
                _ => {
                    let values = self.mapping.entry(key_sym).or_default();
                    values.insert(value_sym);

                    entry_cache_key_sym = Some(key_sym);
                    Some(values)
                }
            };
        }

        Ok(())
    }

    #[apply(Gat!)]
    pub fn try_add_many_grouped<'a, VS: HKT, E, Iter>(&'a mut self, mut iter: Iter) -> Result<(), E>
    where
        Iter: for<'n> LendingIterator<Item<'n> = Result<(InternedValue<'n, KeyInterner>, Feed<'n, VS>), E>>
            + 'a,
        for<'n> Feed<'n, VS>: IntoValues<'n, ValueInterner>,
    {
        let mut entry_cache_key_sym = None;
        let mut entry_cache_values: Option<&mut HashSet<ValueInterner::Symbol>> = None;

        while let Some(result) = iter.next() {
            let (key, values) = result?;

            let key_sym = self.key_interner.get_or_intern(key);

            let value_syms = values
                .into_iter()
                .map(|value| self.value_interner.get_or_intern(value));

            entry_cache_values = match entry_cache_values {
                Some(values) if entry_cache_key_sym == Some(key_sym) => {
                    values.extend(value_syms);
                    Some(values)
                },
                _ => {
                    let values = self.mapping.entry(key_sym).or_default();
                    values.extend(value_syms);

                    entry_cache_key_sym = Some(key_sym);
                    Some(values)
                }
            };
        }

        Ok(())
    }

    pub fn extend_from_tsv_with<'a, R, F>(
        &'a mut self,
        reader: R,
        header: bool,
        record_parser: F,
    ) -> anyhow::Result<()>
    where
        R: io::Read,
        F: for<'n> Fn(
            &'n StringRecord,
        ) -> anyhow::Result<(
            InternedValue<'n, KeyInterner>,
            InternedValue<'n, ValueInterner>,
        )>,
    {
        let csv_reader = csv::ReaderBuilder::new()
            .comment(Some(b'#'))
            .delimiter(b'\t')
            .has_headers(header)
            .from_reader(reader);

        let csv_iter = CsvReaderIter::new(csv_reader).map::<HKT!(
            anyhow::Result<(
                InternedValue<'_, KeyInterner>,
                InternedValue<'_, ValueInterner>
            )>
        ), _>(|[], result| {
                    result.map_err(|e| e.into()).and_then(&record_parser)
                });

        self.try_add_many(csv_iter)?;

        Ok(())
    }

    pub fn extend_from_grouped_tsv_with<'a, VS: HKT, R, F>(
        &'a mut self,
        reader: R,
        header: bool,
        record_parser: F,
    ) -> anyhow::Result<()>
    where
        R: io::Read,
        F: Fn(&'_ StringRecord) -> anyhow::Result<(InternedValue<'_, KeyInterner>, Feed<'_, VS>)>,
        for<'n> Feed<'n, VS>: IntoValues<'n, ValueInterner>,
    {
        let csv_reader = csv::ReaderBuilder::new()
            .comment(Some(b'#'))
            .delimiter(b'\t')
            .has_headers(header)
            .from_reader(reader);

        let csv_iter = CsvReaderIter::new(csv_reader).map::<HKT!(
            anyhow::Result<(InternedValue<'_, KeyInterner>, Feed<'_, VS>)>
        ), _>(|[], result| {
            result.map_err(|e| e.into()).and_then(&record_parser)
        });

        self.try_add_many_grouped(csv_iter)?;

        Ok(())
    }
}

impl<KeyInterner, ValueInterner> InternedMapping<KeyInterner, ValueInterner>
where
    KeyInterner: Default + Interner<ValueHKT = HKTRef<str>>,
    ValueInterner: Default + Interner,
{
    pub fn read_tsv_with<R, F>(reader: R, header: bool, record_parser: F) -> anyhow::Result<Self>
    where
        R: io::Read,
        F: Fn(&StringRecord) -> anyhow::Result<(&str, InternedValue<'_, ValueInterner>)>,
    {
        let mut mapping = Self::default();

        mapping.extend_from_tsv_with(reader, header, record_parser)?;

        Ok(mapping)
    }

    pub fn read_grouped_tsv_with<VS: HKT, R, F>(reader: R, header: bool, record_parser: F) -> anyhow::Result<Self>
    where
        R: io::Read,
        F: Fn(&StringRecord) -> anyhow::Result<(&str, Feed<'_, VS>)>,
        for<'a> Feed<'a, VS>: IntoValues<'a, ValueInterner>,
    {
        let mut mapping = Self::default();

        mapping.extend_from_grouped_tsv_with(reader, header, record_parser)?;

        Ok(mapping)
    }
}

#[apply(Gat!)]
impl<KeyInterner, ValueInterner> InternedMapping<KeyInterner, ValueInterner>
where
    KeyInterner: Resolver,
    ValueInterner: Resolver,
{
    pub fn write_tsv_with<W, F>(
        &self,
        writer: W,
        header: &[&str],
        serializer_fn: F,
    ) -> csv::Result<()>
    where
        W: io::Write,
        F: for<'a> Fn(&mut csv::Writer<W>, InternedValue<'a, KeyInterner>, InternedValue<'a, ValueInterner>) -> csv::Result<()>,
    {
        let mut csv_writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_writer(writer);

        csv_writer.write_record(header)?;

        for (&key_sym, value_syms) in &self.mapping {
            let key = self.key_interner.resolve(key_sym).unwrap();

            for &value_sym in value_syms {
                let value = self.value_interner.resolve(value_sym).unwrap();

                serializer_fn(&mut csv_writer, key, value)?;
            }
        }

        csv_writer.flush()?;
        Ok(())
    }
}

impl<KeyInterner, ValueInterner> InternedMapping<KeyInterner, ValueInterner>
where
    KeyInterner: Interner<ValueHKT = HKTRef<str>>,
    ValueInterner: Interner<ValueHKT = HKTRef<str>>,
{
    pub fn extend_from_tsv<R: io::Read>(&mut self, reader: R, header: bool) -> anyhow::Result<()> {
        self.extend_from_tsv_with(reader, header, |record| Ok((&record[0], &record[1])))
    }

    pub fn write_tsv<W: io::Write>(&self, writer: W, header: &[&str]) -> csv::Result<()> {
        self.write_tsv_with(writer, header, |csv_writer, key, value| {
            csv_writer.write_record([key, value])?;
            Ok(())
        })
    }
}

impl<KeyInterner, ValueInterner> InternedMapping<KeyInterner, ValueInterner>
where
    KeyInterner: Default,
    ValueInterner: Default,
    KeyInterner: for<'a> Interner<ValueHKT = HKTRef<str>>,
    ValueInterner: for<'a> Interner<ValueHKT = HKTRef<str>>,
{
    pub fn read_tsv<R: io::Read>(reader: R, header: bool) -> anyhow::Result<Self> {
        let mut mapping = Self::default();

        mapping.extend_from_tsv(reader, header)?;

        Ok(mapping)
    }
}
