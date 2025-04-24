use core::slice;
use std::io::BufRead;
use std::{collections::HashMap, fs::File, io::BufReader, path::Path};

use serde::{Deserialize, Serialize};
use string_interner::{symbol::SymbolU32, DefaultStringInterner, StringInterner, Symbol};

use crate::taxonomy::tree::node_id::{Entry, NodeIdMap};
use crate::taxonomy::NodeId;

use super::dmp::{parse_fields, DmpError};
use super::TaxId;

#[derive(Copy, Clone, Debug, Hash, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[repr(transparent)]
#[serde(transparent)]
pub struct NameClassSymbol(u32);

impl From<SymbolU32> for NameClassSymbol {
    fn from(sym: SymbolU32) -> Self {
        NameClassSymbol(sym.to_usize() as u32)
    }
}

impl From<NameClassSymbol> for SymbolU32 {
    fn from(sym: NameClassSymbol) -> Self {
        <SymbolU32 as Symbol>::try_from_usize(sym.0 as usize).unwrap()
    }
}

#[derive(Default, Serialize, Deserialize)]
pub struct AllNames {
    pub name_classes: DefaultStringInterner,
    pub taxid_to_names: NodeIdMap<Vec<(NameClassSymbol, String)>>,
    pub name_to_taxids: HashMap<String, Vec<(NameClassSymbol, TaxId)>>,
}

impl AllNames {
    pub const FORMAT_VERSION: u32 = 1;

    pub fn load_filtered_names_dmp<S, P>(classes: &[S], names_dmp: P) -> Result<Self, DmpError>
    where
        S: AsRef<str>,
        P: AsRef<Path>,
    {
        let mut this = AllNames {
            name_classes: StringInterner::new(),
            taxid_to_names: NodeIdMap::new(),
            name_to_taxids: HashMap::new(),
        };

        this.fill_from(names_dmp, |class| {
            classes.iter().any(|c| c.as_ref() == class)
        })?;

        let order: HashMap<NameClassSymbol, usize> = classes
            .iter()
            .enumerate()
            .map(|(i, s)| (this.name_classes.get_or_intern(s.as_ref()).into(), i))
            .collect();

        for names in this.taxid_to_names.values_mut() {
            names.sort_by_key(|(class, _)| order.get(class).copied().unwrap_or(classes.len()));
        }

        for taxids in this.name_to_taxids.values_mut() {
            taxids.sort_by_key(|(class, _)| order.get(class).copied().unwrap_or(classes.len()));
        }

        Ok(this)
    }

    #[allow(clippy::result_large_err)]
    pub fn only_of_class(self, name_class: &str) -> Result<SingleClassNames, Self> {
        let single_class = if let Some(sym) = self.name_classes.get(name_class) {
            sym.into()
        } else {
            return Err(self);
        };

        let taxid_to_name = self
            .taxid_to_names
            .into_iter()
            .filter_map(|(taxid, names)| {
                names
                    .into_iter()
                    .find(|(class, _)| *class == single_class)
                    .map(|(_, single_name)| (taxid, single_name))
            })
            .collect();

        let name_to_taxids = self
            .name_to_taxids
            .into_iter()
            .filter_map(|(name, taxids)| {
                let new_taxids: Vec<TaxId> = taxids
                    .into_iter()
                    .filter(|(class, _)| *class == single_class)
                    .map(|(_, taxid)| taxid)
                    .collect();

                if new_taxids.is_empty() {
                    None
                } else {
                    Some((name, new_taxids))
                }
            })
            .collect();

        Ok(SingleClassNames {
            name_class: name_class.to_owned(),
            taxid_to_name,
            name_to_taxids,
        })
    }
}

#[derive(Debug, Deserialize, Serialize)]
pub struct SingleClassNames {
    pub name_class: String,
    pub taxid_to_name: NodeIdMap<String>,
    pub name_to_taxids: HashMap<String, Vec<TaxId>>,
}

impl SingleClassNames {
    pub const FORMAT_VERSION: u32 = 1;

    pub fn load_names_dmp<P: AsRef<Path>>(
        name_class: String,
        names_dmp: P,
    ) -> Result<Self, DmpError> {
        let mut this = SingleClassNames {
            name_class,
            taxid_to_name: NodeIdMap::new(),
            name_to_taxids: HashMap::new(),
        };

        this.fill_from(names_dmp, |_| true)?;
        Ok(this)
    }
}

#[derive(Copy, Clone, Default, Deserialize, Serialize)]
pub struct NoNames {}

impl NoNames {
    pub const FORMAT_VERSION: u32 = 1;

    pub fn new() -> Self {
        NoNames {}
    }
}

pub trait NamesAssoc {
    fn insert(&mut self, taxid: TaxId, name: &str, unique_name: &str, name_class: &str);

    fn forget_taxids<Keep>(&mut self, keep: Keep)
    where
        Keep: FnMut(TaxId) -> bool;

    type NameLookup;
    fn lookup_names(&self, taxid: TaxId) -> Option<&Self::NameLookup>;

    type NamesLookupIter<'a>: Iterator<Item = &'a str> + 'a;
    fn iter_lookup_names(lookup: &Self::NameLookup) -> Self::NamesLookupIter<'_>;

    type TaxIdsLookup;
    fn lookup_taxids(&self, name: &str) -> Option<&Self::TaxIdsLookup>;

    type TaxIdsLookupIter<'a>: Iterator<Item = TaxId> + 'a;
    fn iter_lookup_taxids(lookup: &Self::TaxIdsLookup) -> Self::TaxIdsLookupIter<'_>;

    fn fill_from<P, F>(&mut self, names_dmp: P, class_filter: F) -> Result<(), DmpError>
    where
        P: AsRef<Path>,
        F: Fn(&str) -> bool,
    {
        for line in BufReader::new(File::open(names_dmp)?).lines() {
            let line = line?;

            parse_fields! {
                &line => {
                    let taxid = next as "tax_id";
                    let name = next as "name_txt";
                    let unique_name = next as "unique name";
                    let name_class = next as "name class";
                };

                let taxid = NodeId(taxid.parse()?);

                if class_filter(name_class) {
                    self.insert(taxid, name, unique_name, name_class);
                }
            }
        }

        Ok(())
    }
}

impl NamesAssoc for AllNames {
    fn insert(&mut self, taxid: TaxId, name: &str, unique_name: &str, name_class: &str) {
        let chosen_name = if !name.is_empty() {
            name
        } else if !unique_name.is_empty() {
            unique_name
        } else {
            return;
        };

        let name_class = self.name_classes.get_or_intern(name_class).into();

        let names = match self.taxid_to_names.entry(taxid) {
            Entry::Occupied(occ) => occ.into_mut(),
            Entry::Vacant(vac) => vac.insert(Vec::new()),
        };

        names.push((name_class, chosen_name.to_owned()));

        self.name_to_taxids
            .entry(chosen_name.to_owned())
            .or_default()
            .push((name_class, taxid));
    }

    fn forget_taxids<Keep: FnMut(TaxId) -> bool>(&mut self, mut keep: Keep) {
        self.taxid_to_names.retain(|taxid, _| keep(taxid));

        self.name_to_taxids.retain(|_, taxids| {
            taxids.retain(|(_, taxid)| keep(*taxid));
            !taxids.is_empty()
        });
    }

    type NameLookup = Vec<(NameClassSymbol, String)>;
    fn lookup_names(&self, taxid: TaxId) -> Option<&Self::NameLookup> {
        self.taxid_to_names.get(taxid)
    }

    type NamesLookupIter<'a> = IterStrSnd<slice::Iter<'a, (NameClassSymbol, String)>>;
    fn iter_lookup_names(lookup: &Self::NameLookup) -> Self::NamesLookupIter<'_> {
        IterStrSnd(lookup.iter())
    }

    type TaxIdsLookup = Vec<(NameClassSymbol, TaxId)>;
    fn lookup_taxids<'a>(&'a self, name: &str) -> Option<&'a Self::TaxIdsLookup> {
        self.name_to_taxids.get(name)
    }

    type TaxIdsLookupIter<'a> = IterCopySnd<slice::Iter<'a, (NameClassSymbol, NodeId)>>;
    fn iter_lookup_taxids(lookup: &Self::TaxIdsLookup) -> Self::TaxIdsLookupIter<'_> {
        IterCopySnd(lookup.iter())
    }
}

pub struct IterStrSnd<I>(I);

impl<'a, T: 'a, I: Iterator<Item = &'a (T, String)>> Iterator for IterStrSnd<I> {
    type Item = &'a str;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next().map(|s| &*s.1)
    }
}

pub struct IterCopySnd<I>(I);

impl<'a, T: 'a, U: 'a + Copy, I: Iterator<Item = &'a (T, U)>> Iterator for IterCopySnd<I> {
    type Item = U;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next().map(|(_, b)| *b)
    }
}

impl NamesAssoc for SingleClassNames {
    fn insert(&mut self, taxid: TaxId, name: &str, unique_name: &str, name_class: &str) {
        if name_class == self.name_class {
            if let Entry::Vacant(vac) = self.taxid_to_name.entry(taxid) {
                let chosen_name = if name.is_empty() { unique_name } else { name };
                vac.insert(chosen_name.to_owned());
            }

            self.name_to_taxids
                .entry(name.to_owned())
                .or_default()
                .push(taxid);
        }
    }

    fn forget_taxids<Keep>(&mut self, mut keep: Keep)
    where
        Keep: FnMut(TaxId) -> bool,
    {
        self.taxid_to_name.retain(|taxid, _| keep(taxid));

        self.name_to_taxids.retain(|_, taxids| {
            taxids.retain(|taxid| keep(*taxid));
            !taxids.is_empty()
        })
    }

    type NameLookup = String;
    fn lookup_names(&self, taxid: TaxId) -> Option<&Self::NameLookup> {
        self.taxid_to_name.get(taxid)
    }

    type NamesLookupIter<'a> = std::iter::Once<&'a str>;
    fn iter_lookup_names(lookup: &Self::NameLookup) -> Self::NamesLookupIter<'_> {
        std::iter::once(lookup)
    }

    type TaxIdsLookup = Vec<TaxId>;
    fn lookup_taxids<'a>(&'a self, name: &str) -> Option<&'a Self::TaxIdsLookup> {
        self.name_to_taxids.get(name)
    }

    type TaxIdsLookupIter<'a> = std::iter::Cloned<std::slice::Iter<'a, TaxId>>;
    fn iter_lookup_taxids(lookup: &Self::TaxIdsLookup) -> Self::TaxIdsLookupIter<'_> {
        lookup.iter().cloned()
    }
}

impl NamesAssoc for NoNames {
    fn insert(&mut self, _taxid: TaxId, _name: &str, _unique_name: &str, _name_class: &str) {}

    fn forget_taxids<Keep>(&mut self, _keep: Keep)
    where
        Keep: FnMut(TaxId) -> bool,
    {
    }

    type NameLookup = !;
    fn lookup_names(&self, _taxid: TaxId) -> Option<&Self::NameLookup> {
        None
    }

    type NamesLookupIter<'a> = std::iter::Empty<&'a str>;
    fn iter_lookup_names(_lookup: &Self::NameLookup) -> Self::NamesLookupIter<'_> {
        std::iter::empty()
    }

    type TaxIdsLookup = !;
    fn lookup_taxids(&self, _name: &str) -> Option<&Self::TaxIdsLookup> {
        None
    }

    type TaxIdsLookupIter<'a> = std::iter::Empty<TaxId>;
    fn iter_lookup_taxids(_lookup: &Self::TaxIdsLookup) -> Self::TaxIdsLookupIter<'_> {
        std::iter::empty()
    }

    fn fill_from<P, F>(&mut self, _names_dmp: P, _class_filter: F) -> Result<(), DmpError>
    where
        F: Fn(&str) -> bool,
        P: AsRef<Path>,
    {
        Ok(())
    }
}
