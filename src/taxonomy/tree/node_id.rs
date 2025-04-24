use std::{fmt::{Debug, Display}, marker::PhantomData, str::FromStr};

use intmap::IntMap;
use serde::{de::{SeqAccess, Visitor}, Deserialize, Deserializer, Serialize, Serializer};


#[derive(Debug, Copy, Clone, Eq, PartialEq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[repr(transparent)]
#[serde(transparent)]
pub struct NodeId(pub u32);

pub fn nodeid_to_u32_vec(v: Vec<NodeId>) -> Vec<u32> {
    // See https://doc.rust-lang.org/std/mem/fn.transmute.html#alternatives
    unsafe {
        // Ensure the original vector is not dropped.
        let mut v = std::mem::ManuallyDrop::new(v);
        Vec::from_raw_parts(v.as_mut_ptr() as *mut u32, v.len(), v.capacity())
    }
}

pub fn u32_to_nodeid_vec(v: Vec<u32>) -> Vec<NodeId> {
    // See https://doc.rust-lang.org/std/mem/fn.transmute.html#alternatives
    unsafe {
        // Ensure the original vector is not dropped.
        let mut v = std::mem::ManuallyDrop::new(v);
        Vec::from_raw_parts(v.as_mut_ptr() as *mut NodeId, v.len(), v.capacity())
    }
}

impl Display for NodeId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl FromStr for NodeId {
    type Err = <u32 as FromStr>::Err;

    #[inline(always)]
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(NodeId(s.parse()?))
    }
}

impl From<NodeId> for u32 {
    #[inline(always)]
    fn from(node: NodeId) -> Self {
        node.0
    }
}

impl From<NodeId> for u64 {
    #[inline(always)]
    fn from(node: NodeId) -> Self {
        node.0 as u64
    }
}

impl From<u32> for NodeId {
    #[inline(always)]
    fn from(node: u32) -> Self {
        NodeId(node)
    }
}

#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[repr(transparent)]
pub struct NodeIdMap<T> {
    inner: IntMap<u32, T>,
}

pub type Entry<'a, T> = intmap::Entry<'a, u32, T>;

impl<T> NodeIdMap<T> {
    #[inline(always)]
    pub fn new() -> Self {
        NodeIdMap {
            inner: IntMap::new(),
        }
    }

    #[inline(always)]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            inner: IntMap::with_capacity(capacity)
        }
    }

    #[inline(always)]
    pub fn reserve(&mut self, additional: usize) {
        self.inner.reserve(additional);
    }

    #[inline(always)]
    pub fn iter(&self) -> impl Iterator<Item = (NodeId, &T)> + '_ {
        self.inner.iter().map(|(k, v)| (NodeId::from(k), v))
    }

    #[inline(always)]
    pub fn values_mut(&mut self) -> impl Iterator<Item = &mut T> + '_ {
        self.inner.values_mut()
    }

    #[inline(always)]
    pub fn contains_key(&self, key: NodeId) -> bool {
        self.inner.contains_key(key.into())
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    #[inline(always)]
    pub fn get(&self, key: NodeId) -> Option<&T> {
        self.inner.get(key.into())
    }

    #[inline(always)]
    pub fn get_mut(&mut self, key: NodeId) -> Option<&mut T> {
        self.inner.get_mut(key.into())
    }

    #[inline(always)]
    pub fn insert(&mut self, key: NodeId, value: T) -> Option<T> {
        self.inner.insert(key.into(), value)
    }

    #[inline(always)]
    pub fn extend<I>(&mut self, iter: I)
    where
        I: IntoIterator<Item = (NodeId, T)>,
    {
        self.inner.extend(iter.into_iter().map(|(k, v)| (k.into(), v)))
    }

    #[inline(always)]
    pub fn remove(&mut self, key: NodeId) -> Option<T> {
        self.inner.remove(key.into())
    }

    #[inline(always)]
    pub fn retain<F>(&mut self, mut f: F)
    where
        F: FnMut(NodeId, &T) -> bool,
    {
        self.inner.retain(move |k, v| f(NodeId::from(k as u32), v))
    }

    #[inline(always)]
    pub fn entry(&mut self, key: NodeId) -> Entry<T> {
        self.inner.entry(key.into())
    }
}

impl<T: Debug> Debug for NodeIdMap<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Debug::fmt(&self.inner, f)
    }
}

impl<T> Default for NodeIdMap<T> {
    fn default() -> Self {
        Self { inner: Default::default() }
    }
}

impl<T> Extend<(NodeId, T)> for NodeIdMap<T> {
    fn extend_one(&mut self, (k, v): (NodeId, T)) {
        self.insert(k, v);
    }

    fn extend_reserve(&mut self, additional: usize) {
        self.inner.extend_reserve(additional);
    }

    fn extend<I: IntoIterator<Item = (NodeId, T)>>(&mut self, iter: I) {
        self.inner.extend(iter.into_iter().map(|(k, v)| (u32::from(k), v)))
    }
}

impl<T> IntoIterator for NodeIdMap<T> {
    type Item = (NodeId, T);

    type IntoIter = impl Iterator<Item = Self::Item>;

    #[inline(always)]
    fn into_iter(self) -> Self::IntoIter {
        self.inner.into_iter().map(|(k, v)| (NodeId::from(k as u32), v))
    }
}

impl<'a, T> IntoIterator for &'a NodeIdMap<T> {
    type Item = (NodeId, &'a T);

    type IntoIter = impl Iterator<Item = Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.inner.iter().map(|(k, v)| (NodeId::from(k as u32), v))
    }
}

impl<T> FromIterator<(NodeId, T)> for NodeIdMap<T> {
    #[inline(always)]
    fn from_iter<I>(iter: I) -> Self
    where
        I: IntoIterator<Item = (NodeId, T)>,
    {
        Self {
            inner: iter.into_iter().map(|(k, v)| (k.into(), v)).collect(),
        }
    }
}

#[derive(Clone, Default, PartialEq, Eq)]
#[repr(transparent)]
pub struct NodeIdSet {
    inner: IntMap<u32, ()>,
}

impl NodeIdSet {
    #[inline(always)]
    pub fn new() -> Self {
        Self {
            inner: IntMap::new(),
        }
    }

    #[inline(always)]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            inner: IntMap::with_capacity(capacity)
        }
    }

    #[inline(always)]
    pub fn reserve(&mut self, additional: usize) {
        self.inner.reserve(additional);
    }

    #[inline(always)]
    pub fn iter(&self) -> impl Iterator<Item = NodeId> + '_ {
        self.inner.iter().map(|(k, _)| NodeId::from(k))
    }

    #[inline(always)]
    pub fn contains(&self, key: NodeId) -> bool {
        self.inner.contains_key(key.into())
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    #[inline(always)]
    pub fn insert(&mut self, key: NodeId) -> bool {
        self.inner.insert(key.into(), ()).is_some()
    }

    #[inline(always)]
    fn clear(&mut self) {
        self.inner.clear()
    }

    #[inline(always)]
    pub fn remove(&mut self, key: NodeId) -> bool {
        self.inner.remove(key.into()).is_some()
    }

    #[inline(always)]
    pub fn retain<F>(&mut self, mut f: F)
    where
        F: FnMut(NodeId) -> bool,
    {
        self.inner.retain(move |k, _| f(NodeId::from(k as u32)))
    }
}

impl Debug for NodeIdSet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_set().entries(self.iter()).finish()
    }
}

impl IntoIterator for NodeIdSet {
    type Item = NodeId;

    type IntoIter = impl Iterator<Item = NodeId>;

    #[inline(always)]
    fn into_iter(self) -> Self::IntoIter {
        self.inner.into_iter().map(|(k, _)| NodeId::from(k as u32))
    }
}

impl<'a> IntoIterator for &'a NodeIdSet {
    type Item = NodeId;

    type IntoIter = impl Iterator<Item = NodeId>;

    fn into_iter(self) -> Self::IntoIter {
        self.inner.iter().map(|(k, _)| NodeId::from(k))
    }
}

impl Extend<NodeId> for NodeIdSet {
    fn extend_one(&mut self, v: NodeId) {
        self.insert(v);
    }

    fn extend_reserve(&mut self, additional: usize) {
        self.inner.extend_reserve(additional);
    }

    fn extend<I: IntoIterator<Item = NodeId>>(&mut self, iter: I) {
        self.inner.extend(iter.into_iter().map(|v| (u32::from(v), ())))
    }
}

impl FromIterator<NodeId> for NodeIdSet {
    #[inline(always)]
    fn from_iter<I>(iter: I) -> Self
    where
        I: IntoIterator<Item = NodeId>,
    {
        Self {
            inner: iter.into_iter().map(|k| (k.into(), ())).collect(),
        }
    }
}

impl Serialize for NodeIdSet {
    #[inline]
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.collect_seq(self)
    }
}

// Copied from serde::de::size_hint (private)
fn size_hint_cautious<Element>(hint: Option<usize>) -> usize {
    const MAX_PREALLOC_BYTES: usize = 1024 * 1024;

    if std::mem::size_of::<Element>() == 0 {
        0
    } else {
        std::cmp::min(
            hint.unwrap_or(0),
            MAX_PREALLOC_BYTES / std::mem::size_of::<Element>(),
        )
    }
}

impl<'de> Deserialize<'de> for NodeIdSet {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct SeqVisitor {
            marker: PhantomData<NodeIdSet>,
        }

        impl<'de> Visitor<'de> for SeqVisitor {
            type Value = NodeIdSet;

            fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
                formatter.write_str("a sequence")
            }

            #[inline]
            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: SeqAccess<'de>,
            {
                let mut values = NodeIdSet::with_capacity(size_hint_cautious::<NodeId>(seq.size_hint()));

                while let Some(value) = seq.next_element()? {
                    values.insert(value);
                }

                Ok(values)
            }
        }

        let visitor = SeqVisitor { marker: PhantomData };
        deserializer.deserialize_seq(visitor)
    }

    fn deserialize_in_place<D>(deserializer: D, place: &mut Self) -> Result<(), D::Error>
    where
        D: Deserializer<'de>,
    {
        struct SeqInPlaceVisitor<'a>(&'a mut NodeIdSet);

        impl<'a, 'de> Visitor<'de> for SeqInPlaceVisitor<'a> {
            type Value = ();

            fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
                formatter.write_str("a sequence")
            }

            #[inline]
            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: SeqAccess<'de>,
            {
                self.0.clear();
                self.0.reserve(size_hint_cautious::<NodeId>(seq.size_hint()));

                // FIXME: try to overwrite old values here? (Vec, VecDeque, LinkedList)
                while let Some(value) = seq.next_element()? {
                    self.0.insert(value);
                }

                Ok(())
            }
        }

        deserializer.deserialize_seq(SeqInPlaceVisitor(place))
    }
}
