use std::{cell::{Ref, RefCell, RefMut}, hash::Hash, marker::PhantomData, ops::{Deref, DerefMut}, sync::{Arc, RwLock}};

use lending_iterator::higher_kinded_types::{Feed, HKTRef, WithLifetime, HKT};
use string_interner::{backend::Backend, DefaultStringInterner, StringInterner};

pub trait Resolver {
    type Symbol: Copy + Eq + Ord + Hash;
    type ValueHKT: for<'a> WithLifetime<'a, T: Copy>;

    fn resolve(&self, sym: Self::Symbol) -> Option<InternedValue<'_, Self>>;
}

pub type InternedValue<'a, I> = Feed<'a, <I as Resolver>::ValueHKT>;

pub trait Interned : Resolver {
    fn get<'a>(&self, value: InternedValue<'a, Self>) -> Option<Self::Symbol>
    where
        Self: 'a;
}

pub trait Interner: Interned {
    fn get_or_intern<'a>(&mut self, value: InternedValue<'a, Self>) -> Self::Symbol
    where
        Self: 'a;
}

impl<B: Backend> Resolver for StringInterner<B>
where
    StringInterner<B>: Default,
    B::Symbol: Ord + Hash,
{
    type Symbol = B::Symbol;
    type ValueHKT = HKTRef<str>;

    fn resolve(&self, symbol: Self::Symbol) -> Option<&str> {
        StringInterner::resolve(self, symbol)
    }
}

impl<B: Backend> Interned for StringInterner<B>
where
    StringInterner<B>: Default,
    B::Symbol: Ord + Hash,
{
    fn get<'a>(&self, value: InternedValue<'a, Self>) -> Option<Self::Symbol>
    where
        Self: 'a,
    {
        StringInterner::get(self, value)
    }
}

impl<B: Backend> Interner for StringInterner<B>
where
    StringInterner<B>: Default,
    B::Symbol: Ord + Hash,
{
    fn get_or_intern<'a>(&mut self, value: InternedValue<'a, Self>) -> Self::Symbol
    where
        Self: 'a
    {
        StringInterner::get_or_intern(self, value)
    }
}

#[derive(Copy)]
pub struct PlaceholderInterner<S, V>(PhantomData<Option<(S, V)>>);

pub type PlaceholderInternerFor<I> =
    PlaceholderInterner<<I as Resolver>::Symbol, <I as Resolver>::ValueHKT>;

impl<S, V> PlaceholderInterner<S, V> {
    pub fn new() -> Self {
        Self(PhantomData)
    }

    pub fn standing_for<I>(_other: I) -> Self
    where
        I: Resolver<Symbol = S, ValueHKT = V>,
    {
        Self(PhantomData)
    }
}

impl<S, V> Clone for PlaceholderInterner<S, V> {
    fn clone(&self) -> Self {
        PlaceholderInterner(PhantomData)
    }
}

impl<S: Copy + Eq + Ord + Hash, V> Resolver for PlaceholderInterner<S, V> {
    type Symbol = S;
    type ValueHKT = HKTRef<V>;

    fn resolve(&self, _sym: Self::Symbol) -> Option<InternedValue<'_, Self>> {
        None
    }
}

impl<S: Copy + Eq + Ord + Hash, V> Interned for PlaceholderInterner<S, V> {
    fn get<'a>(&self, _value: InternedValue<'a, Self>) -> Option<Self::Symbol>
    where
        Self: 'a,
    {
        None
    }
}

pub struct TrivialInterner<T: Copy + Eq + Ord + Hash>(PhantomData<T>);

impl<T: Copy + Eq + Ord + Hash> Default for TrivialInterner<T> {
    fn default() -> Self {
        Self(PhantomData)
    }
}

impl<T: Copy + Eq + Ord + Hash> Resolver for TrivialInterner<T> {
    type Symbol = T;
    type ValueHKT = HKT!(T);

    fn resolve(&self, sym: Self::Symbol) -> Option<InternedValue<'_, Self>> {
        Some(sym)
    }
}

impl<T: Copy + Eq + Ord + Hash + 'static> Interned for TrivialInterner<T> {
    fn get<'a>(&self, value: InternedValue<'a, Self>) -> Option<Self::Symbol>
    where
        Self: 'a,
    {
        Some(value)
    }
}

impl<T: Copy + Eq + Ord + Hash + 'static> Interner for TrivialInterner<T> {
    fn get_or_intern<'a>(&mut self, value: InternedValue<'a, Self>) -> Self::Symbol
    where
        Self: 'a
    {
        value
    }
}

#[derive(Copy, Clone)]
pub struct BorrowedInterner<IR>(IR);

impl<IR> BorrowedInterner<IR> {
    pub fn new(interner: IR) -> Self {
        BorrowedInterner(interner)
    }

    pub fn into_inner(self) -> IR {
        self.0
    }
}

impl<IR: Deref<Target: Resolver>> Resolver for BorrowedInterner<IR> {
    type Symbol = <IR::Target as Resolver>::Symbol;
    type ValueHKT = <IR::Target as Resolver>::ValueHKT;

    fn resolve(&self, sym: Self::Symbol) -> Option<InternedValue<'_, Self>> {
        self.0.resolve(sym)
    }
}

impl<IR: Deref<Target: Interned>> Interned for BorrowedInterner<IR> {
    fn get<'a>(&self, value: InternedValue<'a, Self>) -> Option<Self::Symbol>
    where
        Self: 'a,
    {
        self.0.get(value)
    }
}

impl<IR: DerefMut<Target: Interner>> Interner for BorrowedInterner<IR> {
    fn get_or_intern<'a>(&mut self, value: InternedValue<'a, Self>) -> Self::Symbol
    where
        Self: 'a,
    {
        self.0.get_or_intern(value)
    }
}

#[derive(Clone, Default)]
pub struct SharedInterner<I>(Arc<RefCell<I>>);

impl<I> SharedInterner<I> {
    pub fn new(interner: I) -> Self {
        SharedInterner(Arc::new(RefCell::new(interner)))
    }

    pub fn borrow(&self) -> BorrowedInterner<Ref<I>> {
        BorrowedInterner::new(self.0.borrow())
    }

    pub fn borrow_mut(&self) -> BorrowedInterner<RefMut<I>> {
        BorrowedInterner::new(self.0.borrow_mut())
    }
}

impl<I: Resolver> Resolver for SharedInterner<I> {
    type Symbol = I::Symbol;
    type ValueHKT = I::ValueHKT;

    fn resolve<'a>(&'a self, sym: Self::Symbol) -> Option<InternedValue<'a, Self>> {
        // SAFETY: Interners never invalidate references
        unsafe {
            let interner_ref = self.0.borrow();
            let result: Option<InternedValue<'_, Self>> = interner_ref.resolve(sym);
            std::mem::transmute::<Option<InternedValue<'_, Self>>, Option<InternedValue<'a, Self>>>(result)
        }
    }
}

impl<I: Interned> Interned for SharedInterner<I> {
    fn get<'a>(&self, value: InternedValue<'a, Self>) -> Option<Self::Symbol>
    where
        Self: 'a,
    {
        self.0.borrow().get(value)
    }
}

impl<I: Interner> Interner for SharedInterner<I> {
    fn get_or_intern<'a>(&mut self, value: InternedValue<'a, Self>) -> Self::Symbol
    where
        Self: 'a
    {
        self.0.borrow_mut().get_or_intern(value)
    }
}

#[derive(Clone, Default)]
pub struct SharedSyncInterner<I>(Arc<RwLock<I>>);

impl<I> SharedSyncInterner<I> {
    pub fn new(interner: I) -> Self {
        SharedSyncInterner(Arc::new(RwLock::new(interner)))
    }
}

impl<I: Resolver> Resolver for SharedSyncInterner<I> {
    type Symbol = I::Symbol;
    type ValueHKT = I::ValueHKT;

    fn resolve<'a>(&'a self, sym: Self::Symbol) -> Option<InternedValue<'a, Self>> {
        // SAFETY: The interner never invalidates references
        unsafe {
            let interner_ref = self.0.read().unwrap();
            let result: Option<InternedValue<'_, Self>> = interner_ref.resolve(sym);
            std::mem::transmute::<Option<InternedValue<'_, Self>>, Option<InternedValue<'a, Self>>>(result)
        }
    }
}

impl<I: Interned> Interned for SharedSyncInterner<I> {
    fn get<'a>(&self, value: InternedValue<'a, Self>) -> Option<Self::Symbol>
    where
        Self: 'a,
    {
        self.0.read().unwrap().get(value)
    }
}

impl<I: Interner> Interner for SharedSyncInterner<I> {
    fn get_or_intern<'a>(&mut self, value: InternedValue<'a, Self>) -> Self::Symbol
    where
        Self: 'a
    {
        self.0.write().unwrap().get_or_intern(value)
    }
}

#[derive(Default)]
pub struct PairInterner<I1 = DefaultStringInterner, I2 = I1>(pub I1, pub I2);

impl<I1: Resolver, I2: Resolver> Resolver for PairInterner<I1, I2> {
    type Symbol = (I1::Symbol, I2::Symbol);
    type ValueHKT = HKT!((Feed<'_, I1::ValueHKT>, Feed<'_, I2::ValueHKT>));

    fn resolve(&self, sym: Self::Symbol) -> Option<InternedValue<'_, Self>> {
        Some((self.0.resolve(sym.0)?, self.1.resolve(sym.1)?))
    }
}

impl<I1: Interned, I2: Interned> Interned for PairInterner<I1, I2> {
    fn get<'a>(&self, value: InternedValue<'a, Self>) -> Option<Self::Symbol>
    where
        Self: 'a,
    {
        Some((self.0.get(value.0)?, self.1.get(value.1)?))
    }
}

impl<I1: Interner, I2: Interner> Interner for PairInterner<I1, I2> {
    fn get_or_intern<'a>(&mut self, value: InternedValue<'a, Self>) -> Self::Symbol
    where
        Self: 'a
    {
        (self.0.get_or_intern(value.0), self.1.get_or_intern(value.1))
    }
}
