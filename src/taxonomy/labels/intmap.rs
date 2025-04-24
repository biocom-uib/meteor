use serde::{Deserialize, Serialize};
use string_interner::{
    backend::Backend, DefaultBackend, DefaultStringInterner, DefaultSymbol, Symbol,
};

use crate::taxonomy::{tree::node_id::NodeIdMap, NodeId};


pub type RankSymbol = <DefaultBackend as Backend>::Symbol;

#[derive(Default, Deserialize, Serialize)]
pub struct IntMapRanks {
    pub ranks: NodeIdMap<u32>, // SymbolU32
    pub interner: DefaultStringInterner,
}

impl IntMapRanks {
    pub fn rank_sym_str(&self, rank_sym: DefaultSymbol) -> Option<&str> {
        self.interner.resolve(rank_sym)
    }

    pub fn lookup_rank_sym(&self, rank: &str) -> Option<DefaultSymbol> {
        self.interner.get(rank)
    }

    pub fn find_rank(&self, node: NodeId) -> Option<DefaultSymbol> {
        self.ranks
            .get(node)
            .copied()
            .map(|s| DefaultSymbol::try_from_usize(s as usize).unwrap())
    }

    pub fn node_ranks(&self) -> impl Iterator<Item = (NodeId, DefaultSymbol)> + '_ {
        self.ranks
            .iter()
            .map(|(a, b)| (a, RankSymbol::try_from_usize(*b as usize).unwrap()))
    }
}
