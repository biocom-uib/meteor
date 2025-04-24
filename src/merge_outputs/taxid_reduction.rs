use std::sync::Arc;

use clap::ValueEnum;
use itertools::Itertools;
use polars::{
    datatypes::{PlSmallStr, UInt32Chunked},
    lazy::dsl::Expr,
    prelude::NewChunkedArray,
    series::Series,
};

use crate::{
    csv::polars::ExprExt,
    taxonomy::{
        tree::{
            walk::{LazyLcaCache, LcaCache, LcaCacheIndex, LcaCacheRange, LcaMapResult},
            NodeId,
        }, Taxonomy
    },
};

pub struct NodeRangeList<R> {
    node_ranges: Vec<R>,
}

impl<'a> IntoIterator for &'a NodeRangeList<LcaCacheRange> {
    type Item = Option<LcaCacheRange>;

    type IntoIter = impl Iterator<Item = Self::Item> + 'a;

    fn into_iter(self) -> Self::IntoIter {
        self.node_ranges.iter().copied().map(Some)
    }
}

impl<'a> IntoIterator for &'a NodeRangeList<Option<LcaCacheRange>> {
    type Item = Option<LcaCacheRange>;

    type IntoIter = impl Iterator<Item = Self::Item> + 'a;

    fn into_iter(self) -> Self::IntoIter {
        self.node_ranges.iter().copied()
    }
}

impl NodeRangeList<LcaCacheRange> {
    pub fn from_node_ids<I: IntoIterator<Item = NodeId>>(
        lca_cache: &LcaCache,
        node_ids: I,
    ) -> LcaMapResult<Self> {
        let node_ranges = node_ids
            .into_iter()
            .map(|node| lca_cache.node_range(node))
            .try_collect()?;

        Ok(NodeRangeList {
            node_ranges
        })
    }
}

impl<R> NodeRangeList<R>
where
    for<'a> &'a Self: IntoIterator<Item = Option<LcaCacheRange>>
{
    pub fn is_maximal(&self, range: LcaCacheRange) -> bool {
        !self.into_iter().flatten().any(|other| range < other)
    }

    pub fn most_specific_index_iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.into_iter()
            .positions(|opt_range| opt_range.is_some_and(|range| self.is_maximal(range)))
    }

    pub fn lca_reduce(&self, lca_cache: &LcaCache) -> Option<LcaMapResult<NodeId>> {
        let lca = self
            .into_iter()
            .flatten()
            .map(LcaCacheIndex::from)
            .reduce(|range1, range2| lca_cache.lca_index(range1, range2))?;

        Some(lca_cache.node_id(lca))
    }

    pub fn most_specific_iter<'a>(&'a self, lca_cache: &'a LcaCache) -> impl Iterator<Item = LcaMapResult<NodeId>> + 'a {
        self.into_iter()
            .flatten()
            .filter(|&range| self.is_maximal(range))
            .map(|range| lca_cache.node_id(range))
    }
}

impl NodeRangeList<Option<LcaCacheRange>> {
    pub fn from_u32chunked(lca_cache: &LcaCache, node_ids: &UInt32Chunked) -> LcaMapResult<Self> {
        Self::from_node_ids_opt(lca_cache, node_ids.iter().map(|opt| opt.map(NodeId)))
    }


    pub fn from_node_ids_opt(lca_cache: &LcaCache, node_ids: impl IntoIterator<Item = Option<NodeId>>) -> LcaMapResult<Self> {
        let node_ranges = node_ids
            .into_iter()
            .map(|node| Some(lca_cache.node_range(node?)))
            .map(Option::transpose)
            .try_collect()?;

        Ok(NodeRangeList {
            node_ranges
        })
    }
}

#[derive(Copy, Clone, Eq, PartialEq, PartialOrd, Ord, Hash, ValueEnum)]
pub enum TaxIdReduction {
    /// Compute the single LCA of all the taxids
    Lca,
    /// Drop any ancestors of other present taxids
    KeepMostSpecific,
}

impl TaxIdReduction {
    /// Expr should resolve to a list[32] column
    pub fn most_specific_indices<Tax>(expr: Expr, lca_cache: Arc<LazyLcaCache<Tax>>) -> Expr
    where
        Tax: Taxonomy + Send + Sync + 'static,
    {
        expr.downcast_map_lists(Series::u32, move |taxids| {
            let lca_cache = &**lca_cache;

            let taxids = NodeRangeList::from_u32chunked(lca_cache, taxids)
                .expect("Error mapping beteen taxids and LCA ranges");

            taxids.most_specific_index_iter()
                .map(|ix| Some(ix as u32))
                .collect::<UInt32Chunked>()
        })
    }

    /// Expr should resolve to a list[32] column
    pub fn apply<Tax>(self, expr: Expr, lca_cache: Arc<LazyLcaCache<Tax>>) -> Expr
    where
        Tax: Taxonomy + Send + Sync + 'static,
    {
        let reducer = self.into_reducer_fn();

        expr.downcast_map_lists(Series::u32, move |taxids| {
            let lca_cache = &**lca_cache;

            let taxids = NodeRangeList::from_node_ids(lca_cache, taxids.iter().flatten().map(NodeId))
                .expect("Error mapping beteen taxids and LCA ranges");

            reducer(lca_cache, taxids)
        })
    }

    fn reduce_lca(lca_cache: &LcaCache, node_id_list: NodeRangeList<LcaCacheRange>) -> UInt32Chunked {
        let lca = if let Some(lca) = node_id_list.lca_reduce(lca_cache) {
            let lca = lca.expect("Error mapping between taxids and LCA ranges");
            Some(lca.into())
        } else {
            None
        };

        UInt32Chunked::from_iter_values(PlSmallStr::EMPTY, lca.into_iter())
    }

    fn reduce_keep_most_specific(lca_cache: &LcaCache, node_id_list: NodeRangeList<LcaCacheRange>) -> UInt32Chunked {
        let reduced = node_id_list
            .most_specific_iter(lca_cache)
            .map(|res| res.expect("Error mapping between taxids and LCA ranges"))
            .map(NodeId::into);

        UInt32Chunked::from_iter_values(PlSmallStr::EMPTY, reduced)
    }

    fn into_reducer_fn(self) -> fn(&LcaCache, NodeRangeList<LcaCacheRange>) -> UInt32Chunked {
        match self {
            TaxIdReduction::Lca => Self::reduce_lca,
            TaxIdReduction::KeepMostSpecific => Self::reduce_keep_most_specific,
        }
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use crate::taxonomy::{
        formats::newick::{self},
        tree::walk::{LcaCache, LcaMapResult},
        NodeId,
    };

    use super::NodeRangeList;

    #[test]
    fn test_keep_most_specific_taxids() {
        let tax = newick::tests::sample_taxonomy();
        let lca_cache = LcaCache::compute(&tax);

        let nodes = [
            NodeId(7),
            NodeId(11),
            NodeId(12),
            NodeId(21),
            NodeId(3),
            NodeId(28),
            NodeId(31),
        ];

        let node_list = NodeRangeList::from_node_ids(&lca_cache, nodes).unwrap();

        assert_eq!(
            node_list
                .most_specific_iter(&lca_cache)
                .collect::<LcaMapResult<Vec<NodeId>>>(),
            Ok(vec![
                NodeId(7),
                NodeId(11),
                NodeId(12),
                NodeId(28),
                NodeId(31)
            ]),
        );
    }

    #[test]
    fn test_keep_most_specific_indices() {
        let tax = newick::tests::sample_taxonomy();
        let lca_cache = LcaCache::compute(&tax);

        let nodes = [
            Some(NodeId(7)),
            Some(NodeId(11)),
            None,
            Some(NodeId(21)),
            Some(NodeId(3)),
            Some(NodeId(28)),
            Some(NodeId(31)),
        ];

        let node_list = NodeRangeList::from_node_ids_opt(&lca_cache, nodes).unwrap();

        assert_eq!(
            node_list.most_specific_index_iter().collect_vec(),
            vec![0, 1, 5, 6],
        );
    }
}
