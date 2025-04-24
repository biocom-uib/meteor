use std::{collections::HashSet, mem, sync::Mutex};

use itertools::Itertools;
use serde::{Deserialize, Serialize};
use string_interner::{DefaultSymbol, Symbol};
use thiserror::Error;

use super::{
    labels::intmap::{IntMapRanks, RankSymbol},
    tree::{
        intmap::{RootedIntMapTree, RootedIntMapTreeBuildError, RootedIntMapTreeBuilder},
        walk::RootedTreeWalk,
        RootedTreeMut,
    },
    ContractedNodes, Contractor, NodeId, RootedTree, Taxonomy, TaxonomyMut, TopologyReplacer,
};

#[derive(Deserialize, Serialize)]
pub struct UnlabeledIntMapTaxonomy {
    pub tree: RootedIntMapTree,
    pub ranks: IntMapRanks,

    has_uniform_depths_cache: Mutex<Option<Option<usize>>>,
}


impl UnlabeledIntMapTaxonomy {
    pub const FORMAT_VERSION: u32 = 3;

    pub fn builder() -> GenericTaxonomyBuilder {
        GenericTaxonomyBuilder::default()
    }

    pub fn from_taxonomy<T: Taxonomy>(taxonomy: T) -> Result<Self, TaxonomyBuildError> {
        let mut builder = Self::builder();

        let root = taxonomy.get_root();
        builder.set_root(root)?;

        for (parent, child) in taxonomy.preorder_edges(root) {
            builder.insert_edge(parent, child)?;
        }

        for (node, rank_sym) in taxonomy.node_ranks() {
            builder.set_rank(node, taxonomy.rank_sym_str(rank_sym).unwrap());
        }

        builder.build()
    }

    pub fn topology_health_check(&self) -> bool {
        let check1 = || {
            self.tree.children_lookup.iter().all(|(parent, children)| {
                children
                    .iter()
                    .all(|&child| self.tree.find_parent(child) == Some(parent))
            })
        };

        let check2 = || {
            self.tree.parent_ids.iter().all(|(child, &parent)| {
                let r = self
                    .tree
                    .children_lookup
                    .get(parent)
                    .map_or(false, |children| children.contains(&child));

                if !r {
                    let parent_rank = self.find_rank(parent).and_then(|s| self.rank_sym_str(s));
                    let child_rank = self.find_rank(child).and_then(|s| self.rank_sym_str(s));
                    dbg!((parent, parent_rank, child, child_rank));
                }

                r
            })
        };

        check1() && check2()
    }
}

impl RootedTree for UnlabeledIntMapTaxonomy {
    fn get_root(&self) -> NodeId {
        self.tree.get_root()
    }

    fn fixup_node(&self, node: u32) -> Option<NodeId> {
        self.tree.fixup_node(node)
    }

    fn node_count(&self) -> usize {
        self.tree.node_count()
    }

    fn find_parent(&self, node: NodeId) -> Option<NodeId> {
        self.tree.find_parent(node)
    }

    type Children<'a> = <RootedIntMapTree as RootedTree>::Children<'a>;

    fn iter_children(&self, node: NodeId) -> Self::Children<'_> {
        self.tree.iter_children(node)
    }

    fn is_leaf(&self, node: NodeId) -> bool {
        self.tree.is_leaf(node)
    }
}

impl Taxonomy for UnlabeledIntMapTaxonomy {
    fn has_uniform_depths(&self) -> Option<usize> {
        let mut cache = self.has_uniform_depths_cache.lock().unwrap();

        if let Some(value) = &*cache {
            *value
        } else {
            let mut leaves_depths = self
                .postorder_descendants(self.get_root())
                .filter(|node| self.is_leaf(*node))
                .map(|node| self.strict_ancestors(node).count())
                .peekable();

            let &depth = leaves_depths.peek().expect("The tree has no leaves");

            let value = if leaves_depths.all_equal() {
                Some(depth)
            } else {
                None
            };
            *cache = Some(value);
            value
        }
    }

    type RankSym = DefaultSymbol;

    fn rank_sym_str(&self, rank_sym: Self::RankSym) -> Option<&str> {
        self.ranks.rank_sym_str(rank_sym)
    }

    fn lookup_rank_sym(&self, rank: &str) -> Option<Self::RankSym> {
        self.ranks.lookup_rank_sym(rank)
    }

    fn find_rank(&self, node: NodeId) -> Option<Self::RankSym> {
        self.ranks.find_rank(node)
    }

    type NodeRanks<'a> = impl Iterator<Item = (NodeId, RankSymbol)> + 'a;

    fn node_ranks(&self) -> Self::NodeRanks<'_> {
        self.ranks.node_ranks()
    }
}

impl RootedTreeMut for UnlabeledIntMapTaxonomy {
    type UnderlyingTopology = RootedIntMapTree;

    fn replace_topology_with<Replacer>(&mut self, replacer: Replacer) -> Replacer::Result
    where
        Replacer: TopologyReplacer<Self::UnderlyingTopology>,
    {
        let r = self.tree.replace_topology_with(replacer);

        *self.has_uniform_depths_cache.get_mut().unwrap() = None;

        let mut ranks = mem::take(&mut self.ranks.ranks);
        ranks.retain(|node, _| self.tree.has_node(node));
        self.ranks.ranks = ranks;

        r
    }
}

impl TaxonomyMut for UnlabeledIntMapTaxonomy {
    fn contract(&mut self, new_ranks: HashSet<Self::RankSym>) -> ContractedNodes {
        self.tree.replace_topology_with(Contractor {
            ranks_syms: new_ranks,
            get_rank: |node| self.ranks.find_rank(node),
        })
    }
}

#[derive(Default)]
pub struct GenericTaxonomyBuilder {
    tree_builder: RootedIntMapTreeBuilder,
    ranks: IntMapRanks,
}

#[derive(Debug, Error)]
pub enum TaxonomyBuildError {
    #[error("Error building tree")]
    TreeBuilderError(#[from] RootedIntMapTreeBuildError)
}

impl GenericTaxonomyBuilder {
    pub fn set_root(&mut self, root: NodeId) -> Result<(), TaxonomyBuildError> {
        Ok(self.tree_builder.set_root(root)?)
    }

    pub fn set_rank(&mut self, node: NodeId, rank: &str) {
        let sym = self.ranks.interner.get_or_intern(rank).to_usize() as u32;

        self.ranks.ranks.insert(node, sym);
    }

    pub fn insert_edge(&mut self, parent: NodeId, child: NodeId) -> Result<(), TaxonomyBuildError> {
        Ok(self.tree_builder.insert_edge(parent, child)?)
    }

    pub fn build(self) -> Result<UnlabeledIntMapTaxonomy, TaxonomyBuildError> {
        Ok(UnlabeledIntMapTaxonomy {
            tree: self.tree_builder.build()?,
            has_uniform_depths_cache: Mutex::new(None),
            ranks: self.ranks,
        })
    }
}
