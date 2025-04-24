use std::{collections::HashSet, fmt::Debug, hash::Hash};

use itertools::Itertools;
use serde::{Deserialize, Serialize};
use tree::{node_id::NodeIdMap, util::Loop, walk::RootedTreeWalk, RootedTreeMut, TopologyReplacer};

pub use tree::{NodeId, RootedTree};

pub trait Taxonomy: RootedTree {
    fn has_uniform_depths(&self) -> Option<usize>;

    type RankSym: Eq + Copy + Debug + Hash + Send + Sync + 'static;

    fn rank_sym_str(&self, rank_sym: Self::RankSym) -> Option<&str>;

    fn lookup_rank_sym(&self, rank: &str) -> Option<Self::RankSym>;

    fn find_rank(&self, node: NodeId) -> Option<Self::RankSym>;

    fn find_rank_str(&self, node: NodeId) -> Option<&str> {
        self.find_rank(node).and_then(|sym| self.rank_sym_str(sym))
    }

    type NodeRanks<'a>: Iterator<Item = (NodeId, Self::RankSym)> + 'a
    where
        Self: 'a;

    fn node_ranks(&self) -> Self::NodeRanks<'_>;

    fn rank_ordering(&self) -> Result<Vec<Self::RankSym>, Loop<Self::RankSym>> {
        let rank_graph = self
            .preorder_edges(self.get_root())
            .filter_map(|(parent, child)| {
                let parent_rank = self.find_rank(parent)?;
                let child_rank = self.find_rank(child)?;

                if parent_rank != child_rank {
                    Some((parent_rank, child_rank))
                } else {
                    None
                }
            })
            .into_grouping_map()
            .collect();

        // turbofish because rustc goes crazy
        tree::util::topsort::<Self::RankSym, HashSet<Self::RankSym>>(&rank_graph)
    }
}

#[derive(Debug, Clone)]
pub struct RankedLookup {
    exact: bool,
    first: NodeId,
    rest: Vec<NodeId>,
}

impl RankedLookup {
    pub fn first(&self) -> NodeId {
        self.first
    }

    pub fn is_exact(&self) -> bool {
        self.exact
    }

    pub fn new(exact: bool, node: NodeId) -> Self {
        RankedLookup {
            exact,
            first: node,
            rest: Vec::new(),
        }
    }

    pub fn len(&self) -> usize {
        if self.rest.is_empty() {
            1
        } else {
            self.rest.len()
        }
    }

    pub fn add(&mut self, exact: bool, node: NodeId) {
        if self.exact == exact {
            if self.rest.is_empty() {
                self.rest.push(self.first);
            }
            self.rest.push(node);
        } else if !self.exact {
            self.exact = true;
            self.first = node;
            self.rest.clear();
        }
    }
}

impl From<RankedLookup> for Vec<NodeId> {
    fn from(value: RankedLookup) -> Self {
        if value.len() > 1 {
            value.rest
        } else {
            vec![value.first]
        }
    }
}

pub trait LabeledTaxonomy: Taxonomy + RootedTree {
    type Labels<'a>: Iterator<Item = &'a str>
    where
        Self: 'a;

    fn labels_of(&self, node: NodeId) -> Self::Labels<'_>;

    fn some_label_of(&self, node: NodeId) -> Option<&str> {
        self.labels_of(node).next()
    }

    type NodesWithLabel<'a>: Iterator<Item = NodeId> + 'a
    where
        Self: 'a;

    fn nodes_with_label<'a>(&'a self, label: &'a str) -> Self::NodesWithLabel<'a>;

    fn some_node_with_label(&self, label: &str) -> Option<NodeId> {
        self.nodes_with_label(label).next()
    }

    fn nodes_with_label_and_rank(
        &self,
        rank_sym: Self::RankSym,
        name: &str,
    ) -> Option<RankedLookup> {
        let mut candidate = None;

        for node in self.nodes_with_label(name) {
            let exact = self.find_rank(node) == Some(rank_sym);

            match &mut candidate {
                None => {
                    candidate = Some(RankedLookup::new(exact, node));
                }
                Some(candidate) => {
                    candidate.add(exact, node);
                }
            }
        }
        candidate
    }
}

pub struct MissingRanks<R>(pub Vec<R>);

pub struct Contractor<RankSym, GetRank>
where
    GetRank: Fn(NodeId) -> Option<RankSym>,
{
    ranks_syms: HashSet<RankSym>,
    get_rank: GetRank,
}

#[derive(Clone, Serialize, Deserialize)]
#[repr(transparent)]
#[serde(transparent)]
pub struct ContractedNodes(pub NodeIdMap<NodeId>);

impl<Tree, RankSym, GetRank> TopologyReplacer<Tree> for Contractor<RankSym, GetRank>
where
    Tree: RootedTree,
    RankSym: Eq + Hash,
    GetRank: Fn(NodeId) -> Option<RankSym>,
{
    type Result = ContractedNodes;

    fn replace_topology_with<AddEdge>(self, tree: &Tree, mut add_edge: AddEdge) -> Self::Result
    where
        AddEdge: FnMut(NodeId, NodeId),
    {
        struct StackFrame<'t, Tree: RootedTree + 't> {
            contraction_parent: NodeId,
            children_iter: Tree::Children<'t>,
        }

        impl<'t, Tree: RootedTree> StackFrame<'t, Tree> {
            fn root_new(tree: &'t Tree) -> Self {
                let root = tree.get_root();
                Self::new(root, tree.iter_children(root))
            }

            fn new(contraction_parent: NodeId, children_iter: Tree::Children<'t>) -> Self {
                StackFrame {
                    contraction_parent,
                    children_iter,
                }
            }
        }

        let mut dropped = NodeIdMap::new();
        let mut stack = vec![StackFrame::root_new(tree)];

        while let Some(frame) = stack.last_mut() {
            if let Some(child) = frame.children_iter.next() {
                let child_valid_rank = (self.get_rank)(child)
                    .map_or(false, |rank_sym| self.ranks_syms.contains(&rank_sym));

                let contraction_parent = frame.contraction_parent;

                let grandchildren = tree.iter_children(child);

                if child_valid_rank {
                    add_edge(frame.contraction_parent, child);
                    stack.push(StackFrame::new(child, grandchildren));
                } else {
                    dropped.insert(child, contraction_parent);
                    stack.push(StackFrame::new(contraction_parent, grandchildren));
                }
            } else {
                stack.pop();
            }
        }

        ContractedNodes(dropped)
    }
}

pub trait TaxonomyMut: Taxonomy + RootedTreeMut {
    fn contract(&mut self, new_ranks: HashSet<Self::RankSym>) -> ContractedNodes;
}

pub mod formats;
pub mod labels;
pub mod tree;
pub mod unlabeled_intmap;
