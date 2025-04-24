use std::{fmt::Display, str::FromStr};

use serde::{Deserialize, Serialize};

pub mod node_id;

pub use node_id::{NodeId, nodeid_to_u32_vec, u32_to_nodeid_vec};

#[derive(Debug, Copy, Clone, Eq, PartialEq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Depth(pub usize);

impl Display for Depth {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl FromStr for Depth {
    type Err = <u32 as FromStr>::Err;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Depth(s.parse()?))
    }
}

impl From<Depth> for usize {
    fn from(node: Depth) -> Self {
        node.0
    }
}

impl Depth {
    pub fn zero() -> Self {
        Depth(0)
    }

    pub fn incr(self) -> Self {
        Depth(self.0 + 1)
    }

    pub fn decr(self) -> Option<Self> {
        if self.0 > 0 {
            Some(Depth(self.0 - 1))
        } else {
            None
        }
    }
}

pub trait RootedTree: Sized {
    fn get_root(&self) -> NodeId;

    fn fixup_node(&self, node: u32) -> Option<NodeId>;

    fn node_count(&self) -> usize;

    fn find_parent(&self, node: NodeId) -> Option<NodeId>;

    type Children<'a>: DoubleEndedIterator<Item = NodeId> + 'a
    where
        Self: 'a;

    fn iter_children(&self, node: NodeId) -> Self::Children<'_>;

    fn is_leaf(&self, node: NodeId) -> bool {
        self.iter_children(node).next().is_none()
    }

    fn leftmost_descendant_leaf(&self, node: NodeId) -> NodeId {
        let mut node = node;

        while let Some(child) = self.iter_children(node).next() {
            node = child;
        }

        node
    }
}

pub trait TopologyReplacer<Tree> {
    type Result = ();

    fn replace_topology_with<AddEdge>(self, tree: &Tree, add_edge: AddEdge) -> Self::Result
    where
        AddEdge: FnMut(NodeId, NodeId);
}

pub trait RootedTreeMut : RootedTree {
    type UnderlyingTopology: RootedTree = Self;

    fn replace_topology_with<Replacer>(&mut self, replacer: Replacer) -> Replacer::Result
    where
        Replacer: TopologyReplacer<Self::UnderlyingTopology>;
}

impl<'t, Tree: RootedTree> RootedTree for &'t Tree {
    fn get_root(&self) -> NodeId {
        (*self).get_root()
    }

    fn fixup_node(&self, node: u32) -> Option<NodeId> {
        (*self).fixup_node(node)
    }

    fn node_count(&self) -> usize {
        (*self).node_count()
    }

    fn find_parent(&self, node: NodeId) -> Option<NodeId> {
        (*self).find_parent(node)
    }

    type Children<'a> = Tree::Children<'a>
    where
        't: 'a;

    fn iter_children(&self, node: NodeId) -> Self::Children<'_> {
        (*self).iter_children(node)
    }
}

pub mod intmap;
pub mod newick;
pub mod vec;
pub mod walk;
pub mod postorder_ann;

pub(crate) mod util;
