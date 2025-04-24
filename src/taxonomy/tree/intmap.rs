use std::collections::HashSet;

use serde::{Deserialize, Serialize};
use thiserror::Error;

use super::{
    node_id::{Entry, NodeIdMap},
    walk::RootedTreeWalk,
    NodeId, RootedTree, RootedTreeMut, TopologyReplacer,
};


#[derive(Deserialize, Serialize)]
pub struct RootedIntMapTree {
    pub root: NodeId,

    // NodeId -> parent NodeId
    pub parent_ids: NodeIdMap<NodeId>,

    // NodeId -> children Nodeid
    pub children_lookup: NodeIdMap<Vec<NodeId>>,
}

impl RootedIntMapTree {
    pub fn builder() -> RootedIntMapTreeBuilder {
        RootedIntMapTreeBuilder::default()
    }

    pub fn has_node(&self, node: NodeId) -> bool {
        node == self.root || self.parent_ids.contains_key(node)
    }

    pub fn from_rooted_tree<T: RootedTree>(tree: T) -> Result<Self, RootedIntMapTreeBuildError> {
        let mut builder = Self::builder();

        let root = tree.get_root();
        builder.set_root(root)?;

        for (parent, child) in tree.preorder_edges(root) {
            builder.insert_edge(parent, child)?;
        }

        builder.build()
    }

    pub fn check_is_rooted_tree(&self) -> bool {
        let mut visited = HashSet::new();

        for desc in self.postorder_descendants(self.get_root()) {
            visited.insert(desc);
        }

        visited.len() == self.node_count() - 1
    }
}

impl RootedTree for RootedIntMapTree {
    fn get_root(&self) -> NodeId {
        self.root
    }

    fn fixup_node(&self, node: u32) -> Option<NodeId> {
        if self.has_node(NodeId(node)) {
            Some(NodeId(node))
        } else {
            None
        }
    }

    fn node_count(&self) -> usize {
        1 /* root */ + self.parent_ids.len()
    }

    fn is_leaf(&self, node: NodeId) -> bool {
        if let Some(ch) = self.children_lookup.get(node) {
            ch.is_empty()
        } else {
            true
        }
    }

    type Children<'a> =
        std::iter::Copied<std::iter::Flatten<std::option::IntoIter<&'a Vec<NodeId>>>>;

    fn iter_children(&self, node: NodeId) -> Self::Children<'_> {
        self.children_lookup
            .get(node)
            .into_iter()
            .flatten()
            .copied()
    }

    fn find_parent(&self, node: NodeId) -> Option<NodeId> {
        if node == self.root {
            None
        } else {
            self.parent_ids.get(node).copied()
        }
    }
}

impl RootedTreeMut for RootedIntMapTree {
    type UnderlyingTopology = Self;

    fn replace_topology_with<Replacer>(&mut self, replacer: Replacer) -> Replacer::Result
    where
        Replacer: TopologyReplacer<Self::UnderlyingTopology>,
    {
        use super::node_id::Entry;

        let mut parent_ids = NodeIdMap::new();
        let mut children_lookup = NodeIdMap::<Vec<_>>::new();

        let r = replacer.replace_topology_with(self, |parent, child| {
            parent_ids.insert(child, parent);

            match children_lookup.entry(parent) {
                Entry::Occupied(occ) => {
                    occ.into_mut().push(child);
                }
                Entry::Vacant(vac) => {
                    vac.insert(vec![child]);
                }
            }
        });

        self.parent_ids = parent_ids;
        self.children_lookup = children_lookup;

        r
    }
}

#[derive(Default)]
pub struct RootedIntMapTreeBuilder {
    root: Option<NodeId>,

    parent_ids: NodeIdMap<NodeId>,
    children_lookup: NodeIdMap<Vec<NodeId>>,
}

#[derive(Debug, Error)]
pub enum RootedIntMapTreeBuildError {
    #[error("Root not found")]
    MissingRoot,

    #[error("Attempted to set multiple roots: {} and {}", .0, .1)]
    MultipleRoots(NodeId, NodeId),

    #[error("multiple parents found for {}: {} and {}", .0, .1, .2)]
    MultipleParents(NodeId, NodeId, NodeId),
}

impl RootedIntMapTreeBuilder {
    pub fn set_root(&mut self, root: NodeId) -> Result<(), RootedIntMapTreeBuildError> {
        if let Some(r) = self.root {
            return Err(RootedIntMapTreeBuildError::MultipleRoots(r, root));
        } else {
            self.root = Some(root);
        }

        Ok(())
    }

    pub fn insert_edge(&mut self, parent: NodeId, child: NodeId) -> Result<(), RootedIntMapTreeBuildError> {
        match self.parent_ids.entry(child) {
            Entry::Vacant(vac) => {
                vac.insert(parent);
            }
            Entry::Occupied(occ) => {
                return Err(RootedIntMapTreeBuildError::MultipleParents(
                    child,
                    parent,
                    *occ.get(),
                ));
            }
        }

        match self.children_lookup.entry(parent) {
            Entry::Vacant(vac) => {
                vac.insert(vec![child]);
            }
            Entry::Occupied(occ) => {
                occ.into_mut().push(child);
            }
        }

        Ok(())
    }

    pub fn build(self) -> Result<RootedIntMapTree, RootedIntMapTreeBuildError> {
        Ok(RootedIntMapTree {
            root: self.root.ok_or(RootedIntMapTreeBuildError::MissingRoot)?,
            parent_ids: self.parent_ids,
            children_lookup: self.children_lookup,
        })
    }
}
