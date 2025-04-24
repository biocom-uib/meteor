use newick_rs::SimpleTree;
use serde::{Deserialize, Serialize};

use crate::taxonomy::labels::vec::RootedVecTreeLabels;

use super::{NodeId, RootedTree};


#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct RootedVecTreeDepths {
    depths: Vec<usize>,
}

impl RootedVecTreeDepths {
    pub fn get(&self, node: NodeId) -> Option<usize> {
        self.depths.get(node.0 as usize).copied()
    }

    pub fn iter(&self) -> impl Iterator<Item = (NodeId, usize)> + '_ {
        self.depths
            .iter()
            .enumerate()
            .map(|(i, &depth)| (NodeId::from(i as u32), depth))
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RootedVecTree {
    root: NodeId,
    parent_ids: Vec<NodeId>,
    children_lookup: Vec<Vec<NodeId>>,
}

impl RootedVecTree {
    pub fn from_simple_tree(value: SimpleTree) -> (Self, RootedVecTreeLabels, RootedVecTreeDepths) {
        let mut labels = RootedVecTreeLabels::default();

        let mut depths = RootedVecTreeDepths::default();

        let mut new_node = |label: String, depth: usize| {
            let node_id = labels.push(label);

            depths.depths.push(depth);

            node_id
        };

        struct StackFrame {
            node_id: NodeId,
            children: Vec<NodeId>,
            orig_children_iter: std::vec::IntoIter<SimpleTree>,
        }

        impl StackFrame {
            fn new(node_id: NodeId, orig_children: Vec<SimpleTree>) -> Self {
                StackFrame {
                    node_id,
                    children: Vec::new(),
                    orig_children_iter: orig_children.into_iter(),
                }
            }
        }

        let root = new_node(value.name, 0);
        let mut stack = vec![StackFrame::new(root, value.children)];

        let mut parent_ids = vec![root];
        let mut children_lookup = vec![vec![]];

        while let (depth, Some(frame)) = (stack.len(), stack.last_mut()) {
            if let Some(child) = frame.orig_children_iter.next() {
                let child_id = new_node(child.name, depth);

                assert!(parent_ids.len() as u32 == u32::from(child_id));
                parent_ids.push(frame.node_id);
                frame.children.push(child_id);

                assert!(children_lookup.len() as u32 == u32::from(child_id));
                children_lookup.push(vec![]);

                stack.push(StackFrame::new(child_id, child.children));
            } else {
                let popped = stack.pop().unwrap();
                children_lookup[popped.node_id.0 as usize] = popped.children;
            }
        }

        let tree = RootedVecTree {
            root,
            parent_ids,
            children_lookup,
        };

        (tree, labels, depths)
    }

    pub fn into_simple_tree(&self, mut labels: RootedVecTreeLabels) -> SimpleTree {
        fn build(tree: &RootedVecTree, node_id: NodeId, labels: &mut RootedVecTreeLabels) -> SimpleTree {
            let node_index = node_id.0 as usize;

            let children = tree.children_lookup[node_index]
                .iter()
                .map(|child_id| build(tree, *child_id, labels))
                .collect();

            SimpleTree {
                name: std::mem::take(labels.get_mut(node_id).unwrap()),
                children,
                length: None,
            }
        }

        build(self, self.root, &mut labels)
    }

    pub fn has_node(&self, node: u32) -> bool {
        (node as usize) < self.parent_ids.len()
    }

    pub fn nodes(&self) -> impl Iterator<Item = NodeId> {
        (0..self.parent_ids.len()).map(|id| NodeId::from(id as u32))
    }
}

impl RootedTree for RootedVecTree {
    fn get_root(&self) -> NodeId {
        self.root
    }

    fn fixup_node(&self, node: u32) -> Option<NodeId> {
        if self.has_node(node) {
            Some(NodeId::from(node))
        } else {
            None
        }
    }

    fn node_count(&self) -> usize {
        self.parent_ids.len()
    }

    fn find_parent(&self, node: NodeId) -> Option<NodeId> {
        if node == self.root {
            None
        } else {
            self.parent_ids.get(node.0 as usize).copied()
        }
    }

    type Children<'a> = std::iter::Copied<std::slice::Iter<'a, NodeId>>;

    fn iter_children(&self, node: NodeId) -> Self::Children<'_> {
        self.children_lookup[node.0 as usize].iter().copied()
    }
}
