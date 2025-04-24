use std::{cmp::Ordering, sync::{Arc, LazyLock}};

use itertools::Itertools;

use range_minimum_query::Rmq;

use super::{node_id::NodeIdMap, Depth, NodeId, RootedTree};

pub trait RootedTreeWalk: RootedTree {
    /// node1 could be node2
    fn is_ancestor(&self, node1: NodeId, node2: NodeId) -> bool {
        node1 == node2 || self.strict_ancestors(node2).contains(&node1)
    }

    fn are_independent(&self, node1: NodeId, node2: NodeId) -> bool {
        !self.is_ancestor(node1, node2) && !self.is_ancestor(node2, node1)
    }

    /// Does not include node
    fn strict_ancestors(&self, node: NodeId) -> impl Iterator<Item = NodeId>;

    /// Same as strict_ancestors, but includes node
    fn ancestors(&self, node: NodeId) -> impl Iterator<Item = NodeId> {
        std::iter::once(node).chain(self.strict_ancestors(node))
    }

    fn preorder_edges(&self, node: NodeId) -> impl Iterator<Item = (NodeId, NodeId)>;

    // Includes node
    fn postorder_descendants(&self, node: NodeId) -> impl Iterator<Item = NodeId>;

    // Includes node, not child. Same as advancing postorder_descendants until child is reached
    fn postorder_descendants_from(
        &self,
        node: NodeId,
        child: NodeId,
    ) -> Option<impl Iterator<Item = NodeId>>;

    fn euler_tour(&self, node: NodeId) -> impl Iterator<Item = (NodeId, Depth)>;
}

#[allow(refining_impl_trait)]
impl<Tree: RootedTree> RootedTreeWalk for Tree {
    fn strict_ancestors(&self, node: NodeId) -> AncestorsIter<Self> {
        AncestorsIter::new(self, node)
    }

    fn preorder_edges(&self, node: NodeId) -> PreOrderEdgesIter<Self> {
        PreOrderEdgesIter::new(self, node)
    }

    fn postorder_descendants(&self, node: NodeId) -> PostOrderIter<Self> {
        PostOrderIter::new(self, node)
    }

    fn postorder_descendants_from(
        &self,
        node: NodeId,
        child: NodeId,
    ) -> Option<PostOrderIter<Self>> {
        PostOrderIter::new_from(self, node, child)
    }

    fn euler_tour(&self, node: NodeId) -> EulerTourIter<Self> {
        let depth = self.strict_ancestors(node).count();
        EulerTourIter::new(self, node, Depth(depth))
    }
}

#[derive(Clone)]
pub struct AncestorsIter<'t, Tree> {
    tree: &'t Tree,
    current: NodeId,
}

impl<'t, Tree> AncestorsIter<'t, Tree> {
    pub fn new(tree: &'t Tree, current: NodeId) -> Self {
        Self { tree, current }
    }
}

impl<'t, Tree: RootedTree> Iterator for AncestorsIter<'t, Tree> {
    type Item = NodeId;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(parent) = self.tree.find_parent(self.current) {
            self.current = parent;
            Some(parent)
        } else {
            // root
            None
        }
    }
}

#[derive(Clone)]
pub struct PreOrderEdgesIter<'t, Tree: RootedTree> {
    tree: &'t Tree,
    stack: Vec<(NodeId, Tree::Children<'t>)>,
}

impl<'t, Tree: RootedTree> PreOrderEdgesIter<'t, Tree> {
    pub fn new(tree: &'t Tree, node: NodeId) -> Self {
        Self {
            tree,
            stack: vec![(node, tree.iter_children(node))],
        }
    }
}

impl<'t, Tree: RootedTree> Iterator for PreOrderEdgesIter<'t, Tree> {
    type Item = (NodeId, NodeId);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((current, children)) = self.stack.last_mut() {
            if let Some(child) = children.next() {
                let current = *current;
                self.stack.push((child, self.tree.iter_children(child)));
                Some((current, child))
            } else {
                self.stack.pop();
                self.next()
            }
        } else {
            None
        }
    }
}

#[derive(Clone)]
pub struct PostOrderIter<'t, Tree: RootedTree> {
    tree: &'t Tree,
    stack: Vec<(NodeId, Tree::Children<'t>)>,
}

impl<'t, Tree: RootedTree> PostOrderIter<'t, Tree> {
    pub fn new(tree: &'t Tree, node: NodeId) -> Self {
        Self {
            tree,
            stack: vec![(node, tree.iter_children(node))],
        }
    }

    pub fn new_from(tree: &'t Tree, node: NodeId, child: NodeId) -> Option<Self> {
        let child_ancestors = tree
            .strict_ancestors(child)
            .take_while(|ancestor| *ancestor != node);

        let mut path_to_child: Vec<NodeId> = std::iter::once(child)
            .chain(child_ancestors)
            .collect();

        if tree.find_parent(*path_to_child.last().unwrap()) != Some(node) {
            return None;
        } else {
            path_to_child.push(node);
        }

        path_to_child.reverse();

        let mut stack = Vec::new();

        for (&intermediate, &next_intermediate) in path_to_child.iter().zip(&path_to_child[1..]) {
            // first is (node, node.iter_children()[0])
            // last is (tree.find_parent(child), child)
            let mut intermediate_children = tree.iter_children(intermediate);

            loop {
                let Some(intermediate_child) = intermediate_children.next() else {
                    panic!("The path from node to child is not correct");
                };

                if intermediate_child == next_intermediate {
                    stack.push((intermediate, intermediate_children));
                    break;
                }
            }
        }

        Some(Self {
            tree,
            stack,
        })
    }
}

impl<'t, Tree: RootedTree> Iterator for PostOrderIter<'t, Tree> {
    type Item = NodeId;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((current, children)) = self.stack.last_mut() {
            if let Some(child) = children.next() {
                let grand_children = self.tree.iter_children(child);
                self.stack.push((child, grand_children));
                self.next()
            } else {
                let current = *current;
                self.stack.pop();
                Some(current)
            }
        } else {
            None
        }
    }
}

struct EulerTourIterStackFrame<'t, Tree: RootedTree + 't> {
    node: NodeId,
    depth: Depth,
    children: Tree::Children<'t>,
}

impl<'t, Tree> Clone for EulerTourIterStackFrame<'t, Tree>
where
    Tree: RootedTree + 't,
    Tree::Children<'t>: Clone,
{
    fn clone(&self) -> Self {
        Self {
            node: self.node,
            depth: self.depth,
            children: self.children.clone(),
        }
    }
}

pub struct EulerTourIter<'t, Tree: RootedTree> {
    tree: &'t Tree,
    stack: Vec<EulerTourIterStackFrame<'t, Tree>>,
}

impl<'t, Tree> Clone for EulerTourIter<'t, Tree>
where
    Tree: RootedTree,
    Tree::Children<'t>: Clone,
{
    fn clone(&self) -> Self {
        Self {
            tree: self.tree,
            stack: self.stack.clone(),
        }
    }
}

impl<'t, Tree: RootedTree> EulerTourIter<'t, Tree> {
    pub fn new(tree: &'t Tree, node_id: NodeId, start_depth: Depth) -> Self {
        Self {
            tree,
            stack: vec![EulerTourIterStackFrame {
                node: node_id,
                depth: start_depth,
                children: tree.iter_children(node_id),
            }],
        }
    }
}

impl<'t, Tree: RootedTree> Iterator for EulerTourIter<'t, Tree> {
    type Item = (NodeId, Depth);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(frame) = self.stack.last_mut() {
            let current = frame.node;
            let depth = frame.depth;

            if let Some(child) = frame.children.next() {
                self.stack.push(EulerTourIterStackFrame {
                    node: child,
                    depth: depth.incr(),
                    children: self.tree.iter_children(child),
                });
            } else {
                self.stack.pop();
            }

            Some((current, depth))
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = self.tree.node_count();

        (n, Some(2*n+1))
    }
}

#[derive(Copy, Clone, Debug)]
pub struct LcaCacheRange {
    first_index: usize,
    last_index: usize,
}

impl LcaCacheRange {
    pub fn is_ancestor(&self, other: &LcaCacheRange) -> bool {
        self.first_index <= other.first_index && other.last_index <= self.last_index
    }

    pub fn is_independent(&self, other: &LcaCacheRange) -> bool {
        !self.is_ancestor(other) && other.is_ancestor(self)
    }
}

impl PartialEq for LcaCacheRange {
    fn eq(&self, other: &Self) -> bool {
        self.first_index == other.first_index
    }
}

impl Eq for LcaCacheRange {}

// reachability order
// u <= v  <==>  u is an ancestor of v
impl PartialOrd for LcaCacheRange {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        assert_eq!(
            self.is_ancestor(other) || other.is_ancestor(self),
            other.first_index <= self.last_index && self.first_index <= other.last_index
        );

        if self == other {
            Some(Ordering::Equal)
        } else if self.is_ancestor(other) {
            Some(Ordering::Less)
        } else if other.is_ancestor(self) {
            Some(Ordering::Greater)
        } else {
            None
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct LcaCacheIndex(usize);

impl From<LcaCacheRange> for LcaCacheIndex {
    fn from(range: LcaCacheRange) -> Self {
        LcaCacheIndex(range.last_index)
    }
}

pub struct LcaCache {
    euler_tour_ranges: NodeIdMap<LcaCacheRange>,
    euler_tour_nodes: Vec<NodeId>,
    euler_tour_depths: Rmq,
}

pub type LazyLcaCache<Tree: RootedTree + Send + Sync + 'static> =
    LazyLock<LcaCache, impl FnOnce() -> LcaCache + Send + Sync + 'static>;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum LcaMapError {
    NodeIdToRange(NodeId),
    IndexToNodeId(LcaCacheIndex),
}

pub type LcaMapResult<T> = Result<T, LcaMapError>;

impl LcaCache {
    pub fn compute<Tree>(tree: &Tree) -> Self
    where
        Tree: RootedTree
    {
        let mut ranges = NodeIdMap::new();
        let mut nodes = Vec::with_capacity(2*tree.node_count() + 1);
        let mut depths = Vec::with_capacity(2*tree.node_count() + 1);

        use super::node_id::Entry;

        for (node, depth) in tree.euler_tour(tree.get_root()) {
            let i = nodes.len();
            nodes.push(node);

            match ranges.entry(node) {
                Entry::Vacant(vac) => {
                    vac.insert(LcaCacheRange { first_index: i, last_index: i });
                },
                Entry::Occupied(mut occ) => {
                    occ.get_mut().last_index = i;
                },
            };

            depths.push(depth.0);
        }

        LcaCache {
            euler_tour_ranges: ranges,
            euler_tour_nodes: nodes,
            euler_tour_depths: Rmq::from_iter(depths),
        }
    }

    #[define_opaque(LazyLcaCache)]
    pub fn compute_lazy<Tree: RootedTree + Send + Sync + 'static>(tree: Arc<Tree>) -> LazyLcaCache<Tree> {
        LazyLock::new(move || Self::compute(&*tree))
    }

    pub fn node_range(&self, node: NodeId) -> Result<LcaCacheRange, LcaMapError> {
        if let Some(range) = self.euler_tour_ranges.get(node) {
            Ok(*range)
        } else {
            Err(LcaMapError::NodeIdToRange(node))
        }
    }

    pub fn node_index(&self, node: NodeId) -> Result<LcaCacheIndex, LcaMapError> {
        match self.node_range(node) {
            Ok(range) => Ok(range.into()),
            Err(e) => Err(e),
        }
    }

    pub fn node_id(&self, index: impl Into<LcaCacheIndex>) -> Result<NodeId, LcaMapError> {
        let index: LcaCacheIndex = index.into();
        self.euler_tour_nodes.get(index.0).copied().ok_or(LcaMapError::IndexToNodeId(index))
    }

    pub fn lca_index(
        &self,
        node1: impl Into<LcaCacheIndex>,
        node2: impl Into<LcaCacheIndex>,
    ) -> LcaCacheIndex {
        let node1_index: LcaCacheIndex = node1.into();
        let node2_index: LcaCacheIndex = node2.into();

        let [index1, index2] = std::cmp::minmax(node1_index.0, node2_index.0);

        let lca_index = self
            .euler_tour_depths
            .range_minimum(index1..=index2)
            .expect("Range should not be empty");

        LcaCacheIndex(lca_index)
    }

    pub fn lca_node(&self, node1: NodeId, node2: NodeId) -> Result<NodeId, LcaMapError> {
        self.node_id(self.lca_index(self.node_range(node1)?, self.node_range(node2)?))
    }

    pub fn lca_reduce(&self, iter: impl IntoIterator<Item = NodeId>) -> Option<Result<NodeId, LcaMapError>> {
        let mut iter = iter.into_iter();

        let mut lhs = if let Some(first) = iter.next() {
            match self.node_index(first) {
                Ok(idx) => idx,
                Err(e) => return Some(Err(e)),
            }
        } else {
            return None;
        };

        for rhs in iter {
            match self.node_index(rhs) {
                Ok(rhs) => lhs = self.lca_index(lhs, rhs),
                Err(e) => return Some(Err(e)),
            }
        }

        Some(self.node_id(lhs))
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use std::cmp::Ordering;

    use itertools::Itertools;
    use newick_rs::SimpleTree;

    use crate::taxonomy::{
        formats::newick::{self, NewickTaxonomy}, tree::{newick::newick_literal, walk::{LcaCache, RootedTreeWalk}}, LabeledTaxonomy, NodeId, RootedTree
    };

    #[test]
    fn ancestors() {
        fn chain(n: u32) -> SimpleTree {
            if n == 0 {
                newick_literal!(0)
            } else {
                let sub = chain(n - 1);
                newick_literal!( [ (sub) ] n )
            }
        }

        let simple_t = chain(5);

        let ranks = "a,b,c,d,e,f".split(',').map(|s| s.to_owned()).collect_vec();

        let t = NewickTaxonomy::from_simple_tree(simple_t, ranks);

        assert_eq!(
            t.strict_ancestors(t.nodes_with_label("0").next().unwrap())
                .map(|node| t.labels_of(node).next().unwrap())
                .collect_vec(),
            vec!["1", "2", "3", "4", "5"]
        );
    }

    fn test_tree_0() -> impl LabeledTaxonomy {
        /*
            labels        ranks
              0             a
            /   \
           1     2          b
               /   \
              3     4       c
        */

        let sub = newick_literal! { [3, 4] 2 };
        let simple_t = newick_literal! { [ 1, (sub) ] 0 };

        let ranks = "a,b,c".split(',').map(|s| s.to_owned()).collect_vec();

        let tax = NewickTaxonomy::from_simple_tree(simple_t, ranks);
        assert!(newick::tests::verify_node_ids_equal_labels(&tax));
        tax

    }

    fn test_tree_1() -> impl LabeledTaxonomy {
        /*
            labels        ranks
              0             a
            /   \
           1     2          b
               / | \
              3  4  6       c
                 |
                 5          d
        */

        let sub4 = newick_literal!{ [5] 4 };
        let sub2 = newick_literal! { [3, (sub4), 6] 2 };
        let simple_t = newick_literal! { [ 1, (sub2) ] 0 };

        let ranks = "a,b,c,d".split(',').map(|s| s.to_owned()).collect_vec();

        let tax = NewickTaxonomy::from_simple_tree(simple_t, ranks);
        assert!(newick::tests::verify_node_ids_equal_labels(&tax));
        tax
    }

    #[test]
    fn preorder_edges() {
        let t = test_tree_0();

        assert_eq!(
            t.preorder_edges(t.get_root())
                .map(|(_parent, node)| t.some_label_of(node).unwrap())
                .collect_vec(),
            vec!["1", "2", "3", "4"]
        );
    }

    #[test]
    fn postorder_descendants() {
        let t = test_tree_0();

        assert_eq!(
            t.postorder_descendants(t.get_root())
                .map(|node| t.some_label_of(node).unwrap())
                .collect_vec(),
            vec!["1", "3", "4", "2", "0"]
        );

        let t = test_tree_1();

        assert_eq!(
            t.postorder_descendants(t.get_root())
                .map(|node| t.some_label_of(node).unwrap())
                .collect_vec(),
            vec!["1", "3", "5", "4", "6", "2", "0"]
        );
    }

    #[test]
    fn postorder_descendants_from() {
        let t = test_tree_0();

        assert_eq!(
            t.postorder_descendants_from(t.get_root(), NodeId::from(1u32))
                .unwrap()
                .map(|node| t.some_label_of(node).unwrap())
                .collect_vec(),
            vec!["3", "4", "2", "0"]
        );

        assert_eq!(
            t.postorder_descendants_from(t.get_root(), NodeId::from(3u32))
                .unwrap()
                .map(|node| t.some_label_of(node).unwrap())
                .collect_vec(),
            vec!["4", "2", "0"]
        );

        assert_eq!(
            t.postorder_descendants_from(t.get_root(), NodeId::from(4u32))
                .unwrap()
                .map(|node| t.some_label_of(node).unwrap())
                .collect_vec(),
            vec!["2", "0"]
        );

        assert_eq!(
            t.postorder_descendants_from(t.get_root(), NodeId::from(2u32))
                .unwrap()
                .map(|node| t.some_label_of(node).unwrap())
                .collect_vec(),
            vec!["0"]
        );

        assert!(t
            .postorder_descendants_from(t.get_root(), NodeId::from(0u32))
            .is_none());

        let t = test_tree_1();

        let postorder_desc = t.postorder_descendants(t.get_root()).collect_vec();

        for (i, &child) in postorder_desc.iter().enumerate() {
            for parent in t.strict_ancestors(child) {
                let parent_index = postorder_desc.iter().position(|n| *n == parent).unwrap();

                dbg!((parent, child));

                if i < postorder_desc.len()-1 {
                    assert_eq!(
                        &t.postorder_descendants_from(parent, child).unwrap().collect_vec(),
                        &postorder_desc[i+1..=parent_index],
                    );
                } else {
                    assert!(
                        t.postorder_descendants_from(parent, child).is_none()
                    );
                }
            }
        }
    }

    #[test]
    fn euler_tour() {
        let t = test_tree_0();

        assert_eq!(
            t.euler_tour(t.get_root())
                .map(|(node, depth)| format!("{}-{}", t.some_label_of(node).unwrap(), depth.0))
                .collect_vec(),
            vec!["0-0", "1-1", "0-0", "2-1", "3-2", "2-1", "4-2", "2-1", "0-0"]
        );
    }

    #[test]
    fn lca() {
        let t = test_tree_0();

        let lca_cache = LcaCache::compute(&t);

        let label_node = |label| {
            t.some_node_with_label(label).unwrap()
        };

        let node0 = label_node("0");
        let node1 = label_node("1");
        let node2 = label_node("2");
        let node3 = label_node("3");
        let node4 = label_node("4");

        assert_eq!(lca_cache.lca_node(node1, node2).ok(), Some(node0));
        assert_eq!(lca_cache.lca_node(node0, node1).ok(), Some(node0));
        assert_eq!(lca_cache.lca_node(node1, node3).ok(), Some(node0));
        assert_eq!(lca_cache.lca_node(node3, node4).ok(), Some(node2));
        assert_eq!(lca_cache.lca_node(node3, node3).ok(), Some(node3));
    }

    #[test]
    fn reachability_cmp() {
        let t = test_tree_0();

        let lca_cache = LcaCache::compute(&t);

        let label_key = |label| {
            lca_cache.node_range(t.some_node_with_label(label).unwrap()).unwrap()
        };

        let node0 = label_key("0");
        let node2 = label_key("2");
        let node3 = label_key("3");
        let node4 = label_key("4");

        assert!(node0 < node2);
        assert!(node0 < node2);
        assert!(node2 < node4);
        assert!(node3 > node2);

        assert!(node3.partial_cmp(&node4).is_none());

        assert_eq!(node3.partial_cmp(&node3), Some(Ordering::Equal));
    }
}
