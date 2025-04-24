use std::{
    cmp::Ordering,
    ops::{self, Range},
};

use crate::taxonomy::{
    tree::{node_id::NodeIdMap, walk::RootedTreeWalk},
    NodeId,
};

use super::{
    walk::{LcaCache, LcaMapError},
    RootedTree,
};

#[derive(Clone)]
pub struct PostorderAnnPool<Ann> {
    /// NodeId -> postorder index
    postorder_index: NodeIdMap<u32>,
    /// indexed by postorder traversal
    postorder: Vec<(NodeId, Ann)>,
    /// length: postoder.len()-1, as root has no parent
    parent_index: Vec<u32>,
}

#[derive(Copy, Clone)]
pub struct PostorderAnnSlice<'a, Ann> {
    /// NodeId -> postorder index
    postorder_index: &'a NodeIdMap<u32>,
    /// indexed by postorder traversal
    postorder: &'a [(NodeId, Ann)],
    /// length: postoder.len()-1, as root has no parent
    /// values: shifted by first_leaf_index
    parent_index: &'a [u32],
    /// index corresponding to the original pool this is taken from
    first_leaf_index: usize,
}

pub struct PostorderAnnSliceMut<'a, Ann> {
    /// NodeId -> postorder index
    postorder_index: &'a NodeIdMap<u32>,
    /// indexed by postorder traversal
    postorder: &'a mut [(NodeId, Ann)],
    /// length: postoder.len()-1, as root has no parent
    /// values: shifted by first_leaf_index
    parent_index: &'a [u32],
    /// index corresponding to the original pool this is taken from
    first_leaf_index: usize,
}

impl<Ann: Default> PostorderAnnPool<Ann> {
    pub fn new<Tree: RootedTree>(tree: &Tree, lca: NodeId) -> Self {
        Self::new_with(tree, lca, |_, _| Ann::default())
    }

    pub fn expand_for<'a, Tree: RootedTree>(
        &'a mut self,
        tree: &Tree,
        lca_cache: &LcaCache,
        new_node: NodeId,
    ) -> Result<PostorderAnnSliceMut<'a, Ann>, LcaMapError> {
        self.expand_for_with(tree, lca_cache, new_node, |_, _| Ann::default())
    }

    pub fn new_or_expand_for<'a, Tree: RootedTree>(
        pool: &'a mut Option<Self>,
        tree: &Tree,
        lca_cache: &LcaCache,
        new_node: NodeId,
    ) -> Result<PostorderAnnSliceMut<'a, Ann>, LcaMapError> {
        Self::new_or_expand_for_with(pool, tree, lca_cache, new_node, |_, _| Ann::default())
    }
}

impl<Ann> PostorderAnnPool<Ann> {
    pub fn new_with<Tree: RootedTree, F: FnMut(u32, NodeId) -> Ann>(
        tree: &Tree,
        lca: NodeId,
        mut new_ann: F,
    ) -> Self {
        assert!(tree.node_count() < u32::MAX as usize);

        let mut postorder_index = NodeIdMap::new();

        let postorder: Vec<(NodeId, Ann)> = tree
            .postorder_descendants(lca)
            .enumerate()
            .map(|(i, node)| {
                let i = i as u32;
                postorder_index.insert(node, i);
                (node, new_ann(i, node))
            })
            .collect();

        let postorder_excluding_root = &postorder[..postorder.len() - 1];

        let parent_index = postorder_excluding_root
            .iter()
            .map(|(node, _)| {
                let parent = tree.find_parent(*node).unwrap();
                *postorder_index.get(parent).unwrap()
            })
            .collect();

        Self {
            postorder_index,
            postorder,
            parent_index,
        }
    }

    pub fn new_or_expand_for_with<'a, Tree: RootedTree>(
        pool: &'a mut Option<Self>,
        tree: &Tree,
        lca_cache: &LcaCache,
        node: NodeId,
        new_ann: impl FnMut(u32, NodeId) -> Ann,
    ) -> Result<PostorderAnnSliceMut<'a, Ann>, LcaMapError> {
        if let Some(pool) = pool {
            pool.expand_for_with(tree, lca_cache, node, new_ann)
        } else {
            *pool = Some(Self::new_with(tree, node, new_ann));
            Ok(pool.as_mut().unwrap().slice_root_mut())
        }
    }

    pub fn expand_for_with<'a, Tree: RootedTree, F: FnMut(u32, NodeId) -> Ann>(
        &'a mut self,
        tree: &Tree,
        lca_cache: &LcaCache,
        new_node: NodeId,
        mut new_ann: F,
    ) -> Result<PostorderAnnSliceMut<'a, Ann>, LcaMapError> {
        let previous_lca = self.lca_id();
        let previous_lca_range = lca_cache.node_range(previous_lca)?;

        let new_node_range = lca_cache.node_range(new_node)?;

        let new_lca = match previous_lca_range.partial_cmp(&new_node_range) {
            Some(Ordering::Less | Ordering::Equal) => {
                // current LCA is already an ancestor of node
                return Ok(self
                    .slice_clade_mut(tree, new_node)
                    .expect("new_node should be in the postorder annotation pool"));
            }
            Some(Ordering::Greater) => {
                // node is an ancestor of the current LCA
                new_node
            }
            None => {
                // current LCA and node are unrelated
                // compute their LCA
                let new_lca_range = lca_cache.lca_index(previous_lca_range, new_node_range);
                lca_cache.node_id(new_lca_range)?
            }
        };

        let (mut postorder, mut postorder_index) =
            self.build_postorder_prefix(tree, new_lca, &mut new_ann);

        let prefix_len = postorder.len() as u32;

        postorder.extend(std::mem::take(&mut self.postorder));

        postorder_index.extend(
            std::mem::take(&mut self.postorder_index)
                .iter()
                .map(|(k, parent_i)| (k, *parent_i + prefix_len)),
        );

        let offset = postorder.len() as u32;

        postorder.extend(
            tree.postorder_descendants_from(new_lca, previous_lca)
                .unwrap()
                .enumerate()
                .map(|(i, node)| {
                    let i = offset + i as u32;
                    postorder_index.insert(node, i);
                    (node, new_ann(i, node))
                }),
        );

        let parent_index = Self::expanded_parent_index(
            tree,
            previous_lca,
            prefix_len,
            &postorder,
            &postorder_index,
            &std::mem::take(&mut self.parent_index),
        );

        self.postorder = postorder;
        self.postorder_index = postorder_index;
        self.parent_index = parent_index;

        Ok(self.slice_root_mut())
    }

    fn build_postorder_prefix<Tree: RootedTree, F: FnMut(u32, NodeId) -> Ann>(
        &self,
        tree: &Tree,
        new_lca: NodeId,
        new_ann: &mut F,
    ) -> (Vec<(NodeId, Ann)>, NodeIdMap<u32>) {
        let previous_leftmost_leaf = self.leftmost_leaf();

        let mut postorder_index = NodeIdMap::with_capacity(self.len());

        let postorder: Vec<(NodeId, Ann)> = tree
            .postorder_descendants(new_lca)
            .take_while(|node| *node != previous_leftmost_leaf)
            .enumerate()
            .map(|(i, node)| {
                let i = i as u32;
                postorder_index.insert(node, i);
                (node, new_ann(i, node))
            })
            .collect();

        (postorder, postorder_index)
    }

    fn expanded_parent_index<Tree: RootedTree>(
        tree: &Tree,
        previous_lca: NodeId,
        prefix_len: u32,
        postorder: &[(NodeId, Ann)],
        postorder_index: &NodeIdMap<u32>,
        previous_parent_index: &[u32],
    ) -> Vec<u32> {
        let mut parent_index = Vec::with_capacity(postorder.len() - 1);

        parent_index.extend(postorder[..prefix_len as usize].iter().map(|(node, _)| {
            let parent = tree.find_parent(*node).unwrap();
            *postorder_index.get(parent).unwrap()
        }));

        let previous_len = previous_parent_index.len() + 1;
        let suffix_offset = prefix_len as usize + previous_len;

        parent_index.extend(previous_parent_index.iter().map(|i| prefix_len + i));

        parent_index.push(
            *postorder_index
                .get(tree.find_parent(previous_lca).unwrap())
                .unwrap(),
        );

        let suffix_excluding_lca = &postorder[suffix_offset..postorder.len()-1];

        parent_index.extend(suffix_excluding_lca.iter().map(|(node, _)| {
            let parent = tree.find_parent(*node).unwrap();
            *postorder_index.get(parent).unwrap()
        }));

        parent_index
    }
}

impl<Ann> PostorderAnnPool<Ann> {
    #[expect(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.postorder.len()
    }

    pub fn leftmost_leaf(&self) -> NodeId {
        self.postorder.first().unwrap().0
    }

    pub fn lca_id(&self) -> NodeId {
        self.postorder.last().unwrap().0
    }

    pub fn slice_root(&self) -> PostorderAnnSlice<Ann> {
        PostorderAnnSlice {
            postorder_index: &self.postorder_index,
            postorder: &self.postorder,
            parent_index: &self.parent_index,
            first_leaf_index: 0,
        }
    }

    pub fn slice_root_mut(&mut self) -> PostorderAnnSliceMut<Ann> {
        PostorderAnnSliceMut {
            postorder_index: &self.postorder_index,
            postorder: &mut self.postorder,
            parent_index: &self.parent_index,
            first_leaf_index: 0,
        }
    }

    fn slice_indices<Tree: RootedTree>(&self, tree: &Tree, lca: NodeId) -> Option<(usize, usize)> {
        let leftmost_leaf = tree.leftmost_descendant_leaf(lca);

        let leftmost_leaf_index = *self.postorder_index.get(leftmost_leaf)?;

        let lca_index = *self.postorder_index.get(lca)?;

        Some((leftmost_leaf_index as usize, lca_index as usize))
    }

    pub fn slice_clade<Tree: RootedTree>(
        &self,
        tree: &Tree,
        lca: NodeId,
    ) -> Option<PostorderAnnSlice<Ann>> {
        let (leftmost_leaf_index, lca_index) = self.slice_indices(tree, lca)?;

        Some(PostorderAnnSlice {
            postorder_index: &self.postorder_index,
            postorder: &self.postorder[leftmost_leaf_index..=lca_index],
            parent_index: &self.parent_index[leftmost_leaf_index..lca_index],
            first_leaf_index: leftmost_leaf_index,
        })
    }

    pub fn slice_clade_mut<Tree: RootedTree>(
        &mut self,
        tree: &Tree,
        lca: NodeId,
    ) -> Option<PostorderAnnSliceMut<Ann>> {
        let (leftmost_leaf_index, lca_index) = self.slice_indices(tree, lca)?;

        Some(PostorderAnnSliceMut {
            postorder_index: &self.postorder_index,
            postorder: &mut self.postorder[leftmost_leaf_index..=lca_index],
            parent_index: &self.parent_index[leftmost_leaf_index..lca_index],
            first_leaf_index: leftmost_leaf_index,
        })
    }
}

impl<'a, Ann> ops::Index<usize> for PostorderAnnSlice<'a, Ann> {
    type Output = (NodeId, Ann);

    fn index(&self, index: usize) -> &Self::Output {
        &self.postorder[index]
    }
}

impl<'a, Ann> ops::Index<usize> for PostorderAnnSliceMut<'a, Ann> {
    type Output = (NodeId, Ann);

    fn index(&self, index: usize) -> &Self::Output {
        &self.postorder[index]
    }
}

impl<'a, Ann> ops::IndexMut<usize> for PostorderAnnSliceMut<'a, Ann> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.postorder[index]
    }
}

impl<'a, Ann> PostorderAnnSlice<'a, Ann> {
    pub fn get(&self, node: NodeId) -> Option<&Ann> {
        Some(&self.postorder[self.local_index_of(node)?].1)
    }

    pub fn indices(&self) -> Range<usize> {
        0..self.postorder.len()
    }

    pub fn iter(&self) -> impl DoubleEndedIterator<Item = (NodeId, &Ann)> + ExactSizeIterator {
        self.postorder.iter().map(|(id, ann)| (*id, ann))
    }

    pub fn parent_index(&self, index: usize) -> Option<usize> {
        Some(self.local_index(*self.parent_index.get(index)? as usize))
    }

    pub fn global_lca_index(&self) -> usize {
        self.global_index(self.postorder.len() - 1)
    }

    pub fn lca_id(&self) -> NodeId {
        self.postorder.last().unwrap().0
    }

    pub fn lca(&self) -> &Ann {
        &self.postorder.last().unwrap().1
    }

    fn local_index_of(&self, node: NodeId) -> Option<usize> {
        Some(self.local_index(*self.postorder_index.get(node)? as usize))
    }

    fn local_index(&self, global_index: usize) -> usize {
        global_index - self.first_leaf_index
    }

    fn global_index(&self, local_index: usize) -> usize {
        local_index + self.first_leaf_index
    }
}

impl<'a, Ann> PostorderAnnSliceMut<'a, Ann> {
    pub fn get(&'_ self, node: NodeId) -> Option<&Ann> {
        Some(&self.postorder[self.local_index_of(node)?].1)
    }

    pub fn get_mut(&mut self, node: NodeId) -> Option<&mut Ann> {
        Some(&mut self.postorder[self.local_index_of(node)?].1)
    }

    pub fn indices(&self) -> Range<usize> {
        0..self.postorder.len()
    }

    pub fn iter(&self) -> impl DoubleEndedIterator<Item = (NodeId, &Ann)> + ExactSizeIterator {
        self.postorder.iter().map(|(id, ann)| (*id, ann))
    }

    pub fn iter_mut(&mut self) -> impl DoubleEndedIterator<Item = (NodeId, &mut Ann)> + ExactSizeIterator {
        self.postorder.iter_mut().map(|(id, ann)| (*id, ann))
    }

    pub fn parent_index(&self, index: usize) -> Option<usize> {
        Some(self.local_index(*self.parent_index.get(index)? as usize))
    }

    pub fn global_lca_index(&self) -> usize {
        self.global_index(self.postorder.len() - 1)
    }

    pub fn lca(&self) -> &Ann {
        &self.postorder.last().unwrap().1
    }

    pub fn lca_mut(&mut self) -> &mut Ann {
        &mut self.postorder.last_mut().unwrap().1
    }

    fn local_index_of(&self, node: NodeId) -> Option<usize> {
        Some(self.local_index(*self.postorder_index.get(node)? as usize))
    }

    fn local_index(&self, global_index: usize) -> usize {
        global_index - self.first_leaf_index
    }

    fn global_index(&self, local_index: usize) -> usize {
        local_index + self.first_leaf_index
    }
}

#[cfg(test)]
mod tests {
    use crate::taxonomy::{formats::newick, tree::walk::LcaCache, NodeId};

    use super::PostorderAnnPool;

    #[test]
    fn test_expand() {
        let t = newick::tests::sample_taxonomy();
        let lca_cache = LcaCache::compute(&t);

        let mut pool = PostorderAnnPool::<()>::new(&t, NodeId::from(13u32));
        pool.expand_for(&t, &lca_cache, NodeId::from(3u32)).ok().unwrap();

        let expected = PostorderAnnPool::<()>::new(&t, NodeId::from(3u32));

        assert_eq!(&pool.postorder_index, &expected.postorder_index);

        assert_eq!(&pool.postorder, &expected.postorder);

        assert_eq!(&pool.parent_index, &expected.parent_index);
    }
}
