use std::{fs::File, path::Path};

use itertools::Itertools;
use newick_rs::SimpleTree;
use serde::{Deserialize, Serialize};

use crate::taxonomy::{
    labels::vec::{RootedVecTreeLabels, RootedVecTreeRanks},
    tree::{
        newick::{read_newick_simple_tree, NewickLoadError},
        vec::RootedVecTree,
        RootedTree,
    },
    LabeledTaxonomy, NodeId, Taxonomy,
};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NewickTaxonomy {
    pub tree: RootedVecTree,
    pub labels: RootedVecTreeLabels,
    pub ranks: RootedVecTreeRanks,
}

impl NewickTaxonomy {
    pub const FORMAT_VERSION: u32 = 1;

    pub fn from_simple_tree(value: SimpleTree, ranks: Vec<String>) -> Self {
        let (tree, labels, depths) = RootedVecTree::from_simple_tree(value);

        let ranks = RootedVecTreeRanks::new(depths, ranks);

        NewickTaxonomy {
            tree,
            labels,
            ranks,
        }
    }

    pub fn load_newick(
        path: impl AsRef<Path>,
        ranks: Vec<String>,
    ) -> Result<Self, NewickLoadError> {
        let simple_tree = read_newick_simple_tree(File::open(path)?)?;

        Ok(Self::from_simple_tree(simple_tree, ranks))
    }

    pub fn has_node(&self, node: u32) -> bool {
        self.tree.has_node(node)
    }

    pub fn nodes(&self) -> impl Iterator<Item = NodeId> {
        self.tree.nodes()
    }

    pub fn relabel(&mut self, f: impl Fn(NodeId, &str, &Self) -> String) {
        self.labels = self
            .labels
            .iter()
            .map(|(node, label)| (node, f(node, label, self)))
            .collect();
    }
}

impl RootedTree for NewickTaxonomy {
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

    type Children<'a> = <RootedVecTree as RootedTree>::Children<'a>;

    fn iter_children(&self, node: NodeId) -> Self::Children<'_> {
        self.tree.iter_children(node)
    }
}

impl Taxonomy for NewickTaxonomy {
    fn has_uniform_depths(&self) -> Option<usize> {
        let mut depths = self
            .nodes()
            .filter(|&node| self.is_leaf(node))
            .filter_map(|node| self.ranks.depths().get(node))
            .peekable();

        let &depth = depths.peek().expect("Tree has no leaves");

        if depths.all_equal() {
            Some(depth)
        } else {
            None
        }
    }

    type RankSym = usize;

    fn rank_sym_str(&self, rank_sym: Self::RankSym) -> Option<&str> {
        self.ranks.rank_sym_str(rank_sym)
    }

    fn lookup_rank_sym(&self, rank: &str) -> Option<Self::RankSym> {
        self.ranks.lookup_rank_sym(rank)
    }

    fn find_rank(&self, node: NodeId) -> Option<Self::RankSym> {
        self.ranks.find_rank(node)
    }

    type NodeRanks<'a> = impl Iterator<Item = (NodeId, Self::RankSym)>;

    fn node_ranks(&self) -> Self::NodeRanks<'_> {
        self.ranks.node_ranks()
    }
}

impl LabeledTaxonomy for NewickTaxonomy {
    type Labels<'a> = std::option::IntoIter<&'a str>;

    fn labels_of(&self, node: NodeId) -> Self::Labels<'_> {
        self.labels.get(node).into_iter()
    }

    type NodesWithLabel<'a> = impl Iterator<Item = NodeId> + use<'a>;

    fn nodes_with_label<'a>(&'a self, label: &'a str) -> Self::NodesWithLabel<'a> {
        self.labels.nodes_with_label(label)
    }
}

impl From<NewickTaxonomy> for SimpleTree {
    fn from(taxonomy: NewickTaxonomy) -> Self {
        taxonomy.tree.into_simple_tree(taxonomy.labels)
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use itertools::Itertools;

    use crate::taxonomy::{
        formats::newick::NewickTaxonomy, tree::newick::newick_literal, LabeledTaxonomy, Taxonomy
    };

    pub fn sample_ranks() -> Vec<String> {
        vec![
            "root".to_owned(),
            "superkingdom".to_owned(),
            "phylum".to_owned(),
            "class".to_owned(),
            "order".to_owned(),
            "family".to_owned(),
            "genus".to_owned(),
            "species".to_owned(),
        ]
    }

    fn rank_labeled_taxonomy(mut tax: NewickTaxonomy) -> NewickTaxonomy {
        tax.relabel(|node, label, tax| {
            let prefix = match tax.find_rank_str(node).unwrap_or("") {
                "root" => "R",
                "" => "?",
                "superkingdom" => "S",
                rank => &rank[0..1],
            };

            format!("{prefix}{label}")
        });

        tax
    }

    pub fn sample_taxonomy() -> NewickTaxonomy {
        //                              R0                           root
        //                              S1                           superkingdom
        //                              p2                           phylum
        //                              /\
        //                           /      \
        //                        /            \
        //                     /                  \
        //                  /                        \
        //               c3                           c20            class
        //                |                            /\
        //                |                          /    \
        //                |                        /        \
        //               o4                     o21         o29      order
        //               / \                     |           |
        //             /     \                   |           |
        //           /         \                 |           |
        //         /             \               |           |
        //       f5             f13             f22         f30      family
        //      /  \            /  \            /  \         |
        //     /    \          /    \          /    \        |
        //   g6     g10     g14     g17     g23     g26     g31      genus
        //   /|\    /  \    /  \    /  \    /  \    /  \    /  \
        // s7s8s9 s11 s12 s15 s16 s18 s19 s24 s25 s27 s28 s32 s33    species

        let genera = [
            newick_literal! { [7, 8, 9] 6 },
            newick_literal! { [11, 12]  10 },
            newick_literal! { [15, 16]  14 },
            newick_literal! { [18, 19]  17 },
            newick_literal! { [24, 25]  23 },
            newick_literal! { [27, 28]  26 },
            newick_literal! { [32, 33]  31 },
        ];

        let families = [
            newick_literal! { [ (genera[0].clone()), (genera[1].clone()) ]  5 },
            newick_literal! { [ (genera[2].clone()), (genera[3].clone()) ] 13 },
            newick_literal! { [ (genera[4].clone()), (genera[5].clone()) ] 22 },
            newick_literal! { [ (genera[6].clone()) ] 30 },
        ];

        let orders = [
            newick_literal! { [ (families[0].clone()), (families[1].clone()) ] 4 },
            newick_literal! { [ (families[2].clone()) ] 21 },
            newick_literal! { [ (families[3].clone()) ] 29 },
        ];

        let classes = [
            newick_literal! { [ (orders[0].clone()) ] 3 },
            newick_literal! { [ (orders[1].clone()), (orders[2].clone()) ] 20 },
        ];

        let phylums = [newick_literal! { [ (classes[0].clone()), (classes[1].clone()) ] 2 }];

        let superkingdoms = [newick_literal! { [ (phylums[0].clone()) ] 1 }];

        let root = newick_literal! { [ (superkingdoms[0].clone()) ] 0 };

        rank_labeled_taxonomy(NewickTaxonomy::from_simple_tree(root, sample_ranks()))
    }

    pub fn verify_node_ids_equal_labels(tax: &NewickTaxonomy) -> bool {
        for node in tax.nodes() {
            let labels = tax
                .labels_of(node)
                .map(|label| {
                    if label.chars().next().is_some_and(|x| !x.is_ascii_digit()) {
                        &label[1..]
                    } else {
                        label
                    }
                })
                .collect_vec();

            if labels != [&node.0.to_string()] {
                return false
            }
        }

        true
    }

    #[test]
    pub fn verify_sample_taxonomy() {
        assert!(verify_node_ids_equal_labels(&sample_taxonomy()));
    }

    pub fn ambiguous_taxonomy() -> NewickTaxonomy {
        //                              0                           root
        //                              1                           superkingdom
        //                              2                           phylum
        //                             /\
        //                          /      \
        //                       /            \
        //                    /                  \
        //                 /                        \
        //               3                            20            class
        //               |                            /\
        //               |                          /    \
        //               |                        /        \
        //               4                      21          29      order
        //              / \                     |           |
        //            /     \                   |           |
        //          /         \                 |           |
        //        /             \               |           |
        //       x              13              22          30      family
        //     /  \            /  \            /  \         |
        //    /    \          /    \          /    \        |
        //   x      10      14      17      23      26      31      genus
        //  /|\    /  \    /  \    /  \    /  \    /  \    /  \
        // 7 8 9  11  12  15  16  18  19  24  25  27  28  32  33    species

        let genera = [
            newick_literal! { [7, 8, 9] "x" },
            newick_literal! { [11, 12]  10 },
            newick_literal! { [15, 16]  14 },
            newick_literal! { [18, 19]  17 },
            newick_literal! { [24, 25]  23 },
            newick_literal! { [27, 28]  26 },
            newick_literal! { [32, 33]  31 },
        ];

        let families = [
            newick_literal! { [ (genera[0].clone()), (genera[1].clone()) ] "x" },
            newick_literal! { [ (genera[2].clone()), (genera[3].clone()) ] 13 },
            newick_literal! { [ (genera[4].clone()), (genera[5].clone()) ] 22 },
            newick_literal! { [ (genera[6].clone()) ] 30 },
        ];

        let orders = [
            newick_literal! { [ (families[0].clone()), (families[1].clone()) ] 4 },
            newick_literal! { [ (families[2].clone()) ] 21 },
            newick_literal! { [ (families[3].clone()) ] 29 },
        ];

        let classes = [
            newick_literal! { [ (orders[0].clone()) ] 3 },
            newick_literal! { [ (orders[1].clone()), (orders[2].clone()) ] 20 },
        ];

        let phylums = [newick_literal! { [ (classes[0].clone()), (classes[1].clone()) ] 2 }];

        let superkingdoms = [newick_literal! { [ (phylums[0].clone()) ] 1 }];

        let root = newick_literal! { [ (superkingdoms[0].clone()) ] 0 };

        NewickTaxonomy::from_simple_tree(root, sample_ranks())
    }
}
