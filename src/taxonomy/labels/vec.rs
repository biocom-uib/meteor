use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::taxonomy::{tree::vec::RootedVecTreeDepths, NodeId};


#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct RootedVecTreeLabels {
    labels: Vec<String>,
    label_lookup: HashMap<String, Vec<NodeId>>,
}

impl RootedVecTreeLabels {
    pub fn get(&self, node: NodeId) -> Option<&str> {
        self.labels.get(node.0 as usize).map(|label| label.as_str())
    }

    pub fn get_mut(&mut self, node: NodeId) -> Option<&mut String> {
        self.labels.get_mut(node.0 as usize)
    }

    pub fn nodes_with_label(&self, label: &str) -> impl Iterator<Item = NodeId> + '_ {
        self.label_lookup.get(label).into_iter().flatten().copied()
    }

    pub fn push(&mut self, label: String) -> NodeId {
        let node_id = NodeId(self.labels.len() as u32);

        self.labels.push(label.clone());

        self.label_lookup
            .entry(label)
            .or_insert_with(Vec::new)
            .push(node_id);

        node_id
    }

    pub fn relabel(&self, f: impl Fn(&str, NodeId) -> String) -> Self {
        let mut new_labels = Vec::with_capacity(self.labels.len());

        let mut new_label_lookup = HashMap::new();

        for (node, label) in self.iter() {
            let new_label = f(label, node);

            new_labels.push(new_label.clone());

            new_label_lookup
                .entry(new_label)
                .or_insert_with(Vec::new)
                .push(node);
        }

        RootedVecTreeLabels {
            labels: new_labels,
            label_lookup: new_label_lookup,
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = (NodeId, &'_ str)> + '_ {
        self.labels
            .iter()
            .enumerate()
            .map(|(i, label)| (NodeId::from(i as u32), label.as_str()))
    }
}

impl IntoIterator for RootedVecTreeLabels {
    type Item = (NodeId, String);

    type IntoIter = impl Iterator<Item = (NodeId, String)>;

    fn into_iter(self) -> Self::IntoIter {
        self.labels
            .into_iter()
            .enumerate()
            .map(|(i, label)| (NodeId::from(i as u32), label))
    }
}

impl FromIterator<(NodeId, String)> for RootedVecTreeLabels {
    fn from_iter<T: IntoIterator<Item = (NodeId, String)>>(iter: T) -> Self {
        let iter = iter.into_iter();

        let mut labels = Vec::with_capacity(iter.size_hint().0);
        let mut label_lookup = HashMap::new();

        for (node, label) in iter {
            labels.push(label.clone());

            label_lookup
                .entry(label)
                .or_insert_with(Vec::new)
                .push(node);
        }

        RootedVecTreeLabels {
            labels,
            label_lookup,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RootedVecTreeRanks {
    depths: RootedVecTreeDepths,
    rank_storage: Vec<String>,
}

impl RootedVecTreeRanks {
    pub fn new(depths: RootedVecTreeDepths, rank_storage: Vec<String>) -> Self {
        assert!(depths.iter().all(|(_, depth)| depth < rank_storage.len()));

        Self { depths, rank_storage }
    }
}

impl RootedVecTreeRanks {
    pub fn depths(&self) -> &RootedVecTreeDepths {
        &self.depths
    }

    pub fn rank_sym_str(&self, rank_sym: usize) -> Option<&str> {
        self.rank_storage.get(rank_sym).map(|s| s.as_str())
    }

    pub fn lookup_rank_sym(&self, rank: &str) -> Option<usize> {
        self.rank_storage.iter().position(|r| r == rank)
    }

    pub fn find_rank(&self, node: NodeId) -> Option<usize> {
        self.depths.get(node)
    }

    pub fn node_ranks(&self) -> impl Iterator<Item = (NodeId, usize)> + '_ {
        self.depths.iter()
    }
}
