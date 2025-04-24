use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::hash::Hash;
use std::iter;

fn dfs_visit<Id, Desc>(
    adjacency: &HashMap<Id, Desc>,
    not_visited: &mut HashSet<Id>,
    current_path: &mut Vec<Id>,
    rev_result: &mut Vec<Id>,
    node: Id,
) -> bool
where
    Id: Copy + Debug + Eq + Hash,
    for<'a> &'a Desc: IntoIterator<Item = &'a Id>,
{
    if !not_visited.contains(&node) {
        return true;
    }

    if !current_path.contains(&node) {
        current_path.push(node);
    } else {
        current_path.push(node);
        return false;
    }

    for &child in adjacency.get(&node).into_iter().flatten() {
        // turbofish because rustc goes crazy
        if !dfs_visit::<Id, Desc>(adjacency, not_visited, current_path, rev_result, child) {
            return false;
        }
    }

    assert_eq!(current_path.pop().unwrap(), node);
    not_visited.remove(&node);
    rev_result.push(node);
    true
}

#[derive(Clone, Debug)]
pub struct Loop<Id>(pub Vec<Id>);

pub fn topsort<Id, Desc>(adjacency: &HashMap<Id, Desc>) -> Result<Vec<Id>, Loop<Id>>
where
    Id: Copy + Debug + Eq + Hash,
    for<'a> &'a Desc: IntoIterator<Item = &'a Id>,
{
    let nodes: HashSet<_> = adjacency
        .iter()
        .flat_map(|(parent, children)| iter::once(parent).chain(children))
        .copied()
        .collect();

    let mut rev_result = Vec::with_capacity(nodes.len());

    let mut not_visited = nodes;
    let mut current_path = Vec::new();

    while let Some(&node) = not_visited.iter().next() {
        // turbofish because rustc goes crazy
        let ok = dfs_visit::<Id, Desc>(
            adjacency,
            &mut not_visited,
            &mut current_path,
            &mut rev_result,
            node,
        );

        if !ok {
            return Err(Loop(current_path));
        }

        assert!(current_path.is_empty());
    }

    rev_result.reverse();
    Ok(rev_result)
}

