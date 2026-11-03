//! Generic biological network graph data structures, traversal, and shortest path algorithms.

use crate::BioError;
use std::collections::{BinaryHeap, VecDeque};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Directed or undirected biological network graph.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct BioGraph<N, E> {
    /// Node data items indexed by internal node ID.
    pub nodes: Vec<N>,
    /// Adjacency list: `adjacency[u]` contains `(target_v, edge_data)`.
    pub adjacency: Vec<Vec<(usize, E)>>,
    /// Whether the graph is directed.
    pub is_directed: bool,
}

impl<N, E: Clone> BioGraph<N, E> {
    /// Creates a new empty biological graph.
    pub fn new(is_directed: bool) -> Self {
        Self {
            nodes: Vec::new(),
            adjacency: Vec::new(),
            is_directed,
        }
    }

    /// Adds a new node containing data and returns its unique node index.
    pub fn add_node(&mut self, data: N) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(data);
        self.adjacency.push(Vec::new());
        idx
    }

    /// Adds an edge between source and target nodes with edge payload data.
    pub fn add_edge(&mut self, source: usize, target: usize, data: E) -> crate::Result<()> {
        if source >= self.nodes.len() || target >= self.nodes.len() {
            return Err(BioError::InvalidInput("Node index out of graph bounds".to_string()));
        }

        self.adjacency[source].push((target, data.clone()));
        if !self.is_directed && source != target {
            self.adjacency[target].push((source, data));
        }
        Ok(())
    }

    /// Number of nodes in the graph.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Total number of edges.
    pub fn edge_count(&self) -> usize {
        let total: usize = self.adjacency.iter().map(|v| v.len()).sum();
        if self.is_directed {
            total
        } else {
            total / 2
        }
    }

    /// Returns the immediate neighbor node indices of a given node.
    pub fn neighbors(&self, node: usize) -> Vec<usize> {
        if node < self.adjacency.len() {
            self.adjacency[node].iter().map(|(v, _)| *v).collect()
        } else {
            Vec::new()
        }
    }

    /// Degree (number of connected edges) of a node.
    pub fn degree(&self, node: usize) -> usize {
        if node < self.adjacency.len() {
            self.adjacency[node].len()
        } else {
            0
        }
    }

    /// Breadth-First Search (BFS) traversal starting from a source node.
    pub fn bfs(&self, start: usize) -> Vec<usize> {
        let mut visited = vec![false; self.nodes.len()];
        let mut order = Vec::new();
        let mut queue = VecDeque::new();

        if start >= self.nodes.len() {
            return order;
        }

        visited[start] = true;
        queue.push_back(start);

        while let Some(u) = queue.pop_front() {
            order.push(u);
            for &(v, _) in &self.adjacency[u] {
                if !visited[v] {
                    visited[v] = true;
                    queue.push_back(v);
                }
            }
        }

        order
    }

    /// Computes the shortest path between start and goal using Dijkstra's algorithm.
    pub fn shortest_path<F>(&self, start: usize, goal: usize, mut weight_fn: F) -> Option<(Vec<usize>, f64)>
    where
        F: FnMut(&E) -> f64,
    {
        if start >= self.nodes.len() || goal >= self.nodes.len() {
            return None;
        }

        #[derive(Copy, Clone, PartialEq)]
        struct State {
            cost: f64,
            position: usize,
        }

        impl Eq for State {}

        impl Ord for State {
            fn cmp(&self, other: &Self) -> std::cmp::Ordering {
                other.cost.partial_cmp(&self.cost).unwrap_or(std::cmp::Ordering::Equal)
            }
        }

        impl PartialOrd for State {
            fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
                Some(self.cmp(other))
            }
        }

        let mut dist = vec![f64::INFINITY; self.nodes.len()];
        let mut prev = vec![None; self.nodes.len()];
        let mut heap = BinaryHeap::new();

        dist[start] = 0.0;
        heap.push(State { cost: 0.0, position: start });

        while let Some(State { cost, position }) = heap.pop() {
            if position == goal {
                let mut path = Vec::new();
                let mut curr = Some(goal);
                while let Some(u) = curr {
                    path.push(u);
                    curr = prev[u];
                }
                path.reverse();
                return Some((path, cost));
            }

            if cost > dist[position] {
                continue;
            }

            for (next_node, edge_data) in &self.adjacency[position] {
                let weight = weight_fn(edge_data).max(0.0);
                let next_cost = cost + weight;
                if next_cost < dist[*next_node] {
                    dist[*next_node] = next_cost;
                    prev[*next_node] = Some(position);
                    heap.push(State { cost: next_cost, position: *next_node });
                }
            }
        }

        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_biograph_construction_and_bfs() {
        let mut graph: BioGraph<String, f64> = BioGraph::new(false);
        let a = graph.add_node("GeneA".to_string());
        let b = graph.add_node("GeneB".to_string());
        let c = graph.add_node("GeneC".to_string());

        graph.add_edge(a, b, 1.0).unwrap();
        graph.add_edge(b, c, 2.0).unwrap();

        assert_eq!(graph.node_count(), 3);
        assert_eq!(graph.edge_count(), 2);
        assert_eq!(graph.degree(b), 2);

        let bfs_res = graph.bfs(a);
        assert_eq!(bfs_res, vec![a, b, c]);

        let (path, cost) = graph.shortest_path(a, c, |w| *w).unwrap();
        assert_eq!(path, vec![a, b, c]);
        assert_eq!(cost, 3.0);
    }
}
