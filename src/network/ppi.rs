//! Protein-Protein Interaction (PPI) networks, hub detection, and subnetwork extraction.

use crate::network::graph::BioGraph;
use std::collections::HashMap;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A Protein-Protein Interaction network with interaction confidence edge weights.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PpiNetwork {
    /// Underlying graph representation.
    pub graph: BioGraph<String, f64>,
    /// Index lookup table mapping HGNC symbol to node index.
    pub symbol_to_index: HashMap<String, usize>,
}

impl PpiNetwork {
    /// Creates a new empty PPI network.
    pub fn new() -> Self {
        Self {
            graph: BioGraph::new(false), // PPIs are undirected physical binding events
            symbol_to_index: HashMap::new(),
        }
    }

    /// Adds a protein node if not already present and returns its node index.
    pub fn add_protein(&mut self, symbol: &str) -> usize {
        let sym_upper = symbol.to_uppercase();
        if let Some(&idx) = self.symbol_to_index.get(&sym_upper) {
            idx
        } else {
            let idx = self.graph.add_node(sym_upper.clone());
            self.symbol_to_index.insert(sym_upper, idx);
            idx
        }
    }

    /// Adds an undirected physical interaction with a confidence score ($0.0$ to $1.0$).
    pub fn add_interaction(&mut self, p1: &str, p2: &str, confidence: f64) -> crate::Result<()> {
        let u = self.add_protein(p1);
        let v = self.add_protein(p2);
        self.graph.add_edge(u, v, confidence.clamp(0.0, 1.0))
    }

    /// Retrieves interaction confidence between two proteins if an interaction exists.
    pub fn interaction_confidence(&self, p1: &str, p2: &str) -> Option<f64> {
        let u = *self.symbol_to_index.get(&p1.to_uppercase())?;
        let v = *self.symbol_to_index.get(&p2.to_uppercase())?;

        for &(target, conf) in &self.graph.adjacency[u] {
            if target == v {
                return Some(conf);
            }
        }
        None
    }

    /// Identifies network hub proteins possessing degree $\ge$ `min_degree`.
    pub fn find_hubs(&self, min_degree: usize) -> Vec<(String, usize)> {
        let mut hubs = Vec::new();
        for (i, symbol) in self.graph.nodes.iter().enumerate() {
            let deg = self.graph.degree(i);
            if deg >= min_degree {
                hubs.push((symbol.clone(), deg));
            }
        }
        hubs.sort_by_key(|a| std::cmp::Reverse(a.1));
        hubs
    }

    /// Extracts an induced subnetwork containing only the specified query proteins.
    pub fn extract_subnetwork(&self, protein_symbols: &[&str]) -> PpiNetwork {
        let mut sub = PpiNetwork::new();
        let query_set: std::collections::HashSet<String> =
            protein_symbols.iter().map(|s| s.to_uppercase()).collect();

        // Add nodes
        for sym in &query_set {
            if self.symbol_to_index.contains_key(sym) {
                sub.add_protein(sym);
            }
        }

        // Add edges between members of the query set
        for sym1 in &query_set {
            if let Some(&u) = self.symbol_to_index.get(sym1) {
                for &(v, conf) in &self.graph.adjacency[u] {
                    if u < v {
                        let sym2 = &self.graph.nodes[v];
                        if query_set.contains(sym2) {
                            let _ = sub.add_interaction(sym1, sym2, conf);
                        }
                    }
                }
            }
        }

        sub
    }
}

impl Default for PpiNetwork {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ppi_hubs_and_subnetwork() {
        let mut ppi = PpiNetwork::new();
        ppi.add_interaction("TP53", "MDM2", 0.95).unwrap();
        ppi.add_interaction("TP53", "ATM", 0.90).unwrap();
        ppi.add_interaction("TP53", "CHEK2", 0.85).unwrap();
        ppi.add_interaction("ATM", "CHEK2", 0.88).unwrap();

        let hubs = ppi.find_hubs(3);
        assert_eq!(hubs.len(), 1);
        assert_eq!(hubs[0].0, "TP53");
        assert_eq!(hubs[0].1, 3);

        let sub = ppi.extract_subnetwork(&["TP53", "MDM2"]);
        assert_eq!(sub.graph.node_count(), 2);
        assert_eq!(sub.graph.edge_count(), 1);
        assert_eq!(sub.interaction_confidence("TP53", "MDM2"), Some(0.95));
    }
}
