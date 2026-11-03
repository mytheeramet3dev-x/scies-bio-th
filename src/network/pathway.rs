//! KEGG-style signaling pathway DAG representation and topological sorting.

use crate::network::graph::BioGraph;
use crate::BioError;
use std::collections::{HashMap, VecDeque};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Type of pathway regulatory interaction.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum PathwayInteractionType {
    Activation,
    Inhibition,
    Phosphorylation,
    Dephosphorylation,
    Binding,
}

/// A biological signaling or metabolic pathway structured as a Directed Acyclic Graph (DAG).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PathwayDag {
    /// Pathway HGNC or KEGG identifier (e.g. "hsa04115_p53_signaling").
    pub name: String,
    /// Underlying directed graph.
    pub graph: BioGraph<String, PathwayInteractionType>,
    /// Symbol to node index mapping.
    pub symbol_to_index: HashMap<String, usize>,
}

impl PathwayDag {
    /// Creates a new named pathway DAG.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            graph: BioGraph::new(true), // Strictly directed
            symbol_to_index: HashMap::new(),
        }
    }

    /// Adds a gene or metabolite node if not already present.
    pub fn add_node(&mut self, symbol: &str) -> usize {
        let sym_upper = symbol.to_uppercase();
        if let Some(&idx) = self.symbol_to_index.get(&sym_upper) {
            idx
        } else {
            let idx = self.graph.add_node(sym_upper.clone());
            self.symbol_to_index.insert(sym_upper, idx);
            idx
        }
    }

    /// Adds a directed regulatory step: `source -> target` with interaction type.
    pub fn add_step(
        &mut self,
        source: &str,
        target: &str,
        interaction: PathwayInteractionType,
    ) -> crate::Result<()> {
        let u = self.add_node(source);
        let v = self.add_node(target);
        self.graph.add_edge(u, v, interaction)
    }

    /// Performs Kahn's topological sort, validating that the pathway cascade is acyclic.
    pub fn topological_sort(&self) -> crate::Result<Vec<String>> {
        let n = self.graph.node_count();
        let mut in_degree = vec![0usize; n];

        for adj in &self.graph.adjacency {
            for &(v, _) in adj {
                in_degree[v] += 1;
            }
        }

        let mut queue = VecDeque::new();
        for (i, &deg) in in_degree.iter().enumerate().take(n) {
            if deg == 0 {
                queue.push_back(i);
            }
        }

        let mut sorted = Vec::with_capacity(n);
        while let Some(u) = queue.pop_front() {
            sorted.push(self.graph.nodes[u].clone());
            for &(v, _) in &self.graph.adjacency[u] {
                in_degree[v] -= 1;
                if in_degree[v] == 0 {
                    queue.push_back(v);
                }
            }
        }

        if sorted.len() != n {
            return Err(BioError::InvalidInput("Pathway graph contains a feedback loop/cycle".to_string()));
        }

        Ok(sorted)
    }

    /// Gathers all downstream effectors reachable from a given upstream gene.
    pub fn downstream_effectors(&self, gene: &str) -> Vec<String> {
        let sym_upper = gene.to_uppercase();
        if let Some(&start) = self.symbol_to_index.get(&sym_upper) {
            let reachable = self.graph.bfs(start);
            reachable
                .into_iter()
                .filter(|&u| u != start)
                .map(|u| self.graph.nodes[u].clone())
                .collect()
        } else {
            Vec::new()
        }
    }

    /// Gathers all immediate upstream regulators directly controlling a given target gene.
    pub fn upstream_regulators(&self, gene: &str) -> Vec<String> {
        let sym_upper = gene.to_uppercase();
        let target_idx = match self.symbol_to_index.get(&sym_upper) {
            Some(&idx) => idx,
            None => return Vec::new(),
        };

        let mut upstream = Vec::new();
        for (u, adj) in self.graph.adjacency.iter().enumerate() {
            if adj.iter().any(|&(v, _)| v == target_idx) {
                upstream.push(self.graph.nodes[u].clone());
            }
        }
        upstream
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pathway_dag_topological_order() {
        let mut pathway = PathwayDag::new("DNA_Damage_Response");
        // ATM -> CHEK2 -> TP53 -> CDKN1A (p21)
        pathway.add_step("ATM", "CHEK2", PathwayInteractionType::Phosphorylation).unwrap();
        pathway.add_step("CHEK2", "TP53", PathwayInteractionType::Phosphorylation).unwrap();
        pathway.add_step("TP53", "CDKN1A", PathwayInteractionType::Activation).unwrap();

        let order = pathway.topological_sort().unwrap();
        assert_eq!(order, vec!["ATM", "CHEK2", "TP53", "CDKN1A"]);

        let effectors = pathway.downstream_effectors("ATM");
        assert_eq!(effectors.len(), 3);
        assert!(effectors.contains(&"TP53".to_string()));

        let regulators = pathway.upstream_regulators("TP53");
        assert_eq!(regulators, vec!["CHEK2"]);
    }
}
