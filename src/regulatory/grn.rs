//! Gene Regulatory Network (`regulatory::grn` & `network::grn`) DAG representation,
//! cycle detection, topological ordering, and single-pass regulatory signal propagation.

use crate::BioError;
use std::collections::{HashMap, HashSet, VecDeque};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Type of transcriptional regulatory interaction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum RegulationType {
    /// Transcriptional activation ($+$).
    Activation,
    /// Transcriptional repression ($-$).
    Repression,
}

/// Directed edge from a transcription factor (`tf`) to a target gene (`target`).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct RegulatoryEdge {
    /// Upstream regulator / transcription factor symbol.
    pub tf: String,
    /// Downstream regulated target gene symbol.
    pub target: String,
    /// Activation or repression mode.
    pub interaction: RegulationType,
    /// Coupling weight in $(0.0, 1.0]$.
    pub weight: f64,
}

/// Directed Acyclic Gene Regulatory Network (GRN).
#[derive(Debug, Clone, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct GeneRegulatoryNetwork {
    nodes: HashSet<String>,
    edges: Vec<RegulatoryEdge>,
}

impl GeneRegulatoryNetwork {
    /// Creates an empty Gene Regulatory Network.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a directed regulatory interaction `tf -> target`.
    pub fn add_regulation(
        &mut self,
        tf: &str,
        target: &str,
        interaction: RegulationType,
        weight: f64,
    ) {
        self.nodes.insert(tf.to_string());
        self.nodes.insert(target.to_string());
        self.edges.push(RegulatoryEdge {
            tf: tf.to_string(),
            target: target.to_string(),
            interaction,
            weight: weight.clamp(0.0, 1.0),
        });
    }

    /// Returns all registered gene/TF nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Returns all regulatory edges.
    pub fn edges(&self) -> &[RegulatoryEdge] {
        &self.edges
    }

    /// Computes a deterministic Kahn topological sort over the GRN.
    ///
    /// Returns `Err(BioError::InvalidInput)` if the network contains a directed feedback cycle.
    pub fn topological_order(&self) -> crate::Result<Vec<String>> {
        let mut in_degree: HashMap<String, usize> =
            self.nodes.iter().map(|n| (n.clone(), 0)).collect();
        let mut adj: HashMap<String, Vec<String>> = HashMap::new();

        for edge in &self.edges {
            adj.entry(edge.tf.clone())
                .or_default()
                .push(edge.target.clone());
            *in_degree.entry(edge.target.clone()).or_insert(0) += 1;
        }

        let mut zero_in: Vec<String> = in_degree
            .iter()
            .filter_map(|(k, &deg)| if deg == 0 { Some(k.clone()) } else { None })
            .collect();
        zero_in.sort();

        let mut queue: VecDeque<String> = zero_in.into_iter().collect();
        let mut order = Vec::with_capacity(self.nodes.len());

        while let Some(u) = queue.pop_front() {
            order.push(u.clone());
            if let Some(neighbors) = adj.get(&u) {
                let mut newly_zero = Vec::new();
                for v in neighbors {
                    if let Some(deg) = in_degree.get_mut(v) {
                        *deg -= 1;
                        if *deg == 0 {
                            newly_zero.push(v.clone());
                        }
                    }
                }
                newly_zero.sort();
                for nz in newly_zero {
                    queue.push_back(nz);
                }
            }
        }

        if order.len() != self.nodes.len() {
            return Err(BioError::InvalidInput(
                "Cycle detected in Gene Regulatory Network DAG".to_string(),
            ));
        }

        Ok(order)
    }

    /// Propagates initial TF expression levels along the topological order of the DAG.
    pub fn propagate_steady_state(
        &self,
        initial_expression: &HashMap<String, f64>,
    ) -> crate::Result<HashMap<String, f64>> {
        let topo = self.topological_order()?;
        let mut state: HashMap<String, f64> = self
            .nodes
            .iter()
            .map(|n| (n.clone(), *initial_expression.get(n).unwrap_or(&1.0)))
            .collect();

        for node in &topo {
            let upstream_val = *state.get(node).unwrap_or(&1.0);
            for edge in self.edges.iter().filter(|e| &e.tf == node) {
                let target_entry = state.entry(edge.target.clone()).or_insert(1.0);
                let delta = match edge.interaction {
                    RegulationType::Activation => edge.weight * (upstream_val - 1.0),
                    RegulationType::Repression => -edge.weight * (upstream_val - 1.0),
                };
                *target_entry = (*target_entry + delta).max(0.0);
            }
        }

        Ok(state)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_grn_topological_sort_and_cycle_rejection() {
        // Validation target: topological sort on valid DAG, and detect + reject directed cycles
        let mut grn = GeneRegulatoryNetwork::new();
        grn.add_regulation("OCT4", "SOX2", RegulationType::Activation, 0.8);
        grn.add_regulation("SOX2", "NANOG", RegulationType::Activation, 0.9);
        grn.add_regulation("NANOG", "CDKN2A", RegulationType::Repression, 0.7);

        let order = grn.topological_order().unwrap();
        assert_eq!(order, vec!["OCT4", "SOX2", "NANOG", "CDKN2A"]);

        // Add cycle: CDKN2A -> OCT4
        grn.add_regulation("CDKN2A", "OCT4", RegulationType::Repression, 0.5);
        assert!(grn.topological_order().is_err());
    }
}
