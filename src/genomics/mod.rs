//! Genomics, k-mer spectrum analysis, and De Bruijn graph assembly primitives.
//!
//! Provides data structures for counting k-mers and building De Bruijn Graphs (DBG)
//! for sequence assembly and unitig / contig generation.

use crate::error::{BioError, Result};
use crate::seq::Dna;
use std::collections::HashMap;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Directed edge in a De Bruijn graph corresponding to an observed k-mer.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct DeBruijnEdge {
    /// Full k-mer sequence.
    pub kmer: Vec<u8>,
    /// Prefix (k-1)-mer (source node).
    pub from: Vec<u8>,
    /// Suffix (k-1)-mer (target node).
    pub to: Vec<u8>,
    /// Read coverage / abundance count.
    pub coverage: usize,
}

/// A De Bruijn Graph composed of (k-1)-mer nodes and directed k-mer edges.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct DeBruijnGraph {
    /// K-mer size (edges are length k, nodes are length k - 1).
    pub k: usize,
    /// Adjacency map: (k-1)-mer -> list of outgoing (k-1)-mer target nodes.
    pub adj: HashMap<Vec<u8>, Vec<Vec<u8>>>,
    /// In-degree map: (k-1)-mer -> count of incoming edges.
    pub in_degrees: HashMap<Vec<u8>, usize>,
    /// Edge coverage map: k-mer -> count.
    pub edge_coverages: HashMap<Vec<u8>, usize>,
}

impl DeBruijnGraph {
    /// Creates an empty De Bruijn graph with specified k-mer size.
    pub fn new(k: usize) -> Result<Self> {
        if k < 2 {
            return Err(BioError::InvalidInput("K-mer size must be at least 2".into()));
        }
        Ok(Self {
            k,
            adj: HashMap::new(),
            in_degrees: HashMap::new(),
            edge_coverages: HashMap::new(),
        })
    }

    /// Builds a De Bruijn Graph from a single DNA sequence.
    pub fn from_dna(dna: &Dna, k: usize) -> Result<Self> {
        let mut dbg = Self::new(k)?;
        dbg.add_sequence(dna)?;
        Ok(dbg)
    }

    /// Adds all k-mers from a DNA sequence into the graph.
    pub fn add_sequence(&mut self, dna: &Dna) -> Result<()> {
        let bytes = &dna.sequence;
        let len = bytes.len();
        if len < self.k {
            return Ok(());
        }

        let k = self.k;
        for i in 0..=(len - k) {
            let kmer = bytes[i..i + k].to_vec();
            let prefix = bytes[i..i + k - 1].to_vec();
            let suffix = bytes[i + 1..i + k].to_vec();

            // Track edge coverage
            *self.edge_coverages.entry(kmer).or_insert(0) += 1;

            // Ensure nodes exist in in_degrees map
            self.in_degrees.entry(prefix.clone()).or_insert(0);
            *self.in_degrees.entry(suffix.clone()).or_insert(0) += 1;

            // Add directed edge
            self.adj.entry(prefix).or_default().push(suffix);
        }

        Ok(())
    }

    /// Total number of unique (k-1)-mer nodes in the graph.
    pub fn node_count(&self) -> usize {
        self.in_degrees.len()
    }

    /// Total number of unique k-mer directed edges in the graph.
    pub fn edge_count(&self) -> usize {
        self.edge_coverages.len()
    }

    /// Assembles simple non-branching linear paths (unitigs) from the graph.
    pub fn assemble_unitigs(&self) -> Vec<Dna> {
        let mut unitigs = Vec::new();
        let mut visited_edges: HashMap<(Vec<u8>, Vec<u8>), usize> = HashMap::new();

        // 1. First pass: prioritize starting from boundary or branching nodes (in_degree != 1 or out_degree != 1)
        let mut start_edges = Vec::new();
        for (u, targets) in &self.adj {
            let in_deg = self.in_degrees.get(u).copied().unwrap_or(0);
            let out_deg = targets.len();
            if in_deg != 1 || out_deg != 1 {
                for v in targets {
                    start_edges.push((u.clone(), v.clone()));
                }
            }
        }

        for (u, v) in start_edges {
            let edge_key = (u.clone(), v.clone());
            let used = visited_edges.get(&edge_key).copied().unwrap_or(0);
            let total = self
                .adj
                .get(&u)
                .map(|t| t.iter().filter(|&tgt| tgt == &v).count())
                .unwrap_or(0);

            if used < total {
                let mut path = u.clone();
                path.push(*v.last().unwrap());
                *visited_edges.entry(edge_key).or_insert(0) += 1;

                let mut curr_node = v.clone();
                while let Some(next_targets) = self.adj.get(&curr_node) {
                    let in_deg = self.in_degrees.get(&curr_node).copied().unwrap_or(0);
                    let out_deg = next_targets.len();

                    // Non-branching extension: in_degree == 1 and out_degree == 1
                    if in_deg == 1 && out_deg == 1 {
                        let next_node = &next_targets[0];
                        let next_edge = (curr_node.clone(), next_node.clone());
                        let next_used = visited_edges.get(&next_edge).copied().unwrap_or(0);
                        if next_used == 0 {
                            path.push(*next_node.last().unwrap());
                            *visited_edges.entry(next_edge).or_insert(0) += 1;
                            curr_node = next_node.clone();
                        } else {
                            break;
                        }
                    } else {
                        break;
                    }
                }

                if let Ok(dna) = Dna::new(&path) {
                    unitigs.push(dna);
                }
            }
        }

        // 2. Second pass: handle any isolated non-branching cycles
        for (u, targets) in &self.adj {
            for v in targets {
                let edge_key = (u.clone(), v.clone());
                let used = visited_edges.get(&edge_key).copied().unwrap_or(0);
                let total = targets.iter().filter(|&tgt| tgt == v).count();

                if used < total {
                    let mut path = u.clone();
                    path.push(*v.last().unwrap());
                    *visited_edges.entry(edge_key).or_insert(0) += 1;

                    let mut curr_node = v.clone();
                    while let Some(next_targets) = self.adj.get(&curr_node) {
                        let in_deg = self.in_degrees.get(&curr_node).copied().unwrap_or(0);
                        let out_deg = next_targets.len();

                        if in_deg == 1 && out_deg == 1 {
                            let next_node = &next_targets[0];
                            let next_edge = (curr_node.clone(), next_node.clone());
                            let next_used = visited_edges.get(&next_edge).copied().unwrap_or(0);
                            if next_used == 0 {
                                path.push(*next_node.last().unwrap());
                                *visited_edges.entry(next_edge).or_insert(0) += 1;
                                curr_node = next_node.clone();
                            } else {
                                break;
                            }
                        } else {
                            break;
                        }
                    }

                    if let Ok(dna) = Dna::new(&path) {
                        unitigs.push(dna);
                    }
                }
            }
        }

        unitigs
    }
}

/// Counts frequencies of all k-mers in a DNA sequence.
pub fn count_kmers(dna: &Dna, k: usize) -> Result<HashMap<Vec<u8>, usize>> {
    if k == 0 {
        return Err(BioError::InvalidInput("K must be > 0".into()));
    }
    let mut counts = HashMap::new();
    let bytes = &dna.sequence;
    if bytes.len() < k {
        return Ok(counts);
    }
    for i in 0..=(bytes.len() - k) {
        let kmer = bytes[i..i + k].to_vec();
        *counts.entry(kmer).or_insert(0) += 1;
    }
    Ok(counts)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_count_kmers() {
        let dna = Dna::new(b"ATATAT").unwrap();
        let counts = count_kmers(&dna, 2).unwrap();
        // 2-mers: AT (3), TA (2)
        assert_eq!(counts.get(&b"AT"[..]), Some(&3));
        assert_eq!(counts.get(&b"TA"[..]), Some(&2));
    }

    #[test]
    fn test_debruijn_graph_and_unitig_assembly() {
        // Linear unambiguous sequence: AAGTCG
        let original_seq = b"AAGTCG";
        let dna = Dna::new(original_seq).unwrap();
        let dbg = DeBruijnGraph::from_dna(&dna, 3).unwrap();

        // 3-mers: AAG, AGT, GTC, TCG (4 edges)
        // 2-mer nodes: AA, AG, GT, TC, CG (5 nodes)
        assert_eq!(dbg.edge_count(), 4);
        assert_eq!(dbg.node_count(), 5);

        let unitigs = dbg.assemble_unitigs();
        // Since there is no branching, it should assemble the complete original sequence
        assert!(!unitigs.is_empty());
        assert_eq!(unitigs[0].sequence, original_seq);
    }
}
