//! Phylogenetics, distance matrices, UPGMA tree reconstruction, and Newick tree export.
//!
//! Provides data structures for phylogenetic trees, symmetric distance matrices,
//! the UPGMA (Unweighted Pair Group Method with Arithmetic Mean) clustering algorithm,
//! and Newick format serialization.

use crate::error::{BioError, Result};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A square, symmetric pairwise distance matrix for taxa/sequences.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct DistanceMatrix {
    pub names: Vec<String>,
    pub matrix: Vec<f64>, // Row-major flat matrix N x N
}

impl DistanceMatrix {
    /// Creates a new `DistanceMatrix` from taxon names and a flat $N \times N$ matrix.
    pub fn new(names: Vec<String>, matrix: Vec<f64>) -> Result<Self> {
        let n = names.len();
        if matrix.len() != n * n {
            return Err(BioError::InvalidInput(format!(
                "Distance matrix size ({}) does not match taxon count squared ({} x {} = {})",
                matrix.len(),
                n,
                n,
                n * n
            )));
        }
        Ok(Self { names, matrix })
    }

    /// Number of taxa in the distance matrix.
    pub fn len(&self) -> usize {
        self.names.len()
    }

    /// Returns true if the distance matrix contains no taxa.
    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    /// Gets pairwise distance between taxon $i$ and taxon $j$.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        let n = self.names.len();
        self.matrix[i * n + j]
    }

    /// Sets pairwise distance between taxon $i$ and taxon $j$.
    pub fn set(&mut self, i: usize, j: usize, val: f64) {
        let n = self.names.len();
        self.matrix[i * n + j] = val;
    }
}

/// A node in a rooted phylogenetic tree.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PhyloNode {
    /// Taxon name (populated for leaf nodes).
    pub name: Option<String>,
    /// Branch length leading to parent node.
    pub branch_length: f64,
    /// Child nodes (empty for leaves).
    pub children: Vec<PhyloNode>,
}

impl PhyloNode {
    /// Creates a leaf node with a name and initial branch length.
    pub fn leaf(name: impl Into<String>, branch_length: f64) -> Self {
        Self {
            name: Some(name.into()),
            branch_length,
            children: Vec::new(),
        }
    }

    /// Creates an internal node joining child subtrees.
    pub fn internal(children: Vec<PhyloNode>, branch_length: f64) -> Self {
        Self {
            name: None,
            branch_length,
            children,
        }
    }

    /// Returns true if this node is a leaf (no children).
    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }

    /// Formats this node and its subtrees into Newick format string.
    pub fn to_newick(&self) -> String {
        let mut s = String::new();
        self.write_newick(&mut s);
        s
    }

    fn write_newick(&self, out: &mut String) {
        if self.is_leaf() {
            if let Some(ref name) = self.name {
                out.push_str(name);
            }
        } else {
            out.push('(');
            for (i, child) in self.children.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                child.write_newick(out);
            }
            out.push(')');
            if let Some(ref name) = self.name {
                out.push_str(name);
            }
        }
        if self.branch_length > 0.0 {
            use std::fmt::Write;
            let _ = write!(out, ":{:.4}", self.branch_length);
        }
    }
}

/// A rooted phylogenetic tree.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PhyloTree {
    pub root: PhyloNode,
}

impl PhyloTree {
    /// Exports the entire tree to a Newick format string (terminated with semicolon `;`).
    pub fn to_newick(&self) -> String {
        format!("{};", self.root.to_newick())
    }
}

/// Reconstructs a rooted phylogenetic tree from a distance matrix using the UPGMA algorithm.
pub fn upgma(dist: &DistanceMatrix) -> Result<PhyloTree> {
    let n = dist.len();
    if n == 0 {
        return Err(BioError::InvalidInput("Distance matrix cannot be empty".into()));
    }
    if n == 1 {
        return Ok(PhyloTree {
            root: PhyloNode::leaf(&dist.names[0], 0.0),
        });
    }

    // Cluster representation during clustering:
    // (Node, cluster_height, number_of_leaves)
    struct Cluster {
        node: PhyloNode,
        height: f64,
        size: usize,
    }

    let mut clusters: Vec<Option<Cluster>> = dist
        .names
        .iter()
        .map(|name| {
            Some(Cluster {
                node: PhyloNode::leaf(name, 0.0),
                height: 0.0,
                size: 1,
            })
        })
        .collect();

    // Working distance matrix
    let mut d = dist.matrix.clone();
    let get_d = |mat: &[f64], i: usize, j: usize| mat[i * n + j];
    let set_d = |mat: &mut [f64], i: usize, j: usize, val: f64| {
        mat[i * n + j] = val;
        mat[j * n + i] = val;
    };

    let mut active_clusters = n;

    while active_clusters > 1 {
        // 1. Find pair (min_i, min_j) with minimal distance among active clusters
        let mut min_dist = f64::INFINITY;
        let mut min_i = 0;
        let mut min_j = 0;

        for i in 0..n {
            if clusters[i].is_none() {
                continue;
            }
            for (j, cluster_j) in clusters.iter().enumerate().skip(i + 1) {
                if cluster_j.is_none() {
                    continue;
                }
                let val = get_d(&d, i, j);
                if val < min_dist {
                    min_dist = val;
                    min_i = i;
                    min_j = j;
                }
            }
        }

        if min_dist.is_infinite() {
            break;
        }

        // 2. Join clusters min_i and min_j
        let cluster_i = clusters[min_i].take().unwrap();
        let cluster_j = clusters[min_j].take().unwrap();

        let new_height = min_dist / 2.0;
        let branch_i = (new_height - cluster_i.height).max(0.0);
        let branch_j = (new_height - cluster_j.height).max(0.0);

        let mut child_i = cluster_i.node;
        child_i.branch_length = branch_i;
        let mut child_j = cluster_j.node;
        child_j.branch_length = branch_j;

        let new_node = PhyloNode::internal(vec![child_i, child_j], 0.0);
        let new_size = cluster_i.size + cluster_j.size;

        // 3. Update distances from the new cluster (stored at min_i) to all other active clusters
        for (k, cluster_k) in clusters.iter().enumerate() {
            if k == min_i || k == min_j || cluster_k.is_none() {
                continue;
            }
            let dist_ik = get_d(&d, min_i, k);
            let dist_jk = get_d(&d, min_j, k);
            let avg_dist = (cluster_i.size as f64 * dist_ik + cluster_j.size as f64 * dist_jk)
                / (new_size as f64);
            set_d(&mut d, min_i, k, avg_dist);
        }

        // Replace min_i with new cluster, and min_j remains None
        clusters[min_i] = Some(Cluster {
            node: new_node,
            height: new_height,
            size: new_size,
        });

        active_clusters -= 1;
    }

    let root_cluster = clusters
        .into_iter()
        .find_map(|c| c)
        .ok_or_else(|| BioError::SimulationError("UPGMA clustering failed to produce root".into()))?;

    Ok(PhyloTree {
        root: root_cluster.node,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_upgma_3_taxa() {
        // Distance matrix:
        //      A    B    C
        // A    0    2    4
        // B    2    0    4
        // C    4    4    0
        let names = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let mat = vec![
            0.0, 2.0, 4.0,
            2.0, 0.0, 4.0,
            4.0, 4.0, 0.0,
        ];
        let dm = DistanceMatrix::new(names, mat).unwrap();
        let tree = upgma(&dm).unwrap();

        let newick = tree.to_newick();
        // A and B should cluster first at height 1.0 (branch length 1.0)
        // Then (A, B) clusters with C at height 2.0 (branch length 2.0 to C, 1.0 to (A,B))
        assert!(newick.contains("A:1.0000"));
        assert!(newick.contains("B:1.0000"));
        assert!(newick.contains("C:2.0000"));
        assert!(newick.ends_with(';'));
    }

    #[test]
    fn test_distance_matrix_bounds() {
        let names = vec!["A".to_string(), "B".to_string()];
        let bad_mat = vec![0.0, 1.0]; // Only 2 elements, needs 4
        assert!(DistanceMatrix::new(names, bad_mat).is_err());
    }
}
