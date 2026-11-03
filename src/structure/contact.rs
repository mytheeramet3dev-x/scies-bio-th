//! Pairwise residue distance matrices, contact maps, and residue interaction topology.

use crate::structure::pdb::ProteinStructure3D;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Pairwise Cα–Cα Euclidean distance matrix (Ångströms).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct DistanceMatrix3D {
    /// Dimension $N \times N$ corresponding to residue count.
    pub size: usize,
    /// 1D flattened row-major distance entries.
    pub matrix: Vec<f64>,
}

impl DistanceMatrix3D {
    /// Retrieves distance between residue index $i$ and $j$.
    pub fn get(&self, i: usize, j: usize) -> Option<f64> {
        if i < self.size && j < self.size {
            Some(self.matrix[i * self.size + j])
        } else {
            None
        }
    }
}

/// Binary residue contact map indicating proximity below a distance threshold.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ContactMap3D {
    /// Dimension $N \times N$.
    pub size: usize,
    /// Proximity cutoff threshold in Ångströms (typically 8.0 Å).
    pub threshold: f64,
    /// 1D flattened row-major boolean contact entries.
    pub contacts: Vec<bool>,
}

impl ContactMap3D {
    /// Checks whether residue $i$ and $j$ are in contact.
    pub fn is_contact(&self, i: usize, j: usize) -> bool {
        if i < self.size && j < self.size {
            self.contacts[i * self.size + j]
        } else {
            false
        }
    }

    /// Evaluates contact density among non-local residues ($|i - j| \ge 4$).
    pub fn non_local_contact_density(&self) -> f64 {
        if self.size < 5 {
            return 0.0;
        }
        let mut total_pairs = 0;
        let mut contact_pairs = 0;

        for i in 0..self.size {
            for j in (i + 4)..self.size {
                total_pairs += 1;
                if self.is_contact(i, j) {
                    contact_pairs += 1;
                }
            }
        }

        if total_pairs == 0 {
            0.0
        } else {
            contact_pairs as f64 / total_pairs as f64
        }
    }
}

/// Computes the pairwise Cα–Cα distance matrix for a protein structure.
pub fn calculate_ca_distance_matrix(structure: &ProteinStructure3D) -> DistanceMatrix3D {
    let cas = structure.ca_atoms();
    let n = cas.len();
    let mut matrix = vec![0.0_f64; n * n];

    for i in 0..n {
        for j in i..n {
            let dist = cas[i].distance(cas[j]);
            matrix[i * n + j] = dist;
            matrix[j * n + i] = dist;
        }
    }

    DistanceMatrix3D { size: n, matrix }
}

/// Computes a binary contact map using a designated distance threshold (default 8.0 Å).
pub fn calculate_ca_contact_map(structure: &ProteinStructure3D, threshold_angstrom: f64) -> ContactMap3D {
    let dist_mat = calculate_ca_distance_matrix(structure);
    let n = dist_mat.size;
    let mut contacts = vec![false; n * n];

    for i in 0..n {
        for j in 0..n {
            let d = dist_mat.get(i, j).unwrap_or(f64::INFINITY);
            if d <= threshold_angstrom {
                contacts[i * n + j] = true;
            }
        }
    }

    ContactMap3D {
        size: n,
        threshold: threshold_angstrom,
        contacts,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::pdb::{Atom3D, Chain3D, Residue3D};

    #[test]
    fn test_contact_map() {
        let mut residues = Vec::new();
        // 5 residues along x axis at x = 0, 3, 6, 9, 12
        for i in 0..5 {
            residues.push(Residue3D {
                name: "GLY".to_string(),
                seq_num: i,
                atoms: vec![Atom3D {
                    serial: i as usize,
                    name: "CA".to_string(),
                    res_name: "GLY".to_string(),
                    chain_id: 'A',
                    res_seq: i,
                    x: (i as f64) * 3.0,
                    y: 0.0,
                    z: 0.0,
                    occupancy: 1.0,
                    b_factor: 0.0,
                    element: "C".to_string(),
                }],
            });
        }

        let structure = ProteinStructure3D {
            id: "TEST".to_string(),
            chains: vec![Chain3D { id: 'A', residues }],
        };

        let map = calculate_ca_contact_map(&structure, 7.0);
        assert_eq!(map.size, 5);
        // Distance 0 to 1 is 3.0 <= 7.0 (contact)
        assert!(map.is_contact(0, 1));
        // Distance 0 to 2 is 6.0 <= 7.0 (contact)
        assert!(map.is_contact(0, 2));
        // Distance 0 to 3 is 9.0 > 7.0 (no contact)
        assert!(!map.is_contact(0, 3));
    }
}
