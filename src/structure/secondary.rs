//! Secondary structure assignment (DSSP-like Alpha-Helix, Beta-Sheet, and Coil).

use crate::structure::pdb::Residue3D;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Three-state protein secondary structure classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum SecondaryStructure {
    /// Alpha-helix ('H', 3.6 residues per turn, ~1.5 Å rise).
    Helix,
    /// Beta-strand / Extended sheet ('E').
    Sheet,
    /// Irregular loop, turn, or coil ('C').
    Coil,
}

impl SecondaryStructure {
    /// Single-letter code: 'H', 'E', or 'C'.
    pub fn to_char(self) -> char {
        match self {
            SecondaryStructure::Helix => 'H',
            SecondaryStructure::Sheet => 'E',
            SecondaryStructure::Coil => 'C',
        }
    }
}

/// Assigns DSSP-like secondary structure elements using backbone Cα geometric criteria.
pub fn assign_secondary_structure(residues: &[Residue3D]) -> Vec<SecondaryStructure> {
    let n = residues.len();
    if n < 4 {
        return vec![SecondaryStructure::Coil; n];
    }

    let mut ca_coords: Vec<Option<[f64; 3]>> = Vec::with_capacity(n);
    for res in residues {
        ca_coords.push(res.ca_atom().map(|a| a.coords()));
    }

    let mut sec = vec![SecondaryStructure::Coil; n];

    // Helper to calculate Euclidean distance between Cα indices
    let dist = |i: usize, j: usize| -> Option<f64> {
        let p = ca_coords[i]?;
        let q = ca_coords[j]?;
        let dx = p[0] - q[0];
        let dy = p[1] - q[1];
        let dz = p[2] - q[2];
        Some((dx * dx + dy * dy + dz * dz).sqrt())
    };

    // 1. Detect Alpha-Helices via i, i+3 and i, i+4 proximity (typical ~5.0 - 6.2 Å)
    for i in 0..n.saturating_sub(4) {
        if let (Some(d3), Some(d4)) = (dist(i, i + 3), dist(i, i + 4)) {
            if d3 < 6.5 && d4 < 7.0 {
                sec[i..=i + 4].fill(SecondaryStructure::Helix);
            }
        }
    }

    // 2. Detect Beta-Sheets via extended Cα(i, i+2) distances (6.5 - 7.5 Å) and non-local strand contacts
    for i in 0..n.saturating_sub(2) {
        if sec[i] == SecondaryStructure::Helix {
            continue;
        }
        if let Some(d2) = dist(i, i + 2) {
            if (6.4..=7.6).contains(&d2) {
                // Check if this residue is in contact with another non-local extended residue (|i - j| >= 5)
                for j in 0..n {
                    if (i as isize - j as isize).abs() >= 5 {
                        if let Some(d_ij) = dist(i, j) {
                            if d_ij <= 5.8 {
                                sec[i] = SecondaryStructure::Sheet;
                                if i + 1 < n {
                                    sec[i + 1] = SecondaryStructure::Sheet;
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    sec
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::pdb::Atom3D;

    fn make_test_residue(seq: i32, x: f64, y: f64, z: f64) -> Residue3D {
        Residue3D {
            name: "ALA".to_string(),
            seq_num: seq,
            atoms: vec![Atom3D {
                serial: seq as usize,
                name: "CA".to_string(),
                res_name: "ALA".to_string(),
                chain_id: 'A',
                res_seq: seq,
                x,
                y,
                z,
                occupancy: 1.0,
                b_factor: 10.0,
                element: "C".to_string(),
            }],
        }
    }

    #[test]
    fn test_alpha_helix_assignment() {
        // Idealized alpha helix: radius ~2.3 Å, pitch ~5.4 Å, 3.6 residues/turn
        // theta = i * 100 degrees = i * (2*pi / 3.6), z = i * 1.5 Å
        let mut residues = Vec::new();
        for i in 0..10 {
            let theta = (i as f64) * 100.0 * (std::f64::consts::PI / 180.0);
            let x = 2.3 * theta.cos();
            let y = 2.3 * theta.sin();
            let z = (i as f64) * 1.5;
            residues.push(make_test_residue(i, x, y, z));
        }

        let ss = assign_secondary_structure(&residues);
        // Middle residues should be assigned Helix
        assert_eq!(ss[4], SecondaryStructure::Helix);
        assert_eq!(ss[5], SecondaryStructure::Helix);
    }
}
