//! Molecular interaction scoring functions: Lennard-Jones 6-12 Van der Waals and Coulombic electrostatics.

use crate::structure::Atom3D;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Van der Waals Lennard-Jones parameters for an element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VdwParam {
    /// Collision diameter $\sigma$ in Ångströms.
    pub sigma: f64,
    /// Potential well depth $\epsilon$ in kcal/mol.
    pub epsilon: f64,
}

impl VdwParam {
    /// Default standard parameters for common biological elements (AMBER/CHARMM-like).
    pub fn for_element(element: &str) -> Self {
        match element.to_uppercase().as_str() {
            "H" => Self { sigma: 2.00, epsilon: 0.02 },
            "C" => Self { sigma: 3.40, epsilon: 0.10 },
            "N" => Self { sigma: 3.25, epsilon: 0.17 },
            "O" => Self { sigma: 3.00, epsilon: 0.20 },
            "P" => Self { sigma: 3.74, epsilon: 0.20 },
            "S" => Self { sigma: 3.56, epsilon: 0.25 },
            "F" => Self { sigma: 3.12, epsilon: 0.06 },
            "CL" => Self { sigma: 3.47, epsilon: 0.26 },
            "BR" => Self { sigma: 3.85, epsilon: 0.32 },
            "I" => Self { sigma: 4.10, epsilon: 0.40 },
            _ => Self { sigma: 3.40, epsilon: 0.10 }, // default carbon-like
        }
    }
}

/// Interaction energy score breakdown.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct DockingScore {
    /// Lennard-Jones 6-12 Van der Waals energy in kcal/mol.
    pub vdw_energy: f64,
    /// Coulombic electrostatic energy in kcal/mol.
    pub electrostatic_energy: f64,
    /// Total combined binding score.
    pub total_score: f64,
}

/// Evaluates molecular docking interaction energy between receptor atoms and ligand atoms.
///
/// - **Van der Waals**: Lennard-Jones 12-6 potential with Lorentz-Berthelot mixing rules:
///   $$\sigma_{ij} = \frac{\sigma_i + \sigma_j}{2}, \quad \epsilon_{ij} = \sqrt{\epsilon_i \epsilon_j}$$
///   $$E_{\text{vdW}} = \sum 4\epsilon_{ij} \left[ \left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^6 \right]$$
/// - **Electrostatics**: Distance-dependent dielectric $\epsilon(r) = 4.0 \cdot r$:
///   $$E_{\text{elec}} = \sum \frac{332.0637 \cdot q_i q_j}{4.0 \cdot r_{ij}^2}$$
pub fn calculate_interaction_energy(
    receptor_atoms: &[Atom3D],
    ligand_atoms: &[Atom3D],
    receptor_charges: Option<&[f64]>,
    ligand_charges: Option<&[f64]>,
    distance_cutoff: f64,
) -> DockingScore {
    let mut total_vdw = 0.0;
    let mut total_elec = 0.0;
    let cutoff_sq = distance_cutoff * distance_cutoff;

    let has_charges = receptor_charges.is_some() && ligand_charges.is_some();
    let rec_q = receptor_charges.unwrap_or(&[]);
    let lig_q = ligand_charges.unwrap_or(&[]);

    for (i, rec) in receptor_atoms.iter().enumerate() {
        let p_rec = VdwParam::for_element(&rec.element);
        let q_i = if has_charges && i < rec_q.len() { rec_q[i] } else { 0.0 };

        for (j, lig) in ligand_atoms.iter().enumerate() {
            let dx = rec.x - lig.x;
            let dy = rec.y - lig.y;
            let dz = rec.z - lig.z;
            let dist_sq = dx * dx + dy * dy + dz * dz;

            if dist_sq > cutoff_sq || dist_sq < 0.01 {
                continue;
            }

            let dist = dist_sq.sqrt().max(0.8); // avoid singularity below 0.8 Å
            let p_lig = VdwParam::for_element(&lig.element);

            // Lorentz-Berthelot mixing
            let sigma = 0.5 * (p_rec.sigma + p_lig.sigma);
            let epsilon = (p_rec.epsilon * p_lig.epsilon).sqrt();

            let sr = sigma / dist;
            let sr3 = sr * sr * sr;
            let sr6 = sr3 * sr3;
            let sr12 = sr6 * sr6;

            let e_vdw = 4.0 * epsilon * (sr12 - sr6);
            // Cap excessive steric clashes
            total_vdw += e_vdw.min(50.0);

            // Electrostatics with distance-dependent dielectric epsilon(r) = 4r
            if has_charges && j < lig_q.len() {
                let q_j = lig_q[j];
                if q_i != 0.0 && q_j != 0.0 {
                    // Coulomb constant 332.0637 kcal*Å / (mol * e^2)
                    let e_elec = (332.0637 * q_i * q_j) / (4.0 * dist * dist);
                    total_elec += e_elec;
                }
            }
        }
    }

    DockingScore {
        vdw_energy: total_vdw,
        electrostatic_energy: total_elec,
        total_score: total_vdw + total_elec,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interaction_energy_vdw() {
        let rec = vec![Atom3D {
            serial: 1,
            name: "CA".to_string(),
            res_name: "ALA".to_string(),
            chain_id: 'A',
            res_seq: 1,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            b_factor: 10.0,
            element: "C".to_string(),
        }];

        // Carbon at ~van der Waals minimum (sigma = 3.4 Å, minimum at r = 2^(1/6)*sigma = 3.816 Å)
        let lig = vec![Atom3D {
            serial: 2,
            name: "C1".to_string(),
            res_name: "LIG".to_string(),
            chain_id: 'L',
            res_seq: 1,
            x: 3.816,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            b_factor: 10.0,
            element: "C".to_string(),
        }];

        let score = calculate_interaction_energy(&rec, &lig, None, None, 10.0);
        assert!(score.vdw_energy < 0.0); // Attractive well at optimal distance
        assert_eq!(score.electrostatic_energy, 0.0);
    }
}
