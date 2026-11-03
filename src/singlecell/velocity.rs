//! RNA Velocity estimation from unspliced ($u$) and spliced ($s$) mRNA counts
//! (La Manno et al., steady-state kinetic model $\frac{ds}{dt} = u - \gamma s$).

use crate::BioError;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Unspliced and spliced mRNA read counts for a single gene across cells.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SplicingCounts {
    /// Gene symbol.
    pub gene_symbol: String,
    /// Unspliced nascent pre-mRNA abundance per cell ($u_i$).
    pub unspliced: Vec<f64>,
    /// Spliced mature mRNA abundance per cell ($s_i$).
    pub spliced: Vec<f64>,
}

/// Fitted steady-state degradation slope $\gamma$ and velocity residuals for a gene.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct GeneVelocityMetrics {
    /// Gene symbol.
    pub gene_symbol: String,
    /// Estimated steady-state degradation / splicing ratio $\gamma = \frac{\sum u_i s_i}{\sum s_i^2}$.
    pub gamma_slope: f64,
    /// Per-cell RNA velocity $v_i = u_i - \gamma s_i$ ($> 0$ transcriptional induction, $< 0$ repression).
    pub cell_velocities: Vec<f64>,
}

/// Aggregate RNA velocity magnitude and projected state for a single cell.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CellVelocityProfile {
    /// Cell index ($0\text{..}n$).
    pub cell_index: usize,
    /// Velocity vector across all evaluated genes.
    pub velocity_vector: Vec<f64>,
    /// L2 norm (speed of transcriptomic state transition).
    pub velocity_magnitude: f64,
}

/// Estimates steady-state degradation slopes $\gamma$ and per-cell RNA velocity vectors
/// across a panel of genes.
pub fn estimate_rna_velocity(
    genes: &[SplicingCounts],
) -> crate::Result<(Vec<GeneVelocityMetrics>, Vec<CellVelocityProfile>)> {
    if genes.is_empty() {
        return Err(BioError::InvalidInput(
            "At least one gene is required for RNA velocity estimation".to_string(),
        ));
    }
    let n_cells = genes[0].spliced.len();
    if n_cells == 0
        || genes
            .iter()
            .any(|g| g.spliced.len() != n_cells || g.unspliced.len() != n_cells)
    {
        return Err(BioError::InvalidInput(
            "All genes must have non-empty spliced and unspliced arrays of equal cell length"
                .to_string(),
        ));
    }

    let mut gene_metrics = Vec::with_capacity(genes.len());
    for g in genes {
        let mut sum_us = 0.0;
        let mut sum_ss = 0.0;
        for (&u, &s) in g.unspliced.iter().zip(g.spliced.iter()) {
            sum_us += u.max(0.0) * s.max(0.0);
            sum_ss += s.max(0.0) * s.max(0.0);
        }
        let gamma = if sum_ss > 1e-12 {
            (sum_us / sum_ss).max(0.01)
        } else {
            1.0
        };
        let velocities: Vec<f64> = g
            .unspliced
            .iter()
            .zip(g.spliced.iter())
            .map(|(&u, &s)| u - gamma * s)
            .collect();

        gene_metrics.push(GeneVelocityMetrics {
            gene_symbol: g.gene_symbol.clone(),
            gamma_slope: gamma,
            cell_velocities: velocities,
        });
    }

    let mut cell_profiles = Vec::with_capacity(n_cells);
    for c in 0..n_cells {
        let vec_c: Vec<f64> = gene_metrics
            .iter()
            .map(|gm| gm.cell_velocities[c])
            .collect();
        let mag = vec_c.iter().map(|v| v * v).sum::<f64>().sqrt();
        cell_profiles.push(CellVelocityProfile {
            cell_index: c,
            velocity_vector: vec_c,
            velocity_magnitude: mag,
        });
    }

    Ok((gene_metrics, cell_profiles))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rna_velocity_induction_and_repression() {
        let counts = vec![SplicingCounts {
            gene_symbol: "TERT".to_string(),
            // Cell 0 is strongly inducing (high u relative to s), Cell 2 is repressing (low u)
            unspliced: vec![12.0, 5.0, 1.0],
            spliced: vec![4.0, 5.0, 8.0],
        }];

        let (gm, cp) = estimate_rna_velocity(&counts).unwrap();
        assert!(gm[0].cell_velocities[0] > 0.0); // Induction
        assert!(gm[0].cell_velocities[2] < 0.0); // Repression
        assert_eq!(cp.len(), 3);
    }
}
