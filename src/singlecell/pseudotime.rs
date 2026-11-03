//! Single-cell Diffusion Pseudotime (DPT-lite) and minimum-spanning trajectory inference.
//!
//! Orders single cells along a developmental or reprogramming trajectory starting from a specified
//! root cell, normalizing geodesic transition distances into $[0.0, 1.0]$ and identifying
//! high-degree branch / bifurcation cells.

use crate::stats::dist::euclidean_distance;
use crate::BioError;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A single cell represented by its multi-omics / transcriptomic state vector.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CellStateNode {
    /// Unique single-cell barcode or ID.
    pub cell_id: String,
    /// Expression or multi-omics embedding coordinates.
    pub expression: Vec<f64>,
}

/// Computed pseudotime coordinate and topological role for a single cell.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PseudotimePoint {
    /// Cell identifier.
    pub cell_id: String,
    /// Normalized geodesic pseudotime in $[0.0, 1.0]$ ($0.0 = \text{root state}$).
    pub pseudotime: f64,
    /// Parent cell ID along the minimum-spanning transition tree (`None` for root).
    pub predecessor_id: Option<String>,
    /// True if this cell is a trajectory bifurcation / branch point ($\ge 2$ downstream branches).
    pub is_bifurcation: bool,
}

/// Complete single-cell trajectory produced by [`compute_diffusion_pseudotime`].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PseudotimeTrajectory {
    /// Root cell ID ($t = 0.0$).
    pub root_cell_id: String,
    /// Ordered cells from earliest ($0.0$) to latest ($1.0$) pseudotime.
    pub points: Vec<PseudotimePoint>,
    /// IDs of detected lineage bifurcation / branch points.
    pub bifurcation_cells: Vec<String>,
}

/// Builds a Prim minimum-spanning cell-transition graph in expression space and computes
/// normalized geodesic pseudotime from `root_cell_id`.
pub fn compute_diffusion_pseudotime(
    cells: &[CellStateNode],
    root_cell_id: &str,
) -> crate::Result<PseudotimeTrajectory> {
    let n = cells.len();
    if n < 2 {
        return Err(BioError::InvalidInput(
            "At least 2 cells are required for pseudotime trajectory inference".to_string(),
        ));
    }
    let dim = cells[0].expression.len();
    if dim == 0 || cells.iter().any(|c| c.expression.len() != dim) {
        return Err(BioError::InvalidInput(
            "All cells must have non-empty expression vectors of identical length".to_string(),
        ));
    }

    let root_idx = cells
        .iter()
        .position(|c| c.cell_id == root_cell_id)
        .ok_or_else(|| {
            BioError::InvalidInput(format!("Root cell '{}' not found in dataset", root_cell_id))
        })?;

    // Prim's algorithm rooted at root_idx to form the cell-state transition tree
    let mut in_tree = vec![false; n];
    let mut min_dist = vec![f64::INFINITY; n];
    let mut parent: Vec<Option<usize>> = vec![None; n];
    let mut geodesic = vec![0.0_f64; n];

    min_dist[root_idx] = 0.0;

    for _ in 0..n {
        let mut u = None;
        let mut best_d = f64::INFINITY;
        for (i, &vis) in in_tree.iter().enumerate() {
            if !vis && min_dist[i] < best_d {
                best_d = min_dist[i];
                u = Some(i);
            }
        }
        let u = match u {
            Some(idx) => idx,
            None => break,
        };
        in_tree[u] = true;
        if let Some(p) = parent[u] {
            geodesic[u] = geodesic[p] + min_dist[u];
        }

        for v in 0..n {
            if !in_tree[v] {
                let d = euclidean_distance(&cells[u].expression, &cells[v].expression)?;
                if d < min_dist[v] {
                    min_dist[v] = d;
                    parent[v] = Some(u);
                }
            }
        }
    }

    // Count children in the MST to detect bifurcation nodes (>= 2 children)
    let mut child_counts = vec![0usize; n];
    for p in parent.iter().flatten() {
        child_counts[*p] += 1;
    }

    let max_geo = geodesic.iter().copied().fold(0.0_f64, f64::max).max(1e-12);

    let mut points = Vec::with_capacity(n);
    let mut bifurcation_cells = Vec::new();

    for i in 0..n {
        let is_bif = child_counts[i] >= 2;
        if is_bif {
            bifurcation_cells.push(cells[i].cell_id.clone());
        }
        points.push(PseudotimePoint {
            cell_id: cells[i].cell_id.clone(),
            pseudotime: (geodesic[i] / max_geo).clamp(0.0, 1.0),
            predecessor_id: parent[i].map(|p| cells[p].cell_id.clone()),
            is_bifurcation: is_bif,
        });
    }

    points.sort_by(|a, b| {
        a.pseudotime
            .partial_cmp(&b.pseudotime)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    Ok(PseudotimeTrajectory {
        root_cell_id: root_cell_id.to_string(),
        points,
        bifurcation_cells,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_diffusion_pseudotime_ordering_and_bifurcation() {
        let cells = vec![
            CellStateNode {
                cell_id: "Somatic_Root".to_string(),
                expression: vec![0.0, 0.0, 0.0],
            },
            CellStateNode {
                cell_id: "Rejuvenated_Intermediate".to_string(),
                expression: vec![1.0, 1.0, 0.0],
            },
            CellStateNode {
                cell_id: "Branch_A_Safe_Somatic".to_string(),
                expression: vec![2.0, 1.2, 0.0],
            },
            CellStateNode {
                cell_id: "Branch_B_Pluripotent".to_string(),
                expression: vec![1.0, 2.5, 2.0],
            },
        ];

        let traj = compute_diffusion_pseudotime(&cells, "Somatic_Root").unwrap();
        assert_eq!(traj.points[0].cell_id, "Somatic_Root");
        assert_eq!(traj.points[0].pseudotime, 0.0);
        assert_eq!(traj.points[1].cell_id, "Rejuvenated_Intermediate");
        assert!(traj
            .bifurcation_cells
            .contains(&"Rejuvenated_Intermediate".to_string()));
    }
}
