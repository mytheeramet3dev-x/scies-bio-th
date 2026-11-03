//! Spatial cell-cell neighborhood graphs and paracrine Ligand-Receptor signaling (`imaging::graph`).

use crate::imaging::spatial::CellCentroid;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Directed or undirected spatial proximity edge between two cells in tissue space.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SpatialNeighborEdge {
    /// Source cell ID.
    pub cell_a: String,
    /// Target neighbor cell ID.
    pub cell_b: String,
    /// Euclidean center-to-center distance in micrometers ($\mu\text{m}$).
    pub distance_um: f64,
    /// Spatial contact weight $\exp(-d / d_0)$ in $(0.0, 1.0]$.
    pub contact_weight: f64,
}

/// Spatial cell-cell interaction graph constructed over 2D tissue coordinates.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SpatialCellGraph {
    /// Maximum interaction radius in micrometers ($\mu\text{m}$).
    pub max_radius_um: f64,
    /// Proximity edges connecting neighboring cells.
    pub edges: Vec<SpatialNeighborEdge>,
}

/// Constructs a radius-based spatial neighborhood graph connecting cells within `max_radius_um`.
pub fn build_spatial_cell_graph(cells: &[CellCentroid], max_radius_um: f64) -> SpatialCellGraph {
    let mut edges = Vec::new();
    let decay_len = (max_radius_um * 0.5).max(1.0);

    for i in 0..cells.len() {
        for j in (i + 1)..cells.len() {
            let dx = cells[i].x - cells[j].x;
            let dy = cells[i].y - cells[j].y;
            let dist = (dx * dx + dy * dy).sqrt();
            if dist <= max_radius_um {
                let weight = (-dist / decay_len).exp().clamp(0.01, 1.0);
                edges.push(SpatialNeighborEdge {
                    cell_a: cells[i].cell_id.clone(),
                    cell_b: cells[j].cell_id.clone(),
                    distance_um: dist,
                    contact_weight: weight,
                });
            }
        }
    }

    SpatialCellGraph {
        max_radius_um,
        edges,
    }
}

/// Quantifies spatial paracrine Ligand-Receptor signaling (e.g. SASP `IL6 -> IL6R`)
/// across all neighbor edges of a [`SpatialCellGraph`].
///
/// Uses the mass-action geometric mean $\sqrt{L_i \cdot R_j} \cdot w_{ij}$.
pub fn score_spatial_ligand_receptor(
    graph: &SpatialCellGraph,
    ligand_by_cell: &std::collections::HashMap<String, f64>,
    receptor_by_cell: &std::collections::HashMap<String, f64>,
) -> f64 {
    if graph.edges.is_empty() {
        return 0.0;
    }

    let mut total_flux = 0.0;
    for edge in &graph.edges {
        let l_a = ligand_by_cell.get(&edge.cell_a).copied().unwrap_or(0.0).max(0.0);
        let r_b = receptor_by_cell.get(&edge.cell_b).copied().unwrap_or(0.0).max(0.0);

        let l_b = ligand_by_cell.get(&edge.cell_b).copied().unwrap_or(0.0).max(0.0);
        let r_a = receptor_by_cell.get(&edge.cell_a).copied().unwrap_or(0.0).max(0.0);

        let bidirectional = (l_a * r_b).sqrt() + (l_b * r_a).sqrt();
        total_flux += bidirectional * edge.contact_weight;
    }

    total_flux / (graph.edges.len() as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spatial_cell_graph_and_sasp_signaling() {
        let cells = vec![
            CellCentroid {
                cell_id: "Senescent_Sender".to_string(),
                x: 0.0,
                y: 0.0,
            },
            CellCentroid {
                cell_id: "Neighbor_Receiver".to_string(),
                x: 15.0,
                y: 0.0,
            },
            CellCentroid {
                cell_id: "Distant_Cell".to_string(),
                x: 200.0,
                y: 200.0,
            },
        ];

        let graph = build_spatial_cell_graph(&cells, 30.0);
        assert_eq!(graph.edges.len(), 1);

        let mut il6 = std::collections::HashMap::new();
        il6.insert("Senescent_Sender".to_string(), 10.0);
        let mut il6r = std::collections::HashMap::new();
        il6r.insert("Neighbor_Receiver".to_string(), 8.0);

        let flux = score_spatial_ligand_receptor(&graph, &il6, &il6r);
        assert!(flux > 1.0);
    }
}
