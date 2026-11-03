//! Spatial transcriptomics single-cell transcript assignment (`imaging::spatial`).

use std::collections::HashMap;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Segmented cell nucleus / soma centroid in 2D tissue coordinates ($\mu\text{m}$).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CellCentroid {
    /// Unique cell identifier.
    pub cell_id: String,
    /// X coordinate in micrometers ($\mu\text{m}$).
    pub x: f64,
    /// Y coordinate in micrometers ($\mu\text{m}$).
    pub y: f64,
}

/// Decoded spatial mRNA transcript molecule spot (e.g. MERFISH / Xenium / CosMx).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TranscriptSpot {
    /// Gene symbol of the transcript.
    pub gene_symbol: String,
    /// X coordinate in micrometers ($\mu\text{m}$).
    pub x: f64,
    /// Y coordinate in micrometers ($\mu\text{m}$).
    pub y: f64,
}

/// Single-cell spatial gene count matrix row assigned by coordinate proximity.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SpatialCellAssignment {
    /// Cell identifier.
    pub cell_id: String,
    /// Assigned gene -> transcript molecule count.
    pub gene_counts: HashMap<String, usize>,
    /// Total transcripts assigned to this cell.
    pub total_transcripts: usize,
}

/// Assigns spatial [`TranscriptSpot`] molecules to the nearest [`CellCentroid`] within `max_radius_um`.
pub fn assign_transcripts_to_cells(
    cells: &[CellCentroid],
    transcripts: &[TranscriptSpot],
    max_radius_um: f64,
) -> Vec<SpatialCellAssignment> {
    let max_r_sq = max_radius_um * max_radius_um;
    let mut assignments: Vec<SpatialCellAssignment> = cells
        .iter()
        .map(|c| SpatialCellAssignment {
            cell_id: c.cell_id.clone(),
            gene_counts: HashMap::new(),
            total_transcripts: 0,
        })
        .collect();

    for spot in transcripts {
        let mut best_idx = None;
        let mut best_dist_sq = max_r_sq;

        for (i, cell) in cells.iter().enumerate() {
            let dx = spot.x - cell.x;
            let dy = spot.y - cell.y;
            let dist_sq = dx * dx + dy * dy;
            if dist_sq <= best_dist_sq {
                best_dist_sq = dist_sq;
                best_idx = Some(i);
            }
        }

        if let Some(idx) = best_idx {
            *assignments[idx]
                .gene_counts
                .entry(spot.gene_symbol.clone())
                .or_insert(0) += 1;
            assignments[idx].total_transcripts += 1;
        }
    }

    assignments
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spatial_transcript_assignment() {
        let cells = vec![
            CellCentroid {
                cell_id: "Cell_1".to_string(),
                x: 10.0,
                y: 10.0,
            },
            CellCentroid {
                cell_id: "Cell_2".to_string(),
                x: 50.0,
                y: 50.0,
            },
        ];
        let spots = vec![
            TranscriptSpot {
                gene_symbol: "TERT".to_string(),
                x: 11.5,
                y: 9.0,
            },
            TranscriptSpot {
                gene_symbol: "TP53".to_string(),
                x: 10.2,
                y: 10.8,
            },
            TranscriptSpot {
                gene_symbol: "CDKN2A".to_string(),
                x: 49.0,
                y: 51.0,
            },
        ];

        let assigned = assign_transcripts_to_cells(&cells, &spots, 10.0);
        assert_eq!(assigned[0].total_transcripts, 2);
        assert_eq!(assigned[0].gene_counts.get("TERT"), Some(&1));
        assert_eq!(assigned[1].total_transcripts, 1);
        assert_eq!(assigned[1].gene_counts.get("CDKN2A"), Some(&1));
    }
}
