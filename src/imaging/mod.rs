//! Microscopy and spatial biology primitives: 2D Gaussian fluorescent spot detection,
//! dual-channel co-localization (Pearson & Manders coefficients), and spatial transcriptomics assignment.

pub mod colocalize;
pub mod graph;
pub mod spatial;
pub mod spot;

pub use colocalize::{calculate_colocalization, ColocalizationMetrics};
pub use graph::{
    build_spatial_cell_graph, score_spatial_ligand_receptor, SpatialCellGraph, SpatialNeighborEdge,
};
pub use spatial::{
    assign_transcripts_to_cells, CellCentroid, SpatialCellAssignment, TranscriptSpot,
};
pub use spot::{detect_fluorescent_spots, FluorescentSpot, ImageChannel2D};
