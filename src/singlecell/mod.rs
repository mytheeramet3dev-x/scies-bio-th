//! Single-cell trajectory inference (Diffusion Pseudotime / minimum-spanning state graphs)
//! and RNA splicing velocity kinetics (`singlecell`).

pub mod pseudotime;
pub mod velocity;

pub use pseudotime::{
    compute_diffusion_pseudotime, CellStateNode, PseudotimePoint, PseudotimeTrajectory,
};
pub use velocity::{
    estimate_rna_velocity, CellVelocityProfile, GeneVelocityMetrics, SplicingCounts,
};
