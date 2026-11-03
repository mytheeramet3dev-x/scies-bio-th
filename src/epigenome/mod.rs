//! Epigenomics primitives: DNA methylation, chromatin accessibility, and histone modification states.

pub mod chromatin;
pub mod histone;
pub mod methylation;

pub use chromatin::{evaluate_chromatin_accessibility, ChromatinAccessibilityProfile, ChromatinPeak};
pub use histone::{
    infer_chromatin_state, quantify_h3k27ac_signal, ChromatinState, H3K27acQuantification,
    HistoneMark, HistoneSegment,
};
pub use methylation::{beta_to_m_value, find_cpg_islands, m_value_to_beta, CpgIsland};
