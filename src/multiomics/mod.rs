//! Cross-layer multi-omics integration: RNA-seq + ATAC-seq + Methylation fusion,
//! molecular signature PCA / cosine similarity, and NMF-lite reference cell-type deconvolution.

pub mod deconvolution;
pub mod fusion;
pub mod signature;

pub use deconvolution::{deconvolve_cell_types_nmf, DeconvolutionResult};
pub use fusion::{fuse_multiomics_locus, MultiOmicsFusionResult, OmicsLocusInput};
pub use signature::{
    compare_molecular_signatures, project_pca_2d, MolecularSignature, SignatureComparison,
};
