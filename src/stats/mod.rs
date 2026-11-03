//! Statistical methods for biological data, non-parametric tests, and FDR correction.

pub mod correction;
pub mod dist;
pub mod test;

pub use correction::{benjamini_hochberg, bonferroni_correction};
pub use dist::{cosine_distance, euclidean_distance, manhattan_distance};
pub use test::{
    chi_squared_test, fishers_exact_2x2, kruskal_wallis, mann_whitney_u, wilcoxon_signed_rank,
    ChiSquaredResult, FishersExactResult, KruskalWallisResult, MannWhitneyResult, WilcoxonResult,
};
