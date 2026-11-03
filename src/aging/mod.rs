//! Longevity biology, telomere dynamics, and synthetic immortal genome design.
//!
//! Provides models and precision engineering tools for:
//! - Chromosomal telomere shortening and replicative lifespan prediction (Hayflick limit).
//! - Telomerase (TERT/TERC) activity and oncogenic risk modeling.
//! - Epigenetic aging clocks and biological age acceleration estimation.
//! - Canonical longevity pathway metadata (SIRT1/6, FOXO3, TERT, p53, CDKN2A).
//! - Synthetic genome design pipeline proposing CRISPR sgRNA edits for healthy lifespan extension.

pub mod clock;
pub mod constraints;
pub mod design;
pub mod fidelity;
pub mod guide;
pub mod pathways;
pub mod real_genome;
pub mod report;
pub mod risk;
pub mod simulation;
pub mod telomerase;
pub mod telomere;
pub mod validation;

pub use clock::{calculate_aging_acceleration, estimate_biological_age, AgingProfile};
pub use constraints::{
    enforce_and_auto_remediate, validate_biological_constraints, ConstraintSeverity,
    ConstraintViolation, ConstraintViolationType,
};
pub use design::{
    design_high_fidelity_longevity_package, design_immortal_genome_package, EditType,
    LongevityDesignReport, LongevityStrategy, ProposedLongevityEdit,
};
pub use fidelity::{
    apply_high_fidelity_telomerase, estimate_high_fidelity_biological_age,
    simulate_high_fidelity_shortening, HighFidelityAgeParams, HighFidelityAgeResult,
    HighFidelityTelomereState, TelomeraseKinetics,
};
pub use guide::{design_real_guides, select_optimal_guide, RealGuide};
pub use pathways::{
    find_longevity_gene, LongevityGeneInfo, LongevityPathway, CANONICAL_LONGEVITY_GENES,
};
pub use real_genome::{
    analyze_real_telomere, load_loci_from_fasta, load_loci_from_file,
    load_telomere_regions_from_fasta, load_telomere_regions_from_file, ChromosomeArm,
    TelomereAnalysis, TelomereRegion, CANONICAL_FORWARD_REPEAT, CANONICAL_REVERSE_REPEAT,
    FORWARD_VARIANT_REPEATS, REVERSE_VARIANT_REPEATS,
};
pub use report::{
    HighFidelityDesignReport, HighFidelityProposedEdit, LongevityMetricsComparison,
};
pub use risk::{
    evaluate_oncogenic_risk, OncogenicRiskAssessment, PathwayRiskBreakdown, SafetyLocks,
};
pub use simulation::{
    run_multi_generation_simulation, CellPopulationState, MultiGenerationSimulationConfig,
    SimulationTrajectory,
};
pub use telomerase::{apply_telomerase, cancer_risk_factor, TelomeraseActivity};
pub use telomere::{
    count_telomeric_repeats, predict_replicative_lifespan, simulate_telomere_shortening,
    TelomereState, VERTEBRATE_TELOMERE_REPEAT, VERTEBRATE_TELOMERE_REVERSE,
};
pub use validation::{validate_longevity_package, ValidationOutcome};
