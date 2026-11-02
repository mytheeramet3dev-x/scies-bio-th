//! Longevity biology, telomere dynamics, and synthetic immortal genome design.
//!
//! Provides models and precision engineering tools for:
//! - Chromosomal telomere shortening and replicative lifespan prediction (Hayflick limit).
//! - Telomerase (TERT/TERC) activity and oncogenic risk modeling.
//! - Epigenetic aging clocks and biological age acceleration estimation.
//! - Canonical longevity pathway metadata (SIRT1/6, FOXO3, TERT, p53, CDKN2A).
//! - Synthetic genome design pipeline proposing CRISPR sgRNA edits for healthy lifespan extension.

pub mod clock;
pub mod design;
pub mod pathways;
pub mod telomerase;
pub mod telomere;

pub use clock::{calculate_aging_acceleration, estimate_biological_age, AgingProfile};
pub use design::{
    design_immortal_genome_package, EditType, LongevityDesignReport, LongevityStrategy,
    ProposedLongevityEdit,
};
pub use pathways::{
    find_longevity_gene, LongevityGeneInfo, LongevityPathway, CANONICAL_LONGEVITY_GENES,
};
pub use telomerase::{apply_telomerase, cancer_risk_factor, TelomeraseActivity};
pub use telomere::{
    count_telomeric_repeats, predict_replicative_lifespan, simulate_telomere_shortening,
    TelomereState, VERTEBRATE_TELOMERE_REPEAT, VERTEBRATE_TELOMERE_REVERSE,
};
