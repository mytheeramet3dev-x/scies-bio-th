//! Personalized Cancer Neoantigen & mRNA Vaccine Design Engine.
//!
//! Provides automated pipelines for:
//! - Extracting patient-specific somatic neoantigens from cancer mutations.
//! - Predicting MHC Class I (HLA) binding affinity & Differential Agretopicity Index (DAI).
//! - Assembling multi-epitope polytope constructs with proteasome cleavage linkers.
//! - Designing clinical-grade mRNA vaccines with human codon optimization and regulatory UTRs.

pub mod mhc;
pub mod mrna_design;
pub mod neoantigen;

pub use mhc::{predict_mhc_binding, HlaAllele, MhcBindingResult, PredictionConfidence};
pub use mrna_design::{human_codon_optimize, MrnaVaccineConstruct};
pub use neoantigen::Neoantigen;
