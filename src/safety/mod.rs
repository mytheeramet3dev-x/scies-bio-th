//! Biological safety formalism: DNA damage checkpoints (p53/Rb/ATM), proto-oncogene risk,
//! innate immune nucleic-acid sensing, and a composable pre-intervention safety audit engine.

pub mod audit;
pub mod checkpoint;
pub mod immune;
pub mod oncogene;

pub use audit::{
    evaluate_intervention_safety, AuditRuleResult, AuditSeverity, SafetyAuditReport, SafetyContext,
    SafetyRule,
};
pub use checkpoint::{
    evaluate_checkpoint_activation, CellFateDecision, CheckpointInput, CheckpointResponse,
};
pub use immune::{estimate_innate_immune_trigger, DeliveryPayloadProfile, ImmuneTriggerAssessment};
pub use oncogene::{
    evaluate_oncogenic_activation_risk, OncogeneInterventionInput, OncogenicRiskAssessment,
};
