//! Pharmacokinetics and Pharmacodynamics (`pharma`): one-compartment dosing kinetics,
//! Hill concentration-response curves, and therapeutic window optimization.

pub mod pkpd;

pub use pkpd::{
    evaluate_hill_response, simulate_pkpd_regimen, simulate_two_compartment_pkpd,
    AdvancedPkPdReport, AgingLeverModulation, DosingEvent, HillPharmacodynamics,
    IntermittentDosingSchedule, PkPdParameters, PkPdRegimenReport, RegimenDosingRecommendation,
    TwoCompartmentPk,
};
