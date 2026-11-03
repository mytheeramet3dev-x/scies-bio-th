//! Innate immune response trigger estimation (TLR9 CpG sensing, RIG-I/MDA5 dsRNA, cGAS-STING).

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Nucleic acid and vector delivery characteristics for innate immune modeling.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct DeliveryPayloadProfile {
    /// Fraction of unmethylated CpG dinucleotides in the construct ($0.0\text{--}1.0$).
    pub unmethylated_cpg_fraction: f64,
    /// Fraction of uridine bases modified with N1-methylpseudouridine ($\text{m}^1\Psi$, $0.0\text{--}1.0$).
    pub pseudouridine_modification_fraction: f64,
    /// Double-stranded RNA (dsRNA) impurity fraction ($0.0\text{--}0.05$).
    pub dsrna_impurity_fraction: f64,
    /// Viral vector capsid dose in multiplicity of infection (MOI, $0.0$ for non-viral LNP).
    pub viral_vector_moi: f64,
}

/// Innate immune pathway activation estimate.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ImmuneTriggerAssessment {
    /// Endosomal TLR7/8/9 activation score in $[0.0, 1.0]$.
    pub tlr_activation: f64,
    /// Cytosolic RIG-I / MDA5 / cGAS-STING interferon response in $[0.0, 1.0]$.
    pub cytosolic_ifn_response: f64,
    /// Composite cytokine storm / pyroptosis trigger risk in $[0.0, 1.0]$.
    pub composite_immune_risk: f64,
    /// Whether the payload passes clinical innate-immunogenicity thresholds ($< 0.25$).
    pub is_immunologically_stealth: bool,
}

/// Estimates innate immune recognition and Type-I Interferon trigger probability.
pub fn estimate_innate_immune_trigger(
    payload: &DeliveryPayloadProfile,
) -> ImmuneTriggerAssessment {
    let unmodified_u = (1.0 - payload.pseudouridine_modification_fraction.clamp(0.0, 1.0)).max(0.0);
    let tlr_activation = (0.55 * payload.unmethylated_cpg_fraction.clamp(0.0, 1.0) * 4.0
        + 0.45 * unmodified_u)
        .clamp(0.0, 1.0);

    let dsrna_signal = (payload.dsrna_impurity_fraction.max(0.0) / 0.01).clamp(0.0, 1.0);
    let viral_signal = (payload.viral_vector_moi.max(0.0) / 100.0).clamp(0.0, 1.0);
    let cytosolic_ifn_response = (0.65 * dsrna_signal + 0.35 * viral_signal).clamp(0.0, 1.0);

    let composite_immune_risk = (0.50 * tlr_activation + 0.50 * cytosolic_ifn_response).clamp(0.0, 1.0);

    ImmuneTriggerAssessment {
        tlr_activation,
        cytosolic_ifn_response,
        composite_immune_risk,
        is_immunologically_stealth: composite_immune_risk < 0.25,
    }
}
