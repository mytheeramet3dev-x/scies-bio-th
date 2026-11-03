//! High-fidelity oncogenic transformation risk modeling and multi-layered safety locks.
//!
//! Evaluates carcinogenic hazards across specific cellular failure modes:
//! 1. Replicative crisis bypass (M1/M2) via unconstrained telomerase activation.
//! 2. Accidental proto-oncogene transactivation via off-target CRISPR cleavage.
//! 3. Checkpoint loss or tumor suppressor exhaustion (TP53 / CDKN2A).
//! 4. Chromosomal instability (dicentric bridges) triggered by shelterin uncapping.
//!
//! Mitigated by active synthetic safeguards (Super-p53 tandem copies, inducible suicide switches).

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Detailed breakdown of neoplastic risk across distinct cellular pathways.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PathwayRiskBreakdown {
    /// Risk of bypassing cellular senescence (M1) and crisis (M2) checkpoints without DNA repair.
    pub crisis_bypass_risk: f64,
    /// Risk of spurious off-target cleavage activating proto-oncogenes (e.g. MYC, KRAS).
    pub proto_oncogene_transactivation_risk: f64,
    /// Risk arising from insufficient p53 / p16 cell-cycle arrest surveillance.
    pub checkpoint_loss_risk: f64,
    /// Risk of telomere-to-telomere end fusions and break-fusion-bridge cycles due to low shelterin.
    pub chromosomal_instability_risk: f64,
    /// Unmitigated aggregate risk score before synthetic safeguard intervention.
    pub total_unmitigated_risk: f64,
}

/// Active genetic and synthetic safeguards guarding against oncogenesis.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SafetyLocks {
    /// Extra dosage / pachyderm Super-p53 genomic copy active.
    pub super_p53_active: bool,
    /// Inducible cell suicide circuit (e.g. iCasp9 / AP1903 or HSV-TK).
    pub inducible_suicide_switch: bool,
    /// High-fidelity engineered nuclease variant (e.g. SpCas9-HF1, HiFi-Cas9) minimizing off-targets.
    pub high_fidelity_nuclease: bool,
    /// Promoter restriction confining expression to post-mitotic or controlled cell states.
    pub cell_cycle_restriction: bool,
}

impl SafetyLocks {
    /// Maximum standard safety configuration for balanced and conservative designs.
    pub fn standard_protection() -> Self {
        Self {
            super_p53_active: true,
            inducible_suicide_switch: true,
            high_fidelity_nuclease: true,
            cell_cycle_restriction: true,
        }
    }

    /// Minimal or absent safeguards corresponding to aggressive unmodified constructs.
    pub fn none() -> Self {
        Self {
            super_p53_active: false,
            inducible_suicide_switch: false,
            high_fidelity_nuclease: false,
            cell_cycle_restriction: false,
        }
    }

    /// Calculates aggregate risk reduction factor ($0.0$ to $0.95$).
    pub fn cumulative_mitigation_efficiency(&self) -> f64 {
        let mut retention = 1.0_f64;
        if self.super_p53_active {
            retention *= 0.55; // 45% reduction in crisis/checkpoint failure
        }
        if self.inducible_suicide_switch {
            retention *= 0.50; // 50% reduction in malignant clone expansion
        }
        if self.high_fidelity_nuclease {
            retention *= 0.65; // 35% reduction in off-target oncogene transactivation
        }
        if self.cell_cycle_restriction {
            retention *= 0.70; // 30% reduction in unscheduled proliferation
        }
        (1.0_f64 - retention).clamp(0.0_f64, 0.95_f64)
    }
}

/// Comprehensive oncogenic transformation risk assessment report.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct OncogenicRiskAssessment {
    /// Unmitigated baseline risk breakdown by pathway.
    pub raw_pathways: PathwayRiskBreakdown,
    /// Active safeguards deployed in this synthetic design.
    pub safeguards: SafetyLocks,
    /// Fractional mitigation achieved by active safeguards (e.g. 0.88 = 88% reduction).
    pub mitigation_factor: f64,
    /// Net transformed risk index after safety lock suppression (0.0 to 1.0).
    pub net_risk_index: f64,
    /// Whether this design falls within clinical safety tolerance (< 0.15).
    pub is_acceptable_safety_profile: bool,
}

/// Evaluates high-fidelity oncogenic transformation hazards across pathways,
/// integrating the formal `safety::oncogene` scoring engine.
pub fn evaluate_oncogenic_risk(
    telomerase_expression: f64,
    sgrna_offtarget_scores: &[f64],
    shelterin_capping: f64,
    safeguards: &SafetyLocks,
) -> OncogenicRiskAssessment {
    use crate::safety::oncogene::{
        evaluate_oncogenic_activation_risk, OncogeneInterventionInput,
    };

    // 1. Crisis bypass: scales non-linearly with unconstrained TERT
    let crisis_bypass = if telomerase_expression <= 0.05 {
        0.01
    } else {
        (telomerase_expression.powf(1.4) * 0.55).min(0.95)
    };

    // 2. Proto-oncogene transactivation from sgRNA off-target potential
    let avg_offtarget: f64 = if sgrna_offtarget_scores.is_empty() {
        0.02
    } else {
        sgrna_offtarget_scores.iter().sum::<f64>() / sgrna_offtarget_scores.len() as f64
    };
    let proto_oncogene_risk = (avg_offtarget * 0.45).clamp(0.01, 0.85);

    // 3. Formal safety::oncogene risk model
    let formal_onco = evaluate_oncogenic_activation_risk(&OncogeneInterventionInput {
        tert_fold: 1.0 + telomerase_expression * 2.5,
        myc_fold: 1.0,
        reprogramming_intensity: 0.15,
        tumor_suppressor_fold: if safeguards.super_p53_active { 1.8 } else { 1.0 },
        offtarget_cleavage_burden: proto_oncogene_risk,
        has_suicide_switch: safeguards.inducible_suicide_switch,
    });

    // 4. Checkpoint loss risk: baseline tissue risk modulated by telomerase stress
    let checkpoint_loss = (0.08 + telomerase_expression * 0.20).min(0.75);

    // 5. Chromosomal instability from defective shelterin capping (dicentric bridge formation)
    let instability_risk = ((1.0 - shelterin_capping) * 0.65).clamp(0.02, 0.90);

    // Weighted raw total combining pathway hazards and formal oncogene drive
    let total_unmitigated = (0.30 * crisis_bypass
        + 0.20 * proto_oncogene_risk
        + 0.15 * formal_onco.proliferative_drive_score
        + 0.15 * checkpoint_loss
        + 0.20 * instability_risk)
        .clamp(0.01, 0.99);

    let raw_breakdown = PathwayRiskBreakdown {
        crisis_bypass_risk: crisis_bypass,
        proto_oncogene_transactivation_risk: proto_oncogene_risk,
        checkpoint_loss_risk: checkpoint_loss,
        chromosomal_instability_risk: instability_risk,
        total_unmitigated_risk: total_unmitigated,
    };

    let mitigation = safeguards.cumulative_mitigation_efficiency();
    let net_risk = (total_unmitigated * (1.0 - mitigation)).clamp(0.005, 0.95);
    let is_acceptable = net_risk < 0.15;

    OncogenicRiskAssessment {
        raw_pathways: raw_breakdown,
        safeguards: safeguards.clone(),
        mitigation_factor: mitigation,
        net_risk_index: net_risk,
        is_acceptable_safety_profile: is_acceptable,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_safety_locks_substantially_reduce_risk() {
        let offtargets = vec![0.15, 0.22];
        let unshielded = evaluate_oncogenic_risk(0.85, &offtargets, 0.70, &SafetyLocks::none());
        assert!(unshielded.net_risk_index > 0.25);
        assert!(!unshielded.is_acceptable_safety_profile);

        let shielded = evaluate_oncogenic_risk(0.85, &offtargets, 0.70, &SafetyLocks::standard_protection());
        assert!(shielded.net_risk_index < 0.10);
        assert!(shielded.is_acceptable_safety_profile);
        assert!(shielded.mitigation_factor > 0.80);
    }
}
