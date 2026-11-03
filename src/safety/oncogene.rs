//! Proto-oncogene activation and neoplastic transformation risk quantification.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Gene-expression and editing inputs for oncogenic risk scoring.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct OncogeneInterventionInput {
    /// TERT expression fold change ($1.0 = \text{baseline somatic}$).
    pub tert_fold: f64,
    /// MYC expression fold change ($1.0 = \text{baseline}$).
    pub myc_fold: f64,
    /// Pluripotency reprogramming duration or intensity ($0.0\text{--}1.0$).
    pub reprogramming_intensity: f64,
    /// Tumor suppressor (TP53 / CDKN2A / PTEN) aggregate fold protection ($1.0 = \text{wild-type}$).
    pub tumor_suppressor_fold: f64,
    /// Cumulative CRISPR off-target cleavage risk ($0.0\text{--}1.0$).
    pub offtarget_cleavage_burden: f64,
    /// Whether an inducible suicide gene (e.g. iCasp9) is integrated.
    pub has_suicide_switch: bool,
}

/// Quantitative oncogenic transformation risk breakdown.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct OncogenicRiskAssessment {
    /// Unmitigated proliferative oncogene drive in $[0.0, 1.0]$.
    pub proliferative_drive_score: f64,
    /// Teratoma / dysplasia dedifferentiation hazard in $[0.0, 1.0]$.
    pub dedifferentiation_hazard: f64,
    /// Insertional / off-target mutagenesis hazard in $[0.0, 1.0]$.
    pub mutagenesis_hazard: f64,
    /// Tumor-suppressor & fail-safe attenuation factor in $[0.05, 1.0]$.
    pub mitigation_multiplier: f64,
    /// Final composite oncogenic risk score in $[0.0, 1.0]$ ($< 0.15$ considered safe).
    pub composite_oncogenic_risk: f64,
}

/// Evaluates proto-oncogene activation and neoplastic transformation risk.
pub fn evaluate_oncogenic_activation_risk(
    input: &OncogeneInterventionInput,
) -> OncogenicRiskAssessment {
    let tert_excess = ((input.tert_fold - 1.0).max(0.0) / 4.0).clamp(0.0, 1.0);
    let myc_excess = ((input.myc_fold - 1.0).max(0.0) / 1.5).clamp(0.0, 1.0);
    let proliferative_drive_score = (0.45 * tert_excess + 0.55 * myc_excess).clamp(0.0, 1.0);

    // Continuous reprogramming (> 0.35) sharply increases teratoma risk
    let dedifferentiation_hazard = if input.reprogramming_intensity > 0.25 {
        ((input.reprogramming_intensity - 0.25) / 0.50).powi(2).clamp(0.0, 1.0)
    } else {
        0.05 * input.reprogramming_intensity.max(0.0)
    };

    let mutagenesis_hazard = input.offtarget_cleavage_burden.clamp(0.0, 1.0);

    let raw_hazard = (0.50 * proliferative_drive_score
        + 0.30 * dedifferentiation_hazard
        + 0.20 * mutagenesis_hazard)
        .clamp(0.0, 1.0);

    // Tumor suppressors (TP53 > 1.0) and iCasp9 suicide switches attenuate composite risk
    let ts_protection = 1.0 / input.tumor_suppressor_fold.max(0.25).powf(1.2);
    let switch_factor = if input.has_suicide_switch { 0.45 } else { 1.0 };
    let mitigation_multiplier = (ts_protection * switch_factor).clamp(0.05, 2.0);

    let composite_oncogenic_risk = (raw_hazard * mitigation_multiplier).clamp(0.0, 1.0);

    OncogenicRiskAssessment {
        proliferative_drive_score,
        dedifferentiation_hazard,
        mutagenesis_hazard,
        mitigation_multiplier,
        composite_oncogenic_risk,
    }
}
