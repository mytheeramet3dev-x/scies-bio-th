//! # Clinical-Grade Reference Anchors & Biomarker Audit Engine (`anatomy::anchors`)
//!
//! Provides human reference intervals for 8 core physiological and longevity biomarkers
//! based on published clinical cohorts (CKD-EPI, Framingham, ACSM VO2max, ADA Diabetes Care,
//! Canela et al. telomeres, and cutaneous biomechanics).

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Standard reference interval for a clinical biomarker.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize))]
pub struct ReferenceRange {
    pub min: f64,
    pub max: f64,
    pub units: &'static str,
    pub description: &'static str,
}

// 1. Renal Glomerular Filtration Rate (eGFR)
pub const REF_EGFR_YOUTHFUL: ReferenceRange = ReferenceRange {
    min: 90.0,
    max: 130.0,
    units: "mL/min/1.73m2",
    description: "CKD-EPI Stage 1 Normal/High Glomerular Filtration Rate",
};

// 2. Resting Systolic Blood Pressure
pub const REF_SYSTOLIC_BP_NORMAL: ReferenceRange = ReferenceRange {
    min: 105.0,
    max: 125.0,
    units: "mmHg",
    description: "AHA/ACC Normotensive Resting Arterial Blood Pressure",
};

// 3. Cardiopulmonary Maximal Oxygen Consumption (VO2max)
pub const REF_VO2MAX_YOUNG_ADULT: ReferenceRange = ReferenceRange {
    min: 40.0,
    max: 55.0,
    units: "mL/kg/min",
    description: "ACSM Fitness Norms for Healthy Young Adults",
};

// 4. Fasting & 2h Post-Challenge Blood Glucose
pub const REF_FASTING_GLUCOSE_NORMAL: ReferenceRange = ReferenceRange {
    min: 70.0,
    max: 99.0,
    units: "mg/dL",
    description: "ADA Clinical Practice Guidelines Euglycemic Fasting Glucose",
};

pub const REF_OGTT_2H_GLUCOSE_NORMAL: ReferenceRange = ReferenceRange {
    min: 70.0,
    max: 139.0,
    units: "mg/dL",
    description: "ADA Normal 2-Hour 75g Oral Glucose Tolerance Test",
};

// 5. Age-Stratified Telomere Length Bands
pub const REF_TELOMERE_YOUNG_ADULT: ReferenceRange = ReferenceRange {
    min: 8_500.0,
    max: 11_000.0,
    units: "bp",
    description: "Flow-FISH Peripheral Blood Mononuclear Cell Telomere Length (Age 18-30)",
};

pub const REF_TELOMERE_MIDDLE_AGE: ReferenceRange = ReferenceRange {
    min: 6_500.0,
    max: 8_500.0,
    units: "bp",
    description: "Flow-FISH PBMC Telomere Length (Age 40-60)",
};

pub const REF_TELOMERE_AGED: ReferenceRange = ReferenceRange {
    min: 4_500.0,
    max: 6_500.0,
    units: "bp",
    description: "Flow-FISH PBMC Telomere Length (Age 70+)",
};

// 6. Systemic Inflammaging / CRP / SASP Index
pub const REF_SASP_YOUTHFUL: ReferenceRange = ReferenceRange {
    min: 0.01,
    max: 0.15,
    units: "index",
    description: "Youthful Low Systemic Circulating Senescence-Associated Secretory Phenotype",
};

// 7. Daily Thymic Naive T-Cell Emigration
pub const REF_THYMIC_NAIVE_T_OUTPUT_YOUTH: ReferenceRange = ReferenceRange {
    min: 8.0,
    max: 15.0,
    units: "10^7 cells/day",
    description: "Young Adult Pre-Involution Thymic Emigration Rate",
};

// 8. Cutaneous Young's Elastic Modulus Stiffness
pub const REF_SKIN_STIFFNESS_YOUTH: ReferenceRange = ReferenceRange {
    min: 8.0,
    max: 20.0,
    units: "kPa",
    description: "Cutometer / Indentation Young's Modulus of Youthful Dermis",
};

/// Set of 8 clinical biomarker values for a subject to be audited against reference ranges.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ClinicalBiomarkerValues {
    pub egfr_ml_min: f64,
    pub systolic_bp_mmhg: f64,
    pub vo2max_ml_kg_min: f64,
    pub fasting_glucose_mg_dl: f64,
    pub ogtt_2h_glucose_mg_dl: f64,
    pub telomere_length_bp: f64,
    pub sasp_crp_proxy: f64,
    pub thymic_naive_t_output_10e7: f64,
    pub skin_stiffness_kpa: f64,
}

/// Evaluation record for an individual biomarker against its clinical anchor range.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ClinicalAnchorCheck {
    pub biomarker_name: String,
    pub observed_value: f64,
    pub reference_min: f64,
    pub reference_max: f64,
    pub units: String,
    pub description: String,
    pub in_range: bool,
}

/// Aggregate report produced by [`audit_clinical_reference_anchors`].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ClinicalReferenceAuditReport {
    pub checks: Vec<ClinicalAnchorCheck>,
    pub passed_count: usize,
    pub total_count: usize,
    pub all_anchors_compliant: bool,
}

/// Audits 8 core physiological and longevity biomarkers against published human clinical reference ranges.
pub fn audit_clinical_reference_anchors(
    values: &ClinicalBiomarkerValues,
    chronological_age: f64,
) -> ClinicalReferenceAuditReport {
    let telo_ref = if chronological_age <= 35.0 {
        &REF_TELOMERE_YOUNG_ADULT
    } else if chronological_age <= 65.0 {
        &REF_TELOMERE_MIDDLE_AGE
    } else {
        &REF_TELOMERE_AGED
    };

    let make_check = |name: &str, val: f64, ref_range: &ReferenceRange| -> ClinicalAnchorCheck {
        ClinicalAnchorCheck {
            biomarker_name: name.to_string(),
            observed_value: val,
            reference_min: ref_range.min,
            reference_max: ref_range.max,
            units: ref_range.units.to_string(),
            description: ref_range.description.to_string(),
            in_range: val >= ref_range.min && val <= ref_range.max,
        }
    };

    let checks = vec![
        make_check("eGFR", values.egfr_ml_min, &REF_EGFR_YOUTHFUL),
        make_check("Systolic BP", values.systolic_bp_mmhg, &REF_SYSTOLIC_BP_NORMAL),
        make_check("VO2max", values.vo2max_ml_kg_min, &REF_VO2MAX_YOUNG_ADULT),
        make_check("Fasting Glucose", values.fasting_glucose_mg_dl, &REF_FASTING_GLUCOSE_NORMAL),
        make_check("2h Glucose (OGTT)", values.ogtt_2h_glucose_mg_dl, &REF_OGTT_2H_GLUCOSE_NORMAL),
        make_check("Telomere Length", values.telomere_length_bp, telo_ref),
        make_check("SASP / CRP Index", values.sasp_crp_proxy, &REF_SASP_YOUTHFUL),
        make_check("Thymic Naive T-Cell Output", values.thymic_naive_t_output_10e7, &REF_THYMIC_NAIVE_T_OUTPUT_YOUTH),
        make_check("Skin Stiffness (E)", values.skin_stiffness_kpa, &REF_SKIN_STIFFNESS_YOUTH),
    ];

    let passed_count = checks.iter().filter(|c| c.in_range).count();
    let total_count = checks.len();

    ClinicalReferenceAuditReport {
        checks,
        passed_count,
        total_count,
        all_anchors_compliant: passed_count == total_count,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_youthful_biomarkers_comply_with_clinical_anchors() {
        let youth = ClinicalBiomarkerValues {
            egfr_ml_min: 110.0,
            systolic_bp_mmhg: 118.0,
            vo2max_ml_kg_min: 48.0,
            fasting_glucose_mg_dl: 85.0,
            ogtt_2h_glucose_mg_dl: 110.0,
            telomere_length_bp: 9_600.0,
            sasp_crp_proxy: 0.05,
            thymic_naive_t_output_10e7: 12.0,
            skin_stiffness_kpa: 14.0,
        };

        let audit = audit_clinical_reference_anchors(&youth, 25.0);
        assert_eq!(audit.passed_count, audit.total_count);
        assert!(audit.all_anchors_compliant);
    }

    #[test]
    fn test_aged_biomarkers_flag_out_of_bounds_clinical_anchors() {
        let aged = ClinicalBiomarkerValues {
            egfr_ml_min: 52.0, // Low eGFR
            systolic_bp_mmhg: 154.0, // Systolic hypertension
            vo2max_ml_kg_min: 22.0, // Low aerobic capacity
            fasting_glucose_mg_dl: 125.0, // Impaired fasting glucose
            ogtt_2h_glucose_mg_dl: 180.0, // Impaired glucose tolerance
            telomere_length_bp: 5_200.0, // Appropriate for age 75
            sasp_crp_proxy: 0.45, // Elevated inflammaging
            thymic_naive_t_output_10e7: 1.5, // Severe thymic involution
            skin_stiffness_kpa: 48.0, // Advanced glycation crosslinked stiff skin
        };

        let audit = audit_clinical_reference_anchors(&aged, 75.0);
        // Only telomere length should pass the age-stratified 70+ band; others must flag as out-of-bounds
        let telomere_check = audit.checks.iter().find(|c| c.biomarker_name == "Telomere Length").unwrap();
        assert!(telomere_check.in_range);
        assert!(audit.passed_count < audit.total_count);
        assert!(!audit.all_anchors_compliant);
    }
}
