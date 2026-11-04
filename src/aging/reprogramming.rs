//! In vivo partial epigenetic reprogramming engine (pulsed OSK without dedifferentiation)
//! coupled with transient *TERT* mRNA telomere restoration and *TP53*/*CDKN2A* safety locks.
//!
//! Models cyclical or transient in vivo induction of pluripotency transcription factors
//! (**Oct4, Sox2, Klf4 — OSK**, omitting oncogenic **c-Myc**) to reset biological age,
//! restore histone marks, and clear cellular senescence while preserving somatic cell identity
//! and preventing teratoma formation. Also models optional transient *TERT* mRNA pulsing
//! gated by [`crate::safety::oncogene::evaluate_oncogenic_activation_risk`].
//!
//! References:
//! - Ocampo et al. (2016) *Cell* 167:1719–1733 (Cyclic partial reprogramming in vivo extends lifespan).
//! - Lu et al. (2020) *Nature* 588:124–129 (TET1/2-dependent vision restoration by viral OSK).
//! - Gill et al. (2022) *eLife* 11:e71624 (Maturation phase transient reprogramming resets human epigenetic clock).
//! - Ramunas et al. (2015) *FASEB J* 29:1930–1939 (Transient delivery of modified TERT mRNA extends human telomeres).
//! - Alle et al. (2024) *Nature Communications* 15:1022 (Continuous vs pulsed OSK dynamics and dedifferentiation barriers).

use crate::safety::{GuardAction, InterventionSafetyCoupling};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Maximum safe composite oncogenic risk threshold defined by `safety::oncogene`.
pub const MAX_SAFE_ONCOGENIC_RISK_THRESHOLD: f64 = 0.15;

/// Delivery vector platform for OSK transcription factors.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DeliveryVector {
    /// Adeno-Associated Virus (AAV) with Doxycycline-regulatable Tet-On/Off promoter.
    AavDoxInducible,
    /// Modified mRNA encapsulated in lipid nanoparticles (non-integrating, self-limiting half-life).
    MrnaLipidNanoparticle,
    /// Lentivirus equipped with synthetic suicide locks / caspase switches.
    LentiviralToxLock,
}

/// Safety verdict for an in vivo partial reprogramming regimen.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SafetyVerdict {
    /// Optimal regimen: robust epigenetic rejuvenation with somatic identity intact (>90%).
    OptimalRejuvenation,
    /// Suboptimal: washout interval is too short (<3 days), allowing factor accumulation.
    SuboptimalInsufficientWashout,
    /// Lethal risk: continuous pulse exceeds point-of-no-return (>4.5 days), inducing teratoma/dedifferentiation.
    SevereTeratomaRiskExceededDedifferentiationThreshold,
    /// Ineffective: TET1/TET2 methylcytosine dioxygenase activity is insufficient (<0.35).
    TetDeficiencyIneffectiveDemethylation,
    /// Carcinogenic hazard: oncogenic c-Myc was included in the reprogramming factor cocktail.
    CarcinogenicMycRisk,
    /// Oncogenic hazard: TERT expression fold exceeded safe threshold without sufficient TP53/CDKN2A protection.
    OncogenicTelomeraseThresholdExceeded,
}

/// Configuration regimen for cyclical/pulsed in vivo OSK partial reprogramming and transient *TERT* mRNA pulsing.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OskPulseRegimen {
    /// Somatic tissue or lineage name (e.g. "Fibroblast", "Cardiomyocyte", "RetinalGanglion").
    pub somatic_lineage: String,
    /// Gene delivery vehicle.
    pub delivery_vector: DeliveryVector,
    /// Duration of active factor expression per cycle in days (recommended: 1.5 - 2.5 days).
    pub pulse_duration_days: f64,
    /// Quiescent washout/recovery interval between pulses in days (recommended: 4.0 - 7.0 days).
    pub washout_interval_days: f64,
    /// Total number of consecutive pulsed cycles.
    pub total_cycles: usize,
    /// Functional competence of TET1/TET2 methylcytosine dioxygenases (0.0 to 1.0).
    pub tet1_tet2_competence: f64,
    /// Whether oncogenic c-Myc is included in the factor cocktail (True = lethal oncogenic hazard).
    pub oncogene_cmyc_included: bool,
    /// Relative factor expression dosage (1.0 = standard physiological reprogramming induction).
    pub factor_dosage: f64,
    /// Transient *TERT* mRNA expression fold pulse (1.0 = baseline inactive; 2.0-2.5 = transient telomerase pulse).
    pub transient_tert_mrna_fold: f64,
    /// Tumor suppressor (*TP53* / *CDKN2A*-iCasp9) safeguard fold expression (1.0 = wild-type; >=1.5 = reinforced).
    pub tp53_cdkn2a_safety_fold: f64,
}

impl Default for OskPulseRegimen {
    fn default() -> Self {
        Self {
            somatic_lineage: "SomaticFibroblast".to_string(),
            delivery_vector: DeliveryVector::AavDoxInducible,
            pulse_duration_days: 2.0,
            washout_interval_days: 5.0,
            total_cycles: 6,
            tet1_tet2_competence: 0.90,
            oncogene_cmyc_included: false,
            factor_dosage: 1.0,
            transient_tert_mrna_fold: 2.2,
            tp53_cdkn2a_safety_fold: 1.8,
        }
    }
}

/// Longitudinal metric captured across an individual reprogramming cycle.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CycleMetric {
    /// 1-based cycle index.
    pub cycle: usize,
    /// Retention of lineage-specific enhancer epigenetic bookmarks (1.0 = 100% somatic identity).
    pub somatic_identity_retention: f64,
    /// Rejuvenation achieved in this specific cycle (years reversed).
    pub cycle_rejuvenation_years: f64,
    /// Cumulative biological years reversed since start of regimen.
    pub cumulative_rejuvenation_years: f64,
    /// Fraction of nuclear DNA damage / gamma-H2AX foci cleared (0.0 to 1.0).
    pub dna_damage_cleared: f64,
    /// Instantaneous teratoma / dedifferentiation hazard score (0.0 to 100.0).
    pub teratoma_hazard_score: f64,
    /// True if identity and hazard metrics remained within safe physiological boundaries.
    pub in_safe_operating_window: bool,
}

/// Comprehensive report output after simulating an in vivo partial reprogramming regimen.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InVivoReprogrammingReport {
    /// Somatic lineage analyzed.
    pub somatic_lineage: String,
    /// Vector platform used.
    pub delivery_vector: DeliveryVector,
    /// Initial baseline epigenetic clock age in years.
    pub initial_epigenetic_age: f64,
    /// Final predicted epigenetic clock age after completing the regimen.
    pub final_epigenetic_age: f64,
    /// Total biological age rejuvenation in years.
    pub total_rejuvenation_years: f64,
    /// Final retention fraction of somatic identity bookmarks ($I_{\text{somatic}} \in [0.0, 1.0]$).
    pub final_somatic_identity_retention: f64,
    /// Overall oncogenic/teratoma hazard index ($0.0$ to $1.0$, where $< 0.10$ is safe).
    pub teratoma_hazard_index: f64,
    /// Composite oncogenic risk score evaluated via `safety::oncogene` ($< 0.15$ required).
    pub composite_oncogenic_risk: f64,
    /// Overall safety and efficacy verdict.
    pub safety_verdict: SafetyVerdict,
    /// Longitudinal telemetry per cycle.
    pub per_cycle_metrics: Vec<CycleMetric>,
    /// Estimated stem cell niche fraction ($0.0$ to $1.0$, where $\le 0.15$ is non-hyperplastic).
    pub stem_pool_fraction: f64,
    /// Recommended autonomous closed-loop safety guard action.
    pub recommended_guard_action: GuardAction,
    /// Bioengineering recommendations and safety alerts.
    pub recommendations: Vec<String>,
}

impl InVivoReprogrammingReport {
    /// Extracts the standardized intervention safety coupling payload.
    pub fn safety_coupling(&self) -> InterventionSafetyCoupling {
        InterventionSafetyCoupling {
            oncogenic_risk: self.composite_oncogenic_risk,
            stem_pool_fraction: self.stem_pool_fraction,
            identity_retention: self.final_somatic_identity_retention,
            recommended_guard_action: self.recommended_guard_action,
        }
    }
}

/// Central biological formula: maps lineage integrity and chromatin stability to identity competence ($I_{\text{comp}} \in [0.0, 1.0]$).
///
/// Under cellular identity theory, somatic lineage bookmarks (enhancer-promoter loops) act as the physical anchor
/// for productive chromatin remodeling. When identity drops below the critical commitment threshold (0.70),
/// cells enter dedifferentiated / neoplastic states where factor binding causes aberrant transcriptional noise,
/// causing identity competence and functional rejuvenation efficacy to plummet toward 0.
pub fn compute_identity_competence(lineage_integrity: f64, chromatin_stability: f64) -> f64 {
    let combined_identity = (lineage_integrity.clamp(0.0, 1.0) * chromatin_stability.clamp(0.0, 1.0)).sqrt();
    let threshold = 0.70;
    if combined_identity < threshold {
        ((combined_identity / threshold).powi(4) * 0.05).max(0.0)
    } else {
        let normalized = (combined_identity - threshold) / (1.0 - threshold);
        normalized.powf(1.5).clamp(0.0, 1.0)
    }
}

/// Computes effective productive reprogramming dosage:
/// `productive_dosage = nominal_dosage * identity_competence`
///
/// Gated by saturation kinetics (dosage above 1.25x saturates catalytic TET demethylation while accelerating identity erosion).
pub fn compute_productive_dosage(nominal_dosage: f64, identity_competence: f64) -> f64 {
    nominal_dosage.clamp(0.0, 1.25) * identity_competence.clamp(0.0, 1.0)
}

/// Simulates an in vivo partial reprogramming regimen starting with an explicit initial somatic identity fraction.
pub fn simulate_in_vivo_osk_pulse_regimen_with_identity(
    initial_epigenetic_age: f64,
    initial_somatic_identity: f64,
    regimen: &OskPulseRegimen,
) -> InVivoReprogrammingReport {
    let mut current_age = initial_epigenetic_age;
    let mut somatic_identity = initial_somatic_identity.clamp(0.0, 1.0);
    let mut cumulative_rejuvenation = 0.0f64;
    let mut cycle_metrics = Vec::with_capacity(regimen.total_cycles);
    let mut recommendations = Vec::new();

    // Thresholds (based on Ocampo 2016, Lu 2020, Gill 2022)
    let safe_pulse_max_days = 4.0;
    let point_of_no_return_days = 5.5;
    let min_washout_days = 3.0;

    let mut peak_hazard = 0.0f64;

    // Evaluate proto-oncogene activation & TERT safety lock via `safety::oncogene`
    #[cfg(feature = "safety")]
    let composite_oncogenic_risk = {
        let has_suicide_switch = regimen.tp53_cdkn2a_safety_fold >= 1.4
            || regimen.delivery_vector == DeliveryVector::LentiviralToxLock;
        let onc_input = crate::safety::oncogene::OncogeneInterventionInput {
            tert_fold: regimen.transient_tert_mrna_fold.max(1.0),
            myc_fold: if regimen.oncogene_cmyc_included { 2.8 } else { 1.0 },
            reprogramming_intensity: (regimen.pulse_duration_days / 10.0).clamp(0.0, 1.0),
            tumor_suppressor_fold: regimen.tp53_cdkn2a_safety_fold.max(0.1),
            offtarget_cleavage_burden: 0.01,
            has_suicide_switch,
        };
        crate::safety::oncogene::evaluate_oncogenic_activation_risk(&onc_input)
            .composite_oncogenic_risk
    };
    #[cfg(not(feature = "safety"))]
    let composite_oncogenic_risk = 0.02;

    for cycle in 1..=regimen.total_cycles {
        // 1. Identity erosion: short pulses (<= 2.5 d) at physiological dosage (<= 1.25x) minimally
        // stress enhancer memory. Above 4.0 d or 1.25x dosage, pioneer factor invasion erodes bookmarks.
        let pulse_overrun = (regimen.pulse_duration_days - safe_pulse_max_days).max(0.0);
        let dose_overrun = (regimen.factor_dosage - 1.25).max(0.0);
        let base_erosion = 0.015 * (regimen.pulse_duration_days / 2.0) * regimen.factor_dosage;
        let overrun_erosion =
            (0.045 * pulse_overrun + 0.08 * pulse_overrun.powi(2) + 0.035 * dose_overrun)
                * regimen.factor_dosage;
        let cycle_erosion = base_erosion + overrun_erosion;

        // Washout recovery: somatic enhancers re-bind lineage factors during resting interval
        let washout_recovery = if regimen.washout_interval_days >= min_washout_days {
            (0.012 * (regimen.washout_interval_days / 5.0)).min(cycle_erosion * 0.90)
        } else {
            0.002 // Insufficient washout -> unrecovered identity stress
        };

        if somatic_identity > 0.0 {
            somatic_identity = (somatic_identity - cycle_erosion + washout_recovery).clamp(0.0, 1.0);
        }

        // 2. Epigenetic rejuvenation: active TET1/2 demethylation of age-hypermethylated CpGs.
        let tet_factor = (regimen.tet1_tet2_competence - 0.20).max(0.0) / 0.80;
        let chromatin_stability = somatic_identity.clamp(0.0, 1.0);
        let identity_competence = compute_identity_competence(somatic_identity, chromatin_stability);
        let productive_dosage = compute_productive_dosage(regimen.factor_dosage, identity_competence);

        let oncogenic_competence = if regimen.oncogene_cmyc_included {
            0.15
        } else if composite_oncogenic_risk > MAX_SAFE_ONCOGENIC_RISK_THRESHOLD {
            (1.0 - (composite_oncogenic_risk - MAX_SAFE_ONCOGENIC_RISK_THRESHOLD) * 1.8)
                .clamp(0.15, 1.0)
        } else {
            1.0
        };
        let base_rejuv_per_cycle = 2.8
            * (regimen.pulse_duration_days.min(3.5) / 2.0)
            * tet_factor
            * productive_dosage
            * oncogenic_competence;

        // Gompertz-Horvath Epigenetic Reversal Floor
        let bio_age_floor = if initial_epigenetic_age <= 30.0 {
            18.0
        } else if initial_epigenetic_age >= 85.0 {
            55.0 + (initial_epigenetic_age - 85.0) * 0.40
        } else {
            20.0 + (initial_epigenetic_age - 30.0) * 0.20
        };

        // Diminishing returns as age approaches the donor-specific biological floor
        let age_headroom = (current_age - bio_age_floor).max(0.0) / 45.0;
        let raw_cycle_rejuv = (base_rejuv_per_cycle * age_headroom.min(1.2)).max(0.0);

        // When somatic identity erodes below the lineage commitment threshold (0.85), previously
        // accumulated functional epigenetic youth degrades as somatic chromatin bookmarks collapse.
        let lineage_integrity_retention = if somatic_identity >= 0.85 {
            1.0
        } else if somatic_identity <= 0.50 {
            0.0
        } else {
            ((somatic_identity - 0.50) / 0.35).clamp(0.0, 1.0).powi(2)
        };
        let target_cumulative =
            (cumulative_rejuvenation + raw_cycle_rejuv) * lineage_integrity_retention;
        let cycle_rejuv = target_cumulative - cumulative_rejuvenation;
        cumulative_rejuvenation = target_cumulative.max(0.0);
        current_age =
            (initial_epigenetic_age - cumulative_rejuvenation).clamp(bio_age_floor, initial_epigenetic_age);

        // 3. Nuclear DNA damage repair
        let dna_damage_cleared = (0.18 * cycle as f64 * tet_factor * identity_competence).min(0.85);

        // 4. Teratoma / dedifferentiation hazard calculation
        let mut cycle_hazard = 0.0f64;
        if regimen.pulse_duration_days > safe_pulse_max_days {
            let severity = (regimen.pulse_duration_days - safe_pulse_max_days) / 1.5;
            cycle_hazard += severity * 40.0;
        }
        if somatic_identity < 0.80 {
            cycle_hazard += (0.80 - somatic_identity) * 120.0;
        }
        if regimen.washout_interval_days < min_washout_days {
            cycle_hazard += (min_washout_days - regimen.washout_interval_days) * 15.0;
        }
        if regimen.oncogene_cmyc_included {
            cycle_hazard += 65.0;
        }
        cycle_hazard = cycle_hazard.clamp(0.0, 100.0);
        if cycle_hazard > peak_hazard {
            peak_hazard = cycle_hazard;
        }

        let in_safe_window = somatic_identity >= 0.75
            && cycle_hazard <= 20.0
            && !regimen.oncogene_cmyc_included
            && composite_oncogenic_risk < MAX_SAFE_ONCOGENIC_RISK_THRESHOLD;

        cycle_metrics.push(CycleMetric {
            cycle,
            somatic_identity_retention: somatic_identity,
            cycle_rejuvenation_years: cycle_rejuv,
            cumulative_rejuvenation_years: cumulative_rejuvenation,
            dna_damage_cleared,
            teratoma_hazard_score: cycle_hazard,
            in_safe_operating_window: in_safe_window,
        });
    }

    // Determine overall safety verdict
    let teratoma_hazard_index = (peak_hazard / 100.0).clamp(0.0, 1.0);
    let safety_verdict = if regimen.oncogene_cmyc_included {
        recommendations.push("CRITICAL ALERT: c-Myc must be strictly excluded from in vivo reprogramming to prevent rapid fibrosarcoma/teratoma induction (use OSK or OSK+Lin28 only).".to_string());
        SafetyVerdict::CarcinogenicMycRisk
    } else if composite_oncogenic_risk >= MAX_SAFE_ONCOGENIC_RISK_THRESHOLD {
        recommendations.push(format!(
            "ONCOGENIC LOCKOUT: Composite oncogenic risk ({:.3}) exceeded safe threshold ({:.2}). Reduce transient TERT fold or reinforce TP53/CDKN2A safeguard.",
            composite_oncogenic_risk, MAX_SAFE_ONCOGENIC_RISK_THRESHOLD
        ));
        SafetyVerdict::OncogenicTelomeraseThresholdExceeded
    } else if regimen.pulse_duration_days >= point_of_no_return_days || somatic_identity < 0.70 {
        recommendations.push(format!(
            "DANGER: Pulse duration of {:.1} days crossed the point-of-no-return threshold ({:.1} days). Somatic identity dropped to {:.1}%. Reduce pulse length to <= 2.5 days.",
            regimen.pulse_duration_days, point_of_no_return_days, somatic_identity * 100.0
        ));
        SafetyVerdict::SevereTeratomaRiskExceededDedifferentiationThreshold
    } else if regimen.tet1_tet2_competence < 0.35 {
        recommendations.push("INEFFECTIVE: Insufficient TET1/2 catalytic competence. Active demethylation is stalled; co-deliver alpha-ketoglutarate or TET1 catalytic domain.".to_string());
        SafetyVerdict::TetDeficiencyIneffectiveDemethylation
    } else if regimen.washout_interval_days < min_washout_days {
        recommendations.push(format!(
            "CAUTION: Washout interval ({:.1} days) is shorter than recommended minimum ({:.1} days). Factors will accumulate across cycles.",
            regimen.washout_interval_days, min_washout_days
        ));
        SafetyVerdict::SuboptimalInsufficientWashout
    } else {
        recommendations.push(format!(
            "OPTIMAL: Regimen safely reversed {:.1} years of epigenetic clock age while preserving {:.1}% somatic identity with composite oncogenic risk {:.3} (< {:.2}).",
            cumulative_rejuvenation, somatic_identity * 100.0, composite_oncogenic_risk, MAX_SAFE_ONCOGENIC_RISK_THRESHOLD
        ));
        SafetyVerdict::OptimalRejuvenation
    };

    let stem_pool_fraction = if somatic_identity < 0.70 {
        (0.20 + (0.70 - somatic_identity) * 0.80).clamp(0.05, 0.75)
    } else {
        (0.05 + 0.008 * (regimen.total_cycles as f64) * (regimen.pulse_duration_days / 2.0)).min(0.14)
    };

    let recommended_guard_action = if regimen.oncogene_cmyc_included
        || composite_oncogenic_risk >= MAX_SAFE_ONCOGENIC_RISK_THRESHOLD
        || somatic_identity < 0.70
    {
        crate::safety::GuardAction::AbortSuicideSwitch
    } else if regimen.washout_interval_days < min_washout_days {
        crate::safety::GuardAction::PauseWashout
    } else if composite_oncogenic_risk > 0.08 || somatic_identity < 0.85 || regimen.factor_dosage > 1.20 {
        crate::safety::GuardAction::TitrateDown
    } else {
        crate::safety::GuardAction::Continue
    };

    InVivoReprogrammingReport {
        somatic_lineage: regimen.somatic_lineage.clone(),
        delivery_vector: regimen.delivery_vector,
        initial_epigenetic_age,
        final_epigenetic_age: current_age,
        total_rejuvenation_years: cumulative_rejuvenation,
        final_somatic_identity_retention: somatic_identity,
        teratoma_hazard_index,
        composite_oncogenic_risk,
        safety_verdict,
        per_cycle_metrics: cycle_metrics,
        stem_pool_fraction,
        recommended_guard_action,
        recommendations,
    }
}

/// Simulates an in vivo partial reprogramming regimen over an aged somatic cell population.
///
/// Note: Nuclear OSK epigenetic reprogramming does **not** directly clear mutated mtDNA or
/// restore mitochondrial copy number (which requires `PINK1`/`Parkin` mitophagy).
pub fn simulate_in_vivo_osk_pulse_regimen(
    initial_epigenetic_age: f64,
    regimen: &OskPulseRegimen,
) -> InVivoReprogrammingReport {
    simulate_in_vivo_osk_pulse_regimen_with_identity(initial_epigenetic_age, 1.0, regimen)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_optimal_pulsed_osk_regimen_preserves_identity_and_rejuvenates() {
        let regimen = OskPulseRegimen {
            somatic_lineage: "HumanDermalFibroblast".to_string(),
            delivery_vector: DeliveryVector::AavDoxInducible,
            pulse_duration_days: 2.0,
            washout_interval_days: 5.0,
            total_cycles: 6,
            tet1_tet2_competence: 0.90,
            oncogene_cmyc_included: false,
            factor_dosage: 1.0,
            transient_tert_mrna_fold: 2.2,
            tp53_cdkn2a_safety_fold: 1.8,
        };

        let report = simulate_in_vivo_osk_pulse_regimen(65.0, &regimen);
        assert_eq!(report.safety_verdict, SafetyVerdict::OptimalRejuvenation);
        assert!(report.total_rejuvenation_years >= 10.0);
        assert!(report.final_epigenetic_age < 55.0);
        assert!(report.final_somatic_identity_retention > 0.90);
        assert!(report.teratoma_hazard_index < 0.10);
        assert!(report.composite_oncogenic_risk < MAX_SAFE_ONCOGENIC_RISK_THRESHOLD);
        assert_eq!(report.per_cycle_metrics.len(), 6);
        for m in &report.per_cycle_metrics {
            assert!(m.in_safe_operating_window);
        }
    }

    #[test]
    fn test_continuous_over_induction_triggers_teratoma_lock() {
        let regimen = OskPulseRegimen {
            somatic_lineage: "AdultCardiomyocyte".to_string(),
            delivery_vector: DeliveryVector::AavDoxInducible,
            pulse_duration_days: 6.5, // Exceeds 5.5 day point-of-no-return
            washout_interval_days: 2.0,
            total_cycles: 4,
            tet1_tet2_competence: 0.90,
            oncogene_cmyc_included: false,
            factor_dosage: 1.0,
            transient_tert_mrna_fold: 1.5,
            tp53_cdkn2a_safety_fold: 1.8,
        };

        let report = simulate_in_vivo_osk_pulse_regimen(70.0, &regimen);
        assert_eq!(
            report.safety_verdict,
            SafetyVerdict::SevereTeratomaRiskExceededDedifferentiationThreshold
        );
        assert!(report.final_somatic_identity_retention < 0.70);
        assert!(report.teratoma_hazard_index > 0.40);
        // Functional epigenetic rejuvenation must collapse when somatic identity collapses
        assert!(report.total_rejuvenation_years < 1.0);
    }

    #[test]
    fn test_identity_loss_directly_penalizes_functional_rejuvenation() {
        let base = OskPulseRegimen {
            somatic_lineage: "Fibroblast".to_string(),
            delivery_vector: DeliveryVector::MrnaLipidNanoparticle,
            pulse_duration_days: 2.0,
            washout_interval_days: 5.0,
            total_cycles: 6,
            tet1_tet2_competence: 0.95,
            oncogene_cmyc_included: false,
            factor_dosage: 1.05,
            transient_tert_mrna_fold: 2.2,
            tp53_cdkn2a_safety_fold: 1.8,
        };
        let mut moderate_overdose = base.clone();
        moderate_overdose.factor_dosage = 1.65;
        moderate_overdose.pulse_duration_days = 3.2;

        let mut severe_overdose = base.clone();
        severe_overdose.factor_dosage = 2.50;
        severe_overdose.pulse_duration_days = 4.8;

        let r_base = simulate_in_vivo_osk_pulse_regimen(72.0, &base);
        let r_mod = simulate_in_vivo_osk_pulse_regimen(72.0, &moderate_overdose);
        let r_sev = simulate_in_vivo_osk_pulse_regimen(72.0, &severe_overdose);

        assert!(r_base.final_somatic_identity_retention > r_mod.final_somatic_identity_retention);
        assert!(r_mod.final_somatic_identity_retention > r_sev.final_somatic_identity_retention);
        assert!(r_base.total_rejuvenation_years > r_mod.total_rejuvenation_years);
        assert!(r_mod.total_rejuvenation_years > r_sev.total_rejuvenation_years);
        assert_eq!(r_sev.final_somatic_identity_retention, 0.0);
        assert_eq!(r_sev.total_rejuvenation_years, 0.0);
    }

    #[test]
    fn test_cmyc_inclusion_triggers_carcinogenic_myc_risk() {
        let regimen = OskPulseRegimen {
            somatic_lineage: "NeuralStemCell".to_string(),
            delivery_vector: DeliveryVector::LentiviralToxLock,
            pulse_duration_days: 2.0,
            washout_interval_days: 5.0,
            total_cycles: 4,
            tet1_tet2_competence: 0.90,
            oncogene_cmyc_included: true, // Lethal oncogene
            factor_dosage: 1.0,
            transient_tert_mrna_fold: 1.5,
            tp53_cdkn2a_safety_fold: 1.5,
        };

        let report = simulate_in_vivo_osk_pulse_regimen(50.0, &regimen);
        assert_eq!(report.safety_verdict, SafetyVerdict::CarcinogenicMycRisk);
        assert!(report.teratoma_hazard_index >= 0.65);
    }

    #[test]
    fn test_unguarded_tert_overexpression_triggers_oncogenic_lockout() {
        let regimen = OskPulseRegimen {
            somatic_lineage: "Fibroblast".to_string(),
            delivery_vector: DeliveryVector::MrnaLipidNanoparticle,
            pulse_duration_days: 2.0,
            washout_interval_days: 5.0,
            total_cycles: 4,
            tet1_tet2_competence: 0.90,
            oncogene_cmyc_included: false,
            factor_dosage: 1.0,
            transient_tert_mrna_fold: 4.8, // Excessive TERT overexpression
            tp53_cdkn2a_safety_fold: 0.5,  // Compromised p53/CDKN2A safeguard
        };

        let report = simulate_in_vivo_osk_pulse_regimen(60.0, &regimen);
        assert_eq!(
            report.safety_verdict,
            SafetyVerdict::OncogenicTelomeraseThresholdExceeded
        );
        assert!(report.composite_oncogenic_risk >= MAX_SAFE_ONCOGENIC_RISK_THRESHOLD);
    }

    #[test]
    fn test_tet_deficiency_fails_epigenetic_reset() {
        let regimen = OskPulseRegimen {
            somatic_lineage: "Hepatocyte".to_string(),
            delivery_vector: DeliveryVector::MrnaLipidNanoparticle,
            pulse_duration_days: 2.0,
            washout_interval_days: 5.0,
            total_cycles: 5,
            tet1_tet2_competence: 0.15, // Severe TET deficiency
            oncogene_cmyc_included: false,
            factor_dosage: 1.0,
            transient_tert_mrna_fold: 1.0,
            tp53_cdkn2a_safety_fold: 1.5,
        };

        let report = simulate_in_vivo_osk_pulse_regimen(60.0, &regimen);
        assert_eq!(
            report.safety_verdict,
            SafetyVerdict::TetDeficiencyIneffectiveDemethylation
        );
        // Epigenetic clock should not move significantly without TET
        assert!(report.total_rejuvenation_years < 2.0);
    }

    #[test]
    fn test_identity_loss_strictly_binds_rejuvenation_efficacy_monotonically() {
        let base = OskPulseRegimen {
            somatic_lineage: "Hepatocyte".to_string(),
            delivery_vector: DeliveryVector::MrnaLipidNanoparticle,
            pulse_duration_days: 2.0,
            washout_interval_days: 5.0,
            total_cycles: 6,
            tet1_tet2_competence: 0.95,
            oncogene_cmyc_included: false,
            factor_dosage: 1.0,
            transient_tert_mrna_fold: 2.0,
            tp53_cdkn2a_safety_fold: 1.8,
        };

        // 1. Identity = 1.0 -> near baseline full efficacy
        let r_100 = simulate_in_vivo_osk_pulse_regimen_with_identity(70.0, 1.0, &base);
        assert!(r_100.total_rejuvenation_years >= 10.0);
        let c_100 = r_100.safety_coupling();
        assert_eq!(c_100.recommended_guard_action, crate::safety::GuardAction::Continue);
        assert!(c_100.identity_retention >= 0.85);
        assert!(c_100.stem_pool_fraction <= 0.15);

        // 2. Identity = 0.0 -> rejuvenation efficacy must be 0
        let r_000 = simulate_in_vivo_osk_pulse_regimen_with_identity(70.0, 0.0, &base);
        assert_eq!(r_000.total_rejuvenation_years, 0.0);
        let c_000 = r_000.safety_coupling();
        assert_eq!(c_000.recommended_guard_action, crate::safety::GuardAction::AbortSuicideSwitch);
        assert_eq!(c_000.identity_retention, 0.0);

        // 3. Monotonically decreasing identity steps -> monotonically decreasing rejuvenation years
        let steps = [1.0, 0.90, 0.80, 0.70, 0.60, 0.40, 0.0];
        let mut prev_rejuv = f64::INFINITY;
        for &id in &steps {
            let rep = simulate_in_vivo_osk_pulse_regimen_with_identity(70.0, id, &base);
            assert!(
                rep.total_rejuvenation_years <= prev_rejuv + 1e-9,
                "Efficacy must monotonically decrease with eroding identity: id={} got {}",
                id,
                rep.total_rejuvenation_years
            );
            prev_rejuv = rep.total_rejuvenation_years;
        }
    }
}
