//! High-fidelity biophysical telomere simulation, dynamic attrition, and multi-factorial biological clocks.
//!
//! Replaces simplified linear heuristics with biophysically grounded models accounting for:
//! - 8-oxodG oxidative lesions at GGG triplets.
//! - G-quadruplex (G4) replication fork stalling on lagging strands.
//! - Shelterin protein complex (TRF1, TRF2, POT1, TIN2, TPP1, RAP1) capping dynamics.
//! - Multi-factorial biological age combining DNA methylation, DSBs ($\gamma$-H2AX), and SASP burden.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Physiological state of chromosomal ends modeled under multi-factor biophysical stress.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HighFidelityTelomereState {
    /// Average terminal telomere length in base pairs.
    pub current_length_bp: f64,
    /// Basal end-replication problem loss rate (typically ~50-60 bp/division).
    pub basal_shortening_rate_bp: f64,
    /// Net effective shortening rate per division under active stress conditions.
    pub effective_shortening_rate_bp: f64,
    /// Cumulative loss attributed specifically to oxidative guanine lesions (8-oxodG).
    pub oxidative_attrition_bp: f64,
    /// Cumulative loss attributed to G-quadruplex (G4) replication fork stalling.
    pub g4_stall_attrition_bp: f64,
    /// Shelterin protein complex end-capping integrity index (0.0 = uncapped, 1.0 = fully protected).
    pub shelterin_capping: f64,
    /// Cumulative cellular divisions elapsed.
    pub divisions_elapsed: u32,
    /// Critical Hayflick limit threshold in bp triggering irreversible senescence (~3,800 bp).
    pub critical_hayflick_limit_bp: f64,
    /// Whether telomeres have breached the critical senescence threshold.
    pub is_critically_short: bool,
    /// Whether DNA damage response (ATM/ATR kinase signaling) has been triggered.
    pub ddr_alarm_triggered: bool,
}

impl HighFidelityTelomereState {
    /// Creates a youthful neonatal baseline state (~12,500 bp).
    pub fn youthful_neonatal() -> Self {
        Self {
            current_length_bp: 12_500.0,
            basal_shortening_rate_bp: 55.0,
            effective_shortening_rate_bp: 55.0,
            oxidative_attrition_bp: 0.0,
            g4_stall_attrition_bp: 0.0,
            shelterin_capping: 0.98,
            divisions_elapsed: 0,
            critical_hayflick_limit_bp: 3_800.0,
            is_critically_short: false,
            ddr_alarm_triggered: false,
        }
    }

    /// Evaluates the remaining replicative capacity (Hayflick limit) until cellular arrest.
    pub fn replicative_lifespan(&self) -> u32 {
        if self.current_length_bp <= self.critical_hayflick_limit_bp || self.effective_shortening_rate_bp <= 0.0 {
            return 0;
        }
        let remaining_bp = self.current_length_bp - self.critical_hayflick_limit_bp;
        (remaining_bp / self.effective_shortening_rate_bp).floor() as u32
    }
}

/// Simulates high-fidelity telomeric attrition over divisions under combined oxidative and replication stress.
pub fn simulate_high_fidelity_shortening(
    initial_length_bp: f64,
    divisions: u32,
    oxidative_stress: f64,
    replication_stress_g4: f64,
    initial_shelterin: f64,
) -> HighFidelityTelomereState {
    let basal_rate = 55.0;
    let critical_limit = 3_800.0;

    // 8-oxodG oxidative acceleration (factor kappa_ox ~ 0.35)
    let ox_multiplier = oxidative_stress.clamp(0.0, 5.0) * 0.35;
    // G4 stall acceleration on G-rich lagging strand (factor kappa_g4 ~ 0.25)
    let g4_multiplier = replication_stress_g4.clamp(0.0, 4.0) * 0.25;

    // Shelterin capping buffer: lower shelterin integrity accelerates end degradation
    let capping_factor = (1.0 / initial_shelterin.clamp(0.20, 1.0)).min(3.0);

    let effective_rate = basal_rate * (1.0 + ox_multiplier + g4_multiplier) * capping_factor;

    let total_loss = effective_rate * divisions as f64;
    let current_length = (initial_length_bp - total_loss).max(800.0);

    let oxidative_loss = basal_rate * ox_multiplier * capping_factor * divisions as f64;
    let g4_loss = basal_rate * g4_multiplier * capping_factor * divisions as f64;

    let is_critically_short = current_length <= critical_limit;

    // Dynamic degradation of shelterin complex capping
    let final_shelterin = if current_length > 6_500.0 {
        initial_shelterin.clamp(0.1, 0.99)
    } else {
        (initial_shelterin * ((current_length - 800.0) / 5_700.0)).clamp(0.05, 0.95)
    };

    let ddr_alarm_triggered = is_critically_short || final_shelterin < 0.30;

    HighFidelityTelomereState {
        current_length_bp: current_length,
        basal_shortening_rate_bp: basal_rate,
        effective_shortening_rate_bp: effective_rate,
        oxidative_attrition_bp: oxidative_loss,
        g4_stall_attrition_bp: g4_loss,
        shelterin_capping: final_shelterin,
        divisions_elapsed: divisions,
        critical_hayflick_limit_bp: critical_limit,
        is_critically_short,
        ddr_alarm_triggered,
    }
}

/// Parameters for multi-factorial high-fidelity biological age calculation.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HighFidelityAgeParams {
    /// Chronological age in years.
    pub chronological_age: f64,
    /// Horvath/Hannum epigenetic CpG methylation drift fraction (0.0 to 1.0).
    pub cpg_methylation_fraction: f64,
    /// Terminal telomere length in base pairs.
    pub telomere_length_bp: f64,
    /// Average persistent double-strand break ($\gamma$-H2AX / 53BP1) damage foci per nucleus.
    pub gamma_h2ax_foci: f64,
    /// Senescence-Associated Secretory Phenotype (SASP) inflammatory burden (0.0 to 1.0).
    pub sasp_burden: f64,
}

/// Output breakdown of high-fidelity biological age estimation.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HighFidelityAgeResult {
    /// Integrated biological age in years.
    pub composite_biological_age: f64,
    /// Aging acceleration ($Age_{bio} - Age_{chrono}$).
    pub aging_acceleration: f64,
    /// Component age contribution from DNA methylation drift.
    pub epigenetic_component_years: f64,
    /// Component age contribution from telomeric attrition.
    pub telomeric_component_years: f64,
    /// Component age contribution from persistent DNA DSB foci.
    pub dna_damage_component_years: f64,
    /// Component age contribution from SASP chronic inflammation.
    pub sasp_component_years: f64,
}

/// Evaluates high-fidelity biological age across epigenetics, telomeres, DNA damage, and SASP.
pub fn estimate_high_fidelity_biological_age(params: &HighFidelityAgeParams) -> HighFidelityAgeResult {
    let chrono = params.chronological_age.max(0.0);

    // 1. Epigenetic component (0.0 - 1.0 scaled to ~0-110 years)
    let epigenetic_age = (params.cpg_methylation_fraction * 105.0).clamp(0.0, 120.0);

    // 2. Telomeric component (14,000 bp at age 0, decreasing ~70 bp/year)
    let telomeric_age = ((14_000.0 - params.telomere_length_bp) / 72.0).clamp(0.0, 120.0);

    // 3. DNA damage foci component (~0.3 foci in youth, ~3.0 in old age, scaled to 0-100 years)
    let dsb_age = (params.gamma_h2ax_foci * 30.0).clamp(0.0, 120.0);

    // 4. SASP burden component (0.0 to 1.0 scaled to 0-100 years)
    let sasp_age = (params.sasp_burden * 95.0).clamp(0.0, 120.0);

    // Multi-factor weighted composite:
    // 40% Epigenetics, 25% Telomeres, 15% DNA damage, 10% SASP, 10% Chronological
    let composite = 0.40 * epigenetic_age
        + 0.25 * telomeric_age
        + 0.15 * dsb_age
        + 0.10 * sasp_age
        + 0.10 * chrono;

    let aging_accel = composite - chrono;

    HighFidelityAgeResult {
        composite_biological_age: composite.max(0.0),
        aging_acceleration: aging_accel,
        epigenetic_component_years: epigenetic_age,
        telomeric_component_years: telomeric_age,
        dna_damage_component_years: dsb_age,
        sasp_component_years: sasp_age,
    }
}

/// Kinetic model for Repeat Addition Processivity (RAP) of synthetic telomerase.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TelomeraseKinetics {
    /// Repeat Addition Processivity (average hexanucleotide repeats synthesized per binding event).
    pub repeat_addition_processivity: f64,
    /// Translocation efficiency across RNA template alignment domain (0.0 to 1.0).
    pub translocation_efficiency: f64,
    /// Catalytic turnover rate ($k_{cat}$ in nucleotide additions per second).
    pub k_cat_nucleotides_per_sec: f64,
}

impl TelomeraseKinetics {
    pub fn standard_human() -> Self {
        Self {
            repeat_addition_processivity: 8.5, // ~51 bp per binding event
            translocation_efficiency: 0.92,
            k_cat_nucleotides_per_sec: 12.0,
        }
    }

    pub fn engineered_hyper_processive() -> Self {
        Self {
            repeat_addition_processivity: 18.0, // ~108 bp per binding event
            translocation_efficiency: 0.98,
            k_cat_nucleotides_per_sec: 24.0,
        }
    }
}

/// Simulates high-fidelity elongation accounting for Repeat Addition Processivity and shelterin accessibility.
pub fn apply_high_fidelity_telomerase(
    state: &mut HighFidelityTelomereState,
    kinetics: &TelomeraseKinetics,
    expression_level: f64,
    cycles: u32,
) -> f64 {
    if expression_level <= 0.01 || cycles == 0 {
        return 0.0;
    }

    // Shelterin steric hindrance: dense shelterin (POT1/TRF1) regulates telomerase recruitment
    let accessibility = (1.1 - state.shelterin_capping).clamp(0.20, 0.95);

    let bp_per_cycle = kinetics.repeat_addition_processivity * 6.0 * kinetics.translocation_efficiency;
    let net_extension = bp_per_cycle * expression_level * accessibility * cycles as f64;

    let max_physiological_length = 18_500.0;
    let new_length = (state.current_length_bp + net_extension).min(max_physiological_length);
    let added = new_length - state.current_length_bp;

    state.current_length_bp = new_length;
    if state.current_length_bp > state.critical_hayflick_limit_bp {
        state.is_critically_short = false;
        state.ddr_alarm_triggered = false;
        state.shelterin_capping = (state.shelterin_capping + 0.18).min(0.98);
    }

    added
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_high_fidelity_shortening_with_stresses() {
        // Neonatal baseline (12,500 bp), 50 divisions under high ROS (2.0) and G4 stress (1.5)
        let state = simulate_high_fidelity_shortening(12_500.0, 50, 2.0, 1.5, 0.95);
        assert!(state.effective_shortening_rate_bp > 55.0);
        assert!(state.oxidative_attrition_bp > 0.0);
        assert!(state.g4_stall_attrition_bp > 0.0);
        assert!(state.current_length_bp < 12_500.0);
        assert!(!state.is_critically_short);
    }

    #[test]
    fn test_high_fidelity_biological_age_calculation() {
        let params = HighFidelityAgeParams {
            chronological_age: 50.0,
            cpg_methylation_fraction: 0.35,
            telomere_length_bp: 9_500.0,
            gamma_h2ax_foci: 0.8,
            sasp_burden: 0.15,
        };

        let res = estimate_high_fidelity_biological_age(&params);
        assert!(res.composite_biological_age > 20.0);
        assert!(res.composite_biological_age < 50.0); // Biologically rejuvenated
        assert!(res.aging_acceleration < 0.0);
    }
}
