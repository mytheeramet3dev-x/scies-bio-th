//! Biological aging clocks, epigenetic clock modeling, and physiological aging profiles.
//!
//! Models DNA methylation drift (Horvath/Hannum multi-tissue clock calibration) coupled
//! with telomere length attrition to estimate biological age and aging acceleration.

use crate::aging::telomere::TelomereState;
use std::collections::HashMap;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Comprehensive biological aging profile combining epigenetic, telomeric, and cellular metrics.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct AgingProfile {
    /// Chronological calendar age in years.
    pub chronological_age_years: f64,
    /// Epigenetically and telomerically estimated biological age in years.
    pub biological_age_years: f64,
    /// Telomeric structural status.
    pub telomere: TelomereState,
    /// Senescent cell burden fraction (0.0 = youthful, 1.0 = severe senescent accumulation).
    pub senescence_burden: f64,
    /// Average persistent double-strand break DNA damage foci (gamma-H2AX / 53BP1) per cell nucleus.
    pub dna_damage_foci: f64,
    /// Pathway-specific aging risk scores (e.g. "telomere_exhaustion", "cardiovascular", "cancer").
    pub risk_scores: HashMap<String, f64>,
}

impl AgingProfile {
    /// Creates a baseline youthful adult profile (age 25).
    pub fn youthful_adult(chronological_age: f64) -> Self {
        let telomere = TelomereState::youthful_baseline();
        let bio_age = estimate_biological_age(chronological_age, 0.25, telomere.current_length_bp);

        let mut risk_scores = HashMap::new();
        risk_scores.insert("telomere_exhaustion".to_string(), 0.05);
        risk_scores.insert("genomic_instability".to_string(), 0.08);
        risk_scores.insert("senescence_secretion".to_string(), 0.04);

        Self {
            chronological_age_years: chronological_age,
            biological_age_years: bio_age,
            telomere,
            senescence_burden: 0.05,
            dna_damage_foci: 0.4,
            risk_scores,
        }
    }

    /// Aging acceleration metric ($\text{Age}_{\text{Bio}} - \text{Age}_{\text{Chrono}}$).
    /// Negative values indicate slowed/rejuvenated biological aging.
    pub fn aging_acceleration(&self) -> f64 {
        calculate_aging_acceleration(self.chronological_age_years, self.biological_age_years)
    }
}

/// Estimates biological age from chronological age, global epigenetic methylation drift, and telomere length.
///
/// Weighting calibrated against Horvath multi-tissue epigenetic clock principles:
/// - 55% Epigenetic methylation status
/// - 30% Telomere attrition
/// - 15% Baseline chronological drift
pub fn estimate_biological_age(
    chronological_age: f64,
    methylation_fraction: f64,
    telomere_length_bp: f64,
) -> f64 {
    let chrono = chronological_age.max(0.0);

    // Baseline telomere reference: 14,000 bp at age 0, decreasing ~70 bp/year to 4,000 bp at age 100+
    let telomere_age_proxy = ((14_000.0 - telomere_length_bp) / 70.0).clamp(0.0, 120.0);

    // Methylation age proxy: non-linear scaling across 0.0 to 1.0
    let methylation_age_proxy = (methylation_fraction * 100.0).clamp(0.0, 120.0);

    let composite_bio_age = 0.55 * methylation_age_proxy + 0.30 * telomere_age_proxy + 0.15 * chrono;
    composite_bio_age.max(0.0)
}

/// Calculates aging acceleration in years ($A = \text{Age}_{\text{Bio}} - \text{Age}_{\text{Chrono}}$).
pub fn calculate_aging_acceleration(chronological_age: f64, biological_age: f64) -> f64 {
    biological_age - chronological_age
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_biological_age_estimation() {
        // Person aged 50 with long telomeres (11,000 bp) and low methylation (0.35)
        let bio_age = estimate_biological_age(50.0, 0.35, 11_000.0);
        // Telomere proxy = (14,000 - 11,000) / 70 = 42.85
        // Methylation proxy = 35.0
        // Chrono = 50.0
        // Composite = 0.55*35 + 0.30*42.85 + 0.15*50 = 19.25 + 12.85 + 7.5 = 39.6
        assert!(bio_age < 50.0); // Younger than chronological age
        assert!((bio_age - 39.6).abs() < 1.0);

        let accel = calculate_aging_acceleration(50.0, bio_age);
        assert!(accel < 0.0); // Slowed aging
    }
}
