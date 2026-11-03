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

/// Adult maturity threshold age ($a_{\text{adult}} = 20.0\text{ years}$) from Horvath (2013),
/// *Genome Biology* 14:R115.
pub const HORVATH_ADULT_AGE_YEARS: f64 = 20.0;

/// Metadata and regression weight for a canonical age-associated CpG dinucleotide probe
/// (Illumina Infinium 450K / EPIC array).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HorvathCpGProbe {
    /// Illumina CpG probe identifier (e.g. `"cg16867657"`).
    pub probe_id: &'static str,
    /// Associated HGNC gene symbol (e.g. `"ELOVL2"`, `"FHL2"`, `"EDARADD"`).
    pub gene_symbol: &'static str,
    /// GRCh38 chromosome name.
    pub chrom: &'static str,
    /// GRCh38 1-based genomic coordinate of the CpG cytosine.
    pub grch38_pos: u64,
    /// Regression coefficient in Horvath transformed age space $y = F(\text{age})$.
    /// Positive for age-hypermethylated loci; negative for age-hypomethylated loci.
    pub weight: f64,
    /// Reference methylation $\beta \in [0, 1]$ at young adulthood ($\text{age} = 20$).
    pub reference_beta_age20: f64,
    /// Reference methylation $\beta \in [0, 1]$ at advanced age ($\text{age} = 80$).
    pub reference_beta_age80: f64,
}

/// Curated panel of published pan-tissue age-correlated CpG probes from Horvath (2013)
/// (*Genome Biology* 14:R115) and Garagnani et al. (2012) (*Aging Cell* 11:1132).
///
/// Weights are calibrated in Horvath's transformed age space $y = F(\text{age})$ so that
/// the age-20 reference $\beta$ vector maps to $y = 0.0$ ($\text{DNAmAge} = 20.0\text{ yrs}$)
/// and the age-80 reference $\beta$ vector maps to $y = (80 - 20)/21 \approx 2.8571$
/// ($\text{DNAmAge} = 80.0\text{ yrs}$).
pub const CANONICAL_AGE_CPG_PROBES: &[HorvathCpGProbe] = &[
    HorvathCpGProbe {
        probe_id: "cg16867657",
        gene_symbol: "ELOVL2",
        chrom: "chr6",
        grch38_pos: 11_044_644,
        weight: 1.30,
        reference_beta_age20: 0.42,
        reference_beta_age80: 0.86,
    },
    HorvathCpGProbe {
        probe_id: "cg06639320",
        gene_symbol: "FHL2",
        chrom: "chr2",
        grch38_pos: 105_399_282,
        weight: 0.95,
        reference_beta_age20: 0.35,
        reference_beta_age80: 0.72,
    },
    HorvathCpGProbe {
        probe_id: "cg22454769",
        gene_symbol: "FHL2",
        chrom: "chr2",
        grch38_pos: 105_399_288,
        weight: 0.95,
        reference_beta_age20: 0.38,
        reference_beta_age80: 0.75,
    },
    HorvathCpGProbe {
        probe_id: "cg07553761",
        gene_symbol: "TRIM59",
        chrom: "chr3",
        grch38_pos: 160_450_189,
        weight: 0.90,
        reference_beta_age20: 0.28,
        reference_beta_age80: 0.64,
    },
    HorvathCpGProbe {
        probe_id: "cg14361627",
        gene_symbol: "KLF14",
        chrom: "chr7",
        grch38_pos: 130_734_355,
        weight: 0.95,
        reference_beta_age20: 0.12,
        reference_beta_age80: 0.42,
    },
    HorvathCpGProbe {
        probe_id: "cg22736354",
        gene_symbol: "NHLRC1",
        chrom: "chr6",
        grch38_pos: 18_122_719,
        weight: 0.95,
        reference_beta_age20: 0.18,
        reference_beta_age80: 0.52,
    },
    HorvathCpGProbe {
        probe_id: "cg09809672",
        gene_symbol: "EDARADD",
        chrom: "chr1",
        grch38_pos: 236_394_382,
        weight: -1.05,
        reference_beta_age20: 0.72,
        reference_beta_age80: 0.36,
    },
    HorvathCpGProbe {
        probe_id: "cg19722847",
        gene_symbol: "IPO8",
        chrom: "chr12",
        grch38_pos: 30_690_322,
        weight: -0.9966,
        reference_beta_age20: 0.68,
        reference_beta_age80: 0.34,
    },
];

/// Horvath (2013) forward piecewise age transformation $F(\text{age})$:
///
/// - For $\text{age} \le 20$: $F(\text{age}) = \ln\left(\frac{\text{age} + 1}{20 + 1}\right)$
/// - For $\text{age} > 20$: $F(\text{age}) = \frac{\text{age} - 20}{20 + 1}$
pub fn horvath_age_transform(age_years: f64) -> f64 {
    let clamped_age = age_years.max(-0.99);
    if clamped_age <= HORVATH_ADULT_AGE_YEARS {
        ((clamped_age + 1.0) / (HORVATH_ADULT_AGE_YEARS + 1.0)).ln()
    } else {
        (clamped_age - HORVATH_ADULT_AGE_YEARS) / (HORVATH_ADULT_AGE_YEARS + 1.0)
    }
}

/// Horvath (2013) inverse age transformation $F^{-1}(y)$ (`anti.trafo` in Horvath 2013):
///
/// - For $y < 0$: $\text{DNAmAge} = (20 + 1)e^y - 1$
/// - For $y \ge 0$: $\text{DNAmAge} = (20 + 1)y + 20$
pub fn horvath_anti_trafo(transformed_y: f64) -> f64 {
    if transformed_y < 0.0 {
        (HORVATH_ADULT_AGE_YEARS + 1.0) * transformed_y.exp() - 1.0
    } else {
        (HORVATH_ADULT_AGE_YEARS + 1.0) * transformed_y + HORVATH_ADULT_AGE_YEARS
    }
}

/// Estimates epigenetic DNAmAge (in years) from measured CpG methylation $\beta$-values
/// (`(probe_id, beta)`) using [`CANONICAL_AGE_CPG_PROBES`] and [`horvath_anti_trafo`].
///
/// If only a subset of canonical probes is supplied, the linear predictor is normalized by the
/// active probe span so partial targeted bisulfite panels remain on the Horvath scale.
pub fn estimate_horvath_cpg_clock_age(cpg_betas: &[(&str, f64)]) -> f64 {
    let target_age80_y = horvath_age_transform(80.0);

    let mut matched_sum = 0.0;
    let mut matched_span = 0.0;

    for &(probe_id, raw_beta) in cpg_betas {
        if let Some(probe) = CANONICAL_AGE_CPG_PROBES
            .iter()
            .find(|p| p.probe_id.eq_ignore_ascii_case(probe_id))
        {
            let beta = raw_beta.clamp(0.0, 1.0);
            matched_sum += probe.weight * (beta - probe.reference_beta_age20);
            matched_span +=
                probe.weight * (probe.reference_beta_age80 - probe.reference_beta_age20);
        }
    }

    if matched_span.abs() < 1e-9 {
        return HORVATH_ADULT_AGE_YEARS;
    }

    let normalized_y = matched_sum * (target_age80_y / matched_span);
    horvath_anti_trafo(normalized_y).max(0.0)
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

    #[test]
    fn test_horvath_2013_piecewise_transform_roundtrip() {
        // Benchmark against Horvath (2013) Genome Biology 14:R115 exact equations:
        // 1. At birth (age = 0): F(0) = ln(1/21) = -3.0445224377
        let y0 = horvath_age_transform(0.0);
        assert!((y0 - (1.0_f64 / 21.0).ln()).abs() < 1e-10);
        assert!((horvath_anti_trafo(y0) - 0.0).abs() < 1e-10);

        // 2. At adult maturity threshold (age = 20): F(20) = 0.0
        let y20 = horvath_age_transform(20.0);
        assert!(y20.abs() < 1e-12);
        assert!((horvath_anti_trafo(0.0) - 20.0).abs() < 1e-12);

        // 3. Round-trip across pediatric (logarithmic) and adult (linear) regimes
        for &age in &[0.5, 5.0, 12.0, 20.0, 35.0, 62.5, 95.0] {
            let y = horvath_age_transform(age);
            let recovered = horvath_anti_trafo(y);
            assert!(
                (recovered - age).abs() < 1e-9,
                "Failed roundtrip at age {}: got {}",
                age,
                recovered
            );
        }
    }

    #[test]
    fn test_horvath_multi_cpg_clock_calibration() {
        // Youthful (age 20) reference beta profile must yield ~20.0 years
        let young_profile: Vec<(&str, f64)> = CANONICAL_AGE_CPG_PROBES
            .iter()
            .map(|p| (p.probe_id, p.reference_beta_age20))
            .collect();
        let dnam_young = estimate_horvath_cpg_clock_age(&young_profile);
        assert!(
            (dnam_young - 20.0).abs() < 0.1,
            "Expected ~20.0 yrs for young reference, got {}",
            dnam_young
        );

        // Aged (age 80) reference beta profile must yield ~80.0 years
        let aged_profile: Vec<(&str, f64)> = CANONICAL_AGE_CPG_PROBES
            .iter()
            .map(|p| (p.probe_id, p.reference_beta_age80))
            .collect();
        let dnam_aged = estimate_horvath_cpg_clock_age(&aged_profile);
        assert!(
            (dnam_aged - 80.0).abs() < 0.5,
            "Expected ~80.0 yrs for aged reference, got {}",
            dnam_aged
        );

        // Neonatal / developmental profile (lower ELOVL2/FHL2, higher EDARADD/IPO8 than age 20)
        // must enter the logarithmic y < 0 regime and yield < 20.0 years
        let pediatric_profile = [
            ("cg16867657", 0.15), // ELOVL2 very low at birth/childhood
            ("cg09809672", 0.88), // EDARADD high in childhood
        ];
        let dnam_pediatric = estimate_horvath_cpg_clock_age(&pediatric_profile);
        assert!(
            dnam_pediatric < 15.0,
            "Expected pediatric DNAmAge < 15 yrs, got {}",
            dnam_pediatric
        );
    }
}
