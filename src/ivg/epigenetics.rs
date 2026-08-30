//! Epigenetic landscape, DNA methylation, histone marks, and imprint reprogramming models.
//!
//! Models the chromatin state transitions essential for somatic cell nuclear transfer (SCNT)
//! and in vitro gametogenesis (IVG).

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Canonical histone post-translational modification marks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum HistoneMark {
    /// Active promoter mark (Trimethylation of Lysine 4 on Histone H3).
    H3K4me3,
    /// Polycomb-mediated repressive mark (Trimethylation of Lysine 27 on Histone H3).
    H3K27me3,
    /// Heterochromatin / somatic barrier mark (Trimethylation of Lysine 9 on Histone H3).
    H3K9me3,
    /// Open chromatin / active enhancer mark (Acetylation of Lysine 9 on Histone H3).
    H3K9ac,
}

/// Genomic imprinting pattern defining parent-of-origin specific gene expression.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ImprintPattern {
    /// Paternal imprinting (e.g. H19 ICR methylated, Igf2 expressed).
    Paternal,
    /// Maternal imprinting (e.g. H19 ICR unmethylated, Igf2r methylated).
    Maternal,
    /// Erased state (Primordial Germ Cells during global epigenetic reprogramming).
    Erased,
    /// Aberrant / loss-of-imprinting (common in improperly reprogrammed SCNT clones).
    Aberrant,
}

/// Epigenetic status of a genome, cell, or gamete.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct EpigeneticProfile {
    /// Global 5-methylcytosine (5mC) DNA methylation level (0.0 to 1.0).
    pub global_methylation: f64,
    /// H3K9me3 heterochromatin repressive mark level (somatic reprogramming barrier, 0.0 to 1.0).
    pub h3k9me3_barrier: f64,
    /// H3K27me3 polycomb repressive mark level (developmental gene silencing, 0.0 to 1.0).
    pub h3k27me3_level: f64,
    /// Imprinted Control Region (ICR) methylation pattern.
    pub imprint_pattern: ImprintPattern,
    /// Biological epigenetic age (estimated in years or normalized developmental units).
    pub epigenetic_age: f64,
}

impl EpigeneticProfile {
    /// Creates a default somatic cell epigenetic profile with donor age.
    pub fn somatic_default(donor_age_years: f64) -> Self {
        Self {
            global_methylation: 0.75,
            h3k9me3_barrier: 0.85,
            h3k27me3_level: 0.70,
            imprint_pattern: ImprintPattern::Paternal, // Arbitrary somatic baseline
            epigenetic_age: donor_age_years.max(0.0),
        }
    }

    /// Creates an epigenetic profile corresponding to naive pluripotent stem cells (ESCs/iPSCs).
    pub fn naive_pluripotent() -> Self {
        Self {
            global_methylation: 0.30,
            h3k9me3_barrier: 0.20,
            h3k27me3_level: 0.40,
            imprint_pattern: ImprintPattern::Paternal,
            epigenetic_age: 0.0,
        }
    }

    /// Creates a mature oocyte (egg) epigenetic profile.
    pub fn mature_oocyte() -> Self {
        Self {
            global_methylation: 0.45,
            h3k9me3_barrier: 0.10,
            h3k27me3_level: 0.35,
            imprint_pattern: ImprintPattern::Maternal,
            epigenetic_age: 0.0,
        }
    }

    /// Creates a mature spermatozoon (sperm) epigenetic profile (protamine-compacted, high 5mC).
    pub fn mature_sperm() -> Self {
        Self {
            global_methylation: 0.90,
            h3k9me3_barrier: 0.15,
            h3k27me3_level: 0.25,
            imprint_pattern: ImprintPattern::Paternal,
            epigenetic_age: 0.0,
        }
    }

    /// Simulates primordial germ cell (PGC) global DNA demethylation.
    ///
    /// `efficiency` (0.0 to 1.0) controls the rate of TET-mediated oxidation and passive dilution.
    pub fn demethylate_pgc(&mut self, efficiency: f64) {
        let eff = efficiency.clamp(0.0, 1.0);
        self.global_methylation = (self.global_methylation * (1.0 - eff * 0.95)).max(0.02);
        self.h3k9me3_barrier = (self.h3k9me3_barrier * (1.0 - eff * 0.90)).max(0.05);
        if self.global_methylation < 0.10 {
            self.imprint_pattern = ImprintPattern::Erased;
            self.epigenetic_age = 0.0;
        } else {
            self.epigenetic_age = (self.epigenetic_age * (1.0 - eff)).max(0.0);
        }
    }

    /// Applies KDM4 histone demethylase treatment (e.g. Kdm4d mRNA injection during SCNT)
    /// to erase the repressive H3K9me3 barrier on the somatic donor nucleus.
    pub fn apply_kdm4_treatment(&mut self, dosage: f64) {
        let dose = dosage.clamp(0.0, 1.0);
        self.h3k9me3_barrier = (self.h3k9me3_barrier * (1.0 - dose * 0.85)).max(0.05);
    }

    /// Computes the somatic reprogramming resistance barrier score (0.0 = completely permissive, 100.0 = impenetrable).
    /// High H3K9me3 and high somatic methylation contribute to SCNT developmental arrest.
    pub fn reprogramming_barrier_score(&self) -> f64 {
        let h3k9_penalty = self.h3k9me3_barrier * 50.0;
        let meth_penalty = self.global_methylation * 35.0;
        let age_penalty = (self.epigenetic_age / 50.0).min(1.0) * 15.0;
        (h3k9_penalty + meth_penalty + age_penalty).min(100.0)
    }

    /// Rejuvenates the epigenetic clock towards zero (blastocyst-stage resetting).
    pub fn reset_epigenetic_age(&mut self) {
        self.epigenetic_age = 0.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_epigenetic_barrier_calculation() {
        let old_somatic = EpigeneticProfile::somatic_default(45.0);
        let barrier = old_somatic.reprogramming_barrier_score();
        assert!(barrier > 70.0); // High barrier for 45yo somatic nucleus

        let mut treated = old_somatic.clone();
        treated.apply_kdm4_treatment(0.95);
        let treated_barrier = treated.reprogramming_barrier_score();
        assert!(treated_barrier < barrier);
        assert!(treated.h3k9me3_barrier < 0.20);
    }

    #[test]
    fn test_pgc_demethylation_and_imprint_erasure() {
        let mut profile = EpigeneticProfile::somatic_default(25.0);
        profile.demethylate_pgc(0.98);
        assert!(profile.global_methylation < 0.10);
        assert_eq!(profile.imprint_pattern, ImprintPattern::Erased);
        assert_eq!(profile.epigenetic_age, 0.0);
    }
}
