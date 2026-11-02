//! In Vitro Gametogenesis (IVG) simulation models.
//!
//! Models the induction of Primordial Germ Cell-Like Cells (PGCLCs) from stem cells,
//! meiotic recombination, epigenetic erasure/re-establishment, and synthetic gamete generation.

use crate::error::{BioError, Result};
use crate::ivg::epigenetics::{EpigeneticProfile, ImprintPattern};
use crate::seq::Dna;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Developmental state along the germline differentiation trajectory.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum StemCellState {
    /// Naive pluripotent stem cell (ESC or fully reprogrammed naive iPSC).
    NaivePluripotent,
    /// Primed pluripotent state (epiblast-like cell ready for PGC induction).
    PrimedPluripotent,
    /// Primordial Germ Cell-Like Cell (PGCLC, Sox17+/Blimp1+/Tfap2c+).
    PrimordialGermCellLike,
    /// Gonocyte (sex-specified embryonic germline precursor).
    Gonocyte,
    /// Mature haploid gamete (Oocyte arrested at MII or mature Spermatozoon).
    MatureGamete,
}

/// Type of gamete being generated in vitro.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum GameteType {
    /// Oocyte (Egg cell, large volume, maternal imprints, ~100k mitochondria).
    Oocyte,
    /// Spermatozoon (Sperm cell, motile, paternal imprints, high protamine compaction).
    Spermatozoon,
}

/// Supported In Vitro Gametogenesis differentiation protocols.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum InductionProtocol {
    /// Hayashi et al. (2016) murine reconstitution protocol (BMP4 + LIF + SCF).
    Hayashi2016,
    /// Saitou et al. (2018) human PGCLC induction protocol (ActA/Wnt signaling).
    Saitou2018,
    /// Feeder-free chemically defined synthetic embryology protocol.
    ChemicallyDefinedFeederFree,
}

/// Quality assessment metrics for an in vitro derived synthetic gamete.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct GameteQuality {
    /// Estimated probability of chromosomal aneuploidy (0.0 = euploid, 1.0 = highly abnormal).
    pub aneuploidy_risk: f64,
    /// Completeness of sex-specific imprint establishment (0.0 to 1.0).
    pub imprint_fidelity: f64,
    /// Functional mitochondrial copy number in cytoplasm.
    pub mitochondrial_count: u64,
    /// Whether first polar body was successfully extruded (for oocytes).
    pub polar_body_extruded: bool,
    /// Overall composite quality score (0.0 to 100.0).
    pub overall_score: f64,
}

/// Synthetic Zona Pellucida glycoprotein matrix enveloping the oocyte.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ZonaPellucida {
    /// ZP1 structural cross-linker integrity (0.0 to 1.0).
    pub zp1_crosslinker_density: f64,
    /// ZP2 cleavage status (true = intact/permissive; false = cleaved/polyspermy block).
    pub zp2_intact: bool,
    /// ZP3 species-specific primary sperm receptor affinity (0.0 to 1.0).
    pub zp3_receptor_affinity: f64,
    /// ZP4 filament stabilization factor (0.0 to 1.0).
    pub zp4_stabilization: f64,
    /// Zona thickness in micrometers (typically 13.0 to 18.0 µm in mammals).
    pub thickness_microns: f64,
    /// Young's Elasticity modulus in kPa (15.0 - 25.0 kPa for microinjection resilience).
    pub elasticity_kpa: f64,
}

impl ZonaPellucida {
    /// Creates an engineered synthetic Zona Pellucida matrix.
    pub fn new(thickness_microns: f64, elasticity_kpa: f64) -> Self {
        Self {
            zp1_crosslinker_density: 0.95,
            zp2_intact: true,
            zp3_receptor_affinity: 0.98,
            zp4_stabilization: 0.92,
            thickness_microns,
            elasticity_kpa,
        }
    }

    /// Default synthetic Zona Pellucida optimized for feline and Cave Lion (Panthera spelaea) oocytes.
    pub fn cave_lion_default() -> Self {
        Self::new(16.5, 18.5)
    }

    /// Simulates cortical granule exocytosis hardening upon fertilization/activation (blocks polyspermy).
    pub fn trigger_cortical_hardening(&mut self) {
        self.zp2_intact = false;
        self.elasticity_kpa *= 2.5;
    }
}

/// An in vitro derived synthetic gamete (Egg or Sperm).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SyntheticGamete {
    /// Gamete sex specification.
    pub gamete_type: GameteType,
    /// Haploid genomic DNA sequence resulting from meiotic reduction.
    pub haploid_genome: Dna,
    /// Gamete epigenetic and imprinting profile.
    pub epigenetic_profile: EpigeneticProfile,
    /// Optional synthetic Zona Pellucida protective shell (for oocytes).
    pub zona_pellucida: Option<ZonaPellucida>,
    /// Gamete quality assessment.
    pub quality: GameteQuality,
    /// NCBI Taxonomy ID of the source species.
    pub tax_id: u32,
}

/// Intermediate PGCLC culture ready for gametogenesis maturation.
#[derive(Debug, Clone)]
pub struct PgclcCulture {
    pub donor_dna: Dna,
    pub tax_id: u32,
    pub induction_efficiency: f64,
    pub epigenetic_profile: EpigeneticProfile,
}

impl PgclcCulture {
    /// Induces PGCLCs from a donor pluripotent stem cell genome using the specified protocol.
    pub fn induce(
        donor_dna: &Dna,
        tax_id: u32,
        protocol: InductionProtocol,
        differentiation_days: u32,
    ) -> Result<Self> {
        if donor_dna.sequence.is_empty() {
            return Err(BioError::ParseError("Donor DNA cannot be empty".to_string()));
        }

        let base_efficiency = match protocol {
            InductionProtocol::Hayashi2016 => 0.82,
            InductionProtocol::Saitou2018 => 0.76,
            InductionProtocol::ChemicallyDefinedFeederFree => 0.68,
        };

        // Induction kinetics peak around day 4-6
        let day_factor = if (4..=8).contains(&differentiation_days) {
            1.0
        } else {
            0.5
        };

        let final_efficiency = base_efficiency * day_factor;
        let mut ep_profile = EpigeneticProfile::naive_pluripotent();
        // Global demethylation occurs during PGCLC stage
        ep_profile.demethylate_pgc(final_efficiency);

        Ok(Self {
            donor_dna: donor_dna.clone(),
            tax_id,
            induction_efficiency: final_efficiency,
            epigenetic_profile: ep_profile,
        })
    }

    /// Matures the PGCLCs into a functional synthetic oocyte or sperm through simulated meiosis.
    pub fn mature_gamete(
        &self,
        gamete_type: GameteType,
        crossover_seed: u64,
    ) -> Result<SyntheticGamete> {
        if self.induction_efficiency < 0.30 {
            return Err(BioError::ReprogrammingFailed(
                "PGCLC induction efficiency too low to support meiotic maturation".to_string(),
            ));
        }

        // Simulate meiotic reduction & crossing over on the donor sequence
        let haploid_seq = simulate_meiotic_reduction(&self.donor_dna.sequence, crossover_seed);
        let haploid_dna = Dna::new(&haploid_seq)?;

        let (ep_profile, mito_count, pb_extruded) = match gamete_type {
            GameteType::Oocyte => {
                let mut profile = EpigeneticProfile::mature_oocyte();
                profile.imprint_pattern = ImprintPattern::Maternal;
                (profile, 100_000, true)
            }
            GameteType::Spermatozoon => {
                let mut profile = EpigeneticProfile::mature_sperm();
                profile.imprint_pattern = ImprintPattern::Paternal;
                (profile, 100, false)
            }
        };

        let imprint_fidelity = (self.induction_efficiency * 0.95).min(1.0);
        let aneuploidy_risk = (1.0 - self.induction_efficiency) * 0.25;
        let overall_score = (imprint_fidelity * 60.0 + (1.0 - aneuploidy_risk) * 40.0).clamp(0.0, 100.0);

        let quality = GameteQuality {
            aneuploidy_risk,
            imprint_fidelity,
            mitochondrial_count: mito_count,
            polar_body_extruded: pb_extruded,
            overall_score,
        };

        let zona = if gamete_type == GameteType::Oocyte {
            Some(ZonaPellucida::cave_lion_default())
        } else {
            None
        };

        Ok(SyntheticGamete {
            gamete_type,
            haploid_genome: haploid_dna,
            epigenetic_profile: ep_profile,
            zona_pellucida: zona,
            quality,
            tax_id: self.tax_id,
        })
    }
}

/// Simulates meiotic crossing-over and reduction from diploid to haploid state.
fn simulate_meiotic_reduction(diploid_bytes: &[u8], seed: u64) -> Vec<u8> {
    if diploid_bytes.is_empty() {
        return Vec::new();
    }
    // Simple deterministic pseudo-random shuffling / crossover selection based on seed
    let len = diploid_bytes.len();
    let mut out = diploid_bytes.to_vec();
    let mut state = seed.wrapping_add(0x9E3779B97F4A7C15);

    let mut i = 0;
    while i + 1 < len {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        let swap = (state >> 32) as u32 % 100 < 5; // 5% crossover perturbation rate
        if swap {
            // Meiotic crossover segment swap between adjacent nucleotides
            out.swap(i, i + 1);
            i += 2;
        } else {
            i += 1;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ivg_pgclc_to_synthetic_oocyte() {
        let donor_seq = Dna::new(b"ATGCGATCGATCGATCGATC").unwrap();
        let culture = PgclcCulture::induce(&donor_seq, 9606, InductionProtocol::Hayashi2016, 5).unwrap();
        assert!(culture.induction_efficiency > 0.70);

        let oocyte = culture.mature_gamete(GameteType::Oocyte, 42).unwrap();
        assert_eq!(oocyte.gamete_type, GameteType::Oocyte);
        assert_eq!(oocyte.quality.mitochondrial_count, 100_000);
        assert!(oocyte.quality.polar_body_extruded);
        assert!(oocyte.quality.overall_score > 70.0);
    }

    #[test]
    fn test_meiotic_reduction_crossover() {
        let original = b"ATGCGATCGATCGATCGATC";
        // Run with seeds to check that meiotic reduction maintains length and valid bases
        for seed in [1, 42, 12345, 99999] {
            let reduced = simulate_meiotic_reduction(original, seed);
            assert_eq!(reduced.len(), original.len());
            assert!(Dna::new(&reduced).is_ok());
        }
    }
}
