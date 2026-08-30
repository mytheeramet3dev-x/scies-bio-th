//! Synthetic embryology and developmental progression simulation engine.
//!
//! Models embryogenesis from Zygote to Blastocyst, lineage specification (ICM vs TE),
//! Embryonic Genome Activation (EGA), and synthetic embryo (Stem Cell-derived blastoid) development.

use crate::error::{BioError, Result};
use crate::ivg::gamete::{GameteType, SyntheticGamete};
use crate::ivg::scnt::ScntConstruct;
use crate::seq::Dna;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Morphological and developmental stages of embryogenesis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum DevelopmentalStage {
    /// Day 0: Fertilized 1-cell zygote.
    Zygote,
    /// Day 1: 2-cell stage, onset of Embryonic Genome Activation (EGA).
    TwoCell,
    /// Day 2: 4-cell cleavage stage.
    FourCell,
    /// Day 2.5: 8-cell stage (pre-compaction).
    EightCell,
    /// Day 3-4: Morula stage (cellular compaction and polarity establishment).
    Morula,
    /// Day 5: Blastocyst with distinct Inner Cell Mass (ICM) and Trophectoderm (TE).
    Blastocyst,
    /// Day 6: Hatched blastocyst breaking through synthetic Zona Pellucida.
    HatchedBlastocyst,
    /// Day 7-14: Gastrulation (formation of Ectoderm, Mesoderm, and Endoderm germ layers).
    Gastrula,
    /// Day 15-45: Early organogenesis (cardiac contraction, neural tube closure, limb bud outgrowth).
    Organogenesis,
    /// Day 45-90: Fetal maturation (subcutaneous fat, dense fur coat, skeletal ossification).
    FetalDevelopment,
    /// Day 100-110+: Full-term mature fetus ready for artificial decannulation & emergence.
    FullTermNeonate,
}

/// The biotechnological origin of the embryo.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum EmbryoOrigin {
    /// Standard In Vitro Fertilization (Gamete + Gamete fusion).
    InVitroFertilization,
    /// In Vitro Gametogenesis (Stem cell-derived synthetic egg and sperm).
    InVitroGametogenesis,
    /// Somatic Cell Nuclear Transfer (SCNT clone).
    SomaticCellNuclearTransfer,
    /// Synthetic Embryoid / Blastoid (direct stem cell self-assembly without gametes).
    StemCellSyntheticBlastoid,
}

/// Artificial Womb (Biobag / Extra-Uterine Life Support Device).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ArtificialWombChamber {
    /// Sterile synthetic amniotic fluid circulation rate in mL/hour.
    pub amniotic_fluid_rate_ml_hr: f64,
    /// Pumpless arteriovenous membrane oxygenator PaO2 in mmHg (target: 28-35 mmHg).
    pub oxygenator_pao2_mmhg: f64,
    /// Umbilical vascular cannulation circuit patency (0.0 to 1.0).
    pub umbilical_cannulation_patency: f64,
    /// Temperature in Celsius (37.5 - 38.5°C physiological feline core).
    pub chamber_temp_celsius: f64,
}

impl Default for ArtificialWombChamber {
    fn default() -> Self {
        Self {
            amniotic_fluid_rate_ml_hr: 50.0,
            oxygenator_pao2_mmhg: 30.0,
            umbilical_cannulation_patency: 0.98,
            chamber_temp_celsius: 38.0,
        }
    }
}

/// Physiological and environmental parameters for in vitro / artificial womb culture.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct EmbryoCultureEnvironment {
    /// Oxygen concentration percentage (optimal physiological is ~5.0% low O2 in pre-implantation).
    pub oxygen_percent: f64,
    /// Presence of dynamic microfluidic laminar flow to clear toxic metabolites.
    pub microfluidic_flow: bool,
    /// Culture medium optimized for embryonic genome activation (e.g. sequential KSOM/G1-G2).
    pub sequential_media_optimized: bool,
    /// Artificial womb (extra-uterine support device) connected post-blastocyst stage.
    pub artificial_womb_supported: bool,
    /// Connected artificial womb chamber parameters.
    pub womb_chamber: Option<ArtificialWombChamber>,
}

impl Default for EmbryoCultureEnvironment {
    fn default() -> Self {
        Self {
            oxygen_percent: 5.0,
            microfluidic_flow: true,
            sequential_media_optimized: true,
            artificial_womb_supported: false,
            womb_chamber: None,
        }
    }
}

/// A biological or synthetic embryo undergoing development.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SyntheticEmbryo {
    /// Current developmental stage.
    pub stage: DevelopmentalStage,
    /// Creation origin mode.
    pub origin: EmbryoOrigin,
    /// Diploid embryonic nuclear DNA.
    pub nuclear_genome: Dna,
    /// NCBI Taxonomy ID of species.
    pub tax_id: u32,
    /// Total cell count.
    pub total_cell_count: usize,
    /// Inner Cell Mass (ICM) cell count (pluripotent epiblast precursors).
    pub icm_cell_count: usize,
    /// Trophectoderm (TE) cell count (placenta/extraembryonic precursors).
    pub trophectoderm_cell_count: usize,
    /// Overall developmental viability score (0.0 to 100.0).
    pub viability_score: f64,
    /// Cumulative culture age in days.
    pub age_days: f64,
    /// Estimated fetal weight in grams (during artificial womb gestation).
    pub fetal_weight_grams: f64,
    /// Whether embryonic development has arrested.
    pub is_arrested: bool,
}

impl SyntheticEmbryo {
    /// Creates a zygote via in vitro gamete fertilization.
    pub fn fertilize(oocyte: SyntheticGamete, sperm: SyntheticGamete) -> Result<Self> {
        if oocyte.gamete_type != GameteType::Oocyte {
            return Err(BioError::IncompatibleGamete("First gamete must be an Oocyte".to_string()));
        }
        if sperm.gamete_type != GameteType::Spermatozoon {
            return Err(BioError::IncompatibleGamete("Second gamete must be a Spermatozoon".to_string()));
        }
        if oocyte.tax_id != sperm.tax_id {
            return Err(BioError::IncompatibleGamete(format!(
                "Species mismatch between gametes (TaxID {} vs {})",
                oocyte.tax_id, sperm.tax_id
            )));
        }

        // Combine haploid sets into diploid zygote sequence
        let mut diploid_seq = oocyte.haploid_genome.sequence.clone();
        diploid_seq.extend_from_slice(&sperm.haploid_genome.sequence);
        let diploid_dna = Dna::new(&diploid_seq)?;

        let combined_quality = (oocyte.quality.overall_score * 0.6 + sperm.quality.overall_score * 0.4).clamp(0.0, 100.0);

        Ok(Self {
            stage: DevelopmentalStage::Zygote,
            origin: EmbryoOrigin::InVitroGametogenesis,
            nuclear_genome: diploid_dna,
            tax_id: oocyte.tax_id,
            total_cell_count: 1,
            icm_cell_count: 0,
            trophectoderm_cell_count: 0,
            viability_score: combined_quality,
            age_days: 0.0,
            fetal_weight_grams: 0.001,
            is_arrested: false,
        })
    }

    /// Creates an embryonic construct from an SCNT nuclear transfer event.
    pub fn from_scnt(scnt: ScntConstruct) -> Result<Self> {
        let viability = scnt.predict_blastocyst_competence() * 100.0;

        Ok(Self {
            stage: DevelopmentalStage::Zygote,
            origin: EmbryoOrigin::SomaticCellNuclearTransfer,
            nuclear_genome: scnt.nuclear_dna,
            tax_id: scnt.donor_tax_id,
            total_cell_count: 1,
            icm_cell_count: 0,
            trophectoderm_cell_count: 0,
            viability_score: viability,
            age_days: 0.0,
            fetal_weight_grams: 0.001,
            is_arrested: false,
        })
    }

    /// Advances embryonic development through time under specified culture conditions.
    pub fn advance_development(
        &mut self,
        additional_hours: f64,
        env: &EmbryoCultureEnvironment,
    ) -> Result<DevelopmentalStage> {
        if self.is_arrested {
            return Err(BioError::DevelopmentArrested(format!(
                "Embryo development arrested at stage {:?}",
                self.stage
            )));
        }

        self.age_days += additional_hours / 24.0;

        // Environmental quality modifier
        let env_modifier = if env.oxygen_percent <= 7.0 && env.sequential_media_optimized {
            1.0
        } else {
            0.75
        };

        self.viability_score = (self.viability_score * env_modifier).clamp(0.0, 100.0);

        // Pre-implantation stage progression
        if self.age_days >= 6.0 && self.stage < DevelopmentalStage::HatchedBlastocyst {
            self.stage = DevelopmentalStage::HatchedBlastocyst;
            self.total_cell_count = 256;
            self.icm_cell_count = 64;
            self.trophectoderm_cell_count = 192;
        } else if self.age_days >= 5.0 && self.stage < DevelopmentalStage::Blastocyst {
            if self.viability_score < 20.0 {
                self.is_arrested = true;
                return Err(BioError::DevelopmentArrested(
                    "Embryo arrested prior to blastocyst expansion due to low viability".to_string(),
                ));
            }
            self.stage = DevelopmentalStage::Blastocyst;
            self.total_cell_count = 128;
            self.icm_cell_count = 32;
            self.trophectoderm_cell_count = 96;
        } else if self.age_days >= 3.0 && self.stage < DevelopmentalStage::Morula {
            self.stage = DevelopmentalStage::Morula;
            self.total_cell_count = 32;
        } else if self.age_days >= 2.0 && self.stage < DevelopmentalStage::FourCell {
            self.stage = DevelopmentalStage::FourCell;
            self.total_cell_count = 4;
        } else if self.age_days >= 1.0 && self.stage < DevelopmentalStage::TwoCell {
            self.stage = DevelopmentalStage::TwoCell;
            self.total_cell_count = 2;
        }

        // Post-implantation / Artificial Womb Gestation Progression
        if self.age_days >= 110.0 {
            if env.artificial_womb_supported {
                self.stage = DevelopmentalStage::FullTermNeonate;
                self.total_cell_count = 2_000_000_000;
                self.fetal_weight_grams = 1850.0; // ~1.85 kg newborn Cave Lion cub!
            } else {
                self.is_arrested = true;
                return Err(BioError::DevelopmentArrested(
                    "Full-term gestation requires artificial womb support".to_string(),
                ));
            }
        } else if self.age_days >= 45.0 {
            if env.artificial_womb_supported {
                self.stage = DevelopmentalStage::FetalDevelopment;
                self.total_cell_count = 500_000_000;
                self.fetal_weight_grams = 350.0;
            }
        } else if self.age_days >= 15.0 {
            if env.artificial_womb_supported {
                self.stage = DevelopmentalStage::Organogenesis;
                self.total_cell_count = 10_000_000;
                self.fetal_weight_grams = 15.0;
            }
        } else if self.age_days >= 7.0 && self.stage <= DevelopmentalStage::HatchedBlastocyst {
            if env.artificial_womb_supported {
                self.stage = DevelopmentalStage::Gastrula;
                self.total_cell_count = 1024;
                self.fetal_weight_grams = 0.05;
            } else {
                self.is_arrested = true;
                return Err(BioError::DevelopmentArrested(
                    "Post-implantation gastrulation requires artificial womb or surrogate uterine support".to_string(),
                ));
            }
        }

        Ok(self.stage)
    }

    /// Evaluates implantation potential or readiness for Artificial Womb biobag cannulation.
    pub fn evaluate_implantation_potential(&self) -> f64 {
        if self.is_arrested || self.stage < DevelopmentalStage::Blastocyst {
            return 0.0;
        }
        let icm_ratio = self.icm_cell_count as f64 / self.total_cell_count.max(1) as f64;
        let quality_factor = self.viability_score / 100.0;
        let icm_score = if (0.15..=0.35).contains(&icm_ratio) { 1.0 } else { 0.5 };

        (quality_factor * icm_score * 100.0).clamp(0.0, 100.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ivg::epigenetics::EpigeneticProfile;
    use crate::ivg::gamete::GameteQuality;

    #[test]
    fn test_synthetic_fertilization_and_blastocyst_growth() {
        let seq_egg = Dna::new(b"ATGCGATCGATC").unwrap();
        let seq_sperm = Dna::new(b"GCTAGCTAGCTA").unwrap();

        let oocyte = SyntheticGamete {
            gamete_type: GameteType::Oocyte,
            haploid_genome: seq_egg,
            epigenetic_profile: EpigeneticProfile::mature_oocyte(),
            zona_pellucida: Some(crate::ivg::gamete::ZonaPellucida::cave_lion_default()),
            quality: GameteQuality {
                aneuploidy_risk: 0.02,
                imprint_fidelity: 0.98,
                mitochondrial_count: 100_000,
                polar_body_extruded: true,
                overall_score: 95.0,
            },
            tax_id: 9606,
        };

        let sperm = SyntheticGamete {
            gamete_type: GameteType::Spermatozoon,
            haploid_genome: seq_sperm,
            epigenetic_profile: EpigeneticProfile::mature_sperm(),
            zona_pellucida: None,
            quality: GameteQuality {
                aneuploidy_risk: 0.03,
                imprint_fidelity: 0.96,
                mitochondrial_count: 100,
                polar_body_extruded: false,
                overall_score: 92.0,
            },
            tax_id: 9606,
        };

        let mut embryo = SyntheticEmbryo::fertilize(oocyte, sperm).unwrap();
        assert_eq!(embryo.stage, DevelopmentalStage::Zygote);

        let env = EmbryoCultureEnvironment::default();
        // Advance 120 hours (5 days)
        let stage = embryo.advance_development(120.0, &env).unwrap();
        assert_eq!(stage, DevelopmentalStage::Blastocyst);
        assert_eq!(embryo.total_cell_count, 128);
        assert!(embryo.evaluate_implantation_potential() > 80.0);
    }
}
