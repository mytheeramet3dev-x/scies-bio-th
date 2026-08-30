//! Somatic Cell Nuclear Transfer (SCNT) cloning and nuclear reprogramming engine.
//!
//! Models nuclear enucleation, electrofusion/microinjection, ooplasm-driven chromatin remodeling,
//! and mitochondrial heteroplasmy.

use crate::error::{BioError, Result};
use crate::ivg::epigenetics::EpigeneticProfile;
use crate::ivg::gamete::{GameteType, SyntheticGamete};
use crate::seq::Dna;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Somatic cell donating the diploid nuclear genome for cloning.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SomaticDonorCell {
    /// Full nuclear genomic sequence.
    pub nuclear_dna: Dna,
    /// NCBI Taxonomy ID of donor organism.
    pub tax_id: u32,
    /// Epigenetic status of donor somatic nucleus (carries somatic H3K9me3/5mC barriers).
    pub epigenetic_profile: EpigeneticProfile,
    /// Cell lineage description (e.g. "Dermal Fibroblast", "Ear Skin Biopsy").
    pub cell_type: String,
}

impl SomaticDonorCell {
    /// Creates a somatic donor cell from a DNA sequence and donor age.
    pub fn new(nuclear_dna: Dna, tax_id: u32, donor_age_years: f64, cell_type: &str) -> Self {
        Self {
            nuclear_dna,
            tax_id,
            epigenetic_profile: EpigeneticProfile::somatic_default(donor_age_years),
            cell_type: cell_type.to_string(),
        }
    }
}

/// Recipient oocyte after physical enucleation (maternal chromosomes removed).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct EnucleatedOocyte {
    /// NCBI Taxonomy ID of the oocyte donor species (surrogate cytoplast provider).
    pub cytoplast_tax_id: u32,
    /// Oocyte maternal mitochondrial genome.
    pub mitochondrial_dna: Dna,
    /// Maternal reprogramming competence factor of the ooplasm (0.0 to 1.0).
    pub ooplasm_competence: f64,
}

impl EnucleatedOocyte {
    /// Enucleates a mature synthetic oocyte, stripping its nuclear genome but retaining ooplasm and mtDNA.
    pub fn from_oocyte(oocyte: &SyntheticGamete) -> Result<Self> {
        if oocyte.gamete_type != GameteType::Oocyte {
            return Err(BioError::IncompatibleGamete(
                "Enucleation recipient must be an oocyte (egg), not a sperm".to_string(),
            ));
        }

        Ok(Self {
            cytoplast_tax_id: oocyte.tax_id,
            mitochondrial_dna: oocyte.haploid_genome.clone(), // Default proxy for mtDNA
            ooplasm_competence: (oocyte.quality.overall_score / 100.0).clamp(0.1, 1.0),
        })
    }
}

/// Mode of artificial oocyte activation after nuclear transfer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ActivationMethod {
    /// Strontium Chloride (SrCl2) + Ionomycin chemical calcium oscillation.
    StrontiumIonomycin,
    /// High-voltage DC electrofusion pulse.
    Electrofusion,
    /// Combined electrical pulse + 6-DMAP protein kinase inhibition.
    ElectrofusionPlusKinaseInhibitor,
}

/// SCNT Reconstructed zygote undergoing epigenetic nuclear reprogramming.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ScntConstruct {
    /// Transferred somatic nuclear genome.
    pub nuclear_dna: Dna,
    /// Donor organism TaxID.
    pub donor_tax_id: u32,
    /// Recipient oocyte cytoplast.
    pub cytoplast: EnucleatedOocyte,
    /// Reprogramming epigenetic profile.
    pub epigenetic_profile: EpigeneticProfile,
    /// Percentage of donor mitochondrial carry-over (Heteroplasmy, typically < 1-3%).
    pub donor_mt_dna_percent: f64,
    /// Fusion and activation efficiency (0.0 to 1.0).
    pub activation_efficiency: f64,
}

impl ScntConstruct {
    /// Fuses a somatic donor nucleus into an enucleated oocyte cytoplast.
    pub fn reconstruct(
        donor: SomaticDonorCell,
        cytoplast: EnucleatedOocyte,
        activation: ActivationMethod,
        apply_kdm4_demethylase_treatment: bool,
    ) -> Result<Self> {
        let is_cross_species = donor.tax_id != cytoplast.cytoplast_tax_id;
        
        let base_activation = match activation {
            ActivationMethod::StrontiumIonomycin => 0.85,
            ActivationMethod::Electrofusion => 0.78,
            ActivationMethod::ElectrofusionPlusKinaseInhibitor => 0.92,
        };

        // Cross-species nuclear transfer incurs an incompatibility penalty on ooplasm reprogramming
        let species_penalty = if is_cross_species { 0.70 } else { 1.0 };
        let final_activation = base_activation * cytoplast.ooplasm_competence * species_penalty;

        let mut ep_profile = donor.epigenetic_profile;
        if apply_kdm4_demethylase_treatment {
            // Apply KDM4A/KDM4D histone demethylase mRNA injection to remove somatic H3K9me3 blocks
            ep_profile.apply_kdm4_treatment(0.90);
        }

        Ok(Self {
            nuclear_dna: donor.nuclear_dna,
            donor_tax_id: donor.tax_id,
            cytoplast,
            epigenetic_profile: ep_profile,
            donor_mt_dna_percent: 1.5, // Standard ~1.5% somatic mtDNA contamination in SCNT
            activation_efficiency: final_activation,
        })
    }

    /// Simulates maternal ooplasmic reprogramming of the transferred somatic chromatin.
    ///
    /// Returns the predicted blastocyst formation competence rate (0.0 to 1.0).
    pub fn predict_blastocyst_competence(&self) -> f64 {
        let barrier = self.epigenetic_profile.reprogramming_barrier_score();
        // Lower barrier = higher blastocyst rate
        let barrier_factor = (100.0 - barrier) / 100.0;
        let competence = self.activation_efficiency * barrier_factor.max(0.05);
        competence.clamp(0.01, 0.95)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ivg::gamete::{GameteQuality, GameteType};

    #[test]
    fn test_scnt_reconstruction_and_kdm4_enhancement() {
        let donor_dna = Dna::new(b"ATGCGATCGATC").unwrap();
        let donor = SomaticDonorCell::new(donor_dna.clone(), 9689, 10.0, "Lion Fibroblast");

        let oocyte_gamete = SyntheticGamete {
            gamete_type: GameteType::Oocyte,
            haploid_genome: donor_dna.clone(),
            epigenetic_profile: EpigeneticProfile::mature_oocyte(),
            zona_pellucida: Some(crate::ivg::gamete::ZonaPellucida::cave_lion_default()),
            quality: GameteQuality {
                aneuploidy_risk: 0.05,
                imprint_fidelity: 0.95,
                mitochondrial_count: 100_000,
                polar_body_extruded: true,
                overall_score: 90.0,
            },
            tax_id: 9689,
        };

        let enucleated = EnucleatedOocyte::from_oocyte(&oocyte_gamete).unwrap();

        // Without KDM4
        let scnt_standard = ScntConstruct::reconstruct(
            donor.clone(),
            enucleated.clone(),
            ActivationMethod::ElectrofusionPlusKinaseInhibitor,
            false,
        ).unwrap();

        // With KDM4 demethylase treatment
        let scnt_enhanced = ScntConstruct::reconstruct(
            donor,
            enucleated,
            ActivationMethod::ElectrofusionPlusKinaseInhibitor,
            true,
        ).unwrap();

        assert!(scnt_enhanced.predict_blastocyst_competence() > scnt_standard.predict_blastocyst_competence());
    }
}
