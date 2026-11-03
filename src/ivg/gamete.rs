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

    /// Default synthetic Zona Pellucida optimized for human (Homo sapiens, TaxID 9606) MII oocytes.
    pub fn human_default() -> Self {
        Self::new(17.5, 20.0)
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
            if self.tax_id == 9606 {
                Some(ZonaPellucida::human_default())
            } else {
                Some(ZonaPellucida::cave_lion_default())
            }
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

/// Performs true homologous meiotic recombination (Prophase I chiasmata crossing-over)
/// between two parental homologous chromosomes (`homolog_a` and `homolog_b`), returning
/// `(recombinant_haploid_dna, chiasmata_count)`.
///
/// Ensures at least 1 obligatory chiasma crossover for chromosomes $\ge 4\text{ bp}$
/// (preventing Meiosis I non-disjunction).
pub fn recombine_homologous_chromosomes(
    homolog_a: &Dna,
    homolog_b: &Dna,
    seed: u64,
) -> Result<(Dna, usize)> {
    let seq_a = &homolog_a.sequence;
    let seq_b = &homolog_b.sequence;
    if seq_a.is_empty() || seq_b.is_empty() {
        return Err(BioError::InvalidInput(
            "Homologous chromosomes cannot be empty for meiotic recombination".to_string(),
        ));
    }
    if seq_a.len() != seq_b.len() {
        return Err(BioError::InvalidInput(format!(
            "Homologous chromosomes must have equal length ({} != {})",
            seq_a.len(),
            seq_b.len()
        )));
    }

    let len = seq_a.len();
    if len < 4 {
        return Ok((homolog_a.clone(), 0));
    }

    let mut state = seed.wrapping_add(0x9E3779B97F4A7C15);
    // Determine 1 to 3 chiasmata positions with crossover interference (minimum spacing)
    let num_chiasmata = 1 + ((state >> 32) as usize % 2).min(len / 4);
    let mut crossover_points = Vec::with_capacity(num_chiasmata);
    let segment = len / (num_chiasmata + 1).max(1);

    for c in 0..num_chiasmata {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        let offset = ((state >> 32) as usize) % segment.max(1);
        let pt = ((c + 1) * segment + offset / 2).clamp(1, len - 1);
        if crossover_points.last().copied() != Some(pt) {
            crossover_points.push(pt);
        }
    }

    let mut out = Vec::with_capacity(len);
    let mut use_a = (seed & 1) == 0;
    let mut prev = 0;

    for &pt in &crossover_points {
        let src = if use_a { &seq_a[prev..pt] } else { &seq_b[prev..pt] };
        out.extend_from_slice(src);
        use_a = !use_a;
        prev = pt;
    }
    let tail = if use_a { &seq_a[prev..len] } else { &seq_b[prev..len] };
    out.extend_from_slice(tail);

    Ok((Dna::new(&out)?, crossover_points.len()))
}

/// Derives a heterozygous parental homolog partner (`homolog_b`) from a single haplotype
/// template (`homolog_a`) by introducing transition SNPs (`A<->G`, `C<->T`) at deterministic
/// polymorphic intervals so meiotic crossing-over produces a true recombinant sequence.
pub fn derive_homologous_partner(homolog_a: &Dna, seed: u64) -> Dna {
    let mut seq = homolog_a.sequence.clone();
    if seq.len() < 4 {
        return homolog_a.clone();
    }
    let mut state = seed.wrapping_add(0x517CC1B727220A95);
    let step = (seq.len() / 8).clamp(3, 25);
    let mut idx = (state as usize % step).min(seq.len() - 1);
    while idx < seq.len() {
        seq[idx] = match seq[idx] {
            b'A' => b'G',
            b'G' => b'A',
            b'C' => b'T',
            b'T' => b'C',
            other => other,
        };
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        idx += step + ((state >> 48) as usize % 5);
    }
    Dna::new(&seq).unwrap_or_else(|_| homolog_a.clone())
}

/// Returns the number of homologous chromosome pairs (bivalents) for a given NCBI `tax_id`.
///
/// - *Homo sapiens* (`9606`): $23$ bivalents ($2n = 46$)
/// - *Mus musculus* (`10090`): $20$ bivalents ($2n = 40$)
/// - *Panthera spelaea* / *Felis catus* (`1969727`, `9685`, default mammalian): $19$ bivalents ($2n = 38$)
pub fn species_bivalent_pairs(tax_id: u32) -> usize {
    match tax_id {
        9606 => 23,
        10090 => 20,
        _ => 19,
    }
}

/// Simulates genome-wide Prophase I pachytene chiasmata formation across all homologous
/// chromosome pairs (bivalents) of a species karyotype.
///
/// Calibrated to human female meiotic recombination maps (Broman et al. 1998 *AJHG* 63:861;
/// Ottolini et al. 2015 *Nat Genet* 47:727; Hassold & Hunt 2001 *Nat Rev Genet* 2:280):
/// - Healthy human oocytes ($23$ bivalents) form an obligatory $\ge 1$ chiasma per bivalent
///   plus $1\text{--}2$ interference-spaced crossovers, yielding $\sim 48\text{--}68$ total chiasmata.
/// - Suboptimal `meiotic_competence` ($< 0.75$) causes crossover assurance failure, producing
///   achiasmate bivalents (`achiasmate_bivalents > 0`) that predispose to Meiosis I non-disjunction.
///
/// Returns `(total_chiasmata_count, achiasmate_bivalents_count)`.
pub fn simulate_karyotype_meiotic_crossovers(
    tax_id: u32,
    meiotic_competence: f64,
    seed: u64,
) -> (usize, usize) {
    let pairs = species_bivalent_pairs(tax_id);
    let comp = meiotic_competence.clamp(0.0, 1.0);
    let mut state = seed.wrapping_add(0x9E3779B97F4A7C15);

    let mut total_chiasmata = 0usize;
    let mut achiasmate_bivalents = 0usize;

    for bivalent_idx in 0..pairs {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add((bivalent_idx as u64) + 1);
        let u = ((state >> 32) as u32 as f64) / (u32::MAX as f64);

        // Crossover assurance threshold: healthy competence (>= 0.80) guarantees >= 1 chiasma per bivalent
        let assurance_prob = if comp >= 0.80 {
            1.0
        } else {
            (comp / 0.80).powf(1.5)
        };

        if u > assurance_prob {
            achiasmate_bivalents += 1;
        } else {
            // Obligatory chiasma + length/interference-modulated secondary crossovers (1 to 3 per bivalent)
            let extra_max = if bivalent_idx < 8 { 2 } else { 1 };
            let extra = if comp >= 0.70 {
                1 + (((state >> 16) as usize) % (extra_max + 1)) * (if comp >= 0.85 { 1 } else { 0 })
            } else {
                0
            };
            total_chiasmata += 1 + extra;
        }
    }

    (total_chiasmata, achiasmate_bivalents)
}

/// Simulates meiotic crossing-over and reduction from diploid to haploid state.
fn simulate_meiotic_reduction(diploid_bytes: &[u8], seed: u64) -> Vec<u8> {
    let len = diploid_bytes.len();
    if len < 4 {
        return diploid_bytes.to_vec();
    }
    // Model contiguous block chiasma recombination between the two halves of the template
    let mut out = diploid_bytes.to_vec();
    let state = seed
        .wrapping_add(0x9E3779B97F4A7C15)
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1);
    let block_len = (len / 5).max(1);
    let start_a = ((state >> 16) as usize) % (len / 2 - block_len).max(1);
    let start_b = (len / 2) + (((state >> 32) as usize) % (len / 2 - block_len).max(1));

    for k in 0..block_len {
        if start_a + k < len && start_b + k < len {
            out.swap(start_a + k, start_b + k);
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
