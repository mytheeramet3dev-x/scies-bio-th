//! Tri-Layer Paleo-Genomic Reconstruction Engine & Chassis Edit Planner.
//!
//! Integrates Layer 1 (Phylogenetic ASR), Layer 2 (Paleoproteomics & Ancient DNA),
//! Layer 2.5 (Pseudogene Resurrection), and Layer 3 (Phenotypic Selection) to produce
//! a [`PaleoGenomeDraft`] complete with ranked [`ChassisGenomicEdit`] instructions
//! against the closest extant relative's genome.

use crate::paleo::phenotype::{PathwaySelectionProfile, PhenotypicEvidence};
use crate::paleo::resurrection::ResurrectedGene;
use crate::paleo::types::{translate_dna_codon, PaleoError, ReconstructedLocus, SiteSource};
use crate::seq::Dna;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Functional impact classification of a nucleotide edit on the extant chassis genome.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum EditEffectType {
    /// Reverses a premature nonsense stop codon (`TAA`/`TAG`/`TGA`) to restore a full-length protein.
    PseudogeneStopReversal,
    /// Missense change that alters the encoded amino acid ($dN$).
    NonSynonymous {
        chassis_aa: u8,
        ancestral_aa: u8,
    },
    /// Silent synonymous substitution that preserves the amino acid ($dS$).
    Synonymous {
        amino_acid: u8,
    },
    /// Non-coding or unaligned nucleotide substitution.
    NonCoding,
}

/// A single targeted nucleotide edit required to convert the extant chassis sequence
/// toward the reconstructed extinct locus.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ChassisGenomicEdit {
    /// Target gene symbol.
    pub gene_symbol: String,
    /// 0-based nucleotide coordinate within the gene locus.
    pub locus_position: usize,
    /// Original nucleotide in the extant relative (chassis).
    pub chassis_base: u8,
    /// Reconstructed nucleotide in the extinct target.
    pub target_base: u8,
    /// Posterior probability of `target_base` (`0.25..=1.0`).
    pub posterior_confidence: f64,
    /// Primary evidence source backing this site.
    pub evidence_source: SiteSource,
    /// Functional coding impact of the edit.
    pub effect: EditEffectType,
    /// Composite priority score combining posterior confidence, coding impact, and phenotypic pathway selection.
    pub priority_score: f64,
}

/// Global accuracy and provenance metrics across all reconstructed loci in a [`PaleoGenomeDraft`].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ReconstructionMetrics {
    pub total_bases_reconstructed: usize,
    pub overall_base_accuracy_expected: f64,
    pub fraction_anchored_by_adna: f64,
    pub fraction_anchored_by_proteomics: f64,
    pub fraction_inferred_by_asr: f64,
    pub fraction_pseudogene_repaired: f64,
    pub mean_shannon_entropy_bits: f64,
    pub total_chassis_edits: usize,
    pub non_synonymous_edits: usize,
    pub stop_reversal_edits: usize,
}

/// Complete paleo-genomic reconstruction draft and de-extinction chassis edit plan.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PaleoGenomeDraft {
    /// Name of the reconstructed extinct species (e.g., `"Mammuthus primigenius"`).
    pub target_species: String,
    /// Closest extant species used as the reference chassis (e.g., `"Elephas maximus"`).
    pub closest_extant_chassis: String,
    /// Reconstructed coding loci with per-site posterior distributions.
    pub reconstructed_loci: Vec<ReconstructedLocus>,
    /// Resurrected pseudogenes with repair audit logs.
    pub resurrected_genes: Vec<ResurrectedGene>,
    /// Pathway-level natural selection profile derived from osteology & isotope evidence.
    pub selection_profile: PathwaySelectionProfile,
    /// Ranked list of genomic edits to transform the extant chassis into the extinct target.
    pub chassis_edit_plan: Vec<ChassisGenomicEdit>,
    /// Summary reconstruction fidelity metrics.
    pub summary_metrics: ReconstructionMetrics,
}

/// Builder and orchestrator for Tri-Layer Paleo-Genomic Reconstruction.
pub struct PaleoReconstructionEngine {
    pub target_species: String,
    pub closest_extant_chassis: String,
    loci_with_chassis: Vec<(ReconstructedLocus, Dna)>,
    resurrected_with_chassis: Vec<(ResurrectedGene, Dna)>,
    phenotypic_evidences: Vec<Box<dyn PhenotypicEvidence>>,
}

impl PaleoReconstructionEngine {
    /// Creates a new reconstruction engine for `target_species` using `closest_extant_chassis`.
    pub fn new(
        target_species: impl Into<String>,
        closest_extant_chassis: impl Into<String>,
    ) -> Self {
        Self {
            target_species: target_species.into(),
            closest_extant_chassis: closest_extant_chassis.into(),
            loci_with_chassis: Vec::new(),
            resurrected_with_chassis: Vec::new(),
            phenotypic_evidences: Vec::new(),
        }
    }

    /// Registers a reconstructed locus along with the corresponding extant chassis sequence.
    pub fn add_reconstructed_locus(
        &mut self,
        locus: ReconstructedLocus,
        chassis_sequence: Dna,
    ) -> Result<(), PaleoError> {
        if locus.sites.len() != chassis_sequence.sequence.len() {
            return Err(PaleoError::LengthMismatch {
                expected: locus.sites.len(),
                found: chassis_sequence.sequence.len(),
            });
        }
        self.loci_with_chassis.push((locus, chassis_sequence));
        Ok(())
    }

    /// Registers a resurrected pseudogene along with its inactivated extant chassis sequence.
    pub fn add_resurrected_gene(
        &mut self,
        resurrected: ResurrectedGene,
        chassis_pseudogene: Dna,
    ) -> Result<(), PaleoError> {
        if resurrected.locus.sites.len() != chassis_pseudogene.sequence.len() {
            return Err(PaleoError::LengthMismatch {
                expected: resurrected.locus.sites.len(),
                found: chassis_pseudogene.sequence.len(),
            });
        }
        self.resurrected_with_chassis
            .push((resurrected, chassis_pseudogene));
        Ok(())
    }

    /// Adds a non-molecular paleontological evidence source (osteology, histology, isotopes, etc.).
    pub fn add_phenotypic_evidence<E: PhenotypicEvidence + 'static>(&mut self, evidence: E) {
        self.phenotypic_evidences.push(Box::new(evidence));
    }

    /// Executes the synthesis and generates the complete [`PaleoGenomeDraft`].
    pub fn build_draft(&self) -> Result<PaleoGenomeDraft, PaleoError> {
        if self.loci_with_chassis.is_empty() && self.resurrected_with_chassis.is_empty() {
            return Err(PaleoError::EmptyTaxonSet);
        }

        let ev_refs: Vec<&dyn PhenotypicEvidence> =
            self.phenotypic_evidences.iter().map(|b| b.as_ref()).collect();
        let selection_profile = PathwaySelectionProfile::from_evidences(&ev_refs);

        let mut all_loci = Vec::with_capacity(self.loci_with_chassis.len());
        let mut all_resurrected = Vec::with_capacity(self.resurrected_with_chassis.len());
        let mut chassis_edits = Vec::new();

        let mut total_bases = 0usize;
        let mut sum_conf = 0.0;
        let mut sum_entropy = 0.0;
        let mut count_adna = 0usize;
        let mut count_prot = 0usize;
        let mut count_asr = 0usize;
        let mut count_pseudo = 0usize;

        // Helper closure to diff a ReconstructedLocus against its chassis Dna
        let mut process_locus = |locus: &ReconstructedLocus, chassis: &Dna| {
            let len = locus.sites.len();
            total_bases += len;
            for (i, site) in locus.sites.iter().enumerate() {
                sum_conf += site.max_prob();
                sum_entropy += site.entropy_bits();
                match locus.provenance.site_sources[i] {
                    SiteSource::AncientDna => count_adna += 1,
                    SiteSource::Paleoproteomics => count_prot += 1,
                    SiteSource::AncestralReconstruction => count_asr += 1,
                    SiteSource::PseudogeneRepair => count_pseudo += 1,
                }
            }

            let pathway_weight = selection_profile.gene_priority_weight(&locus.gene_symbol);
            let consensus = locus.to_consensus_dna();
            let is_coding = len % 3 == 0;

            for pos in 0..len {
                let ch_base = chassis.sequence[pos].to_ascii_uppercase();
                let tg_base = consensus.sequence[pos].to_ascii_uppercase();
                if ch_base == tg_base {
                    continue;
                }

                let effect = if is_coding {
                    let codon_start = (pos / 3) * 3;
                    let ch_codon = [
                        chassis.sequence[codon_start].to_ascii_uppercase(),
                        chassis.sequence[codon_start + 1].to_ascii_uppercase(),
                        chassis.sequence[codon_start + 2].to_ascii_uppercase(),
                    ];
                    let tg_codon = [
                        consensus.sequence[codon_start].to_ascii_uppercase(),
                        consensus.sequence[codon_start + 1].to_ascii_uppercase(),
                        consensus.sequence[codon_start + 2].to_ascii_uppercase(),
                    ];
                    let ch_aa = translate_dna_codon(ch_codon);
                    let tg_aa = translate_dna_codon(tg_codon);
                    if ch_aa == b'*' && tg_aa != b'*' {
                        EditEffectType::PseudogeneStopReversal
                    } else if ch_aa != tg_aa {
                        EditEffectType::NonSynonymous {
                            chassis_aa: ch_aa,
                            ancestral_aa: tg_aa,
                        }
                    } else {
                        EditEffectType::Synonymous { amino_acid: tg_aa }
                    }
                } else {
                    EditEffectType::NonCoding
                };

                let impact_mult = match effect {
                    EditEffectType::PseudogeneStopReversal => 3.0,
                    EditEffectType::NonSynonymous { .. } => 2.0,
                    EditEffectType::NonCoding => 1.0,
                    EditEffectType::Synonymous { .. } => 0.4,
                };

                let site_conf = locus.sites[pos].max_prob();
                let priority_score = site_conf * impact_mult * pathway_weight;

                chassis_edits.push(ChassisGenomicEdit {
                    gene_symbol: locus.gene_symbol.clone(),
                    locus_position: pos,
                    chassis_base: ch_base,
                    target_base: tg_base,
                    posterior_confidence: site_conf,
                    evidence_source: locus.provenance.site_sources[pos],
                    effect,
                    priority_score,
                });
            }
        };

        for (locus, chassis) in &self.loci_with_chassis {
            process_locus(locus, chassis);
            all_loci.push(locus.clone());
        }

        for (res, chassis) in &self.resurrected_with_chassis {
            process_locus(&res.locus, chassis);
            all_resurrected.push(res.clone());
        }

        // Sort edits descending by priority_score, then ascending by gene_symbol and locus_position
        chassis_edits.sort_by(|a, b| {
            b.priority_score
                .partial_cmp(&a.priority_score)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.gene_symbol.cmp(&b.gene_symbol))
                .then_with(|| a.locus_position.cmp(&b.locus_position))
        });

        let non_synonymous_edits = chassis_edits
            .iter()
            .filter(|e| matches!(e.effect, EditEffectType::NonSynonymous { .. }))
            .count();
        let stop_reversal_edits = chassis_edits
            .iter()
            .filter(|e| matches!(e.effect, EditEffectType::PseudogeneStopReversal))
            .count();

        let denom = (total_bases.max(1)) as f64;
        let summary_metrics = ReconstructionMetrics {
            total_bases_reconstructed: total_bases,
            overall_base_accuracy_expected: sum_conf / denom,
            fraction_anchored_by_adna: count_adna as f64 / denom,
            fraction_anchored_by_proteomics: count_prot as f64 / denom,
            fraction_inferred_by_asr: count_asr as f64 / denom,
            fraction_pseudogene_repaired: count_pseudo as f64 / denom,
            mean_shannon_entropy_bits: sum_entropy / denom,
            total_chassis_edits: chassis_edits.len(),
            non_synonymous_edits,
            stop_reversal_edits,
        };

        Ok(PaleoGenomeDraft {
            target_species: self.target_species.clone(),
            closest_extant_chassis: self.closest_extant_chassis.clone(),
            reconstructed_loci: all_loci,
            resurrected_genes: all_resurrected,
            selection_profile,
            chassis_edit_plan: chassis_edits,
            summary_metrics,
        })
    }
}
