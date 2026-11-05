//! Layer 2.5: Pseudogene Resurrection Engine.
//!
//! Modern descendants frequently retain inactivated remnants (pseudogenes) of genes that were
//! functional in their extinct ancestors — for example, birds (`Gallus gallus`) carry mutated
//! tooth-enamel genes (`ENAM`, `AMELX`, `MMP20`) interrupted by premature stop codons (`TGA`/`TAA`/`TAG`).
//!
//! [`PseudogeneResurrector`] combines the high-proximity sequence of the inactivated ingroup
//! with functional outgroup sequences (e.g., Crocodilians and Testudines) to repair premature
//! stop codons and ambiguous reading-frame defects back to their ancestral functional coding state.

use crate::paleo::asr::{ExtantSequence, SubstitutionModel};
use crate::paleo::types::{
    base_to_index, index_to_base, translate_dna_codon, BasePosterior, LocusProvenance, PaleoError,
    ReconstructedLocus, SiteSource,
};
use crate::seq::Dna;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Category of molecular repair applied to resurrect an inactivated pseudogene.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum RepairType {
    /// Reversal of an internal nonsense (premature stop `TAA`/`TAG`/`TGA`) codon to an ancestral sense codon.
    PrematureStopReversal,
    /// Restoration of an ambiguous or frameshift-damaged triplet (`N`) using functional outgroups.
    FrameshiftDefectRepair,
}

/// Audit record of a single codon repair performed during pseudogene resurrection.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct RepairEvent {
    /// 0-based codon index in the coding sequence.
    pub codon_index: usize,
    /// 0-based nucleotide offset (`codon_index * 3`).
    pub base_offset: usize,
    /// Original broken triplet in the extant pseudogene.
    pub pseudogene_codon: [u8; 3],
    /// Repaired ancestral sense codon.
    pub resurrected_codon: [u8; 3],
    /// Restored single-letter amino acid.
    pub restored_amino_acid: u8,
    /// Classification of the repair.
    pub repair_type: RepairType,
    /// Posterior probability of the restored sense codon.
    pub posterior_confidence: f64,
}

/// Output of resurrecting an inactivated pseudogene into a functional ancestral coding sequence.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ResurrectedGene {
    /// Reconstructed functional coding locus with per-site posteriors.
    pub locus: ReconstructedLocus,
    /// Ordered list of molecular repair events applied to restore open reading frame integrity.
    pub repairs: Vec<RepairEvent>,
}

/// Pseudogene resurrection engine using Extant Phylogenetic Bracketing.
pub struct PseudogeneResurrector;

impl PseudogeneResurrector {
    /// Repairs an inactivated pseudogene from a close extant relative (`inactivated_ingroup`)
    /// using one or more functional outgroup sequences (`functional_outgroups`, paired with their
    /// evolutionary divergence time in Ma).
    pub fn resurrect_inactivated_gene(
        gene_symbol: &str,
        target_taxon: &str,
        inactivated_ingroup: &Dna,
        ingroup_branch_ma: f64,
        functional_outgroups: &[(ExtantSequence, f64)],
        model: &SubstitutionModel,
    ) -> Result<ResurrectedGene, PaleoError> {
        if functional_outgroups.is_empty() {
            return Err(PaleoError::EmptyTaxonSet);
        }
        let seq_len = inactivated_ingroup.sequence.len();
        if seq_len == 0 || seq_len % 3 != 0 {
            return Err(PaleoError::FrameMismatch(seq_len));
        }
        for (out_seq, _) in functional_outgroups {
            if out_seq.sequence.sequence.len() != seq_len {
                return Err(PaleoError::LengthMismatch {
                    expected: seq_len,
                    found: out_seq.sequence.sequence.len(),
                });
            }
        }

        let num_codons = seq_len / 3;
        let ingroup_p = model.transition_matrix(ingroup_branch_ma.max(1.0));
        let outgroup_mats: Vec<[[f64; 4]; 4]> = functional_outgroups
            .iter()
            .map(|(_, t)| model.transition_matrix(t.max(1.0)))
            .collect();

        let mut sites = Vec::with_capacity(seq_len);
        let mut provenance =
            LocusProvenance::new_uniform(seq_len, SiteSource::AncestralReconstruction);
        let mut repairs = Vec::new();

        for codon_idx in 0..num_codons {
            let base_i = codon_idx * 3;
            let in_c = [
                inactivated_ingroup.sequence[base_i].to_ascii_uppercase(),
                inactivated_ingroup.sequence[base_i + 1].to_ascii_uppercase(),
                inactivated_ingroup.sequence[base_i + 2].to_ascii_uppercase(),
            ];
            let is_terminal = codon_idx + 1 == num_codons;
            let in_aa = translate_dna_codon(in_c);
            let is_premature_stop = !is_terminal && in_aa == b'*';
            let is_damaged_triplet = in_aa == b'X';

            // Compute single-site likelihoods from ingroup + outgroups
            let mut site_lh = [[1.0f64; 4]; 3];
            for offset in 0..3 {
                let pos = base_i + offset;
                // Ingroup contribution (down-weighted if damaged/stop)
                if let Some(in_b_idx) = base_to_index(in_c[offset]) {
                    for anc_b in 0..4 {
                        site_lh[offset][anc_b] *= ingroup_p[anc_b][in_b_idx];
                    }
                }
                // Outgroups contribution
                for (out_i, (out_ext, _)) in functional_outgroups.iter().enumerate() {
                    if let Some(out_b_idx) = base_to_index(out_ext.sequence.sequence[pos]) {
                        for anc_b in 0..4 {
                            site_lh[offset][anc_b] *= outgroup_mats[out_i][anc_b][out_b_idx];
                        }
                    }
                }
            }

            if is_premature_stop || is_damaged_triplet {
                // Condition over the 61 sense codons (excluding TAA, TAG, TGA)
                let mut m1 = [0.0; 4];
                let mut m2 = [0.0; 4];
                let mut m3 = [0.0; 4];
                let mut best_codon = [b'A', b'A', b'A'];
                let mut best_prob = -1.0;
                let mut total_sense_prob = 0.0;

                for i1 in 0..4 {
                    for i2 in 0..4 {
                        for i3 in 0..4 {
                            let cand = [index_to_base(i1), index_to_base(i2), index_to_base(i3)];
                            let aa = translate_dna_codon(cand);
                            if aa == b'*' {
                                continue; // Exclude nonsense codons from functional ancestral state
                            }
                            let p = site_lh[0][i1] * site_lh[1][i2] * site_lh[2][i3];
                            total_sense_prob += p;
                            m1[i1] += p;
                            m2[i2] += p;
                            m3[i3] += p;
                            if p > best_prob {
                                best_prob = p;
                                best_codon = cand;
                            }
                        }
                    }
                }

                let codon_conf = if total_sense_prob > 0.0 {
                    (best_prob / total_sense_prob).clamp(0.0, 1.0)
                } else {
                    0.25
                };

                sites.push(BasePosterior::from_weights(m1));
                sites.push(BasePosterior::from_weights(m2));
                sites.push(BasePosterior::from_weights(m3));

                provenance.site_sources[base_i] = SiteSource::PseudogeneRepair;
                provenance.site_sources[base_i + 1] = SiteSource::PseudogeneRepair;
                provenance.site_sources[base_i + 2] = SiteSource::PseudogeneRepair;

                repairs.push(RepairEvent {
                    codon_index: codon_idx,
                    base_offset: base_i,
                    pseudogene_codon: in_c,
                    resurrected_codon: best_codon,
                    restored_amino_acid: translate_dna_codon(best_codon),
                    repair_type: if is_premature_stop {
                        RepairType::PrematureStopReversal
                    } else {
                        RepairType::FrameshiftDefectRepair
                    },
                    posterior_confidence: codon_conf,
                });
            } else {
                sites.push(BasePosterior::from_weights(site_lh[0]));
                sites.push(BasePosterior::from_weights(site_lh[1]));
                sites.push(BasePosterior::from_weights(site_lh[2]));
            }
        }

        let locus = ReconstructedLocus::new(gene_symbol, target_taxon, sites, provenance);
        Ok(ResurrectedGene { locus, repairs })
    }
}
