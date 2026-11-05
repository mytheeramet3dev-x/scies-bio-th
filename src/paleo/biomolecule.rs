//! Layer 2: Paleobiomolecular Anchoring (Paleoproteomics & Ancient DNA).
//!
//! Provides two sequence-level Bayesian anchoring engines:
//! 1. [`PaleoproteomeTranslator`]: Reverse-translates ancient peptide fragments (e.g., fossil
//!    bone collagen `COL1A1` or tooth enamelin `ENAM`) into nucleotide posteriors by combining
//!    Mass Spectrometry (MS/MS) amino acid likelihoods with phylogenetic priors and extant codon usage bias.
//! 2. [`AncientDnaAligner`]: Integrates degraded short ancient DNA (aDNA) reads using the Briggs
//!    post-mortem cytosine deamination model ($5'\ \text{C}\to\text{T}$ and $3'\ \text{G}\to\text{A}$ overhangs).

use crate::paleo::types::{
    base_to_index, index_to_base, translate_dna_codon, BasePosterior, PaleoError,
    ReconstructedLocus, SiteSource,
};
use crate::seq::Dna;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A fossil peptide fragment recovered via paleoproteomic mass spectrometry (LC-MS/MS).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PeptideFragment {
    /// 0-based codon index where this peptide aligns within the coding sequence.
    pub codon_offset: usize,
    /// Single-letter IUPAC amino acid sequence (e.g. `b"GPPGPQGPR"`).
    /// Note: `b'O'` (hydroxyproline, common in fossil collagen) is automatically mapped to canonical Proline (`b'P'`).
    pub amino_acids: Vec<u8>,
    /// Mass-spectrometry peptide-spectrum match (PSM) confidence $q \in (0.0, 1.0]$.
    pub ms_confidence: f64,
}

impl PeptideFragment {
    /// Creates a new `PeptideFragment`.
    pub fn new(codon_offset: usize, amino_acids: &[u8], ms_confidence: f64) -> Self {
        Self {
            codon_offset,
            amino_acids: amino_acids.to_vec(),
            ms_confidence: ms_confidence.clamp(0.01, 0.9999),
        }
    }
}

/// Bayesian paleoproteomic reverse-translation engine.
pub struct PaleoproteomeTranslator;

impl PaleoproteomeTranslator {
    /// Anchors a [`ReconstructedLocus`] using recovered fossil peptide fragments and optional
    /// extant relative reference coding sequence for codon-usage bias estimation.
    ///
    /// Returns the number of codons anchored by proteomic evidence.
    pub fn anchor_with_peptides(
        locus: &mut ReconstructedLocus,
        fragments: &[PeptideFragment],
        codon_usage_ref: Option<&Dna>,
    ) -> Result<usize, PaleoError> {
        if locus.sites.len() % 3 != 0 {
            return Err(PaleoError::FrameMismatch(locus.sites.len()));
        }
        let num_codons = locus.sites.len() / 3;
        let codon_weights = Self::compute_codon_bias_weights(codon_usage_ref);
        let mut anchored_codons = 0usize;

        for frag in fragments {
            let q = frag.ms_confidence.clamp(0.01, 0.9999);
            let err_bg = (1.0 - q) / 20.0;
            let hit_lh = q + err_bg;

            for (k, &raw_aa) in frag.amino_acids.iter().enumerate() {
                let codon_idx = frag.codon_offset + k;
                if codon_idx >= num_codons {
                    break;
                }
                // Map hydroxyproline 'O' -> Proline 'P'
                let obs_aa = match raw_aa.to_ascii_uppercase() {
                    b'O' => b'P',
                    other => other,
                };
                if obs_aa == b'X' || obs_aa == b'?' {
                    continue;
                }

                let base_i = codon_idx * 3;
                let p1 = locus.sites[base_i].probs;
                let p2 = locus.sites[base_i + 1].probs;
                let p3 = locus.sites[base_i + 2].probs;

                let mut m1 = [0.0; 4];
                let mut m2 = [0.0; 4];
                let mut m3 = [0.0; 4];

                // Enumerate all 64 codons
                for i1 in 0..4 {
                    for i2 in 0..4 {
                        for i3 in 0..4 {
                            let codon = [index_to_base(i1), index_to_base(i2), index_to_base(i3)];
                            let translated = translate_dna_codon(codon);
                            if translated == b'*' {
                                continue;
                            }
                            let flat_idx = (i1 << 4) | (i2 << 2) | i3;
                            let prior_c = p1[i1] * p2[i2] * p3[i3] * codon_weights[flat_idx];
                            let lh = if translated == obs_aa { hit_lh } else { err_bg };
                            let post_c = prior_c * lh;

                            m1[i1] += post_c;
                            m2[i2] += post_c;
                            m3[i3] += post_c;
                        }
                    }
                }

                locus.sites[base_i] = BasePosterior::from_weights(m1);
                locus.sites[base_i + 1] = BasePosterior::from_weights(m2);
                locus.sites[base_i + 2] = BasePosterior::from_weights(m3);

                locus.provenance.site_sources[base_i] = SiteSource::Paleoproteomics;
                locus.provenance.site_sources[base_i + 1] = SiteSource::Paleoproteomics;
                locus.provenance.site_sources[base_i + 2] = SiteSource::Paleoproteomics;
                anchored_codons += 1;
            }
        }

        locus.refresh_metrics();
        Ok(anchored_codons)
    }

    /// Computes relative synonymous codon usage weights (64 entries) with Laplace pseudocounts.
    fn compute_codon_bias_weights(reference: Option<&Dna>) -> [f64; 64] {
        let mut counts = [1.0f64; 64]; // Laplace smoothing pseudocount = 1.0
        if let Some(dna) = reference {
            for chunk in dna.sequence.chunks_exact(3) {
                if let (Some(i1), Some(i2), Some(i3)) = (
                    base_to_index(chunk[0]),
                    base_to_index(chunk[1]),
                    base_to_index(chunk[2]),
                ) {
                    counts[(i1 << 4) | (i2 << 2) | i3] += 1.0;
                }
            }
        }
        // Normalize within each synonymous amino acid family so codon usage only disambiguates synonymous choices
        let mut weights = [1.0f64; 64];
        for i1 in 0..4 {
            for i2 in 0..4 {
                for i3 in 0..4 {
                    let c_idx = (i1 << 4) | (i2 << 2) | i3;
                    let aa =
                        translate_dna_codon([index_to_base(i1), index_to_base(i2), index_to_base(i3)]);
                    let mut family_sum = 0.0;
                    let mut family_size = 0.0;
                    for j1 in 0..4 {
                        for j2 in 0..4 {
                            for j3 in 0..4 {
                                let aa2 = translate_dna_codon([
                                    index_to_base(j1),
                                    index_to_base(j2),
                                    index_to_base(j3),
                                ]);
                                if aa2 == aa {
                                    family_sum += counts[(j1 << 4) | (j2 << 2) | j3];
                                    family_size += 1.0;
                                }
                            }
                        }
                    }
                    if family_sum > 0.0 {
                        weights[c_idx] = (counts[c_idx] / family_sum) * family_size;
                    }
                }
            }
        }
        weights
    }
}

/// A short, degraded ancient DNA (aDNA) sequencing read aligned to a locus.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct AncientDnaRead {
    /// 0-based nucleotide start coordinate in the locus.
    pub base_offset: usize,
    /// Observed read bases (e.g. `30..100` bp).
    pub sequence: Vec<u8>,
    /// Sequencing base-call accuracy ($1 - \epsilon$, typically `0.99` to `0.999`).
    pub base_accuracy: f64,
    /// 5'-overhang cytosine deamination rate ($\text{C}\to\text{T}$ probability at position 0, e.g. `0.25`).
    pub deamination_5prime_max: f64,
    /// 3'-overhang guanine deamination complement rate ($\text{G}\to\text{A}$ probability at 3' end, e.g. `0.20`).
    pub deamination_3prime_max: f64,
}

impl AncientDnaRead {
    /// Creates a new `AncientDnaRead` with standard Pleistocene permafrost damage parameters.
    pub fn new_pleistocene(base_offset: usize, sequence: &[u8]) -> Self {
        Self {
            base_offset,
            sequence: sequence.to_vec(),
            base_accuracy: 0.995,
            deamination_5prime_max: 0.25,
            deamination_3prime_max: 0.20,
        }
    }
}

/// Deamination-aware Ancient DNA alignment and Bayesian integration engine.
pub struct AncientDnaAligner;

impl AncientDnaAligner {
    /// Decay constant $\lambda$ for post-mortem single-stranded overhang damage from read ends.
    const OVERHANG_DECAY_LAMBDA: f64 = 0.35;

    /// Integrates aligned `AncientDnaRead` fragments into `locus` using the Briggs post-mortem
    /// deamination likelihood model.
    ///
    /// Returns the number of nucleotide read positions integrated.
    pub fn integrate_reads(
        locus: &mut ReconstructedLocus,
        reads: &[AncientDnaRead],
    ) -> Result<usize, PaleoError> {
        let locus_len = locus.sites.len();
        let mut integrated_bases = 0usize;

        for read in reads {
            let read_len = read.sequence.len();
            if read_len == 0 {
                continue;
            }
            let acc = read.base_accuracy.clamp(0.50, 0.9999);
            let seq_err = (1.0 - acc) / 3.0;

            for (read_pos, &obs_b) in read.sequence.iter().enumerate() {
                let locus_pos = read.base_offset + read_pos;
                if locus_pos >= locus_len {
                    break;
                }
                let Some(obs_idx) = base_to_index(obs_b) else {
                    continue;
                };

                let dist_5p = read_pos as f64;
                let dist_3p = (read_len - 1 - read_pos) as f64;

                let delta_5p = (read.deamination_5prime_max.clamp(0.0, 0.80))
                    * (-Self::OVERHANG_DECAY_LAMBDA * dist_5p).exp();
                let delta_3p = (read.deamination_3prime_max.clamp(0.0, 0.80))
                    * (-Self::OVERHANG_DECAY_LAMBDA * dist_3p).exp();

                // Compute likelihood P(obs_idx | true_idx) for each true base in [A=0, C=1, G=2, T=3]
                let mut updated_weights = [0.0; 4];
                for true_idx in 0..4 {
                    let lh = Self::observation_likelihood(
                        true_idx, obs_idx, acc, seq_err, delta_5p, delta_3p,
                    );
                    updated_weights[true_idx] = locus.sites[locus_pos].probs[true_idx] * lh;
                }

                locus.sites[locus_pos] = BasePosterior::from_weights(updated_weights);
                locus.provenance.site_sources[locus_pos] = SiteSource::AncientDna;
                integrated_bases += 1;
            }
        }

        locus.refresh_metrics();
        Ok(integrated_bases)
    }

    /// Computes $P(\text{observed base} \mid \text{true ancestral base})$ accounting for both
    /// post-mortem cytosine deamination ($C \to T$ at 5', $G \to A$ at 3') and sequencing error.
    #[inline]
    fn observation_likelihood(
        true_idx: usize,
        obs_idx: usize,
        acc: f64,
        seq_err: f64,
        delta_5p: f64,
        delta_3p: f64,
    ) -> f64 {
        // Post-mortem chemical transition matrix D[true_idx][post_damage_idx]
        let mut d_row = [0.0; 4];
        match true_idx {
            // C (1) -> T (3) with probability delta_5p
            1 => {
                d_row[1] = 1.0 - delta_5p;
                d_row[3] = delta_5p;
            }
            // G (2) -> A (0) with probability delta_3p
            2 => {
                d_row[2] = 1.0 - delta_3p;
                d_row[0] = delta_3p;
            }
            other => {
                d_row[other] = 1.0;
            }
        }

        // Convolve with sequencing error channel
        let mut total_lh = 0.0;
        for dmg_idx in 0..4 {
            let p_dmg = d_row[dmg_idx];
            if p_dmg > 0.0 {
                let p_seq = if dmg_idx == obs_idx { acc } else { seq_err };
                total_lh += p_dmg * p_seq;
            }
        }
        total_lh
    }
}
