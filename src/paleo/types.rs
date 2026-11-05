//! Core types, probabilistic nucleotide distributions, and error handling for paleo-genomics.

use crate::seq::{Dna, Protein};
use std::fmt;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Canonical nucleotide order used across all posterior arrays: `[A, C, G, T]`.
pub const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Maps an ASCII DNA base (`A`, `C`, `G`, `T`, case-insensitive) to index `0..4`.
#[inline]
pub fn base_to_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

/// Maps index `0..4` to canonical uppercase DNA base (`A`, `C`, `G`, `T`).
#[inline]
pub fn index_to_base(idx: usize) -> u8 {
    BASES[idx & 3]
}

/// Translates a 3-base DNA sense/stop codon (`[b1, b2, b3]`) into its single-letter IUPAC amino acid
/// (or `b'*'` for stop codons, `b'X'` for ambiguous bases).
pub fn translate_dna_codon(codon: [u8; 3]) -> u8 {
    let c0 = codon[0].to_ascii_uppercase();
    let c1 = codon[1].to_ascii_uppercase();
    let c2 = codon[2].to_ascii_uppercase();
    match &[c0, c1, c2] {
        b"TTT" | b"TTC" => b'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => b'L',
        b"ATT" | b"ATC" | b"ATA" => b'I',
        b"ATG" => b'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => b'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => b'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => b'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => b'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => b'A',
        b"TAT" | b"TAC" => b'Y',
        b"TAA" | b"TAG" | b"TGA" => b'*',
        b"CAT" | b"CAC" => b'H',
        b"CAA" | b"CAG" => b'Q',
        b"AAT" | b"AAC" => b'N',
        b"AAA" | b"AAG" => b'K',
        b"GAT" | b"GAC" => b'D',
        b"GAA" | b"GAG" => b'E',
        b"TGT" | b"TGC" => b'C',
        b"TGG" => b'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => b'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => b'G',
        _ => b'X',
    }
}

/// Errors arising during paleo-genomic and paleoproteomic reconstruction.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum PaleoError {
    /// No extant reference sequences were provided.
    EmptyTaxonSet,
    /// Aligned sequences have mismatched lengths.
    LengthMismatch { expected: usize, found: usize },
    /// A leaf taxon in the phylogenetic tree was not found in the sequence set.
    TaxonNotFound(String),
    /// Invalid model or evidence parameter.
    InvalidParameter(String),
    /// Coding sequence length is not a multiple of 3 (invalid reading frame).
    FrameMismatch(usize),
}

impl fmt::Display for PaleoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyTaxonSet => write!(f, "PaleoError: empty extant taxon set"),
            Self::LengthMismatch { expected, found } => write!(
                f,
                "PaleoError: sequence length mismatch (expected {expected} bp, found {found} bp)"
            ),
            Self::TaxonNotFound(name) => {
                write!(f, "PaleoError: leaf taxon '{name}' not found in sequence alignment")
            }
            Self::InvalidParameter(msg) => write!(f, "PaleoError: invalid parameter ({msg})"),
            Self::FrameMismatch(len) => write!(
                f,
                "PaleoError: coding sequence length ({len} bp) is not divisible by 3"
            ),
        }
    }
}

impl std::error::Error for PaleoError {}

/// Posterior probability distribution over `[A, C, G, T]` at a single nucleotide site.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct BasePosterior {
    /// Normalized probabilities `[P(A), P(C), P(G), P(T)]`.
    pub probs: [f64; 4],
}

impl BasePosterior {
    /// Creates a flat uniform distribution `[0.25, 0.25, 0.25, 0.25]`.
    pub fn uniform() -> Self {
        Self { probs: [0.25; 4] }
    }

    /// Creates a normalized `BasePosterior` from raw unnormalized weights.
    pub fn from_weights(weights: [f64; 4]) -> Self {
        let mut clean = [0.0; 4];
        let mut sum = 0.0;
        for i in 0..4 {
            let w = if weights[i].is_finite() && weights[i] > 0.0 {
                weights[i]
            } else {
                0.0
            };
            clean[i] = w;
            sum += w;
        }
        if sum <= 1e-15 {
            Self::uniform()
        } else {
            Self {
                probs: [
                    clean[0] / sum,
                    clean[1] / sum,
                    clean[2] / sum,
                    clean[3] / sum,
                ],
            }
        }
    }

    /// Creates a posterior peaked on a known base `b` with probability `confidence` in `[0.25, 1.0]`.
    pub fn from_base(b: u8, confidence: f64) -> Self {
        let conf = confidence.clamp(0.25, 1.0);
        let rem = (1.0 - conf) / 3.0;
        match base_to_index(b) {
            Some(idx) => {
                let mut probs = [rem; 4];
                probs[idx] = conf;
                Self { probs }
            }
            None => Self::uniform(),
        }
    }

    /// Returns the index `0..4` with the maximum posterior probability.
    pub fn max_index(&self) -> usize {
        let mut best_idx = 0;
        let mut best_p = self.probs[0];
        for i in 1..4 {
            if self.probs[i] > best_p {
                best_p = self.probs[i];
                best_idx = i;
            }
        }
        best_idx
    }

    /// Returns the Maximum A Posteriori (MAP) canonical nucleotide (`A`, `C`, `G`, or `T`).
    pub fn most_likely_base(&self) -> u8 {
        index_to_base(self.max_index())
    }

    /// Returns the posterior probability of the most likely base.
    pub fn max_prob(&self) -> f64 {
        self.probs[self.max_index()]
    }

    /// Returns the probability for a specific base (`A`, `C`, `G`, `T`).
    pub fn prob_of(&self, base: u8) -> f64 {
        base_to_index(base).map(|i| self.probs[i]).unwrap_or(0.0)
    }

    /// Converts the posterior to an IUPAC nucleotide symbol:
    /// - If a single base has probability >= `min_confidence`, returns `A`, `C`, `G`, or `T`.
    /// - Otherwise combines the top bases whose cumulative probability reaches `min_confidence`
    ///   into standard IUPAC ambiguity codes (`R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`, `N`).
    pub fn to_iupac(&self, min_confidence: f64) -> u8 {
        let thresh = min_confidence.clamp(0.25, 0.9999);
        if self.max_prob() >= thresh {
            return self.most_likely_base();
        }
        // Sort indices descending by probability
        let mut indices = [0usize, 1, 2, 3];
        indices.sort_by(|&a, &b| {
            self.probs[b]
                .partial_cmp(&self.probs[a])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let mut mask = 0u8;
        let mut cum = 0.0;
        for &idx in &indices {
            if self.probs[idx] <= 0.05 {
                break;
            }
            mask |= 1 << idx; // bit 0=A, 1=C, 2=G, 3=T
            cum += self.probs[idx];
            if cum >= thresh {
                break;
            }
        }

        match mask {
            0b0001 => b'A',
            0b0010 => b'C',
            0b0100 => b'G',
            0b1000 => b'T',
            0b0101 => b'R', // A or G
            0b1010 => b'Y', // C or T
            0b0110 => b'S', // C or G
            0b1001 => b'W', // A or T
            0b1100 => b'K', // G or T
            0b0011 => b'M', // A or C
            0b1110 => b'B', // C, G, T
            0b1101 => b'D', // A, G, T
            0b1011 => b'H', // A, C, T
            0b0111 => b'V', // A, C, G
            _ => b'N',
        }
    }

    /// Computes the Shannon entropy $H = -\sum_b p_b \log_2(p_b)$ in bits (`0.0..=2.0`).
    pub fn entropy_bits(&self) -> f64 {
        let mut h = 0.0;
        for &p in &self.probs {
            if p > 1e-15 {
                h -= p * p.log2();
            }
        }
        h.clamp(0.0, 2.0)
    }
}

/// Primary evidence source that anchored a given nucleotide site.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum SiteSource {
    /// Inferred purely via Phylogenetic Ancestral Sequence Reconstruction (ASR).
    AncestralReconstruction,
    /// Anchored by ancient peptide fragments via paleoproteomic back-translation.
    Paleoproteomics,
    /// Anchored by ancient DNA reads with post-mortem deamination modeling.
    AncientDna,
    /// Repaired from an inactivated pseudogene using functional outgroup bracketing.
    PseudogeneRepair,
}

/// Breakdown of evidence sources across a reconstructed locus.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct LocusProvenance {
    pub site_sources: Vec<SiteSource>,
    pub asr_fraction: f64,
    pub proteomic_fraction: f64,
    pub adna_fraction: f64,
    pub pseudogene_repair_fraction: f64,
}

impl LocusProvenance {
    /// Initializes provenance where all `len` sites originate from `source`.
    pub fn new_uniform(len: usize, source: SiteSource) -> Self {
        let mut prov = Self {
            site_sources: vec![source; len],
            asr_fraction: 0.0,
            proteomic_fraction: 0.0,
            adna_fraction: 0.0,
            pseudogene_repair_fraction: 0.0,
        };
        prov.recompute_fractions();
        prov
    }

    /// Recomputes summary fractions from `site_sources`.
    pub fn recompute_fractions(&mut self) {
        let n = self.site_sources.len();
        if n == 0 {
            self.asr_fraction = 0.0;
            self.proteomic_fraction = 0.0;
            self.adna_fraction = 0.0;
            self.pseudogene_repair_fraction = 0.0;
            return;
        }
        let mut asr = 0usize;
        let mut prot = 0usize;
        let mut adna = 0usize;
        let mut pseudo = 0usize;
        for &s in &self.site_sources {
            match s {
                SiteSource::AncestralReconstruction => asr += 1,
                SiteSource::Paleoproteomics => prot += 1,
                SiteSource::AncientDna => adna += 1,
                SiteSource::PseudogeneRepair => pseudo += 1,
            }
        }
        let nf = n as f64;
        self.asr_fraction = asr as f64 / nf;
        self.proteomic_fraction = prot as f64 / nf;
        self.adna_fraction = adna as f64 / nf;
        self.pseudogene_repair_fraction = pseudo as f64 / nf;
    }
}

/// A reconstructed genomic locus with per-site posterior distributions and uncertainty tracking.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ReconstructedLocus {
    /// Standard gene symbol (e.g., `"COL1A1"`, `"ENAM"`, `"UCP1"`).
    pub gene_symbol: String,
    /// Target extinct taxon name (e.g., `"Tyrannosaurus rex"`, `"Mammuthus primigenius"`).
    pub target_taxon: String,
    /// Per-base posterior distributions `[P(A), P(C), P(G), P(T)]`.
    pub sites: Vec<BasePosterior>,
    /// Mean MAP probability across all sites (`0.25..=1.0`).
    pub mean_posterior_confidence: f64,
    /// Mean Shannon entropy in bits across all sites (`0.0..=2.0`).
    pub mean_entropy_bits: f64,
    /// Contiguous coordinate intervals `[start, end)` where site confidence is `< 0.70`.
    pub high_uncertainty_regions: Vec<(usize, usize)>,
    /// Site-level provenance tracking.
    pub provenance: LocusProvenance,
}

impl ReconstructedLocus {
    /// Constructs a `ReconstructedLocus` and computes all summary uncertainty metrics.
    pub fn new(
        gene_symbol: impl Into<String>,
        target_taxon: impl Into<String>,
        sites: Vec<BasePosterior>,
        provenance: LocusProvenance,
    ) -> Self {
        let mut locus = Self {
            gene_symbol: gene_symbol.into(),
            target_taxon: target_taxon.into(),
            sites,
            mean_posterior_confidence: 0.0,
            mean_entropy_bits: 0.0,
            high_uncertainty_regions: Vec::new(),
            provenance,
        };
        locus.refresh_metrics();
        locus
    }

    /// Recomputes `mean_posterior_confidence`, `mean_entropy_bits`, `high_uncertainty_regions`,
    /// and `provenance` fractions after any Bayesian update.
    pub fn refresh_metrics(&mut self) {
        self.provenance.recompute_fractions();
        if self.sites.is_empty() {
            self.mean_posterior_confidence = 0.0;
            self.mean_entropy_bits = 0.0;
            self.high_uncertainty_regions.clear();
            return;
        }
        let n = self.sites.len() as f64;
        let sum_conf: f64 = self.sites.iter().map(|s| s.max_prob()).sum();
        let sum_ent: f64 = self.sites.iter().map(|s| s.entropy_bits()).sum();
        self.mean_posterior_confidence = sum_conf / n;
        self.mean_entropy_bits = sum_ent / n;

        self.high_uncertainty_regions.clear();
        let mut in_region = false;
        let mut start = 0usize;
        for (i, site) in self.sites.iter().enumerate() {
            let uncertain = site.max_prob() < 0.70;
            if uncertain && !in_region {
                in_region = true;
                start = i;
            } else if !uncertain && in_region {
                in_region = false;
                self.high_uncertainty_regions.push((start, i));
            }
        }
        if in_region {
            self.high_uncertainty_regions.push((start, self.sites.len()));
        }
    }

    /// Returns the Maximum A Posteriori (MAP) consensus [`Dna`] sequence.
    pub fn to_consensus_dna(&self) -> Dna {
        let bytes: Vec<u8> = self.sites.iter().map(|s| s.most_likely_base()).collect();
        Dna::new(&bytes).unwrap_or(Dna { sequence: bytes })
    }

    /// Returns an IUPAC-encoded sequence where ambiguous positions below `min_confidence`
    /// use standard degeneracy symbols (`R`, `Y`, `S`, `W`, `K`, `M`, `N`, etc.).
    pub fn to_iupac_bytes(&self, min_confidence: f64) -> Vec<u8> {
        self.sites
            .iter()
            .map(|s| s.to_iupac(min_confidence))
            .collect()
    }

    /// Translates the MAP consensus coding sequence into a [`Protein`] (ignoring terminal stop `*`).
    pub fn translate_consensus(&self) -> Result<Protein, PaleoError> {
        if self.sites.len() % 3 != 0 {
            return Err(PaleoError::FrameMismatch(self.sites.len()));
        }
        let dna = self.to_consensus_dna();
        let mut aa_bytes = Vec::with_capacity(self.sites.len() / 3);
        for chunk in dna.sequence.chunks_exact(3) {
            let aa = translate_dna_codon([chunk[0], chunk[1], chunk[2]]);
            if aa != b'*' {
                aa_bytes.push(aa);
            }
        }
        Protein::new(&aa_bytes).map_err(|e| PaleoError::InvalidParameter(e.to_string()))
    }
}
