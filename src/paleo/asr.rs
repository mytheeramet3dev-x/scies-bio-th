//! Layer 1: Continuous-Time Markov Chain (CTMC) Ancestral Sequence Reconstruction (ASR).
//!
//! Implements Felsenstein's post-order pruning algorithm on rooted phylogenetic trees
//! (`PhyloNode`) under Jukes-Cantor (JC69) and Kimura 2-Parameter (K80) nucleotide
//! substitution models, along with forward Markov propagation along extinct branches.

use crate::paleo::types::{
    base_to_index, BasePosterior, LocusProvenance, PaleoError, ReconstructedLocus, SiteSource,
};
use crate::phylo::PhyloNode;
use crate::seq::Dna;
use std::collections::HashMap;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Continuous-time Markov chain nucleotide substitution model.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum SubstitutionModel {
    /// Jukes-Cantor (1969): equal base frequencies and equal substitution rates ($\mu$ per Ma).
    JukesCantor { mu: f64 },
    /// Kimura 2-Parameter (1980): transition/transversion rate ratio $\kappa$ and overall rate $\mu$ per Ma.
    Kimura2Param { kappa: f64, mu: f64 },
}

impl SubstitutionModel {
    /// Default vertebrate nuclear coding substitution model (K80 with $\kappa = 2.5, \mu = 0.0022$ per Ma).
    pub fn default_vertebrate_coding() -> Self {
        Self::Kimura2Param {
            kappa: 2.5,
            mu: 0.0022,
        }
    }

    /// Returns true if `(from_idx, to_idx)` is a purine ($A \leftrightarrow G$, `0 <-> 2`)
    /// or pyrimidine ($C \leftrightarrow T$, `1 <-> 3`) transition.
    #[inline]
    pub fn is_transition(from_idx: usize, to_idx: usize) -> bool {
        (from_idx == 0 && to_idx == 2)
            || (from_idx == 2 && to_idx == 0)
            || (from_idx == 1 && to_idx == 3)
            || (from_idx == 3 && to_idx == 1)
    }

    /// Computes exact transition probability $P(X(t) = \text{to} \mid X(0) = \text{from})$ over branch time $t \ge 0$.
    pub fn transition_prob(&self, from_idx: usize, to_idx: usize, branch_time: f64) -> f64 {
        let t = branch_time.max(0.0);
        match *self {
            Self::JukesCantor { mu } => {
                let rate = mu.max(1e-12);
                let exp_term = (-4.0 * rate * t / 3.0).exp();
                if from_idx == to_idx {
                    0.25 + 0.75 * exp_term
                } else {
                    0.25 - 0.25 * exp_term
                }
            }
            Self::Kimura2Param { kappa, mu } => {
                let k = kappa.max(0.1);
                let d = mu.max(1e-12) * t;
                let exp_tv = (-4.0 * d / (k + 2.0)).exp();
                let exp_ts = (-2.0 * d * (k + 1.0) / (k + 2.0)).exp();

                let p_tv = 0.25 - 0.25 * exp_tv;
                let p_ts = 0.25 + 0.25 * exp_tv - 0.50 * exp_ts;
                if from_idx == to_idx {
                    (1.0 - p_ts - 2.0 * p_tv).clamp(0.0, 1.0)
                } else if Self::is_transition(from_idx, to_idx) {
                    p_ts.clamp(0.0, 1.0)
                } else {
                    p_tv.clamp(0.0, 1.0)
                }
            }
        }
    }

    /// Returns the full $4 \times 4$ Markov transition probability matrix $P(t)$ for branch time $t$.
    pub fn transition_matrix(&self, branch_time: f64) -> [[f64; 4]; 4] {
        let mut mat = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                mat[i][j] = self.transition_prob(i, j, branch_time);
            }
        }
        mat
    }
}

/// An aligned extant taxon sequence used as a leaf in Phylogenetic Bracketing.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ExtantSequence {
    pub taxon_name: String,
    pub sequence: Dna,
}

impl ExtantSequence {
    /// Creates a new `ExtantSequence`.
    pub fn new(taxon_name: impl Into<String>, sequence: Dna) -> Self {
        Self {
            taxon_name: taxon_name.into(),
            sequence,
        }
    }
}

/// Felsenstein Ancestral Sequence Reconstruction engine.
pub struct FelsensteinAsr;

impl FelsensteinAsr {
    /// Reconstructs the ancestral nucleotide posterior distribution at the root of `tree`
    /// and propagates it forward along `extinct_branch_ma` (in millions of years) to the target extinct taxon.
    ///
    /// - If `extinct_branch_ma == 0.0`, returns the exact marginal posterior at the phylogenetic bracket ancestor.
    /// - If `extinct_branch_ma > 0.0`, accounts for uninferred lineage-specific evolution along the extinct branch.
    pub fn reconstruct_locus(
        tree: &PhyloNode,
        extant_seqs: &[ExtantSequence],
        model: &SubstitutionModel,
        target_taxon: &str,
        gene_symbol: &str,
        extinct_branch_ma: f64,
    ) -> Result<ReconstructedLocus, PaleoError> {
        if extant_seqs.is_empty() {
            return Err(PaleoError::EmptyTaxonSet);
        }
        let seq_len = extant_seqs[0].sequence.sequence.len();
        let mut seq_map: HashMap<&str, &[u8]> = HashMap::with_capacity(extant_seqs.len());
        for ext in extant_seqs {
            if ext.sequence.sequence.len() != seq_len {
                return Err(PaleoError::LengthMismatch {
                    expected: seq_len,
                    found: ext.sequence.sequence.len(),
                });
            }
            seq_map.insert(ext.taxon_name.as_str(), &ext.sequence.sequence);
        }

        let extinct_trans = model.transition_matrix(extinct_branch_ma.max(0.0));
        let mut sites = Vec::with_capacity(seq_len);

        for pos in 0..seq_len {
            let root_likelihood = Self::prune_site(tree, &seq_map, pos, model)?;
            // Stationary prior pi = [0.25; 4] cancels out during normalization
            let ancestor_post = BasePosterior::from_weights(root_likelihood);

            if extinct_branch_ma <= 1e-9 {
                sites.push(ancestor_post);
            } else {
                // Propagate distribution along the extinct terminal branch: p_ext(j) = sum_i p_anc(i) P_{i->j}(dt)
                let mut ext_weights = [0.0; 4];
                for i in 0..4 {
                    for j in 0..4 {
                        ext_weights[j] += ancestor_post.probs[i] * extinct_trans[i][j];
                    }
                }
                sites.push(BasePosterior::from_weights(ext_weights));
            }
        }

        let provenance = LocusProvenance::new_uniform(seq_len, SiteSource::AncestralReconstruction);
        Ok(ReconstructedLocus::new(
            gene_symbol,
            target_taxon,
            sites,
            provenance,
        ))
    }

    /// Recursive post-order Felsenstein pruning for a single site `pos`.
    fn prune_site(
        node: &PhyloNode,
        seq_map: &HashMap<&str, &[u8]>,
        pos: usize,
        model: &SubstitutionModel,
    ) -> Result<[f64; 4], PaleoError> {
        if node.is_leaf() {
            let name = node.name.as_deref().unwrap_or("");
            let seq = seq_map
                .get(name)
                .ok_or_else(|| PaleoError::TaxonNotFound(name.to_string()))?;
            let obs_base = seq[pos];
            Ok(match base_to_index(obs_base) {
                Some(idx) => {
                    let mut cond = [0.0; 4];
                    cond[idx] = 1.0;
                    cond
                }
                // Ambiguous or gap character 'N' treats all 4 bases as equally plausible
                None => [1.0; 4],
            })
        } else {
            let mut likelihood = [1.0; 4];
            for child in &node.children {
                let child_l = Self::prune_site(child, seq_map, pos, model)?;
                let p_mat = model.transition_matrix(child.branch_length);
                for parent_b in 0..4 {
                    let mut sum_child = 0.0;
                    for child_b in 0..4 {
                        sum_child += p_mat[parent_b][child_b] * child_l[child_b];
                    }
                    likelihood[parent_b] *= sum_child;
                }
            }
            Ok(likelihood)
        }
    }
}
