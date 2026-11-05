//! Layer 1: Continuous-Time Markov Chain (CTMC) Ancestral Sequence Reconstruction (ASR).
//!
//! Implements Felsenstein's post-order pruning algorithm on rooted phylogenetic trees (`PhyloNode`)
//! under HKY85, GTR, Kimura 2-Parameter (K80), and Jukes-Cantor (JC69) nucleotide substitution
//! models, incorporating discrete $\Gamma_4$ rate heterogeneity across sites and dynamic
//! node log-space scaling to eliminate numerical underflow on deep or dinosaur-scale extinct branches.

use crate::paleo::types::{
    base_to_index, BasePosterior, LocusProvenance, PaleoError, ReconstructedLocus, SiteSource,
};
use crate::phylo::PhyloNode;
use crate::seq::Dna;
use std::collections::HashMap;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Continuous-time Markov chain nucleotide substitution models.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum SubstitutionModel {
    /// Jukes-Cantor (1969): equal base frequencies ($\pi_i = 0.25$) and equal substitution rates ($\mu$ per Ma).
    JukesCantor { mu: f64 },

    /// Kimura 2-Parameter (1980): transition/transversion rate ratio $\kappa$ and overall rate $\mu$ per Ma.
    Kimura2Param { kappa: f64, mu: f64 },

    /// Hasegawa-Kishino-Yano (1985): unequal base frequencies $\pi$, transition/transversion ratio $\kappa$,
    /// overall rate $\mu$ per Ma, and optional discrete $\Gamma_4$ shape parameter $\alpha$.
    Hky85 {
        /// Stationary nucleotide frequencies `[pi_A, pi_C, pi_G, pi_T]`.
        pi: [f64; 4],
        /// Transition/transversion rate ratio $\kappa$ (typically $2.0$ to $5.0$).
        kappa: f64,
        /// Overall substitution rate $\mu$ per million years (Ma).
        mu: f64,
        /// Shape parameter $\alpha$ for 4-category discrete Gamma rate heterogeneity ($+\Gamma_4$).
        gamma_shape: Option<f64>,
    },

    /// General Time Reversible (GTR): unequal base frequencies $\pi$, 6 symmetric exchangeabilities
    /// `[r_AC, r_AG, r_AT, r_CG, r_CT, r_GT]`, rate $\mu$, and optional $\Gamma_4$ shape $\alpha$.
    Gtr {
        /// Stationary nucleotide frequencies `[pi_A, pi_C, pi_G, pi_T]`.
        pi: [f64; 4],
        /// Symmetric exchangeability rate multipliers `[r_AC, r_AG, r_AT, r_CG, r_CT, r_GT]`.
        rates: [f64; 6],
        /// Overall substitution rate $\mu$ per million years (Ma).
        mu: f64,
        /// Shape parameter $\alpha$ for 4-category discrete Gamma rate heterogeneity ($+\Gamma_4$).
        gamma_shape: Option<f64>,
    },
}

impl SubstitutionModel {
    /// Default vertebrate nuclear coding substitution model:
    /// Uses HKY85 with empirical vertebrate coding base frequencies ($\text{GC} \approx 46\%$),
    /// transition/transversion ratio $\kappa = 3.0$, rate $\mu = 0.0022$ per Ma,
    /// and 4-category discrete Gamma shape $\alpha = 0.60$.
    pub fn default_vertebrate_coding() -> Self {
        Self::Hky85 {
            pi: [0.27, 0.23, 0.23, 0.27],
            kappa: 3.0,
            mu: 0.0022,
            gamma_shape: Some(0.60),
        }
    }

    /// Stationary nucleotide distribution $\pi = [\pi_A, \pi_C, \pi_G, \pi_T]$.
    pub fn stationary_frequencies(&self) -> [f64; 4] {
        match self {
            Self::JukesCantor { .. } | Self::Kimura2Param { .. } => [0.25; 4],
            Self::Hky85 { pi, .. } | Self::Gtr { pi, .. } => {
                let sum: f64 = pi.iter().copied().filter(|&x| x > 0.0).sum();
                if sum > 1e-9 {
                    [pi[0] / sum, pi[1] / sum, pi[2] / sum, pi[3] / sum]
                } else {
                    [0.25; 4]
                }
            }
        }
    }

    /// Returns true if `(from_idx, to_idx)` is a purine ($A \leftrightarrow G$)
    /// or pyrimidine ($C \leftrightarrow T$) transition.
    #[inline]
    pub fn is_transition(from_idx: usize, to_idx: usize) -> bool {
        (from_idx == 0 && to_idx == 2)
            || (from_idx == 2 && to_idx == 0)
            || (from_idx == 1 && to_idx == 3)
            || (from_idx == 3 && to_idx == 1)
    }

    /// Computes the 4 representative rates for discrete $\Gamma_4$ rate heterogeneity under shape $\alpha$.
    ///
    /// Uses the Wilson-Hilferty transformation at quartile midpoints $p \in \{0.125, 0.375, 0.625, 0.875\}$,
    /// normalized so that $\frac{1}{4}\sum_{k=1}^4 r_k = 1.0000000$.
    pub fn gamma_rates_4(alpha: f64) -> [f64; 4] {
        let a = alpha.max(0.05);
        // Standard normal quantiles for 4 equal-probability bins
        let z = [-1.150349, -0.318639, 0.318639, 1.150349];
        let factor = 1.0 / (9.0 * a);
        let sqrt_factor = factor.sqrt();

        let mut raw = [0.0; 4];
        for i in 0..4 {
            let term = (1.0 - factor + z[i] * sqrt_factor).max(0.01);
            raw[i] = term.powi(3);
        }

        let mean = (raw[0] + raw[1] + raw[2] + raw[3]) * 0.25;
        if mean > 1e-9 {
            [
                raw[0] / mean,
                raw[1] / mean,
                raw[2] / mean,
                raw[3] / mean,
            ]
        } else {
            [1.0; 4]
        }
    }

    /// Computes the 4 discrete gamma rates and weights for this model:
    /// Returns `[(rate_multiplier, category_weight)]`.
    pub fn discrete_rate_categories(&self) -> Vec<(f64, f64)> {
        let shape_opt = match self {
            Self::Hky85 { gamma_shape, .. } | Self::Gtr { gamma_shape, .. } => *gamma_shape,
            _ => None,
        };

        match shape_opt {
            Some(alpha) => {
                let rates = Self::gamma_rates_4(alpha);
                vec![
                    (rates[0], 0.25),
                    (rates[1], 0.25),
                    (rates[2], 0.25),
                    (rates[3], 0.25),
                ]
            }
            None => vec![(1.0, 1.0)],
        }
    }

    /// Computes exact transition probability matrix $P(t) = \exp(Q t)$ for branch time $t \ge 0$.
    pub fn transition_matrix(&self, branch_time: f64) -> [[f64; 4]; 4] {
        let t = branch_time.max(0.0);
        match self {
            Self::JukesCantor { mu } => {
                let rate = mu.max(1e-12);
                let exp_term = (-4.0 * rate * t / 3.0).exp();
                let mut mat = [[0.0; 4]; 4];
                for i in 0..4 {
                    for j in 0..4 {
                        mat[i][j] = if i == j {
                            0.25 + 0.75 * exp_term
                        } else {
                            0.25 - 0.25 * exp_term
                        };
                    }
                }
                mat
            }
            Self::Kimura2Param { kappa, mu } => {
                let k = kappa.max(0.1);
                let d = mu.max(1e-12) * t;
                let exp_tv = (-4.0 * d / (k + 2.0)).exp();
                let exp_ts = (-2.0 * d * (k + 1.0) / (k + 2.0)).exp();

                let p_tv = 0.25 - 0.25 * exp_tv;
                let p_ts = 0.25 + 0.25 * exp_tv - 0.50 * exp_ts;

                let mut mat = [[0.0; 4]; 4];
                for i in 0..4 {
                    for j in 0..4 {
                        mat[i][j] = if i == j {
                            (1.0 - p_ts - 2.0 * p_tv).clamp(0.0, 1.0)
                        } else if Self::is_transition(i, j) {
                            p_ts.clamp(0.0, 1.0)
                        } else {
                            p_tv.clamp(0.0, 1.0)
                        };
                    }
                }
                mat
            }
            Self::Hky85 { kappa, mu, .. } => {
                let freqs = self.stationary_frequencies();
                let pi_a = freqs[0];
                let pi_c = freqs[1];
                let pi_g = freqs[2];
                let pi_t = freqs[3];

                let pi_r = (pi_a + pi_g).max(1e-9);
                let pi_y = (pi_c + pi_t).max(1e-9);
                let k = kappa.max(0.05);

                // Normalization constant beta_0 so that average substitution rate equals mu:
                // rate = 2 * beta_0 * (pi_R * pi_Y + kappa * (pi_A * pi_G + pi_C * pi_T)) = mu
                let denom = 2.0 * (pi_r * pi_y + k * (pi_a * pi_g + pi_c * pi_t));
                let beta_0 = if denom > 1e-12 { mu / denom } else { *mu };

                let lambda_1 = beta_0 * t;
                let lambda_r = beta_0 * t * (1.0 + (k - 1.0) * pi_r);
                let lambda_y = beta_0 * t * (1.0 + (k - 1.0) * pi_y);

                let exp_1 = (-lambda_1).exp();
                let exp_r = (-lambda_r).exp();
                let exp_y = (-lambda_y).exp();

                let is_purine = [true, false, true, false];
                let mut mat = [[0.0; 4]; 4];

                for i in 0..4 {
                    let from_purine = is_purine[i];
                    for j in 0..4 {
                        let to_purine = is_purine[j];
                        let pi_j = freqs[j];
                        let delta_ij = if i == j { 1.0 } else { 0.0 };

                        if from_purine != to_purine {
                            // Transversion between purine and pyrimidine
                            mat[i][j] = (pi_j * (1.0 - exp_1)).max(0.0);
                        } else if from_purine {
                            // Purine -> Purine
                            let p = pi_j
                                + pi_j * (1.0 / pi_r - 1.0) * exp_1
                                + ((delta_ij * pi_r - pi_j) / pi_r) * exp_r;
                            mat[i][j] = p.max(0.0);
                        } else {
                            // Pyrimidine -> Pyrimidine
                            let p = pi_j
                                + pi_j * (1.0 / pi_y - 1.0) * exp_1
                                + ((delta_ij * pi_y - pi_j) / pi_y) * exp_y;
                            mat[i][j] = p.max(0.0);
                        }
                    }

                    // Enforce exact row normalization to eliminate floating point drift
                    let row_sum: f64 = mat[i].iter().sum();
                    if row_sum > 1e-12 {
                        for j in 0..4 {
                            mat[i][j] /= row_sum;
                        }
                    }
                }
                mat
            }
            Self::Gtr { rates, mu, .. } => {
                let freqs = self.stationary_frequencies();
                // Symmetrical rate matrix entries r_ij
                let r_ac = rates[0].max(0.0);
                let r_ag = rates[1].max(0.0);
                let r_at = rates[2].max(0.0);
                let r_cg = rates[3].max(0.0);
                let r_ct = rates[4].max(0.0);
                let r_gt = rates[5].max(0.0);

                let r_mat = [
                    [0.0, r_ac, r_ag, r_at],
                    [r_ac, 0.0, r_cg, r_ct],
                    [r_ag, r_cg, 0.0, r_gt],
                    [r_at, r_ct, r_gt, 0.0],
                ];

                let mut q = [[0.0; 4]; 4];
                for i in 0..4 {
                    for j in 0..4 {
                        if i != j {
                            q[i][j] = r_mat[i][j] * freqs[j];
                        }
                    }
                    let off_sum: f64 = q[i].iter().sum();
                    q[i][i] = -off_sum;
                }

                // Average substitution rate = -sum_i pi_i Q_ii
                let mut avg_rate = 0.0;
                for i in 0..4 {
                    avg_rate -= freqs[i] * q[i][i];
                }

                let scale = if avg_rate > 1e-12 {
                    *mu / avg_rate
                } else {
                    1.0
                };

                let mut q_scaled = [[0.0; 4]; 4];
                for i in 0..4 {
                    for j in 0..4 {
                        q_scaled[i][j] = q[i][j] * scale * t;
                    }
                }

                // Compute matrix exponential exp(Q * t) via scaling and squaring
                Self::matrix_exponential_4x4(&q_scaled)
            }
        }
    }

    /// Computes matrix exponential $\exp(A)$ for a $4 \times 4$ rate matrix using scaling and squaring with 8th-order Taylor polynomial.
    fn matrix_exponential_4x4(a: &[[f64; 4]; 4]) -> [[f64; 4]; 4] {
        let mut max_entry = 0.0f64;
        for i in 0..4 {
            for j in 0..4 {
                max_entry = max_entry.max(a[i][j].abs());
            }
        }

        // Determine squaring count s such that ||A / 2^s|| <= 0.5
        let mut s = 0usize;
        let mut cur = max_entry;
        while cur > 0.5 && s < 16 {
            cur *= 0.5;
            s += 1;
        }
        let divisor = (2.0f64).powi(s as i32);

        let mut a_scaled = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                a_scaled[i][j] = a[i][j] / divisor;
            }
        }

        // Taylor series: E = I + A + A^2/2 + A^3/6 + ... + A^8/40320
        let mut term = [[0.0; 4]; 4];
        for i in 0..4 {
            term[i][i] = 1.0; // I
        }
        let mut sum = term;

        for k in 1..=8 {
            let mut next_term = [[0.0; 4]; 4];
            for i in 0..4 {
                for j in 0..4 {
                    for m in 0..4 {
                        next_term[i][j] += term[i][m] * a_scaled[m][j];
                    }
                    next_term[i][j] /= k as f64;
                }
            }
            term = next_term;
            for i in 0..4 {
                for j in 0..4 {
                    sum[i][j] += term[i][j];
                }
            }
        }

        // Repeated squaring: sum = sum^(2^s)
        for _ in 0..s {
            let mut squared = [[0.0; 4]; 4];
            for i in 0..4 {
                for j in 0..4 {
                    for m in 0..4 {
                        squared[i][j] += sum[i][m] * sum[m][j];
                    }
                }
            }
            sum = squared;
        }

        // Row normalization to enforce stochastic matrix property
        for i in 0..4 {
            let row_sum: f64 = sum[i].iter().map(|&x| x.max(0.0)).sum();
            if row_sum > 1e-12 {
                for j in 0..4 {
                    sum[i][j] = (sum[i][j].max(0.0)) / row_sum;
                }
            } else {
                sum[i] = [0.25; 4];
            }
        }

        sum
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
    pub fn new(taxon_name: impl Into<String>, sequence: Dna) -> Self {
        Self {
            taxon_name: taxon_name.into(),
            sequence,
        }
    }
}

/// Felsenstein Ancestral Sequence Reconstruction engine with node log-space scaling.
pub struct FelsensteinAsr;

impl FelsensteinAsr {
    /// Reconstructs ancestral nucleotide posterior distribution at the root of `tree`
    /// and propagates it forward along `extinct_branch_ma` (in millions of years) to the target extinct taxon.
    ///
    /// Applies discrete $\Gamma_4$ rate categories and dynamic log-space node scaling to guarantee
    /// zero numerical underflow across arbitrary deep phylogenetic trees and dinosaur-scale extinct branches.
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

        let rate_cats = model.discrete_rate_categories();
        let stationary_pi = model.stationary_frequencies();

        // Precompute transition matrices for the extinct terminal branch across all rate categories
        let ext_trans_matrices: Vec<[[f64; 4]; 4]> = rate_cats
            .iter()
            .map(|&(rate_multiplier, _)| {
                model.transition_matrix(extinct_branch_ma.max(0.0) * rate_multiplier)
            })
            .collect();

        let mut sites = Vec::with_capacity(seq_len);

        for pos in 0..seq_len {
            // Marginalize over discrete gamma rate categories with log-scaling protection
            let mut category_root_posteriors = Vec::with_capacity(rate_cats.len());
            let mut category_log_scales = Vec::with_capacity(rate_cats.len());

            for (_cat_idx, &(rate_multiplier, _)) in rate_cats.iter().enumerate() {
                let (scaled_l, log_scale) =
                    Self::prune_site_scaled(tree, &seq_map, pos, model, rate_multiplier)?;
                // Apply stationary prior: posterior_b proportional to pi_b * L(b)
                let mut unnorm = [0.0; 4];
                for b in 0..4 {
                    unnorm[b] = stationary_pi[b] * scaled_l[b];
                }
                let norm = unnorm[0] + unnorm[1] + unnorm[2] + unnorm[3];
                let post = if norm > 1e-15 {
                    [
                        unnorm[0] / norm,
                        unnorm[1] / norm,
                        unnorm[2] / norm,
                        unnorm[3] / norm,
                    ]
                } else {
                    stationary_pi
                };
                category_root_posteriors.push(post);
                category_log_scales.push(log_scale + norm.max(1e-300).ln());
            }

            // Combine rate categories using log-sum-exp normalization
            let max_log = category_log_scales
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);

            let mut final_ancestor_post = [0.0; 4];
            let mut total_cat_weight = 0.0;

            for (cat_idx, &(_, cat_weight)) in rate_cats.iter().enumerate() {
                let relative_weight = cat_weight * (category_log_scales[cat_idx] - max_log).exp();
                total_cat_weight += relative_weight;
                for b in 0..4 {
                    final_ancestor_post[b] += relative_weight * category_root_posteriors[cat_idx][b];
                }
            }

            if total_cat_weight > 1e-15 {
                for b in 0..4 {
                    final_ancestor_post[b] /= total_cat_weight;
                }
            } else {
                final_ancestor_post = stationary_pi;
            }

            let ancestor_post = BasePosterior::from_weights(final_ancestor_post);

            if extinct_branch_ma <= 1e-9 {
                sites.push(ancestor_post);
            } else {
                // Propagate along extinct terminal branch: p_ext(j) = sum_cat w_cat * (sum_i p_anc(i) P_{ij}(t * r_cat))
                let mut ext_weights = [0.0; 4];
                for (cat_idx, &(_, cat_weight)) in rate_cats.iter().enumerate() {
                    let p_mat = &ext_trans_matrices[cat_idx];
                    for i in 0..4 {
                        for j in 0..4 {
                            ext_weights[j] += cat_weight * ancestor_post.probs[i] * p_mat[i][j];
                        }
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

    /// Recursive post-order Felsenstein pruning for a single site `pos` with dynamic node log-scaling.
    ///
    /// Returns `(scaled_likelihood_at_node, accumulated_log_scale_factor)`.
    /// At each internal node, dividing by $\max_b L(b)$ guarantees that the vector is bounded in $[0, 1]$,
    /// completely preventing numerical underflow on deep trees.
    fn prune_site_scaled(
        node: &PhyloNode,
        seq_map: &HashMap<&str, &[u8]>,
        pos: usize,
        model: &SubstitutionModel,
        rate_multiplier: f64,
    ) -> Result<([f64; 4], f64), PaleoError> {
        if node.is_leaf() {
            let name = node.name.as_deref().unwrap_or("");
            let seq = seq_map
                .get(name)
                .ok_or_else(|| PaleoError::TaxonNotFound(name.to_string()))?;
            let obs_base = seq[pos];
            let cond = match base_to_index(obs_base) {
                Some(idx) => {
                    let mut c = [0.0; 4];
                    c[idx] = 1.0;
                    c
                }
                None => [1.0; 4],
            };
            Ok((cond, 0.0))
        } else {
            let mut likelihood = [1.0; 4];
            let mut total_log_scale = 0.0;

            for child in &node.children {
                let (child_l, child_log_scale) =
                    Self::prune_site_scaled(child, seq_map, pos, model, rate_multiplier)?;
                total_log_scale += child_log_scale;

                let p_mat = model.transition_matrix(child.branch_length * rate_multiplier);
                for parent_b in 0..4 {
                    let mut sum_child = 0.0;
                    for child_b in 0..4 {
                        sum_child += p_mat[parent_b][child_b] * child_l[child_b];
                    }
                    likelihood[parent_b] *= sum_child;
                }
            }

            // Dynamic scaling step: divide vector by max element to prevent underflow
            let max_val = likelihood
                .iter()
                .copied()
                .fold(0.0f64, f64::max);

            if max_val > 1e-300 {
                let inv = 1.0 / max_val;
                for b in 0..4 {
                    likelihood[b] *= inv;
                }
                total_log_scale += max_val.ln();
            }

            Ok((likelihood, total_log_scale))
        }
    }
}
