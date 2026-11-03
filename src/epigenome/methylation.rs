//! DNA methylation primitives, CpG island detection (Gardiner-Garden & Frommer), and Beta/M-value conversions.

use crate::seq::Dna;
use crate::BioError;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A detected CpG island genomic region.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CpgIsland {
    /// 0-based start coordinate.
    pub start: usize,
    /// 0-based end coordinate (exclusive).
    pub end: usize,
    /// Length of the island in base pairs.
    pub length: usize,
    /// Fraction of GC content (0.0 to 1.0).
    pub gc_content: f64,
    /// Observed-to-expected CpG ratio.
    pub obs_exp_ratio: f64,
    /// Number of CpG dinucleotides within the island.
    pub cpg_count: usize,
}

/// Detects CpG islands across a DNA sequence according to Gardiner-Garden & Frommer (1987) criteria:
/// 1. Window length $W \ge \text{min\_len}$ (default 200 bp).
/// 2. GC content $\ge \text{min\_gc}$ (default 0.50).
/// 3. Observed/Expected CpG ratio $\ge \text{min\_obs\_exp}$ (default 0.60):
///    $$\text{Obs/Exp} = \frac{\text{Count}(\text{CpG}) \times L}{\text{Count}(C) \times \text{Count}(G)}$$
pub fn find_cpg_islands(
    dna: &Dna,
    window_size: usize,
    step_size: usize,
    min_gc: f64,
    min_obs_exp: f64,
) -> Vec<CpgIsland> {
    let seq = dna.as_bytes();
    let n = seq.len();
    if n < window_size || window_size < 200 {
        return Vec::new();
    }

    let mut candidate_windows = Vec::new();

    for start in (0..=n - window_size).step_by(step_size.max(1)) {
        let end = start + window_size;
        let window = &seq[start..end];

        let mut c_count = 0usize;
        let mut g_count = 0usize;
        let mut cpg_count = 0usize;

        for i in 0..window.len() {
            match window[i] {
                b'C' | b'c' => {
                    c_count += 1;
                    if i + 1 < window.len() && (window[i + 1] == b'G' || window[i + 1] == b'g') {
                        cpg_count += 1;
                    }
                }
                b'G' | b'g' => {
                    g_count += 1;
                }
                _ => {}
            }
        }

        let gc = (c_count + g_count) as f64 / window_size as f64;
        let denom = (c_count * g_count) as f64;
        let obs_exp = if denom > 0.0 {
            (cpg_count as f64 * window_size as f64) / denom
        } else {
            0.0
        };

        if gc >= min_gc && obs_exp >= min_obs_exp {
            candidate_windows.push((start, end, gc, obs_exp, cpg_count));
        }
    }

    if candidate_windows.is_empty() {
        return Vec::new();
    }

    // Merge overlapping windows
    let mut merged = Vec::new();
    let (mut cur_start, mut cur_end, _, _, _) = candidate_windows[0];

    for &(w_start, w_end, _, _, _) in &candidate_windows[1..] {
        if w_start <= cur_end {
            cur_end = cur_end.max(w_end);
        } else {
            merged.push((cur_start, cur_end));
            cur_start = w_start;
            cur_end = w_end;
        }
    }
    merged.push((cur_start, cur_end));

    // Compute exact metrics for merged islands
    merged
        .into_iter()
        .map(|(start, end)| {
            let slice = &seq[start..end];
            let len = end - start;
            let mut c_count = 0usize;
            let mut g_count = 0usize;
            let mut cpg_count = 0usize;

            for i in 0..slice.len() {
                match slice[i] {
                    b'C' | b'c' => {
                        c_count += 1;
                        if i + 1 < slice.len() && (slice[i + 1] == b'G' || slice[i + 1] == b'g') {
                            cpg_count += 1;
                        }
                    }
                    b'G' | b'g' => {
                        g_count += 1;
                    }
                    _ => {}
                }
            }

            let gc = (c_count + g_count) as f64 / len as f64;
            let denom = (c_count * g_count) as f64;
            let obs_exp = if denom > 0.0 {
                (cpg_count as f64 * len as f64) / denom
            } else {
                0.0
            };

            CpgIsland {
                start,
                end,
                length: len,
                gc_content: gc,
                obs_exp_ratio: obs_exp,
                cpg_count,
            }
        })
        .collect()
}

/// Converts methylation Beta-value ($0.0 \le \beta \le 1.0$) to M-value using Du et al. (2010):
/// $$M = \log_2\left(\frac{\beta + \text{offset}}{1 - \beta + \text{offset}}\right)$$
pub fn beta_to_m_value(beta: f64, offset: f64) -> crate::Result<f64> {
    if !(0.0..=1.0).contains(&beta) {
        return Err(BioError::InvalidInput(format!("Beta value must be in [0, 1], got {}", beta)));
    }
    let num = (beta + offset).max(1e-12);
    let den = (1.0 - beta + offset).max(1e-12);
    Ok((num / den).log2())
}

/// Converts methylation M-value back to Beta-value:
/// $$\beta = \frac{2^M}{2^M + 1}$$
pub fn m_value_to_beta(m_val: f64) -> f64 {
    let two_pow_m = 2.0f64.powf(m_val);
    two_pow_m / (two_pow_m + 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_beta_m_roundtrip() {
        let beta = 0.75;
        let m = beta_to_m_value(beta, 0.0).unwrap();
        let rec = m_value_to_beta(m);
        assert!((rec - beta).abs() < 1e-6);
    }

    #[test]
    fn test_cpg_island_detection() {
        // Construct 300 bp synthetic CpG island: "CGCG..." repeated
        let rep = "CGCG".repeat(75);
        let dna = Dna::new(rep.as_bytes()).unwrap();
        let islands = find_cpg_islands(&dna, 200, 50, 0.50, 0.60);
        assert!(!islands.is_empty());
        assert!(islands[0].gc_content >= 0.99);
        assert!(islands[0].obs_exp_ratio > 0.60);
    }
}
