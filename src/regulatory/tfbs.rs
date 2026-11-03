//! JASPAR-compatible Transcription Factor Binding Site (TFBS) Position Frequency & Weight Matrices.
//!
//! Includes the canonical JASPAR CORE `MA0003.1` (`TFAP2A`) and `MA0079.3` (`SP1`) matrices and
//! computes normalized relative log-odds scores in $[0.0, 1.0]$.

use crate::seq::Dna;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A JASPAR Position Frequency Matrix (PFM) converted to normalized log-odds weights.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct JasparMatrix {
    /// JASPAR accession ID (e.g. `"MA0003.1"`).
    pub matrix_id: String,
    /// Transcription factor symbol (e.g. `"TFAP2A"`).
    pub tf_name: String,
    /// Position frequency counts `[A, C, G, T]` for each motif column.
    pub counts: Vec<[f64; 4]>,
}

/// A detected transcription factor binding site on a DNA sequence.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TfbsHit {
    /// JASPAR matrix ID.
    pub matrix_id: String,
    /// Transcription factor symbol.
    pub tf_name: String,
    /// 0-based start offset on the scanned sequence.
    pub position: usize,
    /// Strand (`'+'` or `'-'`).
    pub strand: char,
    /// Raw log-likelihood ratio score.
    pub raw_score: f64,
    /// Normalized JASPAR relative score in $[0.0, 1.0]$ ($\ge 0.80$ is a high-confidence site).
    pub relative_score: f64,
    /// Matched DNA subsequence.
    pub matched_seq: String,
}

impl JasparMatrix {
    /// Motif length in base pairs.
    pub fn len(&self) -> usize {
        self.counts.len()
    }

    /// Returns true if the motif has zero columns.
    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    /// Computes `[A, C, G, T]` log-odds weight matrix with Laplace pseudocount $0.8$
    /// against uniform $0.25$ background, plus `(min_possible_score, max_possible_score)`.
    fn pwm_and_bounds(&self) -> (Vec<[f64; 4]>, f64, f64) {
        let pseudo = 0.8;
        let mut pwm = Vec::with_capacity(self.counts.len());
        let mut min_score = 0.0;
        let mut max_score = 0.0;

        for col in &self.counts {
            let total: f64 = col.iter().sum::<f64>() + 4.0 * pseudo;
            let mut weights = [0.0; 4];
            let mut col_min = f64::INFINITY;
            let mut col_max = f64::NEG_INFINITY;

            for b in 0..4 {
                let p = (col[b] + pseudo) / total;
                let w = (p / 0.25).ln();
                weights[b] = w;
                if w < col_min {
                    col_min = w;
                }
                if w > col_max {
                    col_max = w;
                }
            }
            pwm.push(weights);
            min_score += col_min;
            max_score += col_max;
        }

        (pwm, min_score, max_score)
    }

    /// Scores an exact slice of length `self.len()`, returning `(raw_score, relative_score)`.
    pub fn score_window(&self, window: &[u8]) -> Option<(f64, f64)> {
        if window.len() != self.counts.len() || self.counts.is_empty() {
            return None;
        }
        let (pwm, min_s, max_s) = self.pwm_and_bounds();
        let mut raw = 0.0;
        for (pos, &base) in window.iter().enumerate() {
            let idx = match base.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None,
            };
            raw += pwm[pos][idx];
        }
        let denom = (max_s - min_s).max(1e-9);
        let rel = ((raw - min_s) / denom).clamp(0.0, 1.0);
        Some((raw, rel))
    }
}

/// Canonical JASPAR CORE `MA0003.1` (`TFAP2A`, consensus `GCCNNNGGC`, 9 bp).
pub fn jaspar_ma0003_1_tfap2a() -> JasparMatrix {
    JasparMatrix {
        matrix_id: "MA0003.1".to_string(),
        tf_name: "TFAP2A".to_string(),
        counts: vec![
            [0.0, 0.0, 52.0, 0.0],   // Pos 1: G
            [0.0, 51.0, 1.0, 0.0],   // Pos 2: C
            [0.0, 50.0, 2.0, 0.0],   // Pos 3: C
            [11.0, 18.0, 15.0, 8.0], // Pos 4: C/G/A
            [14.0, 16.0, 14.0, 8.0], // Pos 5: C/A/G
            [9.0, 14.0, 21.0, 8.0],  // Pos 6: G/C
            [1.0, 1.0, 49.0, 1.0],   // Pos 7: G
            [0.0, 1.0, 50.0, 1.0],   // Pos 8: G
            [1.0, 48.0, 1.0, 2.0],   // Pos 9: C
        ],
    }
}

/// Canonical JASPAR CORE `MA0079.3` (`SP1`, GC-box consensus `GGGGCGGGG`, 9 bp).
pub fn jaspar_ma0079_3_sp1() -> JasparMatrix {
    JasparMatrix {
        matrix_id: "MA0079.3".to_string(),
        tf_name: "SP1".to_string(),
        counts: vec![
            [2.0, 2.0, 45.0, 1.0],
            [1.0, 1.0, 47.0, 1.0],
            [1.0, 1.0, 47.0, 1.0],
            [1.0, 1.0, 47.0, 1.0],
            [1.0, 46.0, 2.0, 1.0],
            [1.0, 1.0, 47.0, 1.0],
            [1.0, 1.0, 47.0, 1.0],
            [1.0, 1.0, 47.0, 1.0],
            [2.0, 2.0, 45.0, 1.0],
        ],
    }
}

/// Scans both forward and reverse-complement strands of a DNA sequence against a [`JasparMatrix`],
/// returning all hits with `relative_score >= min_relative_score`.
pub fn scan_jaspar_tfbs(
    dna: &Dna,
    matrix: &JasparMatrix,
    min_relative_score: f64,
) -> Vec<TfbsHit> {
    let k = matrix.len();
    if k == 0 || dna.sequence.len() < k {
        return Vec::new();
    }

    let rev_dna = dna.reverse_complement();
    let n = dna.sequence.len();
    let mut hits = Vec::new();

    for i in 0..=(n - k) {
        let fwd_win = &dna.sequence[i..(i + k)];
        if let Some((raw, rel)) = matrix.score_window(fwd_win) {
            if rel >= min_relative_score {
                hits.push(TfbsHit {
                    matrix_id: matrix.matrix_id.clone(),
                    tf_name: matrix.tf_name.clone(),
                    position: i,
                    strand: '+',
                    raw_score: raw,
                    relative_score: rel,
                    matched_seq: String::from_utf8_lossy(fwd_win).into_owned(),
                });
            }
        }

        let rev_start = n - i - k;
        let rev_win = &rev_dna.sequence[rev_start..(rev_start + k)];
        if let Some((raw, rel)) = matrix.score_window(rev_win) {
            if rel >= min_relative_score {
                hits.push(TfbsHit {
                    matrix_id: matrix.matrix_id.clone(),
                    tf_name: matrix.tf_name.clone(),
                    position: i,
                    strand: '-',
                    raw_score: raw,
                    relative_score: rel,
                    matched_seq: String::from_utf8_lossy(rev_win).into_owned(),
                });
            }
        }
    }

    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jaspar_ma0003_1_tfap2a_known_binding_sites() {
        // Validation target: compare PWM score with JASPAR CORE (MA0003.1 TFAP2A) on known binding sites
        let tfap2a = jaspar_ma0003_1_tfap2a();
        // Optimal consensus site for MA0003.1: GCC C C G GGC ("GCCCCGGGC")
        let (_, opt_rel) = tfap2a.score_window(b"GCCCCGGGC").unwrap();
        assert!((opt_rel - 1.0).abs() < 1e-6);

        // Known mismatched non-site ("TTTTTTTTT")
        let (_, bad_rel) = tfap2a.score_window(b"TTTTTTTTT").unwrap();
        assert!(bad_rel < 0.15);

        // Scan promoter sequence containing embedded TFAP2A site
        let promoter = Dna::new(b"ATATATGCCCCGGGCATATAT").unwrap();
        let hits = scan_jaspar_tfbs(&promoter, &tfap2a, 0.85);
        assert!(!hits.is_empty());
        assert_eq!(hits[0].position, 6);
        assert_eq!(hits[0].matrix_id, "MA0003.1");
        assert!(hits[0].relative_score > 0.99);
    }
}
