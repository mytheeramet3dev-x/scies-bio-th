//! Core promoter cis-regulatory element detection (`TATA-box`, `Inr`, `BREu`, `DPE`).

use crate::seq::Dna;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Type of core RNA Polymerase II promoter motif.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum CorePromoterMotifType {
    /// Goldberg-Hogness TATA-box (`TATAWAWR`, typically $-31$ to $-24$ upstream of TSS).
    TataBox,
    /// Initiator element (`YYANWYY`, spanning the Transcription Start Site $+1$).
    InitiatorInr,
    /// TFIIB Recognition Element upstream (`SSRCGCC`, immediately upstream of TATA).
    BreUpstream,
    /// Downstream Promoter Element (`RGWYV`, $+28$ to $+32$ downstream of TSS).
    DownstreamPromoterDpe,
}

/// A detected core promoter element.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CorePromoterElement {
    /// Element category.
    pub motif_type: CorePromoterMotifType,
    /// 0-based start offset on the sequence.
    pub start: usize,
    /// 0-based end offset (exclusive).
    pub end: usize,
    /// Sequence string matched.
    pub sequence: String,
    /// Match quality score in $[0.0, 1.0]$.
    pub confidence_score: f64,
}

fn matches_iupac(base: u8, code: u8) -> bool {
    let b = base.to_ascii_uppercase();
    match code {
        b'A' => b == b'A',
        b'C' => b == b'C',
        b'G' => b == b'G',
        b'T' => b == b'T',
        b'W' => b == b'A' || b == b'T',
        b'S' => b == b'C' || b == b'G',
        b'R' => b == b'A' || b == b'G',
        b'Y' => b == b'C' || b == b'T',
        b'V' => b == b'A' || b == b'C' || b == b'G',
        b'N' => matches!(b, b'A' | b'C' | b'G' | b'T'),
        _ => false,
    }
}

fn score_iupac_pattern(window: &[u8], pattern: &[u8]) -> f64 {
    if window.len() != pattern.len() {
        return 0.0;
    }
    let matches = window
        .iter()
        .zip(pattern.iter())
        .filter(|(&b, &p)| matches_iupac(b, p))
        .count();
    matches as f64 / pattern.len() as f64
}

/// Scans a DNA sequence for canonical Pol II core promoter elements (`TATA-box`, `Inr`, `BREu`, `DPE`).
pub fn detect_core_promoter_elements(dna: &Dna, min_score: f64) -> Vec<CorePromoterElement> {
    let seq = &dna.sequence;
    let n = seq.len();
    let mut elements = Vec::new();

    let patterns: &[(CorePromoterMotifType, &[u8])] = &[
        (CorePromoterMotifType::TataBox, b"TATAWAWR"),
        (CorePromoterMotifType::InitiatorInr, b"YYANWYY"),
        (CorePromoterMotifType::BreUpstream, b"SSRCGCC"),
        (CorePromoterMotifType::DownstreamPromoterDpe, b"RGWYV"),
    ];

    for &(motif_type, pat) in patterns {
        let k = pat.len();
        if n < k {
            continue;
        }
        for i in 0..=(n - k) {
            let win = &seq[i..(i + k)];
            let score = score_iupac_pattern(win, pat);
            if score >= min_score {
                elements.push(CorePromoterElement {
                    motif_type,
                    start: i,
                    end: i + k,
                    sequence: String::from_utf8_lossy(win).into_owned(),
                    confidence_score: score,
                });
            }
        }
    }

    elements.sort_by_key(|e| e.start);
    elements
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_tata_box_and_inr() {
        // Contains BREu (GGGCGCC), TATA-box (TATAAAAG), and Inr (TCAGTTT)
        let seq = Dna::new(b"GGGCGCCTATAAAAGGGGGGGGGTCAGTTT").unwrap();
        let elems = detect_core_promoter_elements(&seq, 0.85);
        assert!(elems.iter().any(|e| e.motif_type == CorePromoterMotifType::BreUpstream));
        assert!(elems.iter().any(|e| e.motif_type == CorePromoterMotifType::TataBox));
        assert!(elems.iter().any(|e| e.motif_type == CorePromoterMotifType::InitiatorInr));
    }
}
