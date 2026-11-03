//! Histone modification marks and chromatin state segmentation (ChromHMM-like).

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Canonical histone post-translational modification marks.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum HistoneMark {
    /// H3K4me3: Tri-methylation of histone H3 lysine 4 (Active promoter mark).
    H3K4me3,
    /// H3K4me1: Mono-methylation of histone H3 lysine 4 (Enhancer mark).
    H3K4me1,
    /// H3K27ac: Acetylation of histone H3 lysine 27 (Active enhancer / promoter mark).
    H3K27ac,
    /// H3K36me3: Tri-methylation of histone H3 lysine 36 (Actively transcribed gene body).
    H3K36me3,
    /// H3K27me3: Tri-methylation of histone H3 lysine 27 (Polycomb group facultative repression).
    H3K27me3,
    /// H3K9me3: Tri-methylation of histone H3 lysine 9 (Constitutive heterochromatin).
    H3K9me3,
}

/// Chromatin functional state classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ChromatinState {
    /// Active Promoter: High H3K4me3 + H3K27ac.
    ActivePromoter,
    /// Poised / Bivalent Promoter: H3K4me3 + H3K27me3.
    BivalentPromoter,
    /// Active Enhancer: H3K4me1 + H3K27ac.
    ActiveEnhancer,
    /// Primed Enhancer: H3K4me1 without H3K27ac.
    PrimedEnhancer,
    /// Actively Transcribed Gene: High H3K36me3.
    Transcribed,
    /// Facultative Heterochromatin / Polycomb Repressed: High H3K27me3.
    PolycombRepressed,
    /// Constitutive Heterochromatin: High H3K9me3.
    Heterochromatin,
    /// Quiescent / Low signal chromatin.
    Quiescent,
}

/// Genomic segment with annotated histone profile and chromatin state.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HistoneSegment {
    /// Chromosome name.
    pub chromosome: String,
    /// 0-based start coordinate.
    pub start: usize,
    /// 0-based end coordinate.
    pub end: usize,
    /// Histone marks detected in this segment.
    pub marks: Vec<HistoneMark>,
    /// Inferred chromatin functional state.
    pub state: ChromatinState,
}

/// Infers the chromatin state from a collection of co-occurring histone marks.
pub fn infer_chromatin_state(marks: &[HistoneMark]) -> ChromatinState {
    let has_k4me3 = marks.contains(&HistoneMark::H3K4me3);
    let has_k4me1 = marks.contains(&HistoneMark::H3K4me1);
    let has_k27ac = marks.contains(&HistoneMark::H3K27ac);
    let has_k36me3 = marks.contains(&HistoneMark::H3K36me3);
    let has_k27me3 = marks.contains(&HistoneMark::H3K27me3);
    let has_k9me3 = marks.contains(&HistoneMark::H3K9me3);

    if has_k9me3 {
        ChromatinState::Heterochromatin
    } else if has_k4me3 && has_k27me3 {
        ChromatinState::BivalentPromoter
    } else if has_k4me3 {
        ChromatinState::ActivePromoter
    } else if has_k4me1 && has_k27ac {
        ChromatinState::ActiveEnhancer
    } else if has_k4me1 {
        ChromatinState::PrimedEnhancer
    } else if has_k36me3 {
        ChromatinState::Transcribed
    } else if has_k27me3 {
        ChromatinState::PolycombRepressed
    } else {
        ChromatinState::Quiescent
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chromatin_state_inference() {
        assert_eq!(
            infer_chromatin_state(&[HistoneMark::H3K4me3, HistoneMark::H3K27ac]),
            ChromatinState::ActivePromoter
        );
        assert_eq!(
            infer_chromatin_state(&[HistoneMark::H3K4me3, HistoneMark::H3K27me3]),
            ChromatinState::BivalentPromoter
        );
        assert_eq!(
            infer_chromatin_state(&[HistoneMark::H3K9me3]),
            ChromatinState::Heterochromatin
        );
        assert_eq!(
            infer_chromatin_state(&[HistoneMark::H3K4me1]),
            ChromatinState::PrimedEnhancer
        );
    }
}
