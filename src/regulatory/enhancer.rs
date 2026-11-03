//! Enhancer activity scoring integrating H3K27ac ChIP-seq acetylation and ATAC-seq accessibility
//! (Fulco et al. Activity-by-Contact / ABC model inspired).

use crate::epigenome::histone::quantify_h3k27ac_signal;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Functional enhancer classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum EnhancerClass {
    /// Super-Enhancer: High H3K27ac + high accessibility driving cell-identity genes.
    SuperEnhancer,
    /// Active typical enhancer.
    ActiveEnhancer,
    /// Poised or weak enhancer.
    PoisedEnhancer,
    /// Inactive / closed chromatin.
    Inactive,
}

/// Quantitative enhancer activity and Activity-by-Contact (ABC) score.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct EnhancerAssessment {
    /// Geometric mean of chromatin accessibility and H3K27ac activity in $[0.0, 1.0]$.
    pub enhancer_activity_score: f64,
    /// Distance-attenuated 3D contact weight to target TSS in $[0.0, 1.0]$.
    pub contact_weight: f64,
    /// Combined Activity-by-Contact (ABC) regulatory score in $[0.0, 1.0]$.
    pub abc_score: f64,
    /// Functional classification.
    pub class: EnhancerClass,
}

/// Scores enhancer regulatory activity from H3K27ac ChIP/input signals, ATAC accessibility ($[0,1]$),
/// and genomic distance (`distance_to_tss_bp`) to the target promoter.
pub fn score_enhancer_activity(
    h3k27ac_chip: f64,
    h3k27ac_input: f64,
    atac_accessibility: f64,
    distance_to_tss_bp: usize,
) -> EnhancerAssessment {
    let k27 = quantify_h3k27ac_signal(h3k27ac_chip, h3k27ac_input);
    let atac = atac_accessibility.clamp(0.0, 1.0);

    // Activity = sqrt(ATAC * H3K27ac)
    let enhancer_activity_score = (atac * k27.normalized_activity).sqrt().clamp(0.0, 1.0);

    // Power-law 3D chromatin contact decay: (5000 / max(dist, 5000))^0.7
    let d_eff = (distance_to_tss_bp.max(2500)) as f64;
    let contact_weight = (2500.0 / d_eff).powf(0.65).clamp(0.05, 1.0);
    let abc_score = (enhancer_activity_score * contact_weight).clamp(0.0, 1.0);

    let class = if k27.is_super_enhancer_like && atac >= 0.65 {
        EnhancerClass::SuperEnhancer
    } else if enhancer_activity_score >= 0.45 {
        EnhancerClass::ActiveEnhancer
    } else if atac >= 0.30 || k27.normalized_activity >= 0.25 {
        EnhancerClass::PoisedEnhancer
    } else {
        EnhancerClass::Inactive
    };

    EnhancerAssessment {
        enhancer_activity_score,
        contact_weight,
        abc_score,
        class,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_score_enhancer_activity_classes() {
        let super_enh = score_enhancer_activity(30.0, 1.5, 0.85, 5_000);
        assert_eq!(super_enh.class, EnhancerClass::SuperEnhancer);
        assert!(super_enh.enhancer_activity_score > 0.80);

        let inactive = score_enhancer_activity(0.2, 1.5, 0.05, 100_000);
        assert_eq!(inactive.class, EnhancerClass::Inactive);
    }
}
