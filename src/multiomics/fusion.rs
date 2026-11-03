//! Multi-omics locus fusion integrating transcriptomics (RNA-seq TPM), chromatin accessibility
//! (ATAC-seq), and DNA methylation ($\beta$-value) into a unified regulatory state score.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Multi-layer omics measurements at a gene locus.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct OmicsLocusInput {
    /// Gene symbol or locus identifier.
    pub gene_symbol: String,
    /// RNA-seq expression in Transcripts Per Million ($\text{TPM} \ge 0.0$).
    pub rna_tpm: f64,
    /// ATAC-seq chromatin accessibility fraction in $[0.0, 1.0]$.
    pub atac_accessibility: f64,
    /// Promoter CpG methylation $\beta$-value in $[0.0, 1.0]$ ($0 = \text{unmethylated}$, $1 = \text{hypermethylated}$).
    pub promoter_methylation_beta: f64,
    /// Optional H3K27ac enhancer ABC score in $[0.0, 1.0]$.
    pub enhancer_abc_score: f64,
}

/// Unified multi-omics regulatory activity score and concordance metrics.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MultiOmicsFusionResult {
    /// Normalized log-scale transcriptional activity in $[0.0, 1.0]$.
    pub expression_component: f64,
    /// Epigenomic permissiveness ($0.6 \cdot \text{ATAC} + 0.4 \cdot (1 - \beta)$) in $[0.0, 1.0]$.
    pub epigenomic_permissiveness: f64,
    /// Unified multi-omics activity score in $[0.0, 1.0]$.
    pub unified_activity_score: f64,
    /// Cross-layer concordance between epigenome and transcriptome in $[0.0, 1.0]$.
    pub layer_concordance: f64,
}

/// Fuses RNA-seq, ATAC-seq, promoter CpG methylation, and enhancer activity into a unified score.
pub fn fuse_multiomics_locus(input: &OmicsLocusInput) -> MultiOmicsFusionResult {
    // Log2(TPM + 1) normalized against saturation ~ 10.0 (TPM ~ 1023)
    let log_tpm = (input.rna_tpm.max(0.0) + 1.0).log2();
    let expression_component = (log_tpm / 10.0).clamp(0.0, 1.0);

    let unmethylated = (1.0 - input.promoter_methylation_beta.clamp(0.0, 1.0)).clamp(0.0, 1.0);
    let atac = input.atac_accessibility.clamp(0.0, 1.0);
    let enh = input.enhancer_abc_score.clamp(0.0, 1.0);

    let epigenomic_permissiveness = (0.50 * atac + 0.35 * unmethylated + 0.15 * enh).clamp(0.0, 1.0);

    let unified_activity_score =
        (0.55 * expression_component + 0.45 * epigenomic_permissiveness).clamp(0.0, 1.0);

    let layer_concordance =
        (1.0 - (expression_component - epigenomic_permissiveness).abs()).clamp(0.0, 1.0);

    MultiOmicsFusionResult {
        expression_component,
        epigenomic_permissiveness,
        unified_activity_score,
        layer_concordance,
    }
}
