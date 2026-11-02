//! High-fidelity Longevity Design reporting, confidence scoring, and execution ordering.
//!
//! Generates comprehensive multi-tiered reports containing exact genomic coordinates,
//! candidate sgRNA guides, CFD off-target evaluations, confidence metrics, and ordered execution workflows.

use crate::aging::design::{EditType, LongevityStrategy};
use crate::aging::fidelity::HighFidelityAgeResult;
use crate::aging::risk::OncogenicRiskAssessment;
use crate::edit::crispr::SgRnaCandidate;
use crate::genome::Strand;
use crate::seq::Dna;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// An actionable high-fidelity genomic edit candidate with execution sequence metadata.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HighFidelityProposedEdit {
    /// Target gene HGNC symbol (e.g. "TERT", "SIRT6", "TP53", "ATM").
    pub target_gene: String,
    /// Nature of precision genetic modification.
    pub edit_type: EditType,
    /// Genomic chromosome locus (e.g. "chr5", "chr19").
    pub chromosome: String,
    /// Exact genomic start coordinate (0-based inclusive).
    pub start_pos: u64,
    /// Exact genomic end coordinate (0-based exclusive).
    pub end_pos: u64,
    /// Targeted strand.
    pub strand: Strand,
    /// Target sequence context.
    pub target_sequence: Dna,
    /// Biological rationale and therapeutic mechanism.
    pub rationale: String,
    /// Therapeutic priority tier (1 = critical / top tier, 5 = adjunct).
    pub priority: u8,
    /// Strict sequential execution order (1 = safety locks first, 2 = telomere renewal, etc.).
    pub execution_step: u8,
    /// Engineered candidate sgRNA with PAM and cleavage coordinates.
    pub recommended_sgrna: Option<SgRnaCandidate>,
    /// CFD off-target penalty score (0.0 = safe, 100.0 = high off-target potential).
    pub offtarget_score: f64,
    /// Composite engineering confidence score (0.0 to 1.0) factoring guide GC, CFD, and accessibility.
    pub confidence_score: f64,
}

/// Baseline versus projected longevity metrics.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct LongevityMetricsComparison {
    /// Chronological age of subject.
    pub chronological_age: f64,
    /// Baseline pre-intervention biological age.
    pub baseline_biological_age: f64,
    /// Projected post-intervention biological age.
    pub projected_biological_age: f64,
    /// Baseline terminal telomere length in bp.
    pub baseline_telomere_length_bp: f64,
    /// Projected terminal telomere length in bp.
    pub projected_telomere_length_bp: f64,
    /// Baseline remaining cell divisions (Hayflick limit).
    pub baseline_hayflick_divisions: u32,
    /// Projected remaining cell divisions after renewal.
    pub projected_hayflick_divisions: u32,
    /// Baseline shelterin capping integrity index.
    pub baseline_shelterin_capping: f64,
    /// Projected shelterin capping integrity index.
    pub projected_shelterin_capping: f64,
    /// Projected percentage extension in healthy lifespan.
    pub projected_lifespan_extension_percent: f64,
}

/// Complete High-Fidelity Longevity Genome Design Report.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HighFidelityDesignReport {
    /// Design strategy paradigm applied.
    pub strategy: LongevityStrategy,
    /// Quantitative comparison of aging parameters.
    pub metrics: LongevityMetricsComparison,
    /// Comprehensive multi-factor epigenetic and biological age breakdown.
    pub age_breakdown: HighFidelityAgeResult,
    /// Prioritized and sequentially ordered genomic modifications.
    pub proposed_edits: Vec<HighFidelityProposedEdit>,
    /// Pathway-specific oncogenic hazard evaluation and safety locks.
    pub oncogenic_risk: OncogenicRiskAssessment,
}

impl HighFidelityDesignReport {
    /// Generates a comprehensive human-readable summary report.
    pub fn export_detailed_summary(&self) -> String {
        let mut out = String::new();
        out.push_str("======================================================================\n");
        out.push_str("       SCIES-BIO-TH HIGH-FIDELITY LONGEVITY GENOME DESIGN REPORT      \n");
        out.push_str("======================================================================\n");
        out.push_str(&format!("Strategy: {:?}\n", self.strategy));
        out.push_str(&format!(
            "Chronological Age: {:.1} yrs | Biological Age: {:.1} -> Projected: {:.1} yrs\n",
            self.metrics.chronological_age,
            self.metrics.baseline_biological_age,
            self.metrics.projected_biological_age
        ));
        out.push_str(&format!(
            "Telomere Length: {:.0} bp -> Projected: {:.0} bp | Hayflick Limit: {} -> {} divisions\n",
            self.metrics.baseline_telomere_length_bp,
            self.metrics.projected_telomere_length_bp,
            self.metrics.baseline_hayflick_divisions,
            self.metrics.projected_hayflick_divisions
        ));
        out.push_str(&format!(
            "Shelterin Capping: {:.1}% -> {:.1}% | Projected Lifespan Extension: +{:.1}%\n",
            self.metrics.baseline_shelterin_capping * 100.0,
            self.metrics.projected_shelterin_capping * 100.0,
            self.metrics.projected_lifespan_extension_percent
        ));

        out.push_str("\n--- Oncogenic Risk Assessment & Safeguards ---\n");
        out.push_str(&format!(
            "Net Risk Index: {:.3} / 1.000 [Status: {}]\n",
            self.oncogenic_risk.net_risk_index,
            if self.oncogenic_risk.is_acceptable_safety_profile { "CLINICALLY SAFE" } else { "ELEVATED HAZARD" }
        ));
        out.push_str(&format!(
            "Active Safeguards: Super-p53: {} | Suicide Switch: {} | HiFi-Cas: {} | Mitigation: {:.1}%\n",
            self.oncogenic_risk.safeguards.super_p53_active,
            self.oncogenic_risk.safeguards.inducible_suicide_switch,
            self.oncogenic_risk.safeguards.high_fidelity_nuclease,
            self.oncogenic_risk.mitigation_factor * 100.0
        ));
        out.push_str(&format!(
            "Raw Hazards: Crisis Bypass: {:.2} | Proto-Oncogene: {:.2} | Checkpoint: {:.2} | Instability: {:.2}\n",
            self.oncogenic_risk.raw_pathways.crisis_bypass_risk,
            self.oncogenic_risk.raw_pathways.proto_oncogene_transactivation_risk,
            self.oncogenic_risk.raw_pathways.checkpoint_loss_risk,
            self.oncogenic_risk.raw_pathways.chromosomal_instability_risk
        ));

        out.push_str("\n--- Sequential Intervention Workflow (Step-by-Step) ---\n");
        let mut sorted_edits = self.proposed_edits.clone();
        sorted_edits.sort_by_key(|e| (e.execution_step, e.priority));

        for edit in &sorted_edits {
            out.push_str(&format!(
                "Step {}. [{}] Target: {} ({}:{}-{}) | Priority: {}\n",
                edit.execution_step,
                edit.target_gene,
                edit.target_gene,
                edit.chromosome,
                edit.start_pos,
                edit.end_pos,
                edit.priority
            ));
            out.push_str(&format!("   Action: {:?} | Strand: {:?}\n", edit.edit_type, edit.strand));
            out.push_str(&format!("   Rationale: {}\n", edit.rationale));
            out.push_str(&format!(
                "   Design Confidence: {:.1}% | Off-Target CFD Risk: {:.2}\n",
                edit.confidence_score * 100.0,
                edit.offtarget_score
            ));
            if let Some(ref guide) = edit.recommended_sgrna {
                out.push_str(&format!(
                    "   sgRNA Guide (20nt): {} | PAM @ {}\n",
                    String::from_utf8_lossy(&guide.sequence.sequence),
                    guide.pam_index
                ));
            }
        }
        out.push_str("======================================================================\n");
        out
    }
}
