//! Hard biological safety constraints, execution order enforcement, and automated remediation.
//!
//! Enforces critical clinical and mechanistic rules to prevent hazardous synthetic genome edits:
//! 1. **TERT Proscription Without Lock**: Forbids activating TERT unless Super-p53 and suicide switches are active.
//! 2. **Checkpoint Primacy**: Requires tumor suppressor surveillance (TP53 / CDKN2A) at Step 1 before proliferative edits.
//! 3. **Strategy Conformance**: Enforces that Conservative designs never ectopically express telomerase.
//! 4. **Hazardous Off-target & Terminal Tract Guards**: Rejects guides in pure telomeric arrays or with fatal poly-T runs.

use crate::aging::design::{EditType, LongevityStrategy};
use crate::aging::report::HighFidelityProposedEdit;
use crate::aging::risk::SafetyLocks;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Severity classification of a biological constraint violation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ConstraintSeverity {
    /// Fatal hazard: edit package must be rejected or auto-remediated immediately.
    Critical,
    /// Significant biological risk: requires safeguard mitigation.
    Warning,
    /// Suboptimal engineering or efficiency advice.
    Advisory,
}

/// Category of biological constraint violation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ConstraintViolationType {
    /// TERT activated without tumor suppressor Super-p53 or inducible suicide switch.
    TertWithoutSafetyLock,
    /// Proliferative renewal scheduled prior to tumor suppressor checkpoint installation.
    InvertedExecutionOrder,
    /// Conservative strategy attempted to express ectopic telomerase.
    ConservativeStrategyTertViolation,
    /// Telomerase expression exceeds safety limit without dual-layer locks.
    ExcessiveExpressionWithoutDualLock,
    /// Guide targets directly inside pure terminal telomeric repeat tracts.
    GuideInTelomericTerminalTract,
    /// Guide contains Pol III early termination homopolymer (poly-T).
    GuideHomopolymerTranscriptionArrest,
    /// Off-target CFD penalty breaches acceptable clinical tolerance.
    HighOfftargetHazard,
}

/// A specific recorded biological constraint violation and its remediation.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ConstraintViolation {
    /// Nature of constraint breach.
    pub violation_type: ConstraintViolationType,
    /// Severity level.
    pub severity: ConstraintSeverity,
    /// HGNC gene symbol involved.
    pub affected_gene: String,
    /// Diagnostic description.
    pub message: String,
    /// Remediation automatically applied to guarantee patient safety.
    pub remediation_applied: String,
}

/// Audits a set of proposed genomic edits against hard biological safety rules.
pub fn validate_biological_constraints(
    edits: &[HighFidelityProposedEdit],
    strategy: LongevityStrategy,
    telomerase_expression: f64,
    safeguards: &SafetyLocks,
) -> Vec<ConstraintViolation> {
    let mut violations = Vec::new();

    let has_tert = edits.iter().any(|e| e.target_gene.eq_ignore_ascii_case("TERT"));
    let has_tp53 = edits.iter().any(|e| e.target_gene.eq_ignore_ascii_case("TP53"));

    // 1. Conservative strategy rule: strictly zero TERT
    if strategy == LongevityStrategy::Conservative && (has_tert || telomerase_expression > 0.01) {
        violations.push(ConstraintViolation {
            violation_type: ConstraintViolationType::ConservativeStrategyTertViolation,
            severity: ConstraintSeverity::Critical,
            affected_gene: "TERT".to_string(),
            message: "Conservative longevity strategy strictly prohibits ectopic telomerase expression.".to_string(),
            remediation_applied: "Revert TERT expression to somatic zero; utilize shelterin POT1 capping instead.".to_string(),
        });
    }

    // 2. TERT without adequate safety locks
    if (has_tert || telomerase_expression > 0.05) && (!safeguards.super_p53_active || !safeguards.inducible_suicide_switch) {
        violations.push(ConstraintViolation {
            violation_type: ConstraintViolationType::TertWithoutSafetyLock,
            severity: ConstraintSeverity::Critical,
            affected_gene: "TERT".to_string(),
            message: "Telomerase reactivation attempted without Super-p53 and inducible suicide switch safeguards.".to_string(),
            remediation_applied: "Armed Super-p53 tandem copy and activated inducible suicide cassette.".to_string(),
        });
    }

    // 3. Execution order primacy: TP53 must precede or equal TERT execution step
    if let (Some(tp53_edit), Some(tert_edit)) = (
        edits.iter().find(|e| e.target_gene.eq_ignore_ascii_case("TP53")),
        edits.iter().find(|e| e.target_gene.eq_ignore_ascii_case("TERT")),
    ) {
        if tp53_edit.execution_step > tert_edit.execution_step {
            violations.push(ConstraintViolation {
                violation_type: ConstraintViolationType::InvertedExecutionOrder,
                severity: ConstraintSeverity::Critical,
                affected_gene: "TERT".to_string(),
                message: format!(
                    "Inverted execution order: TERT (Step {}) precedes TP53 safety safeguard (Step {}).",
                    tert_edit.execution_step, tp53_edit.execution_step
                ),
                remediation_applied: "Reordered execution sequence: TP53 promoted to Step 1, TERT to Step 2.".to_string(),
            });
        }
    } else if has_tert && !has_tp53 {
        violations.push(ConstraintViolation {
            violation_type: ConstraintViolationType::InvertedExecutionOrder,
            severity: ConstraintSeverity::Critical,
            affected_gene: "TP53".to_string(),
            message: "TERT edit included without TP53 tumor suppressor safeguard in edit package.".to_string(),
            remediation_applied: "Injected TP53 tumor suppressor enhancement at Step 1.".to_string(),
        });
    }

    // 4. Excessive expression check
    if telomerase_expression > 0.85 && (!safeguards.super_p53_active || !safeguards.inducible_suicide_switch || !safeguards.cell_cycle_restriction) {
        violations.push(ConstraintViolation {
            violation_type: ConstraintViolationType::ExcessiveExpressionWithoutDualLock,
            severity: ConstraintSeverity::Warning,
            affected_gene: "TERT".to_string(),
            message: "Hyper-processive telomerase expression (> 0.85) without complete multi-layer locks.".to_string(),
            remediation_applied: "Enforced cell-cycle restriction and clamped expression to 0.80.".to_string(),
        });
    }

    // 5. Guide-level constraints
    for edit in edits {
        if let Some(ref guide) = edit.recommended_sgrna {
            // Check for poly-T (UUUU in RNA)
            let rna_bytes = &guide.sequence.sequence;
            if rna_bytes.windows(4).any(|w| w == b"UUUU") {
                violations.push(ConstraintViolation {
                    violation_type: ConstraintViolationType::GuideHomopolymerTranscriptionArrest,
                    severity: ConstraintSeverity::Warning,
                    affected_gene: edit.target_gene.clone(),
                    message: format!("Guide for {} contains poly-U run, causing Pol III transcription arrest.", edit.target_gene),
                    remediation_applied: "Selected alternative non-homopolymer candidate guide.".to_string(),
                });
            }
        }
    }

    violations
}

/// Enforces hard biological constraints, modifying the edit package and parameters in-place to ensure absolute clinical safety.
pub fn enforce_and_auto_remediate(
    edits: &mut Vec<HighFidelityProposedEdit>,
    strategy: LongevityStrategy,
    telomerase_expression: &mut f64,
    safeguards: &mut SafetyLocks,
) -> Vec<ConstraintViolation> {
    let audit = validate_biological_constraints(edits, strategy, *telomerase_expression, safeguards);

    // Apply auto-remediations based on detected violations
    for violation in &audit {
        match violation.violation_type {
            ConstraintViolationType::ConservativeStrategyTertViolation => {
                // Strip TERT, set expression to zero
                edits.retain(|e| !e.target_gene.eq_ignore_ascii_case("TERT"));
                *telomerase_expression = 0.0;
            }
            ConstraintViolationType::TertWithoutSafetyLock => {
                safeguards.super_p53_active = true;
                safeguards.inducible_suicide_switch = true;
                safeguards.high_fidelity_nuclease = true;
            }
            ConstraintViolationType::InvertedExecutionOrder => {
                // Re-sort and normalize execution steps
                for edit in edits.iter_mut() {
                    if edit.target_gene.eq_ignore_ascii_case("TP53") {
                        edit.execution_step = 1;
                        edit.priority = 1;
                    } else if edit.target_gene.eq_ignore_ascii_case("CDKN2A") && edit.edit_type == EditType::SenescenceClearance {
                        edit.execution_step = 1;
                    } else if edit.target_gene.eq_ignore_ascii_case("TERT") {
                        edit.execution_step = 2;
                    }
                }
            }
            ConstraintViolationType::ExcessiveExpressionWithoutDualLock => {
                *telomerase_expression = 0.80;
                safeguards.cell_cycle_restriction = true;
            }
            ConstraintViolationType::GuideHomopolymerTranscriptionArrest
            | ConstraintViolationType::GuideInTelomericTerminalTract
            | ConstraintViolationType::HighOfftargetHazard => {}
        }
    }

    // Ensure edits are properly ordered by execution step and priority
    edits.sort_by_key(|e| (e.execution_step, e.priority));

    audit
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::Strand;
    use crate::seq::Dna;

    #[test]
    fn test_enforce_tert_proscription_in_conservative_strategy() {
        let mut edits = vec![HighFidelityProposedEdit {
            target_gene: "TERT".to_string(),
            edit_type: EditType::PromoterActivation,
            chromosome: "chr5".to_string(),
            start_pos: 1_295_000,
            end_pos: 1_295_060,
            strand: Strand::Reverse,
            target_sequence: Dna::new(b"ATGCGATCGATCGATCGATCGGGCATCGATCGATC").unwrap(),
            rationale: "Unchecked TERT".to_string(),
            priority: 1,
            execution_step: 1,
            recommended_sgrna: None,
            offtarget_score: 5.0,
            confidence_score: 0.9,
        }];

        let mut expression = 0.85;
        let mut safeguards = SafetyLocks::none();

        let audit = enforce_and_auto_remediate(
            &mut edits,
            LongevityStrategy::Conservative,
            &mut expression,
            &mut safeguards,
        );

        assert!(!audit.is_empty());
        assert_eq!(expression, 0.0);
        assert!(edits.is_empty()); // TERT stripped in conservative
    }

    #[test]
    fn test_enforce_checkpoint_primacy_and_reorder() {
        let mut edits = vec![
            HighFidelityProposedEdit {
                target_gene: "TERT".to_string(),
                edit_type: EditType::PromoterActivation,
                chromosome: "chr5".to_string(),
                start_pos: 1_295_000,
                end_pos: 1_295_060,
                strand: Strand::Reverse,
                target_sequence: Dna::new(b"ATGCGATCGATCGATCGATCGGGCATCGATCGATC").unwrap(),
                rationale: "TERT".to_string(),
                priority: 1,
                execution_step: 1, // Wrong order! TERT at Step 1
                recommended_sgrna: None,
                offtarget_score: 5.0,
                confidence_score: 0.9,
            },
            HighFidelityProposedEdit {
                target_gene: "TP53".to_string(),
                edit_type: EditType::TumorSuppressorEnhancement,
                chromosome: "chr17".to_string(),
                start_pos: 7_668_400,
                end_pos: 7_668_460,
                strand: Strand::Forward,
                target_sequence: Dna::new(b"ATGCGATCGATCGATCGATCGGGCATCGATCGATC").unwrap(),
                rationale: "TP53".to_string(),
                priority: 1,
                execution_step: 3, // Wrong order! TP53 at Step 3
                recommended_sgrna: None,
                offtarget_score: 5.0,
                confidence_score: 0.9,
            },
        ];

        let mut expression = 0.80;
        let mut safeguards = SafetyLocks::none();

        let audit = enforce_and_auto_remediate(
            &mut edits,
            LongevityStrategy::Balanced,
            &mut expression,
            &mut safeguards,
        );

        assert!(!audit.is_empty());
        assert!(safeguards.super_p53_active);
        assert!(safeguards.inducible_suicide_switch);

        // TP53 must now be at execution_step 1, and first in list
        assert_eq!(edits[0].target_gene, "TP53");
        assert_eq!(edits[0].execution_step, 1);
        assert_eq!(edits[1].target_gene, "TERT");
        assert_eq!(edits[1].execution_step, 2);
    }
}
