//! Composable pre-intervention safety rule engine (`safety::audit`).

use crate::safety::checkpoint::{evaluate_checkpoint_activation, CellFateDecision, CheckpointInput};
use crate::safety::immune::{estimate_innate_immune_trigger, DeliveryPayloadProfile};
use crate::safety::oncogene::{evaluate_oncogenic_activation_risk, OncogeneInterventionInput};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Severity of a safety audit rule violation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum AuditSeverity {
    /// Informational or advisory notice.
    Advisory,
    /// Warning requiring monitoring or dose adjustment.
    Warning,
    /// Critical biological safety violation blocking intervention.
    CriticalBlock,
}

/// Built-in formal safety rules for genetic and epigenetic interventions.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum SafetyRule {
    /// Requires intact p53/Rb checkpoint activation before allowing telomerase/mitogenic boost.
    CheckpointPrimacy,
    /// Caps composite oncogenic risk below the specified threshold (e.g. $0.15$).
    MaxOncogenicRisk(f64),
    /// Caps innate immune interferon trigger risk below the specified threshold (e.g. $0.25$).
    MaxImmuneTrigger(f64),
    /// Forbids proto-oncogene `MYC` overexpression above a specified fold limit (e.g. $1.05$).
    ForbidMycOverexpression(f64),
    /// Requires minimum CRISPR off-target specificity (maximum cleavage burden).
    MaxOffTargetBurden(f64),
}

/// Complete pre-intervention biological context evaluated by the audit engine.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SafetyContext {
    /// Checkpoint pathway state.
    pub checkpoint: CheckpointInput,
    /// Proto-oncogene and tumor-suppressor state.
    pub oncogene: OncogeneInterventionInput,
    /// Nucleic-acid / vector delivery immunogenicity profile.
    pub delivery: DeliveryPayloadProfile,
}

/// Evaluation outcome of a single [`SafetyRule`].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct AuditRuleResult {
    /// Rule name identifier.
    pub rule_name: String,
    /// Whether the condition passed.
    pub passed: bool,
    /// Violation severity if failed.
    pub severity: AuditSeverity,
    /// Quantitative metric observed.
    pub observed_value: f64,
    /// Diagnostic explanation.
    pub message: String,
}

/// Aggregate report produced by [`evaluate_intervention_safety`].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SafetyAuditReport {
    /// True if zero `CriticalBlock` violations occurred.
    pub approved: bool,
    /// Individual rule evaluations.
    pub results: Vec<AuditRuleResult>,
    /// Number of critical blocking violations.
    pub critical_violations: usize,
    /// Number of warnings.
    pub warnings: usize,
}

/// Evaluates a slice of composable [`SafetyRule`] assertions against a [`SafetyContext`].
pub fn evaluate_intervention_safety(
    ctx: &SafetyContext,
    rules: &[SafetyRule],
) -> SafetyAuditReport {
    let cp_resp = evaluate_checkpoint_activation(&ctx.checkpoint);
    let onco_resp = evaluate_oncogenic_activation_risk(&ctx.oncogene);
    let imm_resp = estimate_innate_immune_trigger(&ctx.delivery);

    let mut results = Vec::with_capacity(rules.len());
    let mut critical_violations = 0;
    let mut warnings = 0;

    for rule in rules {
        let res = match rule {
            SafetyRule::CheckpointPrimacy => {
                let needs_guard = ctx.oncogene.tert_fold > 1.1 || ctx.checkpoint.mitogenic_drive_fold > 1.2;
                let has_guard = ctx.checkpoint.tp53_expression_fold >= 1.2
                    && cp_resp.fate != CellFateDecision::GenomicInstabilityBypass;
                let passed = !needs_guard || has_guard;
                AuditRuleResult {
                    rule_name: "CheckpointPrimacy".to_string(),
                    passed,
                    severity: AuditSeverity::CriticalBlock,
                    observed_value: ctx.checkpoint.tp53_expression_fold,
                    message: if passed {
                        "p53/Rb checkpoint primacy verified prior to proliferative intervention".to_string()
                    } else {
                        "Telomerase or mitogenic intervention attempted without TP53 >= 1.2x checkpoint guard".to_string()
                    },
                }
            }
            SafetyRule::MaxOncogenicRisk(max_risk) => {
                let passed = onco_resp.composite_oncogenic_risk <= *max_risk;
                AuditRuleResult {
                    rule_name: "MaxOncogenicRisk".to_string(),
                    passed,
                    severity: AuditSeverity::CriticalBlock,
                    observed_value: onco_resp.composite_oncogenic_risk,
                    message: format!(
                        "Composite oncogenic risk {:.4} (limit <= {:.4})",
                        onco_resp.composite_oncogenic_risk, max_risk
                    ),
                }
            }
            SafetyRule::MaxImmuneTrigger(max_imm) => {
                let passed = imm_resp.composite_immune_risk <= *max_imm;
                AuditRuleResult {
                    rule_name: "MaxImmuneTrigger".to_string(),
                    passed,
                    severity: AuditSeverity::Warning,
                    observed_value: imm_resp.composite_immune_risk,
                    message: format!(
                        "Innate immune trigger score {:.4} (limit <= {:.4})",
                        imm_resp.composite_immune_risk, max_imm
                    ),
                }
            }
            SafetyRule::ForbidMycOverexpression(max_myc) => {
                let passed = ctx.oncogene.myc_fold <= *max_myc;
                AuditRuleResult {
                    rule_name: "ForbidMycOverexpression".to_string(),
                    passed,
                    severity: AuditSeverity::CriticalBlock,
                    observed_value: ctx.oncogene.myc_fold,
                    message: format!(
                        "MYC proto-oncogene fold {:.2}x (limit <= {:.2}x; use OSK without MYC)",
                        ctx.oncogene.myc_fold, max_myc
                    ),
                }
            }
            SafetyRule::MaxOffTargetBurden(max_ot) => {
                let passed = ctx.oncogene.offtarget_cleavage_burden <= *max_ot;
                AuditRuleResult {
                    rule_name: "MaxOffTargetBurden".to_string(),
                    passed,
                    severity: AuditSeverity::CriticalBlock,
                    observed_value: ctx.oncogene.offtarget_cleavage_burden,
                    message: format!(
                        "CRISPR off-target cleavage burden {:.4} (limit <= {:.4})",
                        ctx.oncogene.offtarget_cleavage_burden, max_ot
                    ),
                }
            }
        };

        if !res.passed {
            match res.severity {
                AuditSeverity::CriticalBlock => critical_violations += 1,
                AuditSeverity::Warning => warnings += 1,
                AuditSeverity::Advisory => {}
            }
        }
        results.push(res);
    }

    SafetyAuditReport {
        approved: critical_violations == 0,
        results,
        critical_violations,
        warnings,
    }
}
