//! Pre-flight validation, coordinate integrity checks, and clinical safety compliance audits.
//!
//! Verifies complete synthetic longevity genome design packages prior to clinical report generation:
//! - Validates genomic intervals ($start < end$, positive coordinates).
//! - Validates guide RNA lengths ($20\text{ nt}$) and SpCas9 PAM conformity.
//! - Verifies acceptable oncogenic risk thresholds per design strategy.
//! - Audits compliance against all hard biological safety constraints.

use crate::aging::constraints::{validate_biological_constraints, ConstraintSeverity, ConstraintViolation};
use crate::aging::design::LongevityStrategy;
use crate::aging::report::HighFidelityDesignReport;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Comprehensive outcome of a pre-flight design validation audit.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ValidationOutcome {
    /// Whether the design package satisfies all critical safety criteria and is ready for clinical synthesis.
    pub is_valid: bool,
    /// Fatal integrity or safety errors.
    pub errors: Vec<String>,
    /// Non-fatal biological advisories or optimization suggestions.
    pub warnings: Vec<String>,
    /// Audited biological constraint violations.
    pub constraint_violations: Vec<ConstraintViolation>,
}

/// Validates a high-fidelity longevity design package before final issuance.
pub fn validate_longevity_package(report: &HighFidelityDesignReport) -> ValidationOutcome {
    let mut errors = Vec::new();
    let mut warnings = Vec::new();

    // 1. Verify genomic coordinates and guide sequences
    for edit in &report.proposed_edits {
        if edit.start_pos >= edit.end_pos {
            errors.push(format!(
                "Invalid genomic interval for {}: start ({}) >= end ({})",
                edit.target_gene, edit.start_pos, edit.end_pos
            ));
        }

        if edit.target_sequence.sequence.is_empty() {
            errors.push(format!("Empty target sequence context for {}", edit.target_gene));
        }

        if let Some(ref guide) = edit.recommended_sgrna {
            if guide.sequence.sequence.len() != 20 {
                errors.push(format!(
                    "Invalid guide length for {} (expected 20nt, found {}nt)",
                    edit.target_gene, guide.sequence.sequence.len()
                ));
            }
            if guide.pam_motif.len() < 3 || !guide.pam_motif.ends_with("GG") {
                warnings.push(format!(
                    "Guide for {} has non-canonical PAM motif '{}'",
                    edit.target_gene, guide.pam_motif
                ));
            }
        } else {
            warnings.push(format!("No SpCas9 guide candidate discovered for {}", edit.target_gene));
        }
    }

    // 2. Verify oncogenic risk thresholds per strategy
    let risk_limit = match report.strategy {
        LongevityStrategy::Conservative => 0.05,
        LongevityStrategy::Balanced => 0.15,
        LongevityStrategy::Aggressive => 0.25,
    };

    if report.oncogenic_risk.net_risk_index > risk_limit {
        errors.push(format!(
            "Oncogenic risk index ({:.3}) exceeds safety ceiling ({:.3}) for {:?} strategy",
            report.oncogenic_risk.net_risk_index, risk_limit, report.strategy
        ));
    }

    // 3. Verify hard biological safety constraints using remediated expression level
    let telomerase_level = match report.strategy {
        LongevityStrategy::Conservative => 0.0,
        LongevityStrategy::Balanced => 0.80,
        LongevityStrategy::Aggressive => {
            if report.oncogenic_risk.safeguards.super_p53_active
                && report.oncogenic_risk.safeguards.inducible_suicide_switch
                && report.oncogenic_risk.safeguards.cell_cycle_restriction
            {
                0.95
            } else {
                0.80
            }
        }
    };

    let constraint_violations = validate_biological_constraints(
        &report.proposed_edits,
        report.strategy,
        telomerase_level,
        &report.oncogenic_risk.safeguards,
    );

    for violation in &constraint_violations {
        if violation.severity == ConstraintSeverity::Critical {
            errors.push(format!("[CRITICAL] {}: {}", violation.affected_gene, violation.message));
        } else {
            warnings.push(format!("[{}] {}: {}", format!("{:?}", violation.severity).to_uppercase(), violation.affected_gene, violation.message));
        }
    }

    let is_valid = errors.is_empty();

    ValidationOutcome {
        is_valid,
        errors,
        warnings,
        constraint_violations,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aging::design::design_high_fidelity_longevity_package;
    use crate::aging::fidelity::HighFidelityAgeParams;

    #[test]
    fn test_validate_healthy_balanced_package() {
        let subject = HighFidelityAgeParams {
            chronological_age: 50.0,
            cpg_methylation_fraction: 0.40,
            telomere_length_bp: 7_000.0,
            gamma_h2ax_foci: 1.2,
            sasp_burden: 0.25,
        };

        let report = design_high_fidelity_longevity_package(&subject, LongevityStrategy::Balanced, None);
        let outcome = validate_longevity_package(&report);

        assert!(outcome.is_valid);
        assert!(outcome.errors.is_empty());
    }
}
