//! Telomerase ribonucleoprotein complex (TERT/TERC) activity and oncogenic risk modeling.
//!
//! Telomerase synthesizes telomeric `TTAGGG` repeats de novo onto chromosome 3' overhangs.
//! While somatic repression acts as an evolutionary tumor suppressor, controlled synthetic reactivation
//! provides indefinite cellular replicative capacity when paired with genomic safety locks.

use crate::aging::telomere::TelomereState;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Catalytic and functional activity profile of the Telomerase complex (TERT + TERC).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TelomeraseActivity {
    /// Relative expression level (0.0 = completely silenced somatic, 1.0 = embryonic stem cell).
    pub expression_level: f64,
    /// Repeat addition processivity (average nucleotides added per primer engagement, ~40-60 nt).
    pub processivity: f64,
    /// Whether the enzyme is enzymatically competent and uninhibited.
    pub active: bool,
    /// Presence of synthetic transcription safety lock (e.g. inducible doxycycline/tamoxifen or cell-cycle restricted).
    pub safety_lock_present: bool,
}

impl TelomeraseActivity {
    /// Repressed baseline corresponding to normal adult somatic cells.
    pub fn somatic_silenced() -> Self {
        Self {
            expression_level: 0.0,
            processivity: 0.0,
            active: false,
            safety_lock_present: true,
        }
    }

    /// Youthful, regulated physiological telomerase (e.g. germline/stem cell).
    pub fn physiological_stem_cell() -> Self {
        Self {
            expression_level: 0.75,
            processivity: 48.0,
            active: true,
            safety_lock_present: true,
        }
    }

    /// Synthetic engineered telomerase construct with enhanced processivity.
    pub fn synthetic_enhanced(safety_locked: bool) -> Self {
        Self {
            expression_level: 0.85,
            processivity: 65.0,
            active: true,
            safety_lock_present: safety_locked,
        }
    }
}

/// Applies active telomerase elongation cycles to maintain or expand telomeres.
///
/// Returns the total net base pairs added across all chromosomes.
pub fn apply_telomerase(
    state: &mut TelomereState,
    activity: &TelomeraseActivity,
    cycles: u32,
) -> f64 {
    if !activity.active || activity.expression_level <= 0.0 {
        return 0.0;
    }

    // Net extension per cycle = expression * processivity * efficiency factor
    let extension_per_cycle = activity.expression_level * activity.processivity;
    let net_extension = extension_per_cycle * cycles as f64;

    // Cap physiological maximum telomere length at 18,000 bp (prevent abnormal hyper-elongation)
    let max_length = 18_000.0;
    let new_length = (state.current_length_bp + net_extension).min(max_length);
    let actual_added = new_length - state.current_length_bp;

    state.current_length_bp = new_length;
    if state.current_length_bp > state.critical_length_bp {
        state.is_critically_short = false;
        state.shelterin_integrity = (state.shelterin_integrity + 0.15).min(0.98);
    }

    actual_added
}

/// Evaluates oncogenic neoplastic transformation risk associated with telomerase activity.
///
/// Unregulated TERT expression without tumor suppressor surveillance (p53) incurs a high risk score.
pub fn cancer_risk_factor(activity: &TelomeraseActivity, p53_intact: bool) -> f64 {
    if !activity.active || activity.expression_level <= 0.05 {
        return 0.01; // Negligible somatic baseline
    }

    let raw_risk = activity.expression_level * 0.70;
    let p53_damping = if p53_intact { 0.25 } else { 1.50 };
    let safety_damping = if activity.safety_lock_present { 0.35 } else { 1.20 };

    (raw_risk * p53_damping * safety_damping).clamp(0.01, 0.99)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_apply_telomerase_elongation() {
        let mut state = TelomereState {
            current_length_bp: 4_000.0,
            shortening_rate_bp: 60.0,
            critical_length_bp: 3_800.0,
            shelterin_integrity: 0.50,
            divisions_elapsed: 80,
            is_critically_short: false,
        };

        let activity = TelomeraseActivity::synthetic_enhanced(true);
        // Extension per cycle: 0.85 * 65.0 = 55.25 bp
        // 10 cycles -> ~552.5 bp added
        let added = apply_telomerase(&mut state, &activity, 10);
        assert!(added > 500.0);
        assert_eq!(state.current_length_bp, 4_000.0 + added);
    }

    #[test]
    fn test_cancer_risk_mitigation() {
        let unlocked = TelomeraseActivity {
            expression_level: 0.95,
            processivity: 60.0,
            active: true,
            safety_lock_present: false,
        };
        let locked = TelomeraseActivity {
            expression_level: 0.95,
            processivity: 60.0,
            active: true,
            safety_lock_present: true,
        };

        let risk_unlocked = cancer_risk_factor(&unlocked, false);
        let risk_locked_p53 = cancer_risk_factor(&locked, true);

        // Safety lock + intact p53 should dramatically lower oncogenic risk
        assert!(risk_locked_p53 < risk_unlocked);
        assert!(risk_locked_p53 < 0.15);
    }
}
