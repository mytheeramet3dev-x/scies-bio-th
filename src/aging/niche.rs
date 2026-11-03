//! Adult tissue-resident stem cell niche dynamics, `Wnt`/`Notch` asymmetric vs symmetric division,
//! transit-amplifying progenitor expansion, and stem-cell hyperplasia safety caps.
//!
//! References:
//! - Morrison & Spradling (2008) *Cell* 132:598–611 (Stem cells and niches: mechanisms that promote stem cell maintenance).
//! - Conboy et al. (2003) *Science* 302:1575–1577 (Notch-mediated restoration of regenerative potential to aged muscle).
//! - Clevers, Loh & Nusse (2014) *Science* 346:1248012 (An integral program for tissue renewal and regeneration: Wnt signaling and stem cell control).

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Maximum permitted fraction of adult stem cells in a tissue colony ($15\%$) to prevent hyperplasia.
pub const MAX_SAFE_STEM_CELL_FRACTION: f64 = 0.15;

/// Microenvironmental signaling state of an adult tissue-resident stem cell niche.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct StemCellNicheState {
    /// Canonical `Wnt3a` / $\beta$-catenin niche ligand activity ($0.0$ to $1.0$; promotes symmetric self-renewal and TAC expansion).
    pub wnt_signal: f64,
    /// `Delta`/`Jagged`–`Notch1` juxtacrine receptor activity ($0.0$ to $1.0$; preserves stem quiescence and asymmetric division polarity).
    pub notch_signal: f64,
    /// Niche regenerative competence score ($0.0$ to $1.0$).
    pub niche_regenerative_competence: f64,
}

impl Default for StemCellNicheState {
    fn default() -> Self {
        Self::youthful()
    }
}

impl StemCellNicheState {
    /// Creates a youthful balanced stem cell niche state (`Wnt` and `Notch` in homeostatic equilibrium).
    pub fn youthful() -> Self {
        let wnt = 0.78;
        let notch = 0.82;
        Self {
            wnt_signal: wnt,
            notch_signal: notch,
            niche_regenerative_competence: evaluate_niche_competence(wnt, notch),
        }
    }

    /// Estimates an aged stem cell niche state where `Notch` Delta ligand and structured `Wnt` polarity decline.
    pub fn from_somatic_age(age_years: f64) -> Self {
        let age = age_years.clamp(0.0, 120.0);
        let wnt = (0.85 - age * 0.0065).clamp(0.15, 0.90);
        let notch = (0.88 - age * 0.0072).clamp(0.12, 0.92);
        Self {
            wnt_signal: wnt,
            notch_signal: notch,
            niche_regenerative_competence: evaluate_niche_competence(wnt, notch),
        }
    }

    /// Rejuvenates the stem cell niche via pulsed `Wnt`/`Notch` agonist modulation,
    /// returning the updated regenerative competence score.
    pub fn restore_wnt_notch_signaling(&mut self, wnt_boost: f64, notch_boost: f64) -> f64 {
        self.wnt_signal = (self.wnt_signal + wnt_boost.max(0.0)).clamp(0.10, 0.95);
        self.notch_signal = (self.notch_signal + notch_boost.max(0.0)).clamp(0.10, 0.95);
        self.niche_regenerative_competence =
            evaluate_niche_competence(self.wnt_signal, self.notch_signal);
        self.niche_regenerative_competence
    }
}

/// Evaluates overall stem cell niche regenerative competence ($0.0$ to $1.0$) from `Wnt` and `Notch` levels.
///
/// Requires balanced `Wnt` (proliferation drive) and `Notch` (asymmetric polarity & differentiation control).
pub fn evaluate_niche_competence(wnt_signal: f64, notch_signal: f64) -> f64 {
    let w = wnt_signal.clamp(0.0, 1.0);
    let n = notch_signal.clamp(0.0, 1.0);
    // Geometric + harmonic synergy: unbridled Wnt without Notch causes stem cell exhaustion/fibrosis
    let balance_penalty = 1.0 - 0.35 * (w - n).abs();
    ((0.48 * w + 0.52 * n) * balance_penalty).clamp(0.05, 0.99)
}

/// Division mode probabilities for an activated `AdultStemCell`.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct StemDivisionOutcome {
    /// Probability of asymmetric division (`1 AdultStemCell + 1 TransitAmplifyingProgenitor`).
    pub asymmetric_probability: f64,
    /// Probability of symmetric self-renewal (`2 AdultStemCells`), capped to zero if stem fraction $\ge 15\%$.
    pub symmetric_renewal_probability: f64,
    /// Number of somatic daughter cells generated per `TransitAmplifyingProgenitor` amplification round.
    pub progenitor_amplification_factor: usize,
    /// True if symmetric self-renewal was blocked by the $15\%$ hyperplasia safety guardrail.
    pub hyperplasia_cap_enforced: bool,
}

/// Computes `AdultStemCell` asymmetric vs symmetric division probabilities from `Wnt`/`Notch` signaling,
/// strictly enforcing [`MAX_SAFE_STEM_CELL_FRACTION`] ($15\%$) to prevent hyperplasia.
pub fn compute_stem_division_dynamics(
    niche: &StemCellNicheState,
    current_stem_fraction: f64,
) -> StemDivisionOutcome {
    let w = niche.wnt_signal.clamp(0.0, 1.0);
    let n = niche.notch_signal.clamp(0.0, 1.0);

    let hyperplasia_cap_enforced = current_stem_fraction >= MAX_SAFE_STEM_CELL_FRACTION;

    // Symmetric renewal is driven by Wnt when stem pool is below the 15% safety ceiling
    let symmetric_renewal_probability = if hyperplasia_cap_enforced {
        0.0
    } else {
        let headroom = ((MAX_SAFE_STEM_CELL_FRACTION - current_stem_fraction)
            / MAX_SAFE_STEM_CELL_FRACTION)
            .clamp(0.0, 1.0);
        (0.35 * w * headroom).clamp(0.0, 0.40)
    };

    // Asymmetric division is driven by Notch polarity (Numb/Par3 segregation)
    let asymmetric_probability =
        (0.55 * n + 0.35 * w).clamp(0.10, 1.0 - symmetric_renewal_probability);

    let progenitor_amplification_factor = if niche.niche_regenerative_competence >= 0.70 {
        3
    } else if niche.niche_regenerative_competence >= 0.45 {
        2
    } else {
        1
    };

    StemDivisionOutcome {
        asymmetric_probability,
        symmetric_renewal_probability,
        progenitor_amplification_factor,
        hyperplasia_cap_enforced,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stem_niche_aging_and_wnt_notch_restoration() {
        let young = StemCellNicheState::youthful();
        let mut aged = StemCellNicheState::from_somatic_age(72.0);
        assert!(young.niche_regenerative_competence > aged.niche_regenerative_competence + 0.30);

        let restored = aged.restore_wnt_notch_signaling(0.45, 0.50);
        assert!(restored > 0.75);
    }

    #[test]
    fn test_hyperplasia_cap_blocks_symmetric_stem_expansion_above_15_percent() {
        let niche = StemCellNicheState::youthful();
        let below_cap = compute_stem_division_dynamics(&niche, 0.06);
        assert!(!below_cap.hyperplasia_cap_enforced);
        assert!(below_cap.symmetric_renewal_probability > 0.05);

        let at_cap = compute_stem_division_dynamics(&niche, 0.15);
        assert!(at_cap.hyperplasia_cap_enforced);
        assert_eq!(at_cap.symmetric_renewal_probability, 0.0);
    }
}
