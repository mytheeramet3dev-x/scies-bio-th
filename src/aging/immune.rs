//! Immunosenescence, thymic involution (`FOXN1` epithelial atrophy), T-cell receptor (`TCR`) repertoire
//! diversity, and endogenous innate immunosurveillance (`NK` cell perforin/granzyme cytotoxicity &
//! macrophage efferocytosis of senescent and pre-malignant cells).
//!
//! Note: This module (`aging::immune`) models adaptive/innate immune system aging and endogenous
//! senescent cell clearance, distinct from `safety::immune` which models acute TLR/RIG-I vector immunogenicity.
//!
//! References:
//! - Ovadya et al. (2018) *Nat Commun* 9:5435 (Impaired immune surveillance accelerates accumulation of senescent cells and aging).
//! - Fahy et al. (2019) *Aging Cell* 18:e13028 (Reversal of epigenetic aging and immunosenescent trends in humans — TRIIM trial).
//! - Palmer et al. (2018) *PNAS* 115:1883–1888 (Thymic involution and rising disease incidence with age).

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Adaptive and innate immune competence state governing immunosenescence and endogenous senescent-cell clearance.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ImmunosenescenceState {
    /// Thymic cortical/medullary epithelial (`FOXN1+`) output of recent thymic emigrant naive T cells ($0.0$ to $1.0$).
    pub thymic_output_naive_t: f64,
    /// Normalized Shannon diversity index of the peripheral $\alpha\beta$ T-cell receptor (`TCR`) repertoire ($0.0$ to $1.0$).
    pub tcr_diversity_index: f64,
    /// Natural Killer (`NK`) cell `NKG2D`/`DNAM-1` perforin/granzyme cytotoxicity against `MICA`/`ULBP2+` senescent & mutated cells ($0.0$ to $1.0$).
    pub nk_cytotoxicity: f64,
    /// Macrophage `MerTK`/`TIM-4` phagocytic efferocytosis capacity for clearing senescent cells and apoptotic debris ($0.0$ to $1.0$).
    pub macrophage_efferocytosis: f64,
}

impl Default for ImmunosenescenceState {
    fn default() -> Self {
        Self::youthful()
    }
}

impl ImmunosenescenceState {
    /// Creates a youthful immunocompetent state with active thymopoiesis and robust NK/macrophage surveillance.
    pub fn youthful() -> Self {
        Self {
            thymic_output_naive_t: 0.92,
            tcr_diversity_index: 0.95,
            nk_cytotoxicity: 0.90,
            macrophage_efferocytosis: 0.92,
        }
    }

    /// Estimates age-associated thymic involution and declining innate/adaptive immunosurveillance.
    ///
    /// Reflects exponential thymic involution after puberty ($\sim 3\%$ per year; Palmer et al. 2018)
    /// and age-related impairment of `NKG2D`-mediated `NK` cytotoxicity and macrophage efferocytosis (Ovadya et al. 2018).
    pub fn from_somatic_age(age_years: f64) -> Self {
        let age = age_years.clamp(0.0, 120.0);
        let thymic_output = (0.98 * (-0.032 * (age - 15.0).max(0.0)).exp()).clamp(0.04, 0.98);
        let tcr_diversity = (0.98 - 0.0075 * (age - 18.0).max(0.0)).clamp(0.15, 0.98);
        let nk_cyto = (0.94 - 0.0072 * age).clamp(0.18, 0.95);
        let mac_eff = (0.95 - 0.0076 * age).clamp(0.15, 0.95);

        Self {
            thymic_output_naive_t: thymic_output,
            tcr_diversity_index: tcr_diversity,
            nk_cytotoxicity: nk_cyto,
            macrophage_efferocytosis: mac_eff,
        }
    }

    /// Composite endogenous immune surveillance efficiency ($0.0$ to $1.0$) combining `NK` cytotoxicity,
    /// macrophage efferocytosis, and adaptive `TCR` repertoire breadth.
    pub fn surveillance_efficiency(&self) -> f64 {
        (0.40 * self.nk_cytotoxicity
            + 0.35 * self.macrophage_efferocytosis
            + 0.15 * self.tcr_diversity_index
            + 0.10 * self.thymic_output_naive_t)
            .clamp(0.0, 1.0)
    }

    /// Probability ($0.0$ to $0.98$) that endogenous `NK` cells and macrophages detect and eliminate
    /// a senescent (`MICA`/`ULBP2+`) or mutated/dysplastic somatic cell without exogenous senolytic drugs.
    pub fn endogenous_senescent_clearance_probability(&self) -> f64 {
        let innate_synergy = 0.55 * self.nk_cytotoxicity + 0.45 * self.macrophage_efferocytosis;
        let adaptive_support = 0.50 * self.tcr_diversity_index + 0.50 * self.thymic_output_naive_t;
        (0.75 * innate_synergy + 0.25 * adaptive_support).clamp(0.05, 0.98)
    }

    /// Simulates thymic regeneration (`FOXN1` / `KGF` / `GH`–`Metformin` TRIIM protocol; Fahy et al. 2019)
    /// coupled with `IL-15` / `CD38`-inhibition restoration of `NK` cytotoxicity and macrophage efferocytosis.
    pub fn apply_immune_rejuvenation(
        &mut self,
        foxn1_thymic_regeneration: f64,
        nk_macrophage_activation: f64,
    ) -> f64 {
        let before = self.surveillance_efficiency();
        let thymic_boost = foxn1_thymic_regeneration.clamp(0.0, 1.5);
        let innate_boost = nk_macrophage_activation.clamp(0.0, 1.5);

        self.thymic_output_naive_t =
            (self.thymic_output_naive_t + 0.58 * thymic_boost).clamp(0.04, 0.95);
        self.tcr_diversity_index =
            (self.tcr_diversity_index + 0.42 * thymic_boost).clamp(0.15, 0.96);
        self.nk_cytotoxicity = (self.nk_cytotoxicity + 0.48 * innate_boost).clamp(0.18, 0.96);
        self.macrophage_efferocytosis =
            (self.macrophage_efferocytosis + 0.50 * innate_boost).clamp(0.15, 0.96);

        (self.surveillance_efficiency() - before).max(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_thymic_involution_and_immunosenescence_progression() {
        let young = ImmunosenescenceState::youthful();
        let aged = ImmunosenescenceState::from_somatic_age(72.0);

        assert!(aged.thymic_output_naive_t < 0.20);
        assert!(young.tcr_diversity_index > aged.tcr_diversity_index + 0.30);
        assert!(young.surveillance_efficiency() > aged.surveillance_efficiency() + 0.40);
    }

    #[test]
    fn test_thymic_and_nk_macrophage_rejuvenation_restores_endogenous_clearance() {
        let mut aged = ImmunosenescenceState::from_somatic_age(72.0);
        let before_prob = aged.endogenous_senescent_clearance_probability();
        let gain = aged.apply_immune_rejuvenation(1.0, 1.0);

        assert!(gain > 0.35);
        assert!(aged.thymic_output_naive_t > 0.65);
        assert!(aged.nk_cytotoxicity > 0.80);
        assert!(aged.macrophage_efferocytosis > 0.80);
        assert!(aged.endogenous_senescent_clearance_probability() > before_prob + 0.35);
    }
}
