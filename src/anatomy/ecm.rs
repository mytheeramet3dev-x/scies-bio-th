//! Extracellular Matrix (`ECM`) biomechanics, collagen/elastin architecture,
//! Advanced Glycation End-product (`AGE`, glucosepane) cross-linking, Young's elastic modulus ($E$, $\text{kPa}$),
//! `YAP`/`TAZ` mechanotransduction stress, and Ogston porous-gel steric hindrance.
//!
//! References:
//! - Humphrey et al. (2014) *Nat Rev Mol Cell Biol* 15:802–812 (Mechanotransduction and extracellular matrix homeostasis).
//! - Sell & Monnier (2012) *Gerontology* 58:227–237 (Molecular basis of arterial stiffening: role of glucosepane AGE cross-links).
//! - Dupont et al. (2011) *Nature* 474:179–183 (Role of YAP/TAZ in mechanotransduction).
//! - Stylianopoulos & Jain (2013) *PNAS* 110:18632–18637 (Interstitial transport in fibrotic extracellular matrices).

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Biomechanical and biochemical state of an organ's Extracellular Matrix (ECM).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ExtracellularMatrixState {
    /// Youthful reference Young's elastic modulus ($E_0$) in $\text{kPa}$ for this tissue type.
    pub youthful_modulus_kpa: f64,
    /// Current effective Young's elastic modulus ($E$) in $\text{kPa}$ (increases with AGE cross-links and fibrosis).
    pub youngs_modulus_kpa: f64,
    /// Interstitial fibrillar Collagen I/III volume fraction ($0.0$ to $1.0$; elevated in fibrosis).
    pub collagen_volume_fraction: f64,
    /// Functional elastin fiber integrity ($0.0$ to $1.0$; degrades irreversibly with age unless restored).
    pub elastin_integrity: f64,
    /// Glucosepane / pentosidine Advanced Glycation End-product (AGE) cross-link burden ($0.0$ to $1.0$).
    pub ages_crosslink_burden: f64,
    /// Interstitial glycosaminoglycan / hyaluronic acid hydration ($0.0$ to $1.0$).
    pub hyaluronic_hydration: f64,
}

impl ExtracellularMatrixState {
    /// Creates a youthful ECM state calibrated to an organ's canonical `youthful_modulus_kpa`.
    pub fn youthful(youthful_modulus_kpa: f64) -> Self {
        let e0 = youthful_modulus_kpa.max(0.5);
        Self {
            youthful_modulus_kpa: e0,
            youngs_modulus_kpa: e0,
            collagen_volume_fraction: 0.14,
            elastin_integrity: 0.94,
            ages_crosslink_burden: 0.03,
            hyaluronic_hydration: 0.92,
        }
    }

    /// Estimates age-associated ECM stiffening, glucosepane AGE cross-linking, elastin fragmentation,
    /// and interstitial collagen fibrosis from chronological `age_years`.
    pub fn from_somatic_age(youthful_modulus_kpa: f64, age_years: f64) -> Self {
        let e0 = youthful_modulus_kpa.max(0.5);
        let age = age_years.clamp(0.0, 120.0);
        let excess_age = (age - 20.0).max(0.0);

        let ages = (0.03 + excess_age * 0.0085).clamp(0.02, 0.92);
        let collagen = (0.14 + excess_age * 0.0042).clamp(0.10, 0.65);
        let elastin = (0.95 - excess_age * 0.0068).clamp(0.18, 0.96);
        let hydration = (0.93 - excess_age * 0.0055).clamp(0.25, 0.95);

        // Non-linear strain-stiffening from glucosepane cross-linking and elastin loss
        let stiffening_multiplier =
            1.0 + 1.85 * ages + 1.40 * (collagen - 0.14).max(0.0) + 0.85 * (1.0 - elastin);
        let current_e = e0 * stiffening_multiplier;

        Self {
            youthful_modulus_kpa: e0,
            youngs_modulus_kpa: current_e,
            collagen_volume_fraction: collagen,
            elastin_integrity: elastin,
            ages_crosslink_burden: ages,
            hyaluronic_hydration: hydration,
        }
    }

    /// Recomputes `youngs_modulus_kpa` from current biochemical matrix components.
    pub fn update_modulus(&mut self) {
        let stiffening_multiplier = 1.0
            + 1.85 * self.ages_crosslink_burden
            + 1.40 * (self.collagen_volume_fraction - 0.14).max(0.0)
            + 0.85 * (1.0 - self.elastin_integrity);
        self.youngs_modulus_kpa = self.youthful_modulus_kpa * stiffening_multiplier.max(0.8);
    }

    /// Computes normalized `YAP`/`TAZ` nuclear mechanotransduction stress ($0.0$ to $1.0$)
    /// triggered when `youngs_modulus_kpa` exceeds `youthful_modulus_kpa` (Dupont et al. 2011).
    ///
    /// High `YAP`/`TAZ` mechanical stress accelerates fibrotic senescence and impairs stem cell niche quiescence.
    pub fn yap_taz_mechanotransduction_stress(&self) -> f64 {
        let ratio = (self.youngs_modulus_kpa / self.youthful_modulus_kpa.max(0.5)).max(1.0);
        let excess_ratio = ratio - 1.0;
        // Hill-like mechanosensing curve saturating at ~3x youthful stiffness
        (excess_ratio.powi(2) / (0.65 + excess_ratio.powi(2))).clamp(0.0, 1.0)
    }

    /// Computes the Ogston porous-gel interstitial diffusion factor ($0.05$ to $1.0$)
    /// for a therapeutic payload with molecular weight `molecular_weight_da` (Daltons).
    ///
    /// Dense collagen and glucosepane cross-links shrink interstitial pore size, hindering
    /// transport of large macromolecules (`LNP`, `AAV`) much more strongly than small molecules.
    pub fn interstitial_diffusion_transmittance(&self, molecular_weight_da: f64) -> f64 {
        let mw = molecular_weight_da.max(50.0);
        // Approximate hydrodynamic radius r_h (nm) from MW (Da): r_h ~ 0.066 * MW^(1/3)
        let r_h_nm = 0.066 * mw.cbrt();
        // Effective ECM pore fiber spacing r_pore (nm) shrinks with collagen fibrosis & AGE cross-linking
        let effective_pore_nm = (48.0
            * (1.0 - 0.72 * self.collagen_volume_fraction)
            * (1.0 - 0.65 * self.ages_crosslink_burden)
            * (0.45 + 0.55 * self.hyaluronic_hydration))
            .clamp(6.0, 55.0);

        // Ogston hindrance exp(-(r_h / r_pore)^1.35)
        let ratio = r_h_nm / effective_pore_nm;
        (-ratio.powf(1.35)).exp().clamp(0.05, 1.0)
    }

    /// Executes an ECM remodeling and de-stiffening intervention:
    /// - `ages_breaker_activity`: Glucosepane cross-link lyase / deglycation enzyme activity ($0.0\text{--}1.5$).
    /// - `antifibrotic_mmp_remodeling`: Balanced `MMP`/`TIMP` & `TGF-beta` inhibition clearing excess fibrotic collagen ($0.0\text{--}1.5$).
    /// - `elastin_hyaluronan_boost`: Tropoelastin / `LOXL1` & `HAS2` hyaluronan restoration ($0.0\text{--}1.5$).
    pub fn apply_matrix_rejuvenation(
        &mut self,
        ages_breaker_activity: f64,
        antifibrotic_mmp_remodeling: f64,
        elastin_hyaluronan_boost: f64,
    ) {
        let ab = ages_breaker_activity.clamp(0.0, 1.5);
        let mmp = antifibrotic_mmp_remodeling.clamp(0.0, 1.5);
        let eh = elastin_hyaluronan_boost.clamp(0.0, 1.5);

        self.ages_crosslink_burden =
            (self.ages_crosslink_burden * (1.0 - 0.72 * ab)).clamp(0.02, 0.95);
        let excess_collagen = (self.collagen_volume_fraction - 0.14).max(0.0);
        self.collagen_volume_fraction =
            (0.14 + excess_collagen * (1.0 - 0.70 * mmp)).clamp(0.12, 0.65);
        self.elastin_integrity = (self.elastin_integrity + 0.42 * eh).clamp(0.18, 0.95);
        self.hyaluronic_hydration = (self.hyaluronic_hydration + 0.38 * eh).clamp(0.25, 0.95);

        self.update_modulus();
    }

    /// Advances chronological matrix aging by `delta_years`, accumulating glucosepane `AGE`
    /// cross-links, interstitial collagen fibrosis, elastin fragmentation, and dehydration.
    pub fn apply_chronological_aging(&mut self, delta_years: f64) {
        let dt = delta_years.max(0.0);
        self.ages_crosslink_burden = (self.ages_crosslink_burden + 0.0085 * dt).clamp(0.02, 0.94);
        self.collagen_volume_fraction =
            (self.collagen_volume_fraction + 0.0042 * dt).clamp(0.10, 0.65);
        self.elastin_integrity = (self.elastin_integrity - 0.0068 * dt).clamp(0.16, 0.96);
        self.hyaluronic_hydration = (self.hyaluronic_hydration - 0.0055 * dt).clamp(0.22, 0.95);
        self.update_modulus();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ecm_aging_stiffening_and_yap_taz_stress() {
        let young = ExtracellularMatrixState::youthful(4.0);
        let aged = ExtracellularMatrixState::from_somatic_age(4.0, 75.0);

        assert!(aged.youngs_modulus_kpa > young.youngs_modulus_kpa * 1.8);
        assert!(aged.yap_taz_mechanotransduction_stress() > 0.45);
        assert!(young.yap_taz_mechanotransduction_stress() < 0.01);
    }

    #[test]
    fn test_ecm_matrix_rejuvenation_restores_compliance_and_diffusion() {
        let mut aged = ExtracellularMatrixState::from_somatic_age(12.0, 75.0);
        let before_aav = aged.interstitial_diffusion_transmittance(3_700_000.0);
        let before_mod = aged.youngs_modulus_kpa;

        aged.apply_matrix_rejuvenation(1.0, 1.0, 1.0);
        assert!(aged.youngs_modulus_kpa < before_mod * 0.68);
        assert!(aged.interstitial_diffusion_transmittance(3_700_000.0) > before_aav + 0.12);
    }
}
