//! Cellular proteostasis, cytotoxic protein aggregate accumulation (`amyloid`, `tau`, `lipofuscin`),
//! macroautophagy/lysosomal clearance (`mTORC1` inhibition + `TFEB`), and molecular chaperone refolding (`HSP70`/`HSP90`).
//!
//! References:
//! - Hartl, Bracher & Hayer-Hartl (2011) *Nature* 475:324–332 (Molecular chaperones in protein folding and proteostasis).
//! - Settembre et al. (2011) *Science* 332:1429–1433 (TFEB links autophagy to lysosomal biogenesis).
//! - Laplante & Sabatini (2012) *Cell* 149:274–293 (mTOR signaling in growth control and autophagy regulation).
//! - Kaushik & Cuervo (2015) *Nat Med* 21:1406–1415 (Proteostasis and aging).

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Intracellular and lysosomal cytotoxic protein aggregate state of a cell.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ProteostasisState {
    /// Cross-$\beta$ sheet misfolded amyloid oligomer/fibril burden ($0.0$ to $1.0$).
    pub amyloid_burden: f64,
    /// Hyperphosphorylated cytoskeletal/neurofibrillary tau aggregate burden ($0.0$ to $1.0$).
    pub tau_burden: f64,
    /// Intralysosomal cross-linked oxidized lipid-protein lipofuscin granule burden ($0.0$ to $1.0$).
    pub lipofuscin_burden: f64,
}

impl Default for ProteostasisState {
    fn default() -> Self {
        Self::youthful()
    }
}

impl ProteostasisState {
    /// Creates a youthful baseline proteostasis state with minimal aggregate load.
    pub fn youthful() -> Self {
        Self {
            amyloid_burden: 0.015,
            tau_burden: 0.012,
            lipofuscin_burden: 0.018,
        }
    }

    /// Estimates somatic proteostasis aggregate accumulation as a function of chronological age.
    ///
    /// Note: Nuclear OSK epigenetic reprogramming does **not** directly clear pre-existing
    /// insoluble cytoplasmic/lysosomal protein aggregates.
    pub fn from_somatic_age(age_years: f64) -> Self {
        let age = age_years.clamp(0.0, 120.0);
        let amyloid = (0.012 + 0.000055 * age.powf(1.92)).clamp(0.01, 0.85);
        let tau = (0.010 + 0.000048 * age.powf(1.90)).clamp(0.01, 0.80);
        // Lipofuscin accumulates inexorably in post-mitotic and senescent lysosomes
        let lipofuscin = (0.015 + 0.000072 * age.powf(1.95)).clamp(0.01, 0.90);

        Self {
            amyloid_burden: amyloid,
            tau_burden: tau,
            lipofuscin_burden: lipofuscin,
        }
    }

    /// Composite proteostasis aggregate burden ($0.0$ to $1.0$).
    pub fn total_burden(&self) -> f64 {
        (0.35 * self.amyloid_burden + 0.30 * self.tau_burden + 0.35 * self.lipofuscin_burden)
            .clamp(0.0, 1.0)
    }

    /// Applies combined `HSP70`/`HSP90` chaperone refolding and `mTORC1`/`TFEB` autophagy clearance in-place.
    ///
    /// Returns the total aggregate burden reduction achieved.
    pub fn apply_proteostasis_clearance(
        &mut self,
        mtorc1_inhibition: f64,
        tfeb_activity: f64,
        hsp70: f64,
        hsp90: f64,
    ) -> f64 {
        let before = self.total_burden();
        let auto_eff = autophagy_clearance(mtorc1_inhibition, tfeb_activity);
        let chap_amyloid = chaperone_refolding(hsp70, hsp90, self.amyloid_burden);
        let chap_tau = chaperone_refolding(hsp70, hsp90, self.tau_burden);

        // Amyloid and tau are cleared synergistically by disaggregase chaperones + aggrephagy
        let amyloid_cleared = (1.0 - (1.0 - chap_amyloid * 0.55) * (1.0 - auto_eff * 0.65)).clamp(0.0, 0.92);
        let tau_cleared = (1.0 - (1.0 - chap_tau * 0.60) * (1.0 - auto_eff * 0.65)).clamp(0.0, 0.92);
        // Lipofuscin is cross-linked and chaperone-resistant; cleared primarily via TFEB lysosomal exocytosis & macroautophagy
        let lipofuscin_cleared = (auto_eff * 0.72).clamp(0.0, 0.85);

        self.amyloid_burden = (self.amyloid_burden * (1.0 - amyloid_cleared)).clamp(0.005, 1.0);
        self.tau_burden = (self.tau_burden * (1.0 - tau_cleared)).clamp(0.005, 1.0);
        self.lipofuscin_burden = (self.lipofuscin_burden * (1.0 - lipofuscin_cleared)).clamp(0.008, 1.0);

        (before - self.total_burden()).max(0.0)
    }
}

/// Computes macroautophagy and lysosomal exocytosis clearance efficiency ($0.0$ to $0.95$)
/// driven by `mTORC1` inhibition (derepressing `ULK1`/`ATG13`) and nuclear `TFEB` transcriptional activation.
pub fn autophagy_clearance(mtorc1_inhibition: f64, tfeb_activity: f64) -> f64 {
    let mtor_inh = mtorc1_inhibition.clamp(0.0, 1.0);
    let tfeb = tfeb_activity.clamp(0.0, 2.0);

    // ULK1 initiation requires mTORC1 inhibition; autophagosome-lysosome fusion requires TFEB-driven V-ATPase/LAMP2A
    let ulk1_drive = mtor_inh.powf(0.85);
    let lysosomal_capacity = (tfeb / (0.35 + tfeb)).clamp(0.0, 1.0);

    let synergy = 0.45 * ulk1_drive + 0.35 * lysosomal_capacity + 0.30 * (ulk1_drive * lysosomal_capacity);
    synergy.clamp(0.0, 0.95)
}

/// Computes ATP-dependent disaggregation and refolding efficiency ($0.0$ to $0.90$) by the
/// `HSP70` (`HSPA1A`/Hsp40/Hsp110 triad) and `HSP90` chaperone machinery, accounting for
/// Michaelis-Menten substrate saturation under high aggregate burden.
pub fn chaperone_refolding(hsp70: f64, hsp90: f64, substrate_burden: f64) -> f64 {
    let h70 = hsp70.clamp(0.0, 5.0);
    let h90 = hsp90.clamp(0.0, 5.0);
    let burden = substrate_burden.clamp(0.0, 1.0);

    if h70 <= 0.01 && h90 <= 0.01 {
        return 0.0;
    }

    // Cooperative foldase/disaggregase capacity (HSP70 provides primary disaggregation, HSP90 stabilizes intermediates)
    let effective_chaperone = 0.65 * h70 + 0.35 * h90;
    let foldase_activation = effective_chaperone / (0.80 + effective_chaperone);

    // High aggregate burden saturates chaperone holdase sites (proteostasis overload)
    let saturation_penalty = 1.0 / (1.0 + (burden / 0.65).powi(2));

    (foldase_activation * saturation_penalty * 1.15).clamp(0.0, 0.90)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_autophagy_clearance_synergy() {
        let basal = autophagy_clearance(0.05, 0.10);
        let induced = autophagy_clearance(0.85, 1.20);
        assert!(induced > basal * 3.0);
        assert!(induced > 0.70);
    }

    #[test]
    fn test_chaperone_refolding_and_substrate_saturation() {
        // Low aggregate load refolds efficiently under elevated HSP70/HSP90
        let low_load = chaperone_refolding(2.5, 2.0, 0.15);
        // Severe aggregate overload partially saturates chaperones
        let overloaded = chaperone_refolding(2.5, 2.0, 0.85);
        assert!(low_load > overloaded);
        assert!(low_load > 0.65);
        assert_eq!(chaperone_refolding(0.0, 0.0, 0.5), 0.0);
    }

    #[test]
    fn test_proteostasis_state_aging_and_clearance() {
        let young = ProteostasisState::from_somatic_age(22.0);
        let mut aged = ProteostasisState::from_somatic_age(72.0);
        assert!(aged.total_burden() > young.total_burden() * 4.0);

        let before = aged.total_burden();
        let delta = aged.apply_proteostasis_clearance(0.85, 1.25, 2.2, 1.8);
        assert!(delta > 0.15);
        assert!(aged.total_burden() < before * 0.40);
    }
}
