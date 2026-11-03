//! Microvascular perfusion, age-associated capillary rarefaction, endothelial `eNOS` activity,
//! tissue oxygenation ($p\text{O}_2$), and molecular-weight-dependent drug/vector extravasation bioavailability.
//!
//! References:
//! - Das et al. (2018) *Cell* 173:74–89 (Impairment of an endothelial NAD+–H2S signaling network is a reversible cause of vascular aging).
//! - Grunewald et al. (2021) *Science* 373:eabc8479 (Counteracting age-related VEGF signaling insufficiency promotes healthy aging and extends life span).
//! - Jain & Stylianopoulos (2010) *Nat Rev Clin Oncol* 7:653–664 (Delivering nanomedicine to solid tissues: microvascular permeability and Stokes-Einstein hindrance).

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Canonical molecular weight of small-molecule senolytic Navitoclax (ABT-263) in Daltons ($\approx 974.6\text{ Da}$).
pub const MW_SMALL_MOLECULE_NAVITOCLAX_DA: f64 = 974.6;

/// Canonical molecular weight of mRNA Lipid Nanoparticle (LNP) payload in Daltons ($\approx 8.0 \times 10^5\text{ Da}$).
pub const MW_MRNA_LNP_VECTOR_DA: f64 = 800_000.0;

/// Canonical molecular weight of an AAV9 capsid vector in Daltons ($\approx 3.7 \times 10^6\text{ Da}$).
pub const MW_AAV_CAPSID_VECTOR_DA: f64 = 3_700_000.0;

/// Microvascular perfusion and endothelial oxygenation state of a somatic tissue.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct VascularPerfusionState {
    /// Functional microvascular capillary density fraction ($0.0$ to $1.0$).
    pub capillary_density: f64,
    /// Interstitial tissue partial pressure of oxygen ($p\text{O}_2$) in $\text{mmHg}$ (normoxic somatic tissue $\approx 42\text{ mmHg}$).
    pub po2_mmhg: f64,
    /// Endothelial Nitric Oxide Synthase (`eNOS` / `NOS3`) vasodilatory catalytic activity ($0.0$ to $1.0$).
    pub enos_activity: f64,
}

impl Default for VascularPerfusionState {
    fn default() -> Self {
        Self::youthful()
    }
}

impl VascularPerfusionState {
    /// Creates a youthful normoxic microvascular state (`capillary_density = 0.94`, `po2_mmhg = 43.0`, `enos_activity = 0.92`).
    pub fn youthful() -> Self {
        Self {
            capillary_density: 0.94,
            po2_mmhg: 43.0,
            enos_activity: 0.92,
        }
    }

    /// Estimates age-dependent capillary rarefaction, endothelial `eNOS` decline, and tissue hypoxia.
    pub fn from_somatic_age(age_years: f64) -> Self {
        let age = age_years.clamp(0.0, 120.0);
        let cap_density = (0.98 - age * 0.0058).clamp(0.25, 0.96);
        let enos = (0.96 - age * 0.0062).clamp(0.20, 0.95);
        // Krogh cylinder oxygenation depends on both capillary density and eNOS-mediated flow
        let po2 = (14.0 + 31.0 * (0.65 * cap_density + 0.35 * enos)).clamp(15.0, 46.0);

        Self {
            capillary_density: cap_density,
            po2_mmhg: po2,
            enos_activity: enos,
        }
    }

    /// Rejuvenates microvascular capillary density and endothelial `eNOS` activity
    /// (e.g., via pulsed `VEGF-A` / `SIRT1`–$\text{NAD}^+$ endothelial restoration; Das et al. 2018; Grunewald et al. 2021).
    pub fn restore_microvascular_angiogenesis(&mut self, vegf_sirt1_boost: f64) {
        let boost = vegf_sirt1_boost.clamp(0.0, 1.5);
        self.capillary_density = (self.capillary_density + 0.38 * boost).clamp(0.20, 0.96);
        self.enos_activity = (self.enos_activity + 0.42 * boost).clamp(0.20, 0.96);
        self.po2_mmhg =
            (14.0 + 31.0 * (0.65 * self.capillary_density + 0.35 * self.enos_activity)).clamp(15.0, 46.0);
    }
}

/// Computes tissue extravasation and microvascular drug/vector bioavailability ($0.05$ to $0.99$)
/// as a function of `capillary_density` ($0.0\text{--}1.0$) and `molecular_weight_da` (Daltons).
///
/// Models Krogh-cylinder vascular surface area scaling coupled with Stokes-Einstein trans-endothelial
/// diffusion resistance for macromolecules (e.g., small molecules $< 1\text{ kDa}$ vs LNPs $\sim 800\text{ kDa}$ vs AAV $\sim 3.7\text{ MDa}$).
pub fn drug_bioavailability(capillary_density: f64, molecular_weight_da: f64) -> f64 {
    let cap = capillary_density.clamp(0.05, 1.0);
    let mw = molecular_weight_da.max(50.0);

    // Log10 molecular weight penalty: small molecules (~300-1000 Da) diffuse freely;
    // viral capsids and LNPs (>10^5 Da) face trans-endothelial steric hindrance, especially in rarefied beds
    let log_mw = mw.log10(); // ~3.0 for 1 kDa, ~5.9 for LNP, ~6.57 for AAV
    let size_resistance = ((log_mw - 2.7).max(0.0) * 0.075).clamp(0.0, 0.42);

    // Non-linear perfusion coupling: rarefied capillary beds disproportionately impair large-vector delivery
    let perfusion_efficiency = cap.powf(0.85);
    let net_bioavailability = perfusion_efficiency * (1.0 - size_resistance * (1.25 - 0.50 * cap));

    net_bioavailability.clamp(0.05, 0.99)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vascular_aging_and_drug_bioavailability_scaling() {
        let young = VascularPerfusionState::youthful();
        let aged = VascularPerfusionState::from_somatic_age(72.0);
        assert!(young.capillary_density > aged.capillary_density + 0.30);
        assert!(young.po2_mmhg > aged.po2_mmhg + 8.0);

        // Small molecule vs AAV capsid in aged vs young capillary beds
        let nav_aged = drug_bioavailability(aged.capillary_density, MW_SMALL_MOLECULE_NAVITOCLAX_DA);
        let nav_young = drug_bioavailability(young.capillary_density, MW_SMALL_MOLECULE_NAVITOCLAX_DA);
        let aav_aged = drug_bioavailability(aged.capillary_density, MW_AAV_CAPSID_VECTOR_DA);
        let aav_young = drug_bioavailability(young.capillary_density, MW_AAV_CAPSID_VECTOR_DA);

        assert!(nav_young > nav_aged);
        assert!(aav_young > aav_aged);
        // Macromolecular AAV faces higher trans-endothelial hindrance than small-molecule Navitoclax
        assert!(nav_aged > aav_aged);
    }

    #[test]
    fn test_microvascular_angiogenesis_restores_po2_and_delivery() {
        let mut aged = VascularPerfusionState::from_somatic_age(75.0);
        let before_aav = drug_bioavailability(aged.capillary_density, MW_AAV_CAPSID_VECTOR_DA);
        aged.restore_microvascular_angiogenesis(1.0);
        let after_aav = drug_bioavailability(aged.capillary_density, MW_AAV_CAPSID_VECTOR_DA);

        assert!(aged.capillary_density > 0.85);
        assert!(aged.po2_mmhg > 38.0);
        assert!(after_aav > before_aav + 0.20);
    }
}
