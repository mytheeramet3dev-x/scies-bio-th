//! Closed-Loop Hemodynamic Biomechanics (Moens-Korteweg Pulse Wave Velocity `PWV`),
//! Whole-Body First-Pass Vector Biodistribution (`PBPK`), and Multi-Organ Endocrine Feedback Axes
//! (`HPA`, `GH/IGF-1`, `RAAS`, and `alpha-Klotho / FGF23`).
//!
//! References:
//! - O'Rourke & Safar (2005) *Hypertension* 46:200–206 (Relationship between aortic stiffening and microvascular disease in brain and kidney: pulsatile barotrauma).
//! - Zhang et al. (2013) *Nature* 497:211–216 (Hypothalamic programming of systemic ageing via IKK-beta/NF-kB and GnRH).
//! - Kuro-o et al. (1997) *Nature* 390:45–51 (Mutation of the mouse klotho gene leads to a syndrome resembling ageing).
//! - Poon et al. (2019) *ACS Nano* 13:5785–5798 (Hepatobiliary & splenic reticuloendothelial Kupffer/macrophage sequestration of nanoparticles).

use super::organ::{AnatomicalOrgan, OrganArchitecture};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Systemic hemodynamic and Moens-Korteweg pulse-wave biomechanics state.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HemodynamicBiomechanicsState {
    /// Aortic / arterial Pulse Wave Velocity ($\text{PWV}$, $\text{m/s}$; youthful $\approx 5.8\text{ m/s}$, aged $> 11.5\text{ m/s}$).
    pub pulse_wave_velocity_m_s: f64,
    /// Central aortic systolic-diastolic pulse pressure ($\Delta P$, $\text{mmHg}$; youthful $\approx 40\text{ mmHg}$).
    pub central_pulse_pressure_mmhg: f64,
    /// Pulsatile microvascular barotrauma index ($0.0$ to $1.0$) transmitted to low-impedance
    /// organs (Brain `BBB` & Renal `Glomeruli`).
    pub microvascular_barotrauma_index: f64,
    /// Left ventricular afterload stress index ($0.0$ to $1.0$) driving myocardial hypertrophy and fibrosis.
    pub left_ventricular_afterload_stress: f64,
}

impl HemodynamicBiomechanicsState {
    /// Computes Moens-Korteweg Pulse Wave Velocity ($\text{PWV} = \sqrt{\frac{E \cdot h}{2 \rho R}}$)
    /// and downstream microvascular barotrauma from the arterial wall's Young's elastic modulus (`arterial_modulus_kpa`).
    pub fn from_arterial_modulus(arterial_modulus_kpa: f64) -> Self {
        let e_kpa = arterial_modulus_kpa.max(10.0);
        // Calibrated Moens-Korteweg relation: at youthful E_0 = 28 kPa, PWV ~ 5.80 m/s
        // h/2R ~ 0.062, rho ~ 1.06 g/cm^3 -> PWV = 1.096 * sqrt(E_kPa)
        let pwv = (1.096 * e_kpa.sqrt()).clamp(4.5, 18.0);
        let pulse_pressure = (24.0 + 2.85 * pwv).clamp(35.0, 95.0);

        // Microvascular barotrauma surges non-linearly when PWV exceeds 7.5 m/s (wave reflection arrives in mid-systole)
        let excess_pwv = (pwv - 6.2).max(0.0);
        let barotrauma = (excess_pwv.powi(2) / (14.0 + excess_pwv.powi(2))).clamp(0.0, 1.0);
        let afterload = ((pulse_pressure - 40.0).max(0.0) / 55.0).clamp(0.0, 1.0);

        Self {
            pulse_wave_velocity_m_s: pwv,
            central_pulse_pressure_mmhg: pulse_pressure,
            microvascular_barotrauma_index: barotrauma,
            left_ventricular_afterload_stress: afterload,
        }
    }
}

/// Multi-organ endocrine and humoral longevity axes coupling all 20 organs.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SystemicEndocrineAxes {
    /// Hypothalamic `NF-kB` neuroinflammation index ($0.0\text{--}1.0$; suppresses `GnRH` and drives systemic aging; Zhang et al. 2013).
    pub hypothalamic_nfkb_inflammation: f64,
    /// Hypothalamic `GnRH` & circadian `SCN` neuroendocrine output ($0.0\text{--}1.0$).
    pub hypothalamic_gnrh_scn_tone: f64,
    /// Adrenal `Cortisol / DHEA` stress-steroid ratio (youthful $\approx 1.0$; aged $> 3.2$ causes thymic involution & telomere loss).
    pub adrenal_cortisol_to_dhea_ratio: f64,
    /// Systemic Renin-Angiotensin-Aldosterone (`Angiotensin II`) tone ($0.0\text{--}1.0$; excess drives `TGF-beta` fibrosis & `NOX` ROS).
    pub raas_angiotensin_ii_tone: f64,
    /// Circulating soluble $\alpha\text{-Klotho}$ longevity hormone ($0.0\text{--}1.5$; secreted by Kidney & Choroid Plexus).
    pub circulating_alpha_klotho: f64,
    /// Hepatic / systemic `IGF-1` somatotropic signaling ($0.0\text{--}1.2$).
    pub hepatic_igf1_level: f64,
}

impl Default for SystemicEndocrineAxes {
    fn default() -> Self {
        Self::youthful()
    }
}

impl SystemicEndocrineAxes {
    /// Youthful homeostatic endocrine profile.
    pub fn youthful() -> Self {
        Self {
            hypothalamic_nfkb_inflammation: 0.06,
            hypothalamic_gnrh_scn_tone: 0.94,
            adrenal_cortisol_to_dhea_ratio: 1.00,
            raas_angiotensin_ii_tone: 0.14,
            circulating_alpha_klotho: 1.00,
            hepatic_igf1_level: 0.85,
        }
    }

    /// Estimates age-associated neuroendocrine dysregulation across the HPA, RAAS, Klotho, and Somatotropic axes.
    pub fn from_somatic_age(age_years: f64) -> Self {
        let excess = (age_years.clamp(0.0, 120.0) - 20.0).max(0.0);
        Self {
            hypothalamic_nfkb_inflammation: (0.06 + excess * 0.011).clamp(0.05, 0.92),
            hypothalamic_gnrh_scn_tone: (0.94 - excess * 0.0105).clamp(0.15, 0.96),
            adrenal_cortisol_to_dhea_ratio: (1.00 + excess * 0.042).clamp(0.85, 4.50),
            raas_angiotensin_ii_tone: (0.14 + excess * 0.010).clamp(0.10, 0.90),
            circulating_alpha_klotho: (1.00 - excess * 0.0115).clamp(0.18, 1.30),
            hepatic_igf1_level: (0.85 - excess * 0.0065).clamp(0.30, 1.10),
        }
    }

    /// Computes the net systemic endocrine pro-aging drift multiplier ($0.0$ to $1.5$).
    /// When rejuvenated ($\alpha\text{-Klotho} \ge 0.95$, low `NF-kB`, low `Ang II`, balanced `Cortisol/DHEA`),
    /// systemic pro-aging drift is minimized.
    pub fn systemic_pro_aging_drift_factor(&self) -> f64 {
        let hpa_stress = ((self.adrenal_cortisol_to_dhea_ratio - 1.0).max(0.0) * 0.22).clamp(0.0, 0.65);
        let nfkb_stress = self.hypothalamic_nfkb_inflammation * 0.45;
        let raas_fibrosis = self.raas_angiotensin_ii_tone * 0.35;
        let klotho_protection = (self.circulating_alpha_klotho * 0.48).clamp(0.05, 0.70);
        let gnrh_protection = (self.hypothalamic_gnrh_scn_tone * 0.25).clamp(0.02, 0.30);

        (0.45 + hpa_stress + nfkb_stress + raas_fibrosis - klotho_protection - gnrh_protection)
            .clamp(0.05, 1.50)
    }

    /// Executes neuroendocrine and renal-klotho axis rejuvenation:
    /// - Suppresses hypothalamic `IKK-beta`/`NF-kB` microglial inflammation & restores `GnRH`/`SCN` output.
    /// - Normalizes adrenal `Cortisol/DHEA` zona reticularis steroidogenesis.
    /// - Inhibits `RAAS` (`Angiotensin II` / `AT1R`) fibrotic signaling.
    /// - Upregulates renal & choroid plexus soluble `alpha-Klotho` transcription.
    pub fn apply_endocrine_rejuvenation(
        &mut self,
        hypothalamic_gnrh_boost: f64,
        klotho_raas_modulation: f64,
    ) {
        let hg = hypothalamic_gnrh_boost.clamp(0.0, 1.5);
        let kr = klotho_raas_modulation.clamp(0.0, 1.5);

        self.hypothalamic_nfkb_inflammation =
            (self.hypothalamic_nfkb_inflammation * (1.0 - 0.68 * hg)).clamp(0.04, 0.95);
        self.hypothalamic_gnrh_scn_tone =
            (self.hypothalamic_gnrh_scn_tone + 0.45 * hg).clamp(0.15, 0.98);
        let excess_cort = (self.adrenal_cortisol_to_dhea_ratio - 1.0).max(0.0);
        self.adrenal_cortisol_to_dhea_ratio =
            (1.0 + excess_cort * (1.0 - 0.65 * hg)).clamp(0.90, 4.50);
        self.raas_angiotensin_ii_tone =
            (self.raas_angiotensin_ii_tone * (1.0 - 0.62 * kr)).clamp(0.10, 0.95);
        self.circulating_alpha_klotho =
            (self.circulating_alpha_klotho + 0.52 * kr).clamp(0.18, 1.25);
        self.hepatic_igf1_level = (0.72 + 0.12 * hg).clamp(0.40, 0.95);
    }
}

/// First-pass Reticuloendothelial System (`RES` — Hepatic Kupffer / Splenic Marginal Macrophage)
/// sequestration and organ-specific PBPK biodistribution profile.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct OrganBiodistributionShare {
    pub organ_ta2: String,
    pub organ_name: String,
    /// Fraction of systemic dose sequestered by reticuloendothelial clearance in this organ ($0.0\text{--}1.0$).
    pub res_first_pass_sequestration: f64,
    /// Effective parenchymal bioavailability reaching target cells after first-pass & barrier transport ($0.0\text{--}1.0$).
    pub effective_parenchymal_delivery: f64,
}

/// Computes whole-body closed-loop PBPK biodistribution across a slice of [`AnatomicalOrgan`]s,
/// accounting for:
/// 1. Hepatic (`TA2:3025`) & Splenic (`TA2:5172`) Reticuloendothelial System (`RES`) first-pass uptake
///    (which traps $55\text{--}75\%$ of unshielded nanoparticles unless stealth CD47/PEG-lipid shielding is high).
/// 2. Organ Cardiac Output fraction ($\%\text{ CO}$) and endothelial barrier permeability (`BBB`, `Sinusoidal`, `Fenestrated`).
pub fn compute_whole_body_pbpk_biodistribution(
    organs: &[AnatomicalOrgan],
    molecular_weight_da: f64,
    stealth_cd47_peg_shielding: f64,
    receptor_shuttle_efficiency: f64,
) -> Vec<OrganBiodistributionShare> {
    let shield = stealth_cd47_peg_shielding.clamp(0.0, 1.0);
    let is_nanoparticle_or_capsid = molecular_weight_da > 50_000.0;

    // Unshielded nanoparticles/viral capsids suffer ~68% hepatic + splenic RES sequestration;
    // CD47 "don't-eat-me" peptide + stealth lipid engineering reduces RES trapping down to ~18%.
    let total_res_trap = if is_nanoparticle_or_capsid {
        (0.68 * (1.0 - 0.74 * shield)).clamp(0.12, 0.75)
    } else {
        0.10 * (1.0 - 0.50 * shield)
    };

    let systemic_escape_fraction = (1.0 - total_res_trap).clamp(0.20, 0.95);

    organs
        .iter()
        .map(|organ| {
            let ta2 = organ.ta2_code();
            let res_seq = if ta2 == "TA2:3025" {
                total_res_trap * 0.78 // Liver Kupffer cells & LSECs
            } else if ta2 == "TA2:5172" {
                total_res_trap * 0.22 // Splenic marginal zone macrophages
            } else {
                0.0
            };

            let base_parenchymal = organ
                .mean_parenchymal_bioavailability(molecular_weight_da, receptor_shuttle_efficiency);
            // Boosted by stealth circulation half-life extension for non-RES organs
            let effective_delivery = if ta2 == "TA2:3025" || ta2 == "TA2:5172" {
                base_parenchymal
            } else {
                (base_parenchymal * (0.45 + 0.65 * systemic_escape_fraction)).clamp(0.02, 0.99)
            };

            OrganBiodistributionShare {
                organ_ta2: ta2.to_string(),
                organ_name: organ.organ_name().to_string(),
                res_first_pass_sequestration: res_seq,
                effective_parenchymal_delivery: effective_delivery,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::anatomy::organ::AnatomicalSpatialAtlas;

    #[test]
    fn test_moens_korteweg_pwv_barotrauma_and_rejuvenation() {
        let young_hemo = HemodynamicBiomechanicsState::from_arterial_modulus(28.0);
        let aged_hemo = HemodynamicBiomechanicsState::from_arterial_modulus(68.0);

        assert!(young_hemo.pulse_wave_velocity_m_s < 6.2);
        assert!(aged_hemo.pulse_wave_velocity_m_s > 8.8);
        assert!(aged_hemo.microvascular_barotrauma_index > young_hemo.microvascular_barotrauma_index + 0.30);
    }

    #[test]
    fn test_endocrine_axes_and_pbpk_res_shielding() {
        let mut axes = SystemicEndocrineAxes::from_somatic_age(75.0);
        let aged_drift = axes.systemic_pro_aging_drift_factor();
        axes.apply_endocrine_rejuvenation(1.0, 1.0);
        assert!(axes.systemic_pro_aging_drift_factor() < aged_drift * 0.45);
        assert!(axes.circulating_alpha_klotho > 0.85);

        let atlas = AnatomicalSpatialAtlas::from_somatic_age(70.0);
        let unshielded =
            compute_whole_body_pbpk_biodistribution(&atlas.organs, 800_000.0, 0.0, 0.0);
        let shielded =
            compute_whole_body_pbpk_biodistribution(&atlas.organs, 800_000.0, 0.92, 0.92);

        // Brain (index 1) delivery increases substantially with CD47 stealth + TfR BBB shuttle
        assert!(
            shielded[1].effective_parenchymal_delivery
                > unshielded[1].effective_parenchymal_delivery * 3.0
        );
        // Liver (index 2) unwanted RES sequestration drops with CD47 stealth shielding
        assert!(shielded[2].res_first_pass_sequestration < unshielded[2].res_first_pass_sequestration * 0.40);
    }
}
