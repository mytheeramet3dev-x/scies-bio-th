//! # Organ-Specific Physiological Behaviors & Biomarker Dynamics (`behavior`)
//!
//! Provides deep mathematical and clinical physiology implementations (*Guyton & Hall*, *Boron & Boulpaep*)
//! for all 20 concrete *Terminologia Anatomica* (`TA2`) organ classes.
//!
//! Each organ models its fundamental biophysical action, coupling intrinsic microarchitecture,
//! extracellular matrix stiffness ($E$), and microvascular oxygenation ($p\text{O}_2$) to real-world
//! functional outputs (e.g. Frank-Starling cardiac output, Fick alveolar gas diffusion, hepatic
//! zonated detoxification, glomerular filtration eGFR, and circadian neuroendocrine pacing).

use super::organ::{
    AdrenalSuprarenalOrgan, AnatomicalOrgan, ArterialVascularOrgan, BoneMarrowNicheOrgan,
    CerebralCortexOrgan, ChoroidPlexusCsfOrgan, DermalCutaneousOrgan, GastrointestinalCryptOrgan,
    HepaticLobuleOrgan, HippocampalDentateGyrusOrgan, HypothalamicNeuroendocrineOrgan,
    MyocardialVentricularOrgan, PancreaticIsletOrgan, PituitaryHypophysealOrgan,
    PulmonaryAlveolarOrgan, RenalNephronOrgan, SkeletalMuscleFascicleOrgan,
    SplenicLymphoidPulpOrgan, ThermogenicBrownAdiposeOrgan, ThymicLymphoidOrgan,
    VisceralWhiteAdiposeOrgan,
};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

// =================================================================================================
// 1. HEART MYOCARDIUM: Frank-Starling Mechanics & Cardiac Cycle
// =================================================================================================

/// Real-time metrics of a single cardiac cycle in the left ventricular myocardium.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CardiacCycleMetrics {
    /// Heart rate in beats per minute ($\text{bpm}$).
    pub heart_rate_bpm: f64,
    /// End-Diastolic Volume ($\text{EDV}$, $\text{mL}$; youthful $\approx 120\text{ mL}$).
    pub end_diastolic_volume_ml: f64,
    /// End-Systolic Volume ($\text{ESV}$, $\text{mL}$; youthful $\approx 45\text{ mL}$).
    pub end_systolic_volume_ml: f64,
    /// Stroke Volume ($\text{SV} = \text{EDV} - \text{ESV}$, $\text{mL}$; youthful $\approx 75\text{ mL}$).
    pub stroke_volume_ml: f64,
    /// Ejection Fraction ($\text{EF} = \frac{\text{SV}}{\text{EDV}} \times 100\%$; normal $55\text{--}70\%$).
    pub ejection_fraction_pct: f64,
    /// Cardiac Output ($\text{CO} = \text{SV} \times \text{HR}$, $\text{L/min}$; normal $4.5\text{--}6.5\text{ L/min}$).
    pub cardiac_output_l_min: f64,
    /// Diastolic Compliance Index ($0.0\text{--}1.0$; severely impaired in aged/stiff ventricles).
    pub diastolic_compliance_index: f64,
    /// Myocardial Oxygen Demand ($\text{MVO}_2$ index, arbitrary units).
    pub myocardial_o2_demand: f64,
    /// B-type Natriuretic Peptide secretion rate ($\text{BNP}$, $\text{pg/mL}$; surges during wall stretch).
    pub bnp_secretion_pg_ml: f64,
}

impl MyocardialVentricularOrgan {
    /// Simulates ventricular cardiac mechanics based on the **Frank-Starling Law** and diastolic stiffness.
    ///
    /// In an aged, fibrotic heart ($E_{\text{myo}} \gg 1.8\text{ kPa}$, high glucosepane `AGEs`), diastolic filling
    /// is impeded ($\text{EDV}$ drops or requires elevated filling pressure, modeling `HFpEF`).
    pub fn simulate_cardiac_cycle(
        &self,
        heart_rate_bpm: f64,
        venous_filling_pressure_mmhg: f64,
        aortic_afterload_impedance: f64,
    ) -> CardiacCycleMetrics {
        let hr = heart_rate_bpm.clamp(40.0, 190.0);
        let e_kpa = self.ecm.youngs_modulus_kpa.max(0.5);
        let ages = self.ecm.ages_crosslink_burden.clamp(0.0, 1.0);

        // Diastolic compliance diminishes inversely with myocardial Young's modulus & AGE cross-links
        // Calibrated to youthful myocardium E_0 = 18.0 kPa -> compliance ~ 1.0
        let compliance = (18.0 / (e_kpa * (1.0 + 0.85 * ages))).clamp(0.42, 1.25);

        // End-Diastolic Volume governed by filling pressure & compliance (normal resting EDV ~120 mL; HFpEF ~98-105 mL)
        let edv = (76.0 + 7.2 * venous_filling_pressure_mmhg.clamp(2.0, 25.0) * compliance)
            .clamp(88.0, 160.0);

        // Frank-Starling inotropic contractility factor (preserved ~54-58% EF at rest in HFpEF, impaired contractile reserve)
        let o2_factor = (self.vascular.po2_mmhg / 40.0).clamp(0.72, 1.18);
        let contractility = (0.64 * o2_factor * (1.0 - 0.10 * self.ecm.yap_taz_mechanotransduction_stress()))
            .clamp(0.52, 0.74);

        // Afterload impedance reduces stroke volume ejection
        let afterload_penalty = (1.0 + 0.12 * (aortic_afterload_impedance - 1.0).max(0.0)).clamp(1.0, 1.6);
        let sv = (edv * contractility / afterload_penalty).clamp(42.0, 115.0);
        let esv = (edv - sv).max(15.0);
        let ef = (sv / edv) * 100.0;
        let co = (sv * hr) / 1000.0;

        // Laplace myocardial wall stress sigma = (P * r) / (2 * h) drives BNP release
        let wall_stress = (edv * aortic_afterload_impedance / compliance).clamp(20.0, 600.0);
        let bnp = (15.0 + 0.45 * (wall_stress - 80.0).max(0.0).powf(1.3)).clamp(10.0, 1200.0);
        let mvo2 = (hr * sv * aortic_afterload_impedance * 0.0018).clamp(0.5, 4.5);

        CardiacCycleMetrics {
            heart_rate_bpm: hr,
            end_diastolic_volume_ml: edv,
            end_systolic_volume_ml: esv,
            stroke_volume_ml: sv,
            ejection_fraction_pct: ef,
            cardiac_output_l_min: co,
            diastolic_compliance_index: compliance,
            myocardial_o2_demand: mvo2,
            bnp_secretion_pg_ml: bnp,
        }
    }
}

// =================================================================================================
// 2. PULMONARY ALVEOLI: Fick Gas Exchange, Surfactant & Acid-Base Ventilation
// =================================================================================================

/// Real-time metrics of pulmonary gas exchange across the alveolar-capillary membrane.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PulmonaryGasExchangeMetrics {
    /// Minute ventilation rate ($\text{L/min}$).
    pub minute_ventilation_l_min: f64,
    /// Alveolar ventilation rate ($\text{L/min}$; excluding physiological anatomical dead space).
    pub alveolar_ventilation_l_min: f64,
    /// Alveolar Oxygen Partial Pressure ($P_A\text{O}_2$, $\text{mmHg}$; normal $\approx 102\text{ mmHg}$).
    pub alveolar_po2_mmhg: f64,
    /// Arterial Oxygen Partial Pressure ($P_a\text{O}_2$, $\text{mmHg}$; normal $90\text{--}100\text{ mmHg}$).
    pub arterial_po2_mmhg: f64,
    /// Arterial Oxygen Saturation ($S_a\text{O}_2$, $\%$; normal $96\text{--}99\%$).
    pub arterial_o2_saturation_pct: f64,
    /// Arterial Carbon Dioxide Partial Pressure ($P_a\text{CO}_2$, $\text{mmHg}$; normal $38\text{--}42\text{ mmHg}$).
    pub arterial_pco2_mmhg: f64,
    /// Arterial Blood $\text{pH}$ calculated from Henderson-Hasselbalch equation (normal $7.38\text{--}7.42$).
    pub arterial_ph: f64,
    /// Diffusing capacity of the lung for oxygen/CO ($DL_{\text{CO}}$, $\text{mL/min/mmHg}$).
    pub diffusing_capacity_dlco: f64,
    /// Alveolar Dynamic Compliance ($C_L$, $\text{mL/cmH}_2\text{O}$).
    pub dynamic_lung_compliance: f64,
}

impl PulmonaryAlveolarOrgan {
    /// Simulates alveolar ventilation, Fick's membrane diffusion, and hemoglobin $\text{O}_2$ saturation.
    pub fn simulate_respiratory_gas_exchange(
        &self,
        tidal_volume_ml: f64,
        respiratory_rate_bpm: f64,
        fraction_inspired_o2: f64,
        cardiac_output_l_min: f64,
    ) -> PulmonaryGasExchangeMetrics {
        let rr = respiratory_rate_bpm.clamp(6.0, 45.0);
        let vt = tidal_volume_ml.clamp(200.0, 1200.0);
        let fio2 = fraction_inspired_o2.clamp(0.21, 1.0);
        let _co = cardiac_output_l_min.clamp(2.0, 20.0);

        let minute_vent = (vt * rr) / 1000.0;
        let e_kpa = self.ecm.youngs_modulus_kpa.max(0.5);
        let ages = self.ecm.ages_crosslink_burden.clamp(0.0, 1.0);
        // Physiological dead space = anatomical (150 mL) + alveolar V/Q mismatch in aged/fibrotic lung
        let vq_mismatch_dead_space = 45.0 * (1.0 - self.vascular.capillary_density).max(0.0)
            + 6.0 * (e_kpa - 2.0).max(0.0);
        let dead_space_ml = (150.0 + vq_mismatch_dead_space).clamp(145.0, 240.0);
        let alv_vent = ((vt - dead_space_ml).max(50.0) * rr) / 1000.0;

        // Alveolar Gas Equation: PAO2 = FiO2*(Patm - PH2O) - PaCO2/R
        // Patm = 760, PH2O = 47 -> 713 mmHg; Respiratory Quotient R ~ 0.8
        // PaCO2 = (VCO2 * 0.863) / VA (with resting VA ~ 4.9 L/min yielding PaCO2 ~ 40 mmHg -> 196.0 / VA)
        let pco2 = (196.0 / alv_vent).clamp(25.0, 75.0);
        let pao2 = (fio2 * 713.0 - (pco2 / 0.8)).clamp(50.0, 650.0);

        // Fick's Law of Diffusion: DLCO drops with interstitial fibrosis & basement membrane thickening
        let dlco = (28.0 / ((e_kpa / 2.0).sqrt() * (1.0 + 0.8 * ages))).clamp(8.0, 38.0);

        // A-a gradient (Alveolar - arterial gradient) increases with pulmonary aging
        let a_a_gradient = (4.0 + 0.35 * (e_kpa - 2.0).max(0.0) + (1.0 / (dlco / 25.0))).clamp(3.0, 45.0);
        let pao2_art = (pao2 - a_a_gradient).clamp(35.0, 600.0);

        // Hill equation for Oxygen-Hemoglobin Saturation: n=2.8, P50=26.6 mmHg
        let p50 = 26.6;
        let n = 2.8;
        let p_ratio = pao2_art / p50;
        let sao2 = (p_ratio.powf(n) / (1.0 + p_ratio.powf(n))) * 100.0;

        // Henderson-Hasselbalch equation: pH = 6.1 + log10([HCO3-] / (0.0307 * PaCO2))
        let hco3 = 24.0; // Standard serum bicarbonate mEq/L
        let ph = (6.1 + ((hco3 / (0.0307 * pco2)).max(0.1)).log10()).clamp(7.10, 7.60);

        // Lung compliance CL = deltaV / deltaP
        let compliance = (85.0 / (e_kpa * (1.0 + 0.5 * ages))).clamp(15.0, 110.0);

        PulmonaryGasExchangeMetrics {
            minute_ventilation_l_min: minute_vent,
            alveolar_ventilation_l_min: alv_vent,
            alveolar_po2_mmhg: pao2,
            arterial_po2_mmhg: pao2_art,
            arterial_o2_saturation_pct: sao2.clamp(50.0, 100.0),
            arterial_pco2_mmhg: pco2,
            arterial_ph: ph,
            diffusing_capacity_dlco: dlco,
            dynamic_lung_compliance: compliance,
        }
    }
}

// =================================================================================================
// 3. ARTERIAL VASCULAR WALL: Windkessel Compliance, Blood Pressure & eNOS Shear Stress
// =================================================================================================

/// Real-time arterial hemodynamics, Windkessel cushioning, and endothelial shear tone.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ArterialHemodynamicsMetrics {
    /// Systolic Blood Pressure ($P_{\text{sys}}$, $\text{mmHg}$; normal $110\text{--}120\text{ mmHg}$).
    pub systolic_bp_mmhg: f64,
    /// Diastolic Blood Pressure ($P_{\text{dia}}$, $\text{mmHg}$; normal $70\text{--}80\text{ mmHg}$).
    pub diastolic_bp_mmhg: f64,
    /// Mean Arterial Pressure ($\text{MAP} = \frac{2}{3}P_{\text{dia}} + \frac{1}{3}P_{\text{sys}}$, $\text{mmHg}$).
    pub mean_arterial_pressure_mmhg: f64,
    /// Total Arterial Windkessel Compliance ($C_{\text{art}}$, $\text{mL/mmHg}$; youthful $\approx 1.8\text{ mL/mmHg}$).
    pub windkessel_compliance_ml_mmhg: f64,
    /// Endothelial Shear Stress ($\tau_{\text{wall}}$, $\text{dyn/cm}^2$; normal $15\text{--}25\text{ dyn/cm}^2$).
    pub endothelial_shear_stress_dyn_cm2: f64,
    /// Endothelial Nitric Oxide release rate ($\text{NO}$, arbitrary fold relative to youthful baseline).
    pub nitric_oxide_release_rate: f64,
    /// Downstream Wave Reflection Augmentation Index ($\text{AIx}$, $\%$; rises with aortic stiffness).
    pub augmentation_index_pct: f64,
}

impl ArterialVascularOrgan {
    /// Simulates aortic Windkessel buffering, pulse pressure amplification, and endothelial eNOS shear stress.
    pub fn simulate_arterial_hemodynamics(
        &self,
        stroke_volume_ml: f64,
        cardiac_output_l_min: f64,
        systemic_vascular_resistance_wood: f64,
    ) -> ArterialHemodynamicsMetrics {
        let sv = stroke_volume_ml.clamp(30.0, 130.0);
        let co = cardiac_output_l_min.clamp(2.0, 15.0);
        let svr = systemic_vascular_resistance_wood.clamp(10.0, 30.0); // Normal ~15-20 Wood units (~80 * Wood = dyn*s/cm5)

        let e_kpa = self.ecm.youngs_modulus_kpa.max(10.0);
        let ages = self.ecm.ages_crosslink_burden.clamp(0.0, 1.0);

        // Windkessel Compliance: Cart = V / (E_aorta * (1 + AGEs))
        // Youthful E_0 = 28 kPa -> Cart ~ 1.85 mL/mmHg; Aged E ~ 75 kPa -> Cart ~ 0.62 mL/mmHg
        let compliance = (54.0 / (e_kpa * (1.0 + 0.85 * ages))).clamp(0.55, 2.10);

        // Pulse Pressure = SV / Compliance (Youthful ~ 40 mmHg; Aged ISH ~ 75-85 mmHg)
        let pulse_pressure = (sv / compliance).clamp(34.0, 96.0);

        // Mean Arterial Pressure = CO * SVR * 0.90 (Youthful ~ 85-92 mmHg; Aged ISH ~ 100-110 mmHg)
        let map = (co * svr * 0.90).clamp(65.0, 155.0);
        // In aortic stiffening, early systolic wave reflection augments systolic pressure while lowering diastolic support
        let wave_reflection_shift = 0.08 * (pulse_pressure - 40.0).max(0.0);
        let diastolic = (map - (pulse_pressure / 3.0) - wave_reflection_shift).clamp(58.0, 105.0);
        let systolic = (diastolic + pulse_pressure).clamp(95.0, 205.0);

        // Poiseuille wall shear stress tau = 4*mu*Q / (pi*r^3)
        let shear = (18.0 * (co / 5.0) * (self.vascular.capillary_density / 0.85)).clamp(5.0, 35.0);

        // eNOS activation driven by shear stress and youthful endothelial integrity
        let enos = (self.vascular.enos_activity * (shear / 18.0)).clamp(0.15, 1.40);

        // Augmentation Index AIx reflects early wave reflection from stiff peripheral muscular arteries
        let aix = (10.0 + 0.65 * (e_kpa - 28.0).max(0.0) + 15.0 * ages).clamp(5.0, 65.0);

        ArterialHemodynamicsMetrics {
            systolic_bp_mmhg: systolic,
            diastolic_bp_mmhg: diastolic,
            mean_arterial_pressure_mmhg: map,
            windkessel_compliance_ml_mmhg: compliance,
            endothelial_shear_stress_dyn_cm2: shear,
            nitric_oxide_release_rate: enos,
            augmentation_index_pct: aix,
        }
    }
}

// =================================================================================================
// 4. RENAL NEPHRONS: Starling Glomerular Filtration (eGFR), Reabsorption & Klotho Secretion
// =================================================================================================

/// Real-time metrics of renal filtration, podocyte integrity, and tubular reabsorption.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct RenalFiltrationMetrics {
    /// Estimated Glomerular Filtration Rate ($\text{eGFR}$, $\text{mL/min/1.73 m}^2$; youthful $\ge 110$).
    pub egfr_ml_min: f64,
    /// Renal Blood Flow ($\text{RBF}$, $\text{L/min}$; normal $\approx 1.1\text{ L/min}$).
    pub renal_blood_flow_l_min: f64,
    /// Filtration Fraction ($\text{FF} = \frac{\text{GFR}}{\text{RPF}}$; normal $18\text{--}22\%$).
    pub filtration_fraction_pct: f64,
    /// Proteinuria / Albuminuria Albumin Leak Index ($0.0\text{--}1.0$; surges with podocyte damage).
    pub albuminuria_leak_index: f64,
    /// Distal Tubule Soluble $\alpha$-Klotho secretion rate ($\text{pg/mL}$).
    pub alpha_klotho_secretion_pg_ml: f64,
    /// Juxtaglomerular Renin Release Rate (arbitrary units; surges when renal perfusion pressure drops).
    pub renin_release_rate: f64,
    /// Sodium Fractional Excretion ($\text{FE}_{\text{Na}}$, $\%$; normal $0.5\text{--}1.5\%$).
    pub fractional_sodium_excretion_pct: f64,
}

impl RenalNephronOrgan {
    /// Simulates Starling glomerular filtration, podocyte slit diaphragm permeability, and tubular endocrine output.
    pub fn simulate_renal_filtration(
        &self,
        mean_arterial_pressure_mmhg: f64,
        cardiac_output_l_min: f64,
        barotrauma_index: f64,
    ) -> RenalFiltrationMetrics {
        let rbf = (cardiac_output_l_min * 0.20 * (self.vascular.capillary_density / 0.85)).clamp(0.40, 1.45);
        let rpf = rbf * 0.55; // Renal plasma flow (hematocrit ~ 45%)

        let e_kpa = self.ecm.youngs_modulus_kpa.max(1.0);
        let ages = self.ecm.ages_crosslink_burden.clamp(0.0, 1.0);

        // Autoregulation of glomerular capillary hydraulic pressure P_GC between 75-180 mmHg MAP
        let p_gc = (50.0 + 0.08 * (mean_arterial_pressure_mmhg.clamp(70.0, 180.0) - 90.0)).clamp(45.0, 58.0);
        let p_bs = 12.0; // Bowman's space hydraulic pressure
        let pi_gc = 28.0; // Glomerular oncotic pressure
        let net_filtration_pressure = (p_gc - p_bs - pi_gc).max(4.0);

        // Starling filtration coefficient Kf damaged by glomerulosclerosis, microvascular barotrauma, and fibrosis
        let podocyte_health = (1.0 - 0.45 * barotrauma_index - 0.35 * ages).clamp(0.15, 1.0);
        let k_f = (11.0 * podocyte_health * (5.0 / e_kpa).sqrt()).clamp(1.5, 13.5);

        let egfr = (k_f * net_filtration_pressure).clamp(15.0, 135.0);
        let ff = (egfr / (rpf * 1000.0)) * 100.0;

        // Podocyte effacement and glomerular barrier breakdown causes albuminuria leak
        let albumin_leak = (0.02 + 0.60 * (1.0 - podocyte_health).powi(2)).clamp(0.01, 0.95);

        // Distal convoluted tubule soluble alpha-Klotho secretion
        let klotho = (850.0 * podocyte_health * (self.vascular.po2_mmhg / 38.0)).clamp(120.0, 1250.0);

        // Juxtaglomerular renin release: inverse with renal perfusion pressure + driven by intrarenal microvascular ischemia (nephrosclerosis)
        let base_renin = if mean_arterial_pressure_mmhg < 85.0 {
            (1.0 + 0.05 * (85.0 - mean_arterial_pressure_mmhg).powf(1.4)).clamp(1.0, 6.0)
        } else {
            (1.0 / (1.0 + 0.02 * (mean_arterial_pressure_mmhg - 85.0))).clamp(0.2, 1.0)
        };
        let intrarenal_ischemia =
            1.8 * (0.85 - self.vascular.capillary_density).max(0.0) + 1.4 * barotrauma_index;
        let renin = (base_renin + intrarenal_ischemia).clamp(0.25, 6.0);

        RenalFiltrationMetrics {
            egfr_ml_min: egfr,
            renal_blood_flow_l_min: rbf,
            filtration_fraction_pct: ff.clamp(5.0, 35.0),
            albuminuria_leak_index: albumin_leak,
            alpha_klotho_secretion_pg_ml: klotho,
            renin_release_rate: renin,
            fractional_sodium_excretion_pct: 1.0,
        }
    }
}

// =================================================================================================
// 5. HEPATIC LOBULE: Zonated Metabolism (Zone 1 Gluconeogenesis/Urea vs Zone 3 CYP Detox)
// =================================================================================================

/// Real-time metrics of hepatic zonated metabolism, detoxification, and bile production.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HepaticMetabolicMetrics {
    /// Hepatic Sinusoidal Blood Flow ($\text{L/min}$; $\approx 25\%$ of cardiac output).
    pub total_hepatic_blood_flow_l_min: f64,
    /// Ammonia $\to$ Urea Conversion Efficiency ($\%$; Zone 1 periportal periportal arginase/CPS-1).
    pub urea_cycle_clearance_efficiency_pct: f64,
    /// Cytochrome P450 Phase I/II Detoxification Rate (Zone 3 pericentral clearance multiplier).
    pub cyp_detoxification_capacity_multiplier: f64,
    /// Hepatic Glucose Output ($\text{mg/min}$; balanced by insulin/glucagon).
    pub hepatic_glucose_output_mg_min: f64,
    /// Serum Albumin Synthesis Rate ($\text{g/day}$; normal $12\text{--}15\text{ g/day}$).
    pub albumin_synthesis_rate_g_day: f64,
    /// Bile Acid Secretion Rate ($\text{mL/day}$; normal $500\text{--}800\text{ mL/day}$).
    pub bile_acid_secretion_ml_day: f64,
    /// Systemic Somatotropic IGF-1 Synthesis Output ($\text{ng/mL}$).
    pub igf1_synthesis_ng_ml: f64,
}

impl HepaticLobuleOrgan {
    /// Simulates zonated hepatic parenchymal functions (Zone 1 Urea Cycle vs Zone 3 CYP450 clearance).
    pub fn simulate_zonated_metabolism(
        &self,
        portal_flow_fraction: f64,
        circulating_insulin_pmol: f64,
        gh_somatotroph_stimulus: f64,
        systemic_ammonia_load: f64,
    ) -> HepaticMetabolicMetrics {
        let total_flow = (1.45 * portal_flow_fraction.clamp(0.6, 1.4)).clamp(0.7, 2.2);

        let e_kpa = self.ecm.youngs_modulus_kpa.max(0.5);
        let fibrosis_penalty = (5.2 / e_kpa).clamp(0.25, 1.15);
        let sin_health = (self.vascular.capillary_density / 0.85).clamp(0.3, 1.2);

        // Zone 1 Periportal Urea Cycle: detoxifies ammonia into urea
        let urea_eff = (98.0 * fibrosis_penalty * sin_health * (1.0 / (1.0 + 0.02 * systemic_ammonia_load)))
            .clamp(25.0, 99.0);

        // Zone 3 Pericentral CYP450 Detoxification (CYP3A4/2E1)
        let cyp_cap = (1.10 * fibrosis_penalty * sin_health * (self.vascular.po2_mmhg / 38.0)).clamp(0.20, 1.40);

        // Hepatic Glucose Output: suppressed by insulin, elevated by fasting/glucagon
        let insulin_suppression = (1.0 / (1.0 + circulating_insulin_pmol / 60.0)).clamp(0.15, 1.0);
        let hgo = (150.0 * insulin_suppression * fibrosis_penalty).clamp(30.0, 280.0);

        // Albumin and IGF-1 protein translation by functional hepatocytes
        let albumin = (14.5 * fibrosis_penalty * sin_health).clamp(3.5, 18.0);
        let igf1 = (185.0 * gh_somatotroph_stimulus.clamp(0.2, 1.8) * fibrosis_penalty).clamp(35.0, 260.0);
        let bile = (650.0 * fibrosis_penalty * sin_health).clamp(150.0, 850.0);

        HepaticMetabolicMetrics {
            total_hepatic_blood_flow_l_min: total_flow,
            urea_cycle_clearance_efficiency_pct: urea_eff,
            cyp_detoxification_capacity_multiplier: cyp_cap,
            hepatic_glucose_output_mg_min: hgo,
            albumin_synthesis_rate_g_day: albumin,
            bile_acid_secretion_ml_day: bile,
            igf1_synthesis_ng_ml: igf1,
        }
    }
}

// =================================================================================================
// 6. CEREBRAL CORTEX: Neurovascular Coupling, Synaptic Plasticity & Sleep Glymphatics
// =================================================================================================

/// Real-time metrics of cerebral cortex computation, BBB integrity, and metabolic coupling.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CorticalCognitiveMetrics {
    /// Cerebral Blood Flow ($\text{CBF}$, $\text{mL/100g/min}$; normal $50\text{--}55$).
    pub cerebral_blood_flow_ml_100g_min: f64,
    /// Cerebral Metabolic Rate of Oxygen ($\text{CMRO}_2$, $\mu\text{mol/100g/min}$; normal $\approx 150$).
    pub cmro2_umol_100g_min: f64,
    /// Synaptic Density & Plasticity Index ($0.0\text{--}1.0$).
    pub synaptic_plasticity_index: f64,
    /// Blood-Brain Barrier Transendothelial Electrical Resistance ($\text{TEER}$, $\Omega\cdot\text{cm}^2$).
    pub bbb_teer_ohm_cm2: f64,
    /// Interstitial Glymphatic Waste Clearance Flux ($\mu\text{L/min}$; surges during slow-wave sleep).
    pub glymphatic_waste_clearance_rate: f64,
    /// Neurovascular Coupling Latency ($\text{seconds}$; delayed in stiff microvessels).
    pub neurovascular_coupling_latency_sec: f64,
}

impl CerebralCortexOrgan {
    /// Simulates cortical neurovascular coupling, BBB tight junction impedance, and glymphatic clearance.
    pub fn simulate_cortical_computation(
        &self,
        mean_arterial_pressure_mmhg: f64,
        sleep_slow_wave_fraction: f64,
        microvascular_barotrauma_index: f64,
    ) -> CorticalCognitiveMetrics {
        // Autoregulation maintains CBF ~52 mL/100g/min between MAP 60-150 mmHg
        let autoreg = if mean_arterial_pressure_mmhg < 60.0 {
            (mean_arterial_pressure_mmhg / 60.0).clamp(0.4, 1.0)
        } else if mean_arterial_pressure_mmhg > 150.0 {
            (1.0 + 0.005 * (mean_arterial_pressure_mmhg - 150.0)).clamp(1.0, 1.4)
        } else {
            1.0
        };

        let cbf = (52.0 * autoreg * (self.vascular.capillary_density / 0.85)).clamp(22.0, 75.0);
        let cmro2 = (150.0 * (self.vascular.po2_mmhg / 38.0) * autoreg).clamp(65.0, 185.0);

        // BBB tight junctions (Claudin-5, Occludin) degraded by barotrauma and collagen fibrosis
        let e_kpa = self.ecm.youngs_modulus_kpa.max(0.5);
        let bbb_health = (1.0 - 0.50 * microvascular_barotrauma_index - 0.25 * (e_kpa - 1.5).max(0.0)).clamp(0.15, 1.0);
        let teer = (1800.0 * bbb_health).clamp(250.0, 2200.0);

        // Glymphatic convective flow surges 4-fold during slow-wave sleep via astrocytic AQP4 polarization
        let sleep_factor = 1.0 + 3.2 * sleep_slow_wave_fraction.clamp(0.0, 1.0);
        let glymphatic = (8.5 * sleep_factor * (1.0 / (e_kpa / 1.5).sqrt())).clamp(3.0, 42.0);

        let latency = (0.85 + 1.2 * (e_kpa - 1.5).max(0.0) + 0.9 * microvascular_barotrauma_index).clamp(0.7, 3.8);
        let plasticity = (0.92 * bbb_health * (self.vascular.po2_mmhg / 38.0)).clamp(0.15, 0.98);

        CorticalCognitiveMetrics {
            cerebral_blood_flow_ml_100g_min: cbf,
            cmro2_umol_100g_min: cmro2,
            synaptic_plasticity_index: plasticity,
            bbb_teer_ohm_cm2: teer,
            glymphatic_waste_clearance_rate: glymphatic,
            neurovascular_coupling_latency_sec: latency,
        }
    }
}

// =================================================================================================
// 7. HIPPOCAMPUS & DENTATE GYRUS: Adult Neurogenesis, BDNF & Episodic Memory Encoding
// =================================================================================================

/// Real-time metrics of hippocampal subgranular neurogenesis and synaptic LTP.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HippocampalFunctionMetrics {
    /// Subgranular Zone Adult Neural Stem Cell Proliferation Rate ($\text{new neurons/day}$).
    pub adult_neurogenesis_rate_cells_day: f64,
    /// Synaptic Long-Term Potentiation (LTP) Fidelity ($0.0\text{--}1.0$).
    pub ltp_synaptic_fidelity: f64,
    /// Dentate Gyrus Granule Cell Excitability Index.
    pub granule_cell_excitability: f64,
    /// Cortisol Neurotoxicity Vulnerability Score ($0.0\text{--}1.0$).
    pub cortisol_vulnerability_score: f64,
}

impl HippocampalDentateGyrusOrgan {
    /// Simulates adult hippocampal neurogenesis in the SGZ and LTP retention modulated by BDNF, Klotho, and Cortisol.
    pub fn simulate_neurogenesis_and_memory(
        &self,
        circulating_alpha_klotho: f64,
        adrenal_cortisol_ratio: f64,
        circulating_bdnf_ng_ml: f64,
    ) -> HippocampalFunctionMetrics {
        let klotho_boost = circulating_alpha_klotho.clamp(0.2, 1.5);
        let bdnf = circulating_bdnf_ng_ml.clamp(5.0, 45.0) / 25.0;
        let stress = (adrenal_cortisol_ratio - 1.0).max(0.0);

        let e_kpa = self.ecm.youngs_modulus_kpa.max(0.4);
        let niche_factor = (1.4 / e_kpa).clamp(0.3, 1.2) * (self.vascular.po2_mmhg / 36.0);

        // Neurogenesis in SGZ: stimulated by Klotho & BDNF, repressed by cortisol stress
        let neuro = (750.0 * niche_factor * klotho_boost * bdnf / (1.0 + 0.45 * stress.powi(2)))
            .clamp(40.0, 1100.0);

        let ltp = (0.88 * klotho_boost * bdnf / (1.0 + 0.35 * stress)).clamp(0.12, 0.98);
        let excitability = (1.0 + 0.15 * (klotho_boost - 1.0)).clamp(0.4, 1.3);
        let vuln = (0.20 + 0.25 * stress).clamp(0.10, 0.95);

        HippocampalFunctionMetrics {
            adult_neurogenesis_rate_cells_day: neuro,
            ltp_synaptic_fidelity: ltp,
            granule_cell_excitability: excitability,
            cortisol_vulnerability_score: vuln,
        }
    }
}

// =================================================================================================
// 8. HYPOTHALAMUS: SCN Master Circadian Clock & Neuroendocrine Releasing Hormones
// =================================================================================================

/// Real-time metrics of hypothalamic circadian gating, GnRH, and neuroinflammation.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HypothalamicRegulationMetrics {
    /// Circadian Master Oscillator Tone ($0.0\text{--}1.0$).
    pub scn_circadian_amplitude: f64,
    /// Pulsatile GnRH Secretion Output (arbitrary units; key systemic longevity driver; Zhang et al. 2013).
    pub gnrh_pulsatile_release: f64,
    /// Corticotropin-Releasing Hormone (`CRH`) Drive into hypophyseal portal vessels.
    pub crh_drive: f64,
    /// Growth Hormone-Releasing Hormone (`GHRH`) Drive.
    pub ghrh_drive: f64,
    /// Microglial IKK-$\beta$/NF-$\kappa$B Neuroinflammation Tone ($0.0\text{--}1.0$).
    pub microglial_nfkb_inflammation: f64,
    /// Preoptic Thermostat Setpoint Core Temperature ($\,^\circ\text{C}$).
    pub preoptic_temperature_setpoint_c: f64,
}

impl HypothalamicNeuroendocrineOrgan {
    /// Simulates master hypothalamic regulation, circadian amplitude, and neuroendocrine outputs.
    pub fn simulate_master_regulation(
        &self,
        circadian_phase_hour: f64,
        systemic_inflammaging_sasp: f64,
    ) -> HypothalamicRegulationMetrics {
        // Circadian oscillation: 24-hour cosine rhythm
        let rad = (circadian_phase_hour / 24.0) * 2.0 * std::f64::consts::PI;
        let base_circ = 0.5 + 0.5 * (rad - 1.0).cos();

        let e_kpa = self.ecm.youngs_modulus_kpa.max(0.5);
        let sasp = systemic_inflammaging_sasp.clamp(0.0, 1.0);

        // IKK-beta / NF-kB surges with systemic SASP and tissue aging
        let nfkb = (0.08 + 0.65 * sasp + 0.15 * (e_kpa - 1.4).max(0.0)).clamp(0.05, 0.95);

        // GnRH is actively repressed by hypothalamic NF-kB
        let gnrh = ((1.0 - 0.75 * nfkb) * (self.vascular.po2_mmhg / 38.0)).clamp(0.08, 1.10);
        let scn_amp = ((1.0 - 0.45 * nfkb) * base_circ).clamp(0.15, 1.0);

        // CRH surges under neuroinflammation / stress
        let crh = (1.0 + 1.2 * nfkb).clamp(0.8, 2.5);
        let ghrh = (1.0 - 0.50 * nfkb).clamp(0.2, 1.2);

        HypothalamicRegulationMetrics {
            scn_circadian_amplitude: scn_amp,
            gnrh_pulsatile_release: gnrh,
            crh_drive: crh,
            ghrh_drive: ghrh,
            microglial_nfkb_inflammation: nfkb,
            preoptic_temperature_setpoint_c: 37.0,
        }
    }
}

// =================================================================================================
// 9. PITUITARY GLAND: Hypophyseal Endocrine Cascade (ACTH, GH, TSH, LH/FSH)
// =================================================================================================

/// Real-time metrics of anterior and posterior pituitary endocrine secretion.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PituitaryEndocrineMetrics {
    /// Adrenocorticotropic Hormone (`ACTH`, $\text{pg/mL}$; normal $10\text{--}60\text{ pg/mL}$).
    pub acth_pg_ml: f64,
    /// Growth Hormone pulsatile release (`GH`, $\text{ng/mL}$).
    pub growth_hormone_ng_ml: f64,
    /// Thyroid-Stimulating Hormone (`TSH`, $\mu\text{IU/mL}$; normal $0.5\text{--}4.5$).
    pub tsh_miu_ml: f64,
    /// Antidiuretic Hormone / Vasopressin (`AVP/ADH`, $\text{pg/mL}$).
    pub vasopressin_pg_ml: f64,
}

impl PituitaryHypophysealOrgan {
    /// Simulates the anterior hypophyseal cascade driven by hypothalamic releasing hormones and feedback inhibition.
    pub fn simulate_hypophyseal_cascade(
        &self,
        crh_stimulus: f64,
        ghrh_stimulus: f64,
        feedback_cortisol_ratio: f64,
        feedback_igf1_level: f64,
    ) -> PituitaryEndocrineMetrics {
        // Corticotroph ACTH secretion: driven by CRH, suppressed by high cortisol negative feedback
        let acth = (25.0 * crh_stimulus / (0.5 + 0.5 * feedback_cortisol_ratio)).clamp(8.0, 95.0);

        // Somatotroph GH secretion: stimulated by GHRH, suppressed by IGF-1 negative feedback
        let gh = (4.5 * ghrh_stimulus / (0.4 + 0.6 * feedback_igf1_level.clamp(0.2, 2.0))).clamp(0.4, 9.5);

        let tsh = (1.8 * (self.vascular.po2_mmhg / 38.0)).clamp(0.4, 5.0);
        let avp = 2.4;

        PituitaryEndocrineMetrics {
            acth_pg_ml: acth,
            growth_hormone_ng_ml: gh,
            tsh_miu_ml: tsh,
            vasopressin_pg_ml: avp,
        }
    }
}

// =================================================================================================
// 10. ADRENAL GLAND: Zona Glomerulosa, Fasciculata & Reticularis Steroidogenesis
// =================================================================================================

/// Real-time metrics of adrenal steroid and catecholamine synthesis.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct AdrenalSteroidMetrics {
    /// Plasma Cortisol ($\mu\text{g/dL}$; normal diurnal peak $10\text{--}20$, trough $3\text{--}8$).
    pub cortisol_ug_dl: f64,
    /// Plasma Dehydroepiandrosterone Sulfate (`DHEA-S`, $\mu\text{g/dL}$; plunges with age).
    pub dhea_s_ug_dl: f64,
    /// Cortisol / DHEA-S Catabolic-to-Anabolic Ratio (youthful $\approx 0.05$; aged $> 0.25$).
    pub cortisol_to_dhea_ratio: f64,
    /// Aldosterone Secretion ($\text{ng/dL}$; driven by Angiotensin II & plasma $\text{K}^+$).
    pub aldosterone_ng_dl: f64,
    /// Medullary Epinephrine / Norepinephrine Secretion Rate (arbitrary units).
    pub catecholamine_secretion_rate: f64,
}

impl AdrenalSuprarenalOrgan {
    /// Simulates adrenal steroidogenesis across the cortex (Aldosterone, Cortisol, DHEA-S).
    pub fn simulate_steroidogenesis(
        &self,
        acth_pg_ml: f64,
        angiotensin_ii_tone: f64,
        sympathetic_stress_tone: f64,
    ) -> AdrenalSteroidMetrics {
        let acth_factor = (acth_pg_ml / 25.0).clamp(0.4, 3.5);

        // Cortisol from Zona Fasciculata: driven by ACTH
        let cortisol = (12.0 * acth_factor).clamp(4.0, 38.0);

        // DHEA-S from Zona Reticularis: severely degenerates with aging (adrenopause)
        let e_kpa = self.ecm.youngs_modulus_kpa.max(0.5);
        let reticularis_integrity = (4.0 / e_kpa).clamp(0.15, 1.25);
        let dhea_s = (280.0 * reticularis_integrity * acth_factor.sqrt()).clamp(15.0, 450.0);

        // Normalized Catabolic-to-Anabolic Ratio relative to youthful setpoint (12.0 ug/dL / 280.0 ug/dL = 1.0)
        let ratio = ((cortisol / 12.0) / (dhea_s / 280.0)).clamp(0.50, 8.50);

        // Aldosterone from Zona Glomerulosa: driven by Angiotensin II
        let aldo = (8.5 * (1.0 + 1.8 * angiotensin_ii_tone.clamp(0.0, 1.5))).clamp(4.0, 42.0);
        let catechol = (1.0 + 1.5 * sympathetic_stress_tone.clamp(0.0, 1.0)).clamp(0.8, 3.5);

        AdrenalSteroidMetrics {
            cortisol_ug_dl: cortisol,
            dhea_s_ug_dl: dhea_s,
            cortisol_to_dhea_ratio: ratio,
            aldosterone_ng_dl: aldo,
            catecholamine_secretion_rate: catechol,
        }
    }
}

// =================================================================================================
// 11. GASTROINTESTINAL TRACT: Epithelial Barrier, Incretin Secretion & Nutrient Absorption
// =================================================================================================

/// Real-time metrics of gut epithelial barrier integrity, incretin release, and endotoxemia.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct GastrointestinalFunctionMetrics {
    /// Epithelial Barrier Integrity Index ($0.0\text{--}1.0$; Claudin-1, Occludin tight junctions).
    pub gut_barrier_integrity_index: f64,
    /// Portal Vein Lipopolysaccharide (`LPS`) Endotoxemia Translocation Rate ($\text{EU/mL}$).
    pub portal_lps_endotoxin_translocation_eu_ml: f64,
    /// Incretin Glucagon-Like Peptide-1 (`GLP-1`) Secretion ($\text{pmol/L}$).
    pub glp1_incretin_secretion_pmol_l: f64,
    /// Crypt Lgr5+ Stem Cell Turnover Rate ($\text{renewals/week}$; normal $\approx 1.8$).
    pub crypt_stem_cell_turnover_rate: f64,
    /// Carbohydrate & Amino Acid Absorption Capacity ($\text{g/hr}$).
    pub nutrient_absorption_rate_g_hr: f64,
}

impl GastrointestinalCryptOrgan {
    /// Simulates intestinal epithelial turnover, tight-junction barrier sealing, and GLP-1 incretin secretion.
    pub fn simulate_digestion_and_barrier(
        &self,
        caloric_intake_kcal: f64,
        dietary_fiber_g: f64,
    ) -> GastrointestinalFunctionMetrics {
        let e_kpa = self.ecm.youngs_modulus_kpa.max(0.5);
        let ages = self.ecm.ages_crosslink_burden.clamp(0.0, 1.0);

        // Gut barrier weakened by mucosal fibrosis and AGE cross-links, protected by fiber short-chain fatty acids (SCFAs)
        let scfa_boost = (dietary_fiber_g / 30.0).clamp(0.3, 1.5);
        let barrier = ((1.0 - 0.40 * ages) * scfa_boost * (2.8 / e_kpa).sqrt()).clamp(0.20, 0.98);

        // LPS leaky gut translocates into portal blood when barrier degrades
        let lps = (0.15 + 2.5 * (1.0 - barrier).powi(2)).clamp(0.12, 4.5);

        // GLP-1 released by enteroendocrine L-cells in response to caloric intake
        let glp1 = (12.0 + 0.025 * caloric_intake_kcal.clamp(0.0, 1200.0) * barrier).clamp(8.0, 48.0);
        let turnover = (1.8 * (self.vascular.po2_mmhg / 36.0) * (2.8 / e_kpa).sqrt()).clamp(0.5, 2.4);
        let absorption = (45.0 * barrier).clamp(12.0, 65.0);

        GastrointestinalFunctionMetrics {
            gut_barrier_integrity_index: barrier,
            portal_lps_endotoxin_translocation_eu_ml: lps,
            glp1_incretin_secretion_pmol_l: glp1,
            crypt_stem_cell_turnover_rate: turnover,
            nutrient_absorption_rate_g_hr: absorption,
        }
    }
}

// =================================================================================================
// 12. PANCREATIC ISLETS: Biphasic Glucose-Stimulated Insulin (GSIS) & Glucagon Homeostasis
// =================================================================================================

/// Real-time metrics of pancreatic $\beta$-cell insulin and $\alpha$-cell glucagon secretion.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PancreaticIsletMetrics {
    /// Circulating Insulin Secretion Rate ($\text{pmol/L}$; baseline fasting $\approx 35\text{--}60$).
    pub insulin_secretion_rate_pmol_l: f64,
    /// Circulating Glucagon Secretion Rate ($\text{pg/mL}$; baseline $\approx 50\text{--}100$).
    pub glucagon_secretion_rate_pg_ml: f64,
    /// $\beta$-Cell Functional Mass Index ($0.0\text{--}1.0$).
    pub beta_cell_functional_mass: f64,
    /// Biphasic GSIS First-Phase Acute Insulin Spike ($\text{fold}$).
    pub gsis_first_phase_spike_fold: f64,
}

impl PancreaticIsletOrgan {
    /// Simulates glucose-stimulated insulin secretion (GSIS) and glucagon counter-regulation.
    pub fn simulate_islet_secretion(
        &self,
        blood_glucose_mg_dl: f64,
        glp1_pmol_l: f64,
    ) -> PancreaticIsletMetrics {
        let e_kpa = self.ecm.youngs_modulus_kpa.max(0.5);
        let ages = self.ecm.ages_crosslink_burden.clamp(0.0, 1.0);

        // Islet amyloid polypeptide (IAPP) and ECM stiffening reduce functional beta cell mass
        let beta_mass = ((1.0 - 0.45 * ages) * (3.8 / e_kpa).sqrt() * (self.vascular.po2_mmhg / 38.0))
            .clamp(0.20, 1.0);

        let bg = blood_glucose_mg_dl.clamp(40.0, 400.0);
        let glp1_potentiator = 1.0 + 0.02 * (glp1_pmol_l - 12.0).max(0.0);

        // Biphasic GSIS: Sigmoidal Hill function above threshold 90 mg/dL
        let glucose_stim = if bg > 90.0 {
            let excess = bg - 90.0;
            (excess.powf(1.8) / (70.0_f64.powf(1.8) + excess.powf(1.8))).clamp(0.0, 1.0)
        } else {
            0.0
        };

        let insulin = (45.0 + 280.0 * glucose_stim * beta_mass * glp1_potentiator).clamp(25.0, 450.0);
        let spike = (1.0 + 3.8 * glucose_stim * beta_mass).clamp(1.0, 5.0);

        // Glucagon secreted during hypoglycemia, suppressed by hyperglycemia and paracrine insulin
        let glucagon = if bg < 90.0 {
            (80.0 + 2.5 * (90.0 - bg)).clamp(75.0, 240.0)
        } else {
            (80.0 / (1.0 + 0.02 * (bg - 90.0))).clamp(20.0, 80.0)
        };

        PancreaticIsletMetrics {
            insulin_secretion_rate_pmol_l: insulin,
            glucagon_secretion_rate_pg_ml: glucagon,
            beta_cell_functional_mass: beta_mass,
            gsis_first_phase_spike_fold: spike,
        }
    }
}

// =================================================================================================
// 13. BONE MARROW: Hematopoiesis, Erythropoiesis & Myeloid vs Lymphoid Lineage Skewing
// =================================================================================================

/// Real-time metrics of bone marrow hematopoietic stem cell niche output and lineage commitment.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HematopoieticNicheMetrics {
    /// Daily Erythrocyte (RBC) Production Rate ($\times 10^{11}\text{ cells/day}$; normal $\approx 2.0$).
    pub daily_rbc_production_10e11: f64,
    /// Daily Platelet Production Rate ($\times 10^{11}\text{ cells/day}$; normal $\approx 1.0$).
    pub daily_platelet_production_10e11: f64,
    /// Myeloid-to-Lymphoid Differentiation Ratio (youthful $\approx 1.2$; aged myeloid skewing $> 3.5$).
    pub myeloid_to_lymphoid_skewing_ratio: f64,
    /// HSC Endosteal Niche Quiescence Retention Tone ($0.0\text{--}1.0$).
    pub hsc_quiescence_retention: f64,
}

impl BoneMarrowNicheOrgan {
    /// Simulates hematopoiesis, erythropoietin response, and age-associated myeloid skewing.
    pub fn simulate_hematopoiesis(
        &self,
        erythropoietin_level: f64,
        sympathetic_circadian_tone: f64,
    ) -> HematopoieticNicheMetrics {
        let epo = erythropoietin_level.clamp(0.4, 3.5);
        let e_kpa = self.ecm.youngs_modulus_kpa.max(1.0);
        let niche_health = (8.5 / e_kpa).clamp(0.3, 1.25) * (self.vascular.po2_mmhg / 38.0);

        let rbc = (2.0 * epo * niche_health).clamp(0.6, 4.2);
        let plt = (1.0 * niche_health).clamp(0.4, 2.2);

        // Myeloid skewing increases as marrow adipose tissue expands and niche stiffens
        let skewing = (1.2 + 0.45 * (e_kpa - 8.5).max(0.0)).clamp(1.0, 5.5);
        let quiescence = (0.85 * (1.0 - 0.25 * sympathetic_circadian_tone.clamp(0.0, 1.0))).clamp(0.3, 0.95);

        HematopoieticNicheMetrics {
            daily_rbc_production_10e11: rbc,
            daily_platelet_production_10e11: plt,
            myeloid_to_lymphoid_skewing_ratio: skewing,
            hsc_quiescence_retention: quiescence,
        }
    }
}

// =================================================================================================
// 14. THYMUS: T-Cell Receptor Selection & Naive T-Cell Repertoire Output
// =================================================================================================

/// Real-time metrics of thymic positive/negative selection and naive T cell emigration.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ThymicSelectionMetrics {
    /// Daily Naive $\text{CD4}^+$ and $\text{CD8}^+$ T-cell Output ($\times 10^7\text{ cells/day}$).
    pub daily_naive_t_cell_output_10e7: f64,
    /// Cortical Positive Selection Efficiency ($0.0\text{--}1.0$).
    pub positive_selection_efficiency: f64,
    /// Medullary AIRE-Dependent Negative Selection Competence ($0.0\text{--}1.0$; prevents autoimmunity).
    pub negative_selection_competence: f64,
    /// Corticomedullary Functional Cellularity ($0.0\text{--}1.0$; collapses with thymic involution).
    pub functional_thymic_cellularity: f64,
}

impl ThymicLymphoidOrgan {
    /// Simulates thymic lymphopoiesis and naive T cell output (modeling age-associated involution).
    pub fn simulate_thymic_selection(
        &self,
        bone_marrow_lymphoid_influx: f64,
    ) -> ThymicSelectionMetrics {
        let e_kpa = self.ecm.youngs_modulus_kpa.max(0.5);
        let ages = self.ecm.ages_crosslink_burden.clamp(0.0, 1.0);

        // Thymic involution: adipose replacement and stromal stiffening
        let cellularity = ((3.8 / e_kpa).sqrt() * (1.0 - 0.65 * ages)).clamp(0.04, 1.0);

        let pos = (0.88 * cellularity).clamp(0.10, 0.95);
        let neg = (0.92 * cellularity).clamp(0.15, 0.98);

        let output = (12.0 * cellularity * bone_marrow_lymphoid_influx.clamp(0.3, 1.5)).clamp(0.2, 18.0);

        ThymicSelectionMetrics {
            daily_naive_t_cell_output_10e7: output,
            positive_selection_efficiency: pos,
            negative_selection_competence: neg,
            functional_thymic_cellularity: cellularity,
        }
    }
}

// =================================================================================================
// 15. SPLEEN: Red Pulp Erythrocyte Culling & White Pulp Immunological Screening
// =================================================================================================

/// Real-time metrics of splenic erythrocyte quality control, iron recycling, and immune screening.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SplenicImmunityMetrics {
    /// Senescent Erythrocyte Mechanical Clearance Rate ($\times 10^9\text{ cells/day}$).
    pub senescent_rbc_clearance_10e9: f64,
    /// Heme Iron Recycling Rate ($\text{mg/day}$; normal $20\text{--}25\text{ mg/day}$).
    pub heme_iron_recycling_mg_day: f64,
    /// Marginal Zone Macrophage Phagocytic Clearance Competence ($0.0\text{--}1.0$).
    pub macrophage_phagocytic_clearance: f64,
}

impl SplenicLymphoidPulpOrgan {
    /// Simulates red pulp mechanical erythrocyte culling and marginal zone immune screening.
    pub fn simulate_splenic_culling(
        &self,
        splenic_blood_flow_l_min: f64,
    ) -> SplenicImmunityMetrics {
        let flow = splenic_blood_flow_l_min.clamp(0.10, 0.60);
        let e_kpa = self.ecm.youngs_modulus_kpa.max(0.5);

        let slit_health = (4.2 / e_kpa).clamp(0.3, 1.25);
        let rbc_cull = (200.0 * (flow / 0.25) * slit_health).clamp(50.0, 320.0);
        let iron = (22.0 * slit_health).clamp(6.0, 28.0);
        let phag = (0.90 * slit_health * (self.vascular.po2_mmhg / 36.0)).clamp(0.20, 0.98);

        SplenicImmunityMetrics {
            senescent_rbc_clearance_10e9: rbc_cull,
            heme_iron_recycling_mg_day: iron,
            macrophage_phagocytic_clearance: phag,
        }
    }
}

// =================================================================================================
// 16. SKELETAL MUSCLE: Contractile Force, GLUT4 Glucose Disposal & Myokine Secretion
// =================================================================================================

/// Real-time metrics of skeletal muscle contraction, metabolic uptake, and myokine release.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MuscleFunctionMetrics {
    /// Specific Contractile Force ($\text{N/cm}^2$; youthful $\approx 25\text{ N/cm}^2$).
    pub specific_contractile_force_n_cm2: f64,
    /// GLUT4-Mediated Glucose Disposal Rate ($\text{mg/min}$; surges during contraction & insulin).
    pub glucose_disposal_rate_mg_min: f64,
    /// Irisin Myokine Secretion Rate ($\text{ng/mL}$; induces adipose browning).
    pub irisin_secretion_ng_ml: f64,
    /// Mechanical Work & Heat Generation Rate ($\text{Watts}$).
    pub mechanical_heat_generation_watts: f64,
}

impl SkeletalMuscleFascicleOrgan {
    /// Simulates cross-bridge force generation, GLUT4 translocation, and irisin myokine output.
    pub fn simulate_locomotor_contraction(
        &self,
        recruitment_fraction: f64,
        circulating_insulin_pmol: f64,
    ) -> MuscleFunctionMetrics {
        let recruit = recruitment_fraction.clamp(0.0, 1.0);
        let e_kpa = self.ecm.youngs_modulus_kpa.max(1.0);
        let ages = self.ecm.ages_crosslink_burden.clamp(0.0, 1.0);

        // Sarcomeric force reduced by intermuscular fibrosis and age cross-links (sarcopenia)
        let sarco_health = ((12.0 / e_kpa).sqrt() * (1.0 - 0.35 * ages)).clamp(0.25, 1.25);
        let force = (25.0 * sarco_health * recruit).clamp(0.0, 32.0);

        // Dual glucose uptake: insulin-stimulated + contraction-induced AMPK/GLUT4 translocation
        let ins_uptake = 45.0 * (circulating_insulin_pmol / 60.0).clamp(0.2, 2.5);
        let cont_uptake = 180.0 * recruit;
        let glucose_disp = ((ins_uptake + cont_uptake) * sarco_health).clamp(10.0, 350.0);

        let irisin = (8.5 * (1.0 + 1.8 * recruit) * sarco_health).clamp(4.0, 32.0);
        let heat = ((15.0 + 85.0 * recruit) * sarco_health).clamp(10.0, 120.0);

        MuscleFunctionMetrics {
            specific_contractile_force_n_cm2: force,
            glucose_disposal_rate_mg_min: glucose_disp,
            irisin_secretion_ng_ml: irisin,
            mechanical_heat_generation_watts: heat,
        }
    }
}

// =================================================================================================
// 17. DERMAL CUTANEOUS TISSUE: Barrier TEWL, Cutometer Elasticity & Thermoregulation
// =================================================================================================

/// Real-time metrics of dermal skin barrier, elastic recoil, and heat dissipation.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CutaneousBarrierMetrics {
    /// Transepidermal Water Loss ($\text{TEWL}$, $\text{g/m}^2\text{/h}$; youthful $5\text{--}8$, compromised $> 18$).
    pub transepidermal_water_loss_g_m2_h: f64,
    /// Cutometer Gross Elastic Recovery Parameter $R_2$ ($U_a / U_f$, $\%$; youthful $> 70\%$).
    pub cutometer_elastic_recovery_pct: f64,
    /// Dermal Blood Flow for Heat Dissipation ($\text{mL/min}$).
    pub cutaneous_heat_dissipation_flow_ml_min: f64,
}

impl DermalCutaneousOrgan {
    /// Simulates skin barrier hydration retention, cutometer recoil, and vascular heat exchange.
    pub fn simulate_barrier_and_thermoregulation(
        &self,
        ambient_temp_c: f64,
        core_temp_c: f64,
    ) -> CutaneousBarrierMetrics {
        let e_kpa = self.ecm.youngs_modulus_kpa.max(2.0);
        let ages = self.ecm.ages_crosslink_burden.clamp(0.0, 1.0);

        // Stratum corneum barrier integrity and TEWL
        let barrier_health = (1.0 - 0.40 * ages).clamp(0.2, 1.0);
        let tewl = (6.0 + 12.0 * (1.0 - barrier_health).powi(2)).clamp(4.5, 24.0);

        // Cutometer elastic recovery R2 (Ua/Uf) reflects elastin integrity vs cross-linked stiffness
        let r2 = (78.0 / (1.0 + 0.025 * (e_kpa - 13.0).max(0.0) + 0.6 * ages)).clamp(25.0, 85.0);

        // Dermal vasodilation for thermoregulation
        let delta_t = (core_temp_c - ambient_temp_c).max(0.0);
        let flow = (120.0 + 18.0 * delta_t).clamp(60.0, 650.0);

        CutaneousBarrierMetrics {
            transepidermal_water_loss_g_m2_h: tewl,
            cutometer_elastic_recovery_pct: r2,
            cutaneous_heat_dissipation_flow_ml_min: flow,
        }
    }
}

// =================================================================================================
// 18. CHOROID PLEXUS: Active CSF Secretion & Neuroprotective Ventricular Influx
// =================================================================================================

/// Real-time metrics of choroid plexus CSF production and endocrine delivery.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ChoroidCsfMetrics {
    /// Daily Cerebrospinal Fluid Production Rate ($\text{mL/day}$; normal $\approx 500\text{ mL/day}$).
    pub daily_csf_production_ml_day: f64,
    /// Ventricular Soluble $\alpha$-Klotho Secretion ($\text{pg/mL}$).
    pub ventricular_alpha_klotho_pg_ml: f64,
    /// Blood-CSF Barrier Impedance ($\text{TEER}$, $\Omega\cdot\text{cm}^2$).
    pub blood_csf_barrier_teer: f64,
}

impl ChoroidPlexusCsfOrgan {
    /// Simulates active CSF secretion via epithelial $\text{Na}^+/\text{K}^+$ ATPase and $\alpha$-Klotho delivery.
    pub fn simulate_csf_secretion(
        &self,
        cerebral_perfusion_pressure_mmhg: f64,
    ) -> ChoroidCsfMetrics {
        let cpp = cerebral_perfusion_pressure_mmhg.clamp(40.0, 110.0);
        let e_kpa = self.ecm.youngs_modulus_kpa.max(0.2);

        let epithelial_health = (2.0 / e_kpa).clamp(0.25, 1.25) * (self.vascular.po2_mmhg / 38.0);
        let csf_vol = (500.0 * epithelial_health * (cpp / 75.0).sqrt()).clamp(150.0, 620.0);
        let klotho = (520.0 * epithelial_health).clamp(80.0, 680.0);
        let teer = (250.0 * epithelial_health).clamp(80.0, 320.0);

        ChoroidCsfMetrics {
            daily_csf_production_ml_day: csf_vol,
            ventricular_alpha_klotho_pg_ml: klotho,
            blood_csf_barrier_teer: teer,
        }
    }
}

// =================================================================================================
// 19. VISCERAL WHITE ADIPOSE TISSUE: Lipolysis, Adiponectin & Inflammaging SASP
// =================================================================================================

/// Real-time metrics of visceral adipose lipid storage, adipokine secretion, and SASP flux.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct VisceralAdiposeMetrics {
    /// Free Fatty Acid (`FFA`) Portal Release Flux ($\mu\text{mol/min}$).
    pub free_fatty_acid_portal_flux_umol_min: f64,
    /// Circulating Adiponectin Level ($\mu\text{g/mL}$; youthful anti-inflammatory insulin-sensitizer $\approx 12\text{--}18$).
    pub circulating_adiponectin_ug_ml: f64,
    /// Circulating Leptin Level ($\text{ng/mL}$; normal $5\text{--}15$).
    pub circulating_leptin_ng_ml: f64,
    /// Visceral Adipose IL-6 / TNF-$\alpha$ SASP Paracrine Secretion Rate ($0.0\text{--}1.0$).
    pub visceral_sasp_secretion_rate: f64,
}

impl VisceralWhiteAdiposeOrgan {
    /// Simulates lipolysis, adipokine secretion, and senescence-associated inflammatory output.
    pub fn simulate_lipolysis_and_adipokines(
        &self,
        sympathetic_adrenergic_tone: f64,
        circulating_insulin_pmol: f64,
    ) -> VisceralAdiposeMetrics {
        let e_kpa = self.ecm.youngs_modulus_kpa.max(0.5);
        let ages = self.ecm.ages_crosslink_burden.clamp(0.0, 1.0);

        // Insulin strongly suppresses HSL lipolysis; adrenergic tone activates it
        let ins_suppress = (1.0 / (1.0 + circulating_insulin_pmol / 45.0)).clamp(0.1, 1.0);
        let ffa = (180.0 * sympathetic_adrenergic_tone.clamp(0.2, 2.5) * ins_suppress).clamp(20.0, 480.0);

        // Adiponectin declines with adipose hypertrophy, fibrosis, and aging
        let adiponectin = (15.0 / (1.0 + 0.35 * (e_kpa - 2.5).max(0.0) + 1.2 * ages)).clamp(2.5, 22.0);
        let leptin = (8.5 * (1.0 + 0.25 * (e_kpa - 2.5).max(0.0))).clamp(3.0, 35.0);

        // Visceral WAT is a major reservoir of senescent cells producing IL-6 / TNF-a
        let sasp = (0.05 + 0.45 * ages + 0.20 * (e_kpa - 2.5).max(0.0)).clamp(0.04, 0.95);

        VisceralAdiposeMetrics {
            free_fatty_acid_portal_flux_umol_min: ffa,
            circulating_adiponectin_ug_ml: adiponectin,
            circulating_leptin_ng_ml: leptin,
            visceral_sasp_secretion_rate: sasp,
        }
    }
}

// =================================================================================================
// 20. THERMOGENIC BROWN ADIPOSE TISSUE: UCP1 Uncoupling & Non-Shivering Thermogenesis
// =================================================================================================

/// Real-time metrics of brown adipose mitochondrial uncoupling and metabolic clearing.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct BrownAdiposeMetrics {
    /// Non-Shivering Thermogenesis Heat Generation Rate ($\text{Watts}$).
    pub non_shivering_thermogenesis_watts: f64,
    /// Mitochondrial UCP1 Proton Leak Uncoupling Index ($0.0\text{--}1.0$).
    pub ucp1_uncoupling_activity: f64,
    /// Plasma Triglyceride & Glucose Clearance Flux ($\text{mg/min}$).
    pub lipid_glucose_substrate_combustion_mg_min: f64,
}

impl ThermogenicBrownAdiposeOrgan {
    /// Simulates $\beta_3$-adrenergic activation of UCP1 proton leak and non-shivering thermogenesis.
    pub fn simulate_thermogenesis(
        &self,
        sympathetic_norepinephrine_stimulus: f64,
        ambient_temp_c: f64,
    ) -> BrownAdiposeMetrics {
        let cold_stimulus = ((24.0 - ambient_temp_c).max(0.0) / 10.0).clamp(0.0, 2.0);
        let total_drive = (sympathetic_norepinephrine_stimulus.clamp(0.1, 2.0) + cold_stimulus).clamp(0.2, 3.5);

        let e_kpa = self.ecm.youngs_modulus_kpa.max(0.5);
        let bat_health = (3.2 / e_kpa).clamp(0.2, 1.25) * (self.vascular.po2_mmhg / 38.0);

        let ucp1 = (0.85 * (total_drive / 1.5).clamp(0.1, 1.2) * bat_health).clamp(0.05, 0.98);
        let heat = (2.5 + 45.0 * ucp1).clamp(2.0, 52.0);
        let combust = (15.0 + 85.0 * ucp1).clamp(10.0, 110.0);

        BrownAdiposeMetrics {
            non_shivering_thermogenesis_watts: heat,
            ucp1_uncoupling_activity: ucp1,
            lipid_glucose_substrate_combustion_mg_min: combust,
        }
    }
}

// =================================================================================================
// OVERARCHING ORGAN SNAPSHOT ENUM
// =================================================================================================

/// Universal behavioral snapshot encapsulating the concrete real-time physiological metrics of any organ.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum OrganBehavioralSnapshot {
    ArterialVascular(ArterialHemodynamicsMetrics),
    CerebralCortex(CorticalCognitiveMetrics),
    HepaticLobule(HepaticMetabolicMetrics),
    BoneMarrowNiche(HematopoieticNicheMetrics),
    ThymicLymphoid(ThymicSelectionMetrics),
    DermalCutaneous(CutaneousBarrierMetrics),
    MyocardialVentricular(CardiacCycleMetrics),
    RenalNephron(RenalFiltrationMetrics),
    PulmonaryAlveolar(PulmonaryGasExchangeMetrics),
    PancreaticIslet(PancreaticIsletMetrics),
    SkeletalMuscleFascicle(MuscleFunctionMetrics),
    GastrointestinalCrypt(GastrointestinalFunctionMetrics),
    HypothalamicNeuroendocrine(HypothalamicRegulationMetrics),
    HippocampalDentateGyrus(HippocampalFunctionMetrics),
    ChoroidPlexusCsf(ChoroidCsfMetrics),
    PituitaryHypophyseal(PituitaryEndocrineMetrics),
    AdrenalSuprarenal(AdrenalSteroidMetrics),
    SplenicLymphoidPulp(SplenicImmunityMetrics),
    VisceralWhiteAdipose(VisceralAdiposeMetrics),
    ThermogenicBrownAdipose(BrownAdiposeMetrics),
}

impl AnatomicalOrgan {
    /// Evaluates the organ's intrinsic physiological behavior under canonical resting physiological drive.
    pub fn simulate_intrinsic_behavior(&self) -> OrganBehavioralSnapshot {
        match self {
            AnatomicalOrgan::ArterialVascular(o) => {
                OrganBehavioralSnapshot::ArterialVascular(o.simulate_arterial_hemodynamics(75.0, 5.4, 16.5))
            }
            AnatomicalOrgan::CerebralCortex(o) => {
                OrganBehavioralSnapshot::CerebralCortex(o.simulate_cortical_computation(90.0, 0.85, 0.05))
            }
            AnatomicalOrgan::HepaticLobule(o) => {
                OrganBehavioralSnapshot::HepaticLobule(o.simulate_zonated_metabolism(1.0, 45.0, 1.0, 25.0))
            }
            AnatomicalOrgan::BoneMarrowNiche(o) => {
                OrganBehavioralSnapshot::BoneMarrowNiche(o.simulate_hematopoiesis(1.0, 0.5))
            }
            AnatomicalOrgan::ThymicLymphoid(o) => {
                OrganBehavioralSnapshot::ThymicLymphoid(o.simulate_thymic_selection(1.0))
            }
            AnatomicalOrgan::DermalCutaneous(o) => {
                OrganBehavioralSnapshot::DermalCutaneous(o.simulate_barrier_and_thermoregulation(22.0, 36.8))
            }
            AnatomicalOrgan::MyocardialVentricular(o) => {
                OrganBehavioralSnapshot::MyocardialVentricular(o.simulate_cardiac_cycle(72.0, 6.0, 1.0))
            }
            AnatomicalOrgan::RenalNephron(o) => {
                OrganBehavioralSnapshot::RenalNephron(o.simulate_renal_filtration(90.0, 5.4, 0.0))
            }
            AnatomicalOrgan::PulmonaryAlveolar(o) => {
                OrganBehavioralSnapshot::PulmonaryAlveolar(o.simulate_respiratory_gas_exchange(500.0, 14.0, 0.21, 5.4))
            }
            AnatomicalOrgan::PancreaticIslet(o) => {
                OrganBehavioralSnapshot::PancreaticIslet(o.simulate_islet_secretion(92.0, 14.0))
            }
            AnatomicalOrgan::SkeletalMuscleFascicle(o) => {
                OrganBehavioralSnapshot::SkeletalMuscleFascicle(o.simulate_locomotor_contraction(0.6, 45.0))
            }
            AnatomicalOrgan::GastrointestinalCrypt(o) => {
                OrganBehavioralSnapshot::GastrointestinalCrypt(o.simulate_digestion_and_barrier(450.0, 28.0))
            }
            AnatomicalOrgan::HypothalamicNeuroendocrine(o) => {
                OrganBehavioralSnapshot::HypothalamicNeuroendocrine(o.simulate_master_regulation(12.0, 0.05))
            }
            AnatomicalOrgan::HippocampalDentateGyrus(o) => {
                OrganBehavioralSnapshot::HippocampalDentateGyrus(o.simulate_neurogenesis_and_memory(1.0, 1.0, 22.0))
            }
            AnatomicalOrgan::ChoroidPlexusCsf(o) => {
                OrganBehavioralSnapshot::ChoroidPlexusCsf(o.simulate_csf_secretion(75.0))
            }
            AnatomicalOrgan::PituitaryHypophyseal(o) => {
                OrganBehavioralSnapshot::PituitaryHypophyseal(o.simulate_hypophyseal_cascade(1.0, 1.0, 1.0, 1.0))
            }
            AnatomicalOrgan::AdrenalSuprarenal(o) => {
                OrganBehavioralSnapshot::AdrenalSuprarenal(o.simulate_steroidogenesis(24.0, 0.25, 0.35))
            }
            AnatomicalOrgan::SplenicLymphoidPulp(o) => {
                OrganBehavioralSnapshot::SplenicLymphoidPulp(o.simulate_splenic_culling(0.22))
            }
            AnatomicalOrgan::VisceralWhiteAdipose(o) => {
                OrganBehavioralSnapshot::VisceralWhiteAdipose(o.simulate_lipolysis_and_adipokines(0.35, 45.0))
            }
            AnatomicalOrgan::ThermogenicBrownAdipose(o) => {
                OrganBehavioralSnapshot::ThermogenicBrownAdipose(o.simulate_thermogenesis(0.35, 22.0))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_all_20_organ_concrete_behaviors() {
        // 1. Heart
        let heart = MyocardialVentricularOrgan::youthful();
        let beat = heart.simulate_cardiac_cycle(72.0, 8.0, 1.0);
        assert!(beat.ejection_fraction_pct > 55.0);
        assert!(beat.cardiac_output_l_min > 4.5);

        // 2. Lung
        let lung = PulmonaryAlveolarOrgan::youthful();
        let gas = lung.simulate_respiratory_gas_exchange(500.0, 14.0, 0.21, beat.cardiac_output_l_min);
        assert!(gas.arterial_po2_mmhg > 90.0);
        assert!(gas.arterial_o2_saturation_pct > 96.0);
        assert!(gas.arterial_ph > 7.35 && gas.arterial_ph < 7.45);

        // 3. Artery
        let artery = ArterialVascularOrgan::youthful();
        let hemo = artery.simulate_arterial_hemodynamics(beat.stroke_volume_ml, beat.cardiac_output_l_min, 18.0);
        assert!(hemo.systolic_bp_mmhg > 100.0 && hemo.systolic_bp_mmhg < 130.0);
        assert!(hemo.diastolic_bp_mmhg > 65.0 && hemo.diastolic_bp_mmhg < 85.0);

        // 4. Kidney
        let kidney = RenalNephronOrgan::youthful();
        let renal = kidney.simulate_renal_filtration(hemo.mean_arterial_pressure_mmhg, beat.cardiac_output_l_min, 0.0);
        assert!(renal.egfr_ml_min > 90.0);
        assert!(renal.alpha_klotho_secretion_pg_ml > 500.0);

        // 5. Liver
        let liver = HepaticLobuleOrgan::youthful();
        let hep = liver.simulate_zonated_metabolism(1.0, 50.0, 1.0, 1.0);
        assert!(hep.urea_cycle_clearance_efficiency_pct > 90.0);
        assert!(hep.albumin_synthesis_rate_g_day > 12.0);

        // 6. Brain
        let brain = CerebralCortexOrgan::youthful();
        let cort = brain.simulate_cortical_computation(hemo.mean_arterial_pressure_mmhg, 0.8, 0.0);
        assert!(cort.bbb_teer_ohm_cm2 > 1200.0);
        assert!(cort.glymphatic_waste_clearance_rate > 15.0);

        // 7. Hippocampus
        let hippo = HippocampalDentateGyrusOrgan::youthful();
        let hip = hippo.simulate_neurogenesis_and_memory(1.0, 1.0, 25.0);
        assert!(hip.adult_neurogenesis_rate_cells_day > 400.0);

        // 8. Hypothalamus
        let hypo = HypothalamicNeuroendocrineOrgan::youthful();
        let reg = hypo.simulate_master_regulation(12.0, 0.05);
        assert!(reg.gnrh_pulsatile_release > 0.8);

        // 9. Pituitary
        let pit = PituitaryHypophysealOrgan::youthful();
        let pit_res = pit.simulate_hypophyseal_cascade(reg.crh_drive, reg.ghrh_drive, 1.0, 1.0);
        assert!(pit_res.acth_pg_ml > 15.0);

        // 10. Adrenal
        let adr = AdrenalSuprarenalOrgan::youthful();
        let adr_res = adr.simulate_steroidogenesis(pit_res.acth_pg_ml, 0.2, 0.2);
        assert!(adr_res.dhea_s_ug_dl > 200.0);

        // 11. Gut
        let gut = GastrointestinalCryptOrgan::youthful();
        let gut_res = gut.simulate_digestion_and_barrier(500.0, 30.0);
        assert!(gut_res.gut_barrier_integrity_index > 0.85);

        // 12. Pancreas
        let panc = PancreaticIsletOrgan::youthful();
        let panc_res = panc.simulate_islet_secretion(105.0, gut_res.glp1_incretin_secretion_pmol_l);
        assert!(panc_res.insulin_secretion_rate_pmol_l > 50.0);

        // 13. Bone Marrow
        let marrow = BoneMarrowNicheOrgan::youthful();
        let marrow_res = marrow.simulate_hematopoiesis(1.0, 0.5);
        assert!(marrow_res.daily_rbc_production_10e11 > 1.5);

        // 14. Thymus
        let thymus = ThymicLymphoidOrgan::youthful();
        let thym_res = thymus.simulate_thymic_selection(1.0);
        assert!(thym_res.daily_naive_t_cell_output_10e7 > 8.0);

        // 15. Spleen
        let spleen = SplenicLymphoidPulpOrgan::youthful();
        let spl_res = spleen.simulate_splenic_culling(0.25);
        assert!(spl_res.senescent_rbc_clearance_10e9 > 100.0);

        // 16. Muscle
        let muscle = SkeletalMuscleFascicleOrgan::youthful();
        let mus_res = muscle.simulate_locomotor_contraction(0.6, panc_res.insulin_secretion_rate_pmol_l);
        assert!(mus_res.specific_contractile_force_n_cm2 > 12.0);

        // 17. Skin
        let skin = DermalCutaneousOrgan::youthful();
        let skin_res = skin.simulate_barrier_and_thermoregulation(25.0, 37.0);
        assert!(skin_res.cutometer_elastic_recovery_pct > 65.0);

        // 18. Choroid
        let choroid = ChoroidPlexusCsfOrgan::youthful();
        let chor_res = choroid.simulate_csf_secretion(75.0);
        assert!(chor_res.daily_csf_production_ml_day > 400.0);

        // 19. Visceral Adipose
        let wat = VisceralWhiteAdiposeOrgan::youthful();
        let wat_res = wat.simulate_lipolysis_and_adipokines(0.8, panc_res.insulin_secretion_rate_pmol_l);
        assert!(wat_res.circulating_adiponectin_ug_ml > 10.0);

        // 20. Brown Adipose
        let bat = ThermogenicBrownAdiposeOrgan::youthful();
        let bat_res = bat.simulate_thermogenesis(1.0, 18.0);
        assert!(bat_res.non_shivering_thermogenesis_watts > 15.0);
    }
}
