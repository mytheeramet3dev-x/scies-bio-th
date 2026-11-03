//! # Closed-Loop Whole-Body Multi-Organ Physiological Cycles (`cycles`)
//!
//! Models the 7 major systemic physiological circuits where multiple anatomical organs
//! execute synchronized actions in continuous, closed feedback loops (*Guyton & Hall*):
//! 1. **Cardiovascular-Hemodynamic Cycle** ([`CardiovascularCycle`]): Heart, Arterial Windkessel, Microcirculation, Venous Return, and Baroreceptors.
//! 2. **Cardiopulmonary Ventilation & Gas Exchange Cycle** ([`CardiopulmonaryCycle`]): Alveoli, Heart, Erythrocytes, Brainstem Chemoreceptors, and Acid-Base balance.
//! 3. **Enterohepatic-Metabolic Cycle** ([`EnterohepaticMetabolicCycle`]): Gut, Pancreas, Liver Zones 1–3, Visceral WAT, Spleen, and Bile Acid recycling.
//! 4. **Renal-Cardiovascular Fluid, Electrolyte & RAAS Cycle** ([`RenalCardiovascularRaasCycle`]): Kidney, Liver, Lung ACE, Adrenal, Pituitary AVP, and $\alpha$-Klotho.
//! 5. **Neuroendocrine-Glymphatic Sleep Cycle** ([`NeuroendocrineGlymphaticCycle`]): SCN, Hypothalamus, Pituitary, Adrenal HPA, Choroid Plexus CSF, Cortical AQP4 Glymphatics, and Hippocampal SGZ.
//! 6. **Immune-Lymphatic Surveillance Cycle** ([`ImmuneLymphaticCycle`]): Bone Marrow HSC, Thymus, Spleen, Peripheral Patrolling, and Thoracic Duct.
//! 7. **Energy & Thermoregulatory Cycle** ([`EnergyThermoregulatoryCycle`]): Preoptic Thermostat, BAT UCP1, Muscle Shivering/Irisin, WAT Lipolysis, and Dermal Heat Dissipation.

use super::organ::*;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

// =================================================================================================
// 1. CARDIOVASCULAR-HEMODYNAMIC CYCLE
// =================================================================================================

/// State and closed-loop telemetry of the Cardiovascular-Hemodynamic Circuit.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CardiovascularCycle {
    pub heart_rate_bpm: f64,
    pub stroke_volume_ml: f64,
    pub cardiac_output_l_min: f64,
    pub ejection_fraction_pct: f64,
    pub systolic_bp_mmhg: f64,
    pub diastolic_bp_mmhg: f64,
    pub mean_arterial_pressure_mmhg: f64,
    pub central_venous_pressure_mmhg: f64,
    pub systemic_vascular_resistance_wood: f64,
    pub baroreceptor_firing_rate_hz: f64,
    pub sympathetic_tone: f64,
    pub parasympathetic_tone: f64,
}

impl Default for CardiovascularCycle {
    fn default() -> Self {
        Self::new()
    }
}

impl CardiovascularCycle {
    pub fn new() -> Self {
        Self {
            heart_rate_bpm: 72.0,
            stroke_volume_ml: 75.0,
            cardiac_output_l_min: 5.4,
            ejection_fraction_pct: 62.5,
            systolic_bp_mmhg: 118.0,
            diastolic_bp_mmhg: 76.0,
            mean_arterial_pressure_mmhg: 90.0,
            central_venous_pressure_mmhg: 6.0,
            systemic_vascular_resistance_wood: 16.5,
            baroreceptor_firing_rate_hz: 50.0,
            sympathetic_tone: 0.35,
            parasympathetic_tone: 0.65,
        }
    }

    /// Advances one step through the closed-loop cardiovascular circuit:
    /// Venous Return $\to$ Ventricular End-Diastolic Filling $\to$ Frank-Starling Stroke Volume $\to$
    /// Arterial Windkessel Expansion $\to$ Systolic/Diastolic BP $\to$ Carotid Baroreceptor Reflex.
    pub fn step(
        &mut self,
        heart: &MyocardialVentricularOrgan,
        artery: &ArterialVascularOrgan,
        blood_volume_expansion_factor: f64,
    ) {
        // 1. Venous filling pressure driven by central venous volume and tone
        let cvp = (5.0 * blood_volume_expansion_factor.clamp(0.7, 1.4) + 2.5 * self.sympathetic_tone)
            .clamp(2.0, 16.0);
        self.central_venous_pressure_mmhg = cvp;

        // 2. Frank-Starling Ventricular Beat
        let afterload = (self.mean_arterial_pressure_mmhg / 90.0).clamp(0.7, 2.0);
        let beat = heart.simulate_cardiac_cycle(self.heart_rate_bpm, cvp, afterload);
        self.stroke_volume_ml = beat.stroke_volume_ml;
        self.ejection_fraction_pct = beat.ejection_fraction_pct;
        self.cardiac_output_l_min = beat.cardiac_output_l_min;

        // 3. Arterial Hemodynamics & Windkessel
        let hemo = artery.simulate_arterial_hemodynamics(
            self.stroke_volume_ml,
            self.cardiac_output_l_min,
            self.systemic_vascular_resistance_wood,
        );
        self.systolic_bp_mmhg = hemo.systolic_bp_mmhg;
        self.diastolic_bp_mmhg = hemo.diastolic_bp_mmhg;
        self.mean_arterial_pressure_mmhg = hemo.mean_arterial_pressure_mmhg;

        // 4. Carotid Sinus / Aortic Arch Baroreceptor Reflex
        // Baroreceptors fire at higher frequency when MAP rises, suppressing sympathetic tone
        let setpoint_map = 93.0;
        let delta_p = self.mean_arterial_pressure_mmhg - setpoint_map;
        self.baroreceptor_firing_rate_hz = (50.0 + 1.8 * delta_p).clamp(10.0, 100.0);

        // Baroreflex adjustment of autonomic tone
        let symp = (0.35 - 0.008 * delta_p).clamp(0.10, 0.85);
        let vagal = (0.65 + 0.008 * delta_p).clamp(0.15, 0.90);
        self.sympathetic_tone = symp;
        self.parasympathetic_tone = vagal;

        // Heart rate and SVR modulated by baroreflex autonomic balance + microvascular endothelial tone
        self.heart_rate_bpm = (72.0 + 35.0 * (symp - 0.35) - 25.0 * (vagal - 0.65)).clamp(48.0, 160.0);
        let enos_deficit = (0.85 - artery.vascular.enos_activity).max(0.0);
        let rarefaction = (0.85 - artery.vascular.capillary_density).max(0.0);
        self.systemic_vascular_resistance_wood =
            (16.5 + 8.0 * (symp - 0.35) + 8.5 * enos_deficit + 6.5 * rarefaction).clamp(10.0, 28.0);
    }
}

// =================================================================================================
// 2. CARDIOPULMONARY VENTILATION & GAS EXCHANGE CYCLE
// =================================================================================================

/// State and closed-loop telemetry of the Cardiopulmonary Gas Exchange Circuit.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CardiopulmonaryCycle {
    pub respiratory_rate_bpm: f64,
    pub tidal_volume_ml: f64,
    pub minute_ventilation_l_min: f64,
    pub alveolar_po2_mmhg: f64,
    pub arterial_po2_mmhg: f64,
    pub arterial_pco2_mmhg: f64,
    pub arterial_o2_saturation_pct: f64,
    pub arterial_ph: f64,
    pub whole_body_vo2_consumption_ml_min: f64,
    pub whole_body_vco2_production_ml_min: f64,
}

impl Default for CardiopulmonaryCycle {
    fn default() -> Self {
        Self::new()
    }
}

impl CardiopulmonaryCycle {
    pub fn new() -> Self {
        Self {
            respiratory_rate_bpm: 14.0,
            tidal_volume_ml: 500.0,
            minute_ventilation_l_min: 7.0,
            alveolar_po2_mmhg: 102.0,
            arterial_po2_mmhg: 96.0,
            arterial_pco2_mmhg: 40.0,
            arterial_o2_saturation_pct: 98.2,
            arterial_ph: 7.40,
            whole_body_vo2_consumption_ml_min: 250.0,
            whole_body_vco2_production_ml_min: 200.0,
        }
    }

    /// Advances one step through the cardiopulmonary gas exchange circuit:
    /// Peripheral Metabolic Oxygen Consumption $\to$ Venous Return to Pulmonary Artery $\to$
    /// Alveolar Membrane Diffusion $\to$ Arterial Saturation & $P_a\text{CO}_2$ Exhalation $\to$
    /// Medullary Central Chemoreceptor sensing of CSF pH to adjust Ventilation.
    pub fn step(
        &mut self,
        lung: &PulmonaryAlveolarOrgan,
        cardiac_output_l_min: f64,
        metabolic_intensity_factor: f64,
    ) {
        // 1. Metabolic VO2 and VCO2 based on whole-body activity
        let vo2 = 250.0 * metabolic_intensity_factor.clamp(0.6, 4.0);
        let vco2 = vo2 * 0.82; // Respiratory Exchange Ratio R ~ 0.82
        self.whole_body_vo2_consumption_ml_min = vo2;
        self.whole_body_vco2_production_ml_min = vco2;

        // 2. Alveolar gas exchange & Henderson-Hasselbalch arterial pH
        let gas = lung.simulate_respiratory_gas_exchange(
            self.tidal_volume_ml,
            self.respiratory_rate_bpm,
            0.21, // Room air FiO2
            cardiac_output_l_min,
        );
        self.minute_ventilation_l_min = gas.minute_ventilation_l_min;
        self.alveolar_po2_mmhg = gas.alveolar_po2_mmhg;
        self.arterial_po2_mmhg = gas.arterial_po2_mmhg;
        self.arterial_pco2_mmhg = gas.arterial_pco2_mmhg;
        self.arterial_o2_saturation_pct = gas.arterial_o2_saturation_pct;
        self.arterial_ph = gas.arterial_ph;

        // 3. Brainstem Medullary Central Chemoreceptors
        // Sensing PaCO2 and CSF [H+]: hypercapnia (PaCO2 > 40) or acidemia stimulates ventilation
        let co2_error = self.arterial_pco2_mmhg - 40.0;
        let ph_error = 7.40 - self.arterial_ph;

        let drive = (1.0 + 0.025 * co2_error + 1.2 * ph_error).clamp(0.85, 1.35);
        let target_rr = 14.0 * drive;
        let target_vt = 500.0 * drive;
        self.respiratory_rate_bpm = (0.75 * self.respiratory_rate_bpm + 0.25 * target_rr).clamp(10.0, 24.0);
        self.tidal_volume_ml = (0.75 * self.tidal_volume_ml + 0.25 * target_vt).clamp(420.0, 650.0);
    }
}

// =================================================================================================
// 3. ENTEROHEPATIC-METABOLIC CYCLE
// =================================================================================================

/// State and closed-loop telemetry of the Enterohepatic-Metabolic Circuit.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct EnterohepaticMetabolicCycle {
    pub blood_glucose_mg_dl: f64,
    pub circulating_insulin_pmol_l: f64,
    pub circulating_glucagon_pg_ml: f64,
    pub glp1_incretin_pmol_l: f64,
    pub portal_blood_flow_l_min: f64,
    pub portal_lps_endotoxin_eu_ml: f64,
    pub hepatic_urea_cycle_efficiency_pct: f64,
    pub cyp_detoxification_capacity: f64,
    pub hepatic_glucose_output_mg_min: f64,
    pub bile_acid_pool_circulation_rate: f64,
    pub systemic_ammonia_umol_l: f64,
}

impl Default for EnterohepaticMetabolicCycle {
    fn default() -> Self {
        Self::new()
    }
}

impl EnterohepaticMetabolicCycle {
    pub fn new() -> Self {
        Self {
            blood_glucose_mg_dl: 92.0,
            circulating_insulin_pmol_l: 45.0,
            circulating_glucagon_pg_ml: 65.0,
            glp1_incretin_pmol_l: 14.0,
            portal_blood_flow_l_min: 1.1,
            portal_lps_endotoxin_eu_ml: 0.18,
            hepatic_urea_cycle_efficiency_pct: 95.0,
            cyp_detoxification_capacity: 1.05,
            hepatic_glucose_output_mg_min: 140.0,
            bile_acid_pool_circulation_rate: 6.0, // Re-circulated ~6 times/day
            systemic_ammonia_umol_l: 25.0,
        }
    }

    /// Advances one step through the enterohepatic metabolic circuit:
    /// Ingestion & Gut Absorption $\to$ Incretin GLP-1 $\to$ Pancreatic Biphasic Insulin $\to$
    /// Portal Vein Inflow into Liver Lobule $\to$ Zone 1-3 Glycogen/Detox $\to$ Bile Acid Duodenal
    /// Secretion & 95% Ileal Reabsorption $\to$ Visceral Adipose FFA modulation.
    pub fn step(
        &mut self,
        gut: &GastrointestinalCryptOrgan,
        pancreas: &PancreaticIsletOrgan,
        liver: &HepaticLobuleOrgan,
        caloric_intake_kcal: f64,
        dietary_fiber_g: f64,
    ) {
        // 1. Gastrointestinal absorption & incretin release
        let gut_res = gut.simulate_digestion_and_barrier(caloric_intake_kcal, dietary_fiber_g);
        self.glp1_incretin_pmol_l = gut_res.glp1_incretin_secretion_pmol_l;
        self.portal_lps_endotoxin_eu_ml = gut_res.portal_lps_endotoxin_translocation_eu_ml;

        // 2. Glucose influx from meal
        let meal_glucose_influx = if caloric_intake_kcal > 0.0 {
            (caloric_intake_kcal * 0.15).clamp(0.0, 120.0)
        } else {
            0.0
        };

        // 3. Pancreatic islet response to glucose and incretin
        let panc_res = pancreas.simulate_islet_secretion(self.blood_glucose_mg_dl, self.glp1_incretin_pmol_l);
        self.circulating_insulin_pmol_l = panc_res.insulin_secretion_rate_pmol_l;
        self.circulating_glucagon_pg_ml = panc_res.glucagon_secretion_rate_pg_ml;

        // 4. Hepatic zonated metabolism & urea cycle
        let hep_res = liver.simulate_zonated_metabolism(
            1.0,
            self.circulating_insulin_pmol_l,
            1.0,
            self.systemic_ammonia_umol_l,
        );
        self.hepatic_urea_cycle_efficiency_pct = hep_res.urea_cycle_clearance_efficiency_pct;
        self.cyp_detoxification_capacity = hep_res.cyp_detoxification_capacity_multiplier;
        self.hepatic_glucose_output_mg_min = hep_res.hepatic_glucose_output_mg_min;

        // 5. Systemic glucose mass balance with euglycemic glucagon/insulin feedback
        let glucose_ratio = (self.blood_glucose_mg_dl / 92.0).clamp(0.6, 2.0);
        let peripheral_disposal =
            (105.0 * (self.circulating_insulin_pmol_l / 50.0) * glucose_ratio).clamp(35.0, 320.0);
        let hgo_effective =
            self.hepatic_glucose_output_mg_min * (1.0 + 1.2 * (1.0 - glucose_ratio).max(0.0));
        let net_glucose_flux = meal_glucose_influx + hgo_effective - peripheral_disposal;
        self.blood_glucose_mg_dl = (self.blood_glucose_mg_dl + net_glucose_flux * 0.05).clamp(72.0, 240.0);

        // Ammonia cleared by functional Zone 1 urea cycle
        self.systemic_ammonia_umol_l = (22.0 * (100.0 / self.hepatic_urea_cycle_efficiency_pct)).clamp(15.0, 110.0);
    }
}

// =================================================================================================
// 4. RENAL-CARDIOVASCULAR FLUID, ELECTROLYTE & RAAS CYCLE
// =================================================================================================

/// State and closed-loop telemetry of the Renal-Cardiovascular RAAS Homeostatic Circuit.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct RenalCardiovascularRaasCycle {
    pub egfr_ml_min: f64,
    pub renal_blood_flow_l_min: f64,
    pub plasma_renin_activity_ng_ml_hr: f64,
    pub circulating_angiotensin_ii_pmol_l: f64,
    pub plasma_aldosterone_ng_dl: f64,
    pub extracellular_fluid_volume_l: f64,
    pub circulating_alpha_klotho_pg_ml: f64,
    pub serum_creatinine_mg_dl: f64,
}

impl Default for RenalCardiovascularRaasCycle {
    fn default() -> Self {
        Self::new()
    }
}

impl RenalCardiovascularRaasCycle {
    pub fn new() -> Self {
        Self {
            egfr_ml_min: 105.0,
            renal_blood_flow_l_min: 1.15,
            plasma_renin_activity_ng_ml_hr: 1.2,
            circulating_angiotensin_ii_pmol_l: 18.0,
            plasma_aldosterone_ng_dl: 8.5,
            extracellular_fluid_volume_l: 15.0,
            circulating_alpha_klotho_pg_ml: 950.0,
            serum_creatinine_mg_dl: 0.85,
        }
    }

    /// Advances one step through the renal-cardiovascular RAAS circuit:
    /// Perfusion Pressure $\to$ Glomerular Filtration (eGFR) $\to$ Juxtaglomerular Renin $\to$
    /// Angiotensin II $\to$ Adrenal Aldosterone $\to$ Tubular Sodium & Water Retention $\to$
    /// ECFV Expansion & BP restoration $\to$ Distal Tubule $\alpha$-Klotho Longevity Secretion.
    pub fn step(
        &mut self,
        kidney: &RenalNephronOrgan,
        adrenal: &AdrenalSuprarenalOrgan,
        mean_arterial_pressure_mmhg: f64,
        cardiac_output_l_min: f64,
        barotrauma_index: f64,
    ) {
        // 1. Renal Glomerular Filtration & Endocrine Klotho/Renin Release
        let renal_res = kidney.simulate_renal_filtration(
            mean_arterial_pressure_mmhg,
            cardiac_output_l_min,
            barotrauma_index,
        );
        self.egfr_ml_min = renal_res.egfr_ml_min;
        self.renal_blood_flow_l_min = renal_res.renal_blood_flow_l_min;
        self.circulating_alpha_klotho_pg_ml = renal_res.alpha_klotho_secretion_pg_ml;

        // 2. RAAS Cascade: Renin converts Angiotensinogen -> Ang I -> Pulmonary ACE -> Ang II
        self.plasma_renin_activity_ng_ml_hr = (1.2 * renal_res.renin_release_rate).clamp(0.2, 8.5);
        self.circulating_angiotensin_ii_pmol_l =
            (15.0 * self.plasma_renin_activity_ng_ml_hr).clamp(4.0, 110.0);

        // 3. Adrenal Aldosterone response to Angiotensin II
        let adr_res = adrenal.simulate_steroidogenesis(
            22.0, // Baseline ACTH
            self.circulating_angiotensin_ii_pmol_l / 25.0,
            0.3,
        );
        self.plasma_aldosterone_ng_dl = adr_res.aldosterone_ng_dl;

        // 4. ECFV regulation: Aldosterone promotes Na+ and H2O reabsorption
        let aldo_factor = (self.plasma_aldosterone_ng_dl / 8.5).clamp(0.4, 3.5);
        let gfr_factor = (self.egfr_ml_min / 100.0).clamp(0.2, 1.4);
        let target_ecfv = 15.0 * (1.0 + 0.12 * (aldo_factor - 1.0) - 0.08 * (gfr_factor - 1.0));
        self.extracellular_fluid_volume_l = (0.92 * self.extracellular_fluid_volume_l + 0.08 * target_ecfv)
            .clamp(11.0, 20.0);

        // Serum creatinine inversely related to GFR: Cr = k / GFR
        self.serum_creatinine_mg_dl = (90.0 / self.egfr_ml_min).clamp(0.55, 4.5);
    }
}

// =================================================================================================
// 5. NEUROENDOCRINE-GLYMPHATIC SLEEP CYCLE
// =================================================================================================

/// State and closed-loop telemetry of the Neuroendocrine-Glymphatic Sleep & Longevity Circuit.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct NeuroendocrineGlymphaticCycle {
    pub circadian_hour: f64,
    pub is_slow_wave_sleep: bool,
    pub hypothalamic_gnrh_level: f64,
    pub plasma_cortisol_ug_dl: f64,
    pub plasma_dhea_s_ug_dl: f64,
    pub catabolic_cortisol_to_dhea_ratio: f64,
    pub daily_csf_turnover_ml_day: f64,
    pub cortical_glymphatic_clearance_rate: f64,
    pub adult_hippocampal_neurogenesis_cells_day: f64,
}

impl Default for NeuroendocrineGlymphaticCycle {
    fn default() -> Self {
        Self::new()
    }
}

impl NeuroendocrineGlymphaticCycle {
    pub fn new() -> Self {
        Self {
            circadian_hour: 8.0,
            is_slow_wave_sleep: false,
            hypothalamic_gnrh_level: 0.95,
            plasma_cortisol_ug_dl: 14.0,
            plasma_dhea_s_ug_dl: 260.0,
            catabolic_cortisol_to_dhea_ratio: 0.054,
            daily_csf_turnover_ml_day: 520.0,
            cortical_glymphatic_clearance_rate: 8.5,
            adult_hippocampal_neurogenesis_cells_day: 750.0,
        }
    }

    /// Advances one step through the neuroendocrine-glymphatic sleep circuit:
    /// SCN Circadian Pacing $\to$ Hypothalamic Releasing Hormones $\to$ Pituitary ACTH $\to$
    /// Adrenal Cortisol/DHEA $\to$ Deep Sleep Phase Shift $\to$ Choroid Plexus CSF Outflow $\to$
    /// Astrocytic AQP4 Glymphatic flushing of Amyloid/Tau $\to$ Hippocampal Neurogenesis Protection.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &mut self,
        hypo: &HypothalamicNeuroendocrineOrgan,
        pit: &PituitaryHypophysealOrgan,
        adrenal: &AdrenalSuprarenalOrgan,
        choroid: &ChoroidPlexusCsfOrgan,
        cortex: &CerebralCortexOrgan,
        hippo: &HippocampalDentateGyrusOrgan,
        dt_hours: f64,
        systemic_sasp: f64,
    ) {
        self.circadian_hour = (self.circadian_hour + dt_hours) % 24.0;
        self.is_slow_wave_sleep = self.circadian_hour >= 23.0 || self.circadian_hour <= 5.0;
        let sleep_fraction = if self.is_slow_wave_sleep { 0.95 } else { 0.05 };

        // 1. Hypothalamic master regulation
        let reg = hypo.simulate_master_regulation(self.circadian_hour, systemic_sasp);
        self.hypothalamic_gnrh_level = reg.gnrh_pulsatile_release;

        // 2. Pituitary ACTH & GH cascade
        let pit_res = pit.simulate_hypophyseal_cascade(
            reg.crh_drive,
            reg.ghrh_drive,
            self.catabolic_cortisol_to_dhea_ratio,
            1.0,
        );

        // 3. Adrenal Cortisol & DHEA-S steroidogenesis
        let adr_res = adrenal.simulate_steroidogenesis(pit_res.acth_pg_ml, 0.25, 0.35);
        self.plasma_cortisol_ug_dl = adr_res.cortisol_ug_dl;
        self.plasma_dhea_s_ug_dl = adr_res.dhea_s_ug_dl;
        self.catabolic_cortisol_to_dhea_ratio = adr_res.cortisol_to_dhea_ratio;

        // 4. Choroid Plexus active CSF secretion
        let chor_res = choroid.simulate_csf_secretion(75.0);
        self.daily_csf_turnover_ml_day = chor_res.daily_csf_production_ml_day;

        // 5. Cortical AQP4 Glymphatic interstitial flushing during sleep
        let cort_res = cortex.simulate_cortical_computation(90.0, sleep_fraction, 0.0);
        self.cortical_glymphatic_clearance_rate = cort_res.glymphatic_waste_clearance_rate;

        // 6. Hippocampal SGZ adult neurogenesis
        let hip_res = hippo.simulate_neurogenesis_and_memory(
            chor_res.ventricular_alpha_klotho_pg_ml / 500.0,
            self.catabolic_cortisol_to_dhea_ratio,
            25.0,
        );
        self.adult_hippocampal_neurogenesis_cells_day = hip_res.adult_neurogenesis_rate_cells_day;
    }
}

// =================================================================================================
// 6. IMMUNE-LYMPHATIC SURVEILLANCE CYCLE
// =================================================================================================

/// State and closed-loop telemetry of the Immune-Lymphatic Surveillance Circuit.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ImmuneLymphaticCycle {
    pub daily_naive_t_cell_emigration_10e7: f64,
    pub myeloid_to_lymphoid_skew_ratio: f64,
    pub splenic_senescent_rbc_clearance_10e9: f64,
    pub thoracic_duct_lymph_return_ml_hr: f64,
    pub peripheral_immune_patrolling_competence: f64,
}

impl Default for ImmuneLymphaticCycle {
    fn default() -> Self {
        Self::new()
    }
}

impl ImmuneLymphaticCycle {
    pub fn new() -> Self {
        Self {
            daily_naive_t_cell_emigration_10e7: 12.0,
            myeloid_to_lymphoid_skew_ratio: 1.25,
            splenic_senescent_rbc_clearance_10e9: 200.0,
            thoracic_duct_lymph_return_ml_hr: 120.0,
            peripheral_immune_patrolling_competence: 0.94,
        }
    }

    /// Advances one step through the immune-lymphatic circuit:
    /// Bone Marrow HSC & Progenitor Output $\to$ Thymic Corticomedullary Selection $\to$
    /// Splenic Red Pulp Filtration & White Pulp Screening $\to$ Peripheral Tissue Patrolling $\to$
    /// Lymphatic Drainage via Thoracic Duct into Venous Circulation.
    pub fn step(
        &mut self,
        marrow: &BoneMarrowNicheOrgan,
        thymus: &ThymicLymphoidOrgan,
        spleen: &SplenicLymphoidPulpOrgan,
    ) {
        // 1. Bone marrow lymphopoiesis vs myelopoiesis
        let marrow_res = marrow.simulate_hematopoiesis(1.0, 0.5);
        self.myeloid_to_lymphoid_skew_ratio = marrow_res.myeloid_to_lymphoid_skewing_ratio;

        // 2. Thymic T-cell receptor selection and naive emigration
        let lymph_influx = 1.0 / self.myeloid_to_lymphoid_skew_ratio.sqrt();
        let thym_res = thymus.simulate_thymic_selection(lymph_influx);
        self.daily_naive_t_cell_emigration_10e7 = thym_res.daily_naive_t_cell_output_10e7;

        // 3. Splenic erythrocyte culling and marginal zone immunity
        let spl_res = spleen.simulate_splenic_culling(0.25);
        self.splenic_senescent_rbc_clearance_10e9 = spl_res.senescent_rbc_clearance_10e9;

        // 4. Lymphatic return and systemic patrolling
        self.peripheral_immune_patrolling_competence =
            (0.50 * thym_res.negative_selection_competence + 0.50 * spl_res.macrophage_phagocytic_clearance)
                .clamp(0.15, 0.98);
        self.thoracic_duct_lymph_return_ml_hr = (120.0 * (marrow.vascular.capillary_density / 0.85))
            .clamp(50.0, 160.0);
    }
}

// =================================================================================================
// 7. ENERGY & THERMOREGULATORY CYCLE
// =================================================================================================

/// State and closed-loop telemetry of the Energy Expenditure & Thermoregulatory Circuit.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct EnergyThermoregulatoryCycle {
    pub core_body_temperature_c: f64,
    pub basal_metabolic_rate_kcal_day: f64,
    pub brown_adipose_heat_watts: f64,
    pub muscle_thermogenesis_watts: f64,
    pub cutaneous_heat_loss_watts: f64,
    pub circulating_irisin_ng_ml: f64,
    pub circulating_adiponectin_ug_ml: f64,
}

impl Default for EnergyThermoregulatoryCycle {
    fn default() -> Self {
        Self::new()
    }
}

impl EnergyThermoregulatoryCycle {
    pub fn new() -> Self {
        Self {
            core_body_temperature_c: 37.0,
            basal_metabolic_rate_kcal_day: 1650.0,
            brown_adipose_heat_watts: 18.0,
            muscle_thermogenesis_watts: 25.0,
            cutaneous_heat_loss_watts: 43.0,
            circulating_irisin_ng_ml: 12.0,
            circulating_adiponectin_ug_ml: 14.5,
        }
    }

    /// Advances one step through the energy and thermoregulatory circuit:
    /// Preoptic Core Temperature Sensing $\to$ Sympathetic BAT UCP1 Thermogenesis $\to$
    /// Skeletal Muscle Heat & Irisin Secretion $\to$ WAT Lipolysis $\to$ Cutaneous Vasomotor Dissipation.
    pub fn step(
        &mut self,
        bat: &ThermogenicBrownAdiposeOrgan,
        muscle: &SkeletalMuscleFascicleOrgan,
        wat: &VisceralWhiteAdiposeOrgan,
        skin: &DermalCutaneousOrgan,
        ambient_temp_c: f64,
        locomotor_work_fraction: f64,
    ) {
        let delta_core = (37.0 - self.core_body_temperature_c).clamp(-1.5, 2.0);
        let cold_stim = (delta_core * 1.5).max(0.0);

        // 1. Brown adipose non-shivering UCP1 thermogenesis
        let bat_res = bat.simulate_thermogenesis(1.0 + cold_stim, ambient_temp_c);
        self.brown_adipose_heat_watts = bat_res.non_shivering_thermogenesis_watts;

        // 2. Skeletal muscle locomotion & shivering heat
        let mus_res = muscle.simulate_locomotor_contraction(locomotor_work_fraction + 0.3 * cold_stim, 50.0);
        self.muscle_thermogenesis_watts = mus_res.mechanical_heat_generation_watts;
        self.circulating_irisin_ng_ml = mus_res.irisin_secretion_ng_ml;

        // 3. Visceral WAT lipolysis and adipokines
        let wat_res = wat.simulate_lipolysis_and_adipokines(1.0 + 0.5 * cold_stim, 50.0);
        self.circulating_adiponectin_ug_ml = wat_res.circulating_adiponectin_ug_ml;

        // 4. Cutaneous dermal vasomotor heat loss modulated by preoptic negative feedback
        let temp_error = self.core_body_temperature_c - 37.0;
        let skin_res = skin.simulate_barrier_and_thermoregulation(ambient_temp_c, self.core_body_temperature_c);
        let baseline_heat_gen = self.brown_adipose_heat_watts + self.muscle_thermogenesis_watts;
        let heat_loss = (baseline_heat_gen
            + 25.0 * temp_error
            + 0.04 * (skin_res.cutaneous_heat_dissipation_flow_ml_min - 120.0))
            .clamp(20.0, 180.0);
        self.cutaneous_heat_loss_watts = heat_loss;

        // Thermal balance: Core temperature negative feedback
        let net_heat = baseline_heat_gen - self.cutaneous_heat_loss_watts;
        self.core_body_temperature_c = (self.core_body_temperature_c + net_heat * 0.02).clamp(36.0, 38.5);

        // Basal Metabolic Rate (BMR) in kcal/day
        self.basal_metabolic_rate_kcal_day =
            (1450.0 + 4.5 * self.brown_adipose_heat_watts + 6.0 * self.muscle_thermogenesis_watts).clamp(1200.0, 2600.0);
    }
}

// =================================================================================================
// WHOLE-BODY SYSTEMIC CIRCUITS CONTAINER
// =================================================================================================

/// Complete Whole-Body Multi-Organ Systemic Circuits orchestrating all 7 physiological cycles.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct WholeBodySystemicCircuits {
    pub cardiovascular: CardiovascularCycle,
    pub cardiopulmonary: CardiopulmonaryCycle,
    pub enterohepatic: EnterohepaticMetabolicCycle,
    pub renal_raas: RenalCardiovascularRaasCycle,
    pub neuroendocrine_glymphatic: NeuroendocrineGlymphaticCycle,
    pub immune_lymphatic: ImmuneLymphaticCycle,
    pub energy_thermoregulatory: EnergyThermoregulatoryCycle,
}

impl Default for WholeBodySystemicCircuits {
    fn default() -> Self {
        Self::new()
    }
}

impl WholeBodySystemicCircuits {
    /// Initializes all 7 physiological cycles to youthful homeostatic baselines.
    pub fn new() -> Self {
        Self {
            cardiovascular: CardiovascularCycle::new(),
            cardiopulmonary: CardiopulmonaryCycle::new(),
            enterohepatic: EnterohepaticMetabolicCycle::new(),
            renal_raas: RenalCardiovascularRaasCycle::new(),
            neuroendocrine_glymphatic: NeuroendocrineGlymphaticCycle::new(),
            immune_lymphatic: ImmuneLymphaticCycle::new(),
            energy_thermoregulatory: EnergyThermoregulatoryCycle::new(),
        }
    }

    /// Simulates a synchronized step of all 7 physiological cycles coupled across all 20 organs.
    pub fn step_synchronized(
        &mut self,
        atlas: &AnatomicalSpatialAtlas,
        dt_hours: f64,
        ambient_temp_c: f64,
        caloric_intake_kcal: f64,
        systemic_sasp: f64,
    ) {
        // Extract organs from atlas (indices 0..20)
        let artery = match &atlas.organs[0] { AnatomicalOrgan::ArterialVascular(o) => o, _ => unreachable!() };
        let cortex = match &atlas.organs[1] { AnatomicalOrgan::CerebralCortex(o) => o, _ => unreachable!() };
        let liver = match &atlas.organs[2] { AnatomicalOrgan::HepaticLobule(o) => o, _ => unreachable!() };
        let marrow = match &atlas.organs[3] { AnatomicalOrgan::BoneMarrowNiche(o) => o, _ => unreachable!() };
        let thymus = match &atlas.organs[4] { AnatomicalOrgan::ThymicLymphoid(o) => o, _ => unreachable!() };
        let skin = match &atlas.organs[5] { AnatomicalOrgan::DermalCutaneous(o) => o, _ => unreachable!() };
        let heart = match &atlas.organs[6] { AnatomicalOrgan::MyocardialVentricular(o) => o, _ => unreachable!() };
        let kidney = match &atlas.organs[7] { AnatomicalOrgan::RenalNephron(o) => o, _ => unreachable!() };
        let lung = match &atlas.organs[8] { AnatomicalOrgan::PulmonaryAlveolar(o) => o, _ => unreachable!() };
        let pancreas = match &atlas.organs[9] { AnatomicalOrgan::PancreaticIslet(o) => o, _ => unreachable!() };
        let muscle = match &atlas.organs[10] { AnatomicalOrgan::SkeletalMuscleFascicle(o) => o, _ => unreachable!() };
        let gut = match &atlas.organs[11] { AnatomicalOrgan::GastrointestinalCrypt(o) => o, _ => unreachable!() };
        let hypo = match &atlas.organs[12] { AnatomicalOrgan::HypothalamicNeuroendocrine(o) => o, _ => unreachable!() };
        let hippo = match &atlas.organs[13] { AnatomicalOrgan::HippocampalDentateGyrus(o) => o, _ => unreachable!() };
        let choroid = match &atlas.organs[14] { AnatomicalOrgan::ChoroidPlexusCsf(o) => o, _ => unreachable!() };
        let pit = match &atlas.organs[15] { AnatomicalOrgan::PituitaryHypophyseal(o) => o, _ => unreachable!() };
        let adrenal = match &atlas.organs[16] { AnatomicalOrgan::AdrenalSuprarenal(o) => o, _ => unreachable!() };
        let spleen = match &atlas.organs[17] { AnatomicalOrgan::SplenicLymphoidPulp(o) => o, _ => unreachable!() };
        let wat = match &atlas.organs[18] { AnatomicalOrgan::VisceralWhiteAdipose(o) => o, _ => unreachable!() };
        let bat = match &atlas.organs[19] { AnatomicalOrgan::ThermogenicBrownAdipose(o) => o, _ => unreachable!() };

        // 1. Cardiovascular cycle
        let ecfv_expansion = self.renal_raas.extracellular_fluid_volume_l / 15.0;
        self.cardiovascular.step(heart, artery, ecfv_expansion);

        // 2. Cardiopulmonary cycle
        self.cardiopulmonary.step(lung, self.cardiovascular.cardiac_output_l_min, 1.0);

        // 3. Enterohepatic metabolic cycle
        self.enterohepatic.step(gut, pancreas, liver, caloric_intake_kcal, 30.0);

        // 4. Renal RAAS cycle
        let barotrauma = (self.cardiovascular.mean_arterial_pressure_mmhg - 100.0).max(0.0) / 40.0;
        self.renal_raas.step(
            kidney,
            adrenal,
            self.cardiovascular.mean_arterial_pressure_mmhg,
            self.cardiovascular.cardiac_output_l_min,
            barotrauma.clamp(0.0, 1.0),
        );

        // 5. Neuroendocrine glymphatic sleep cycle
        self.neuroendocrine_glymphatic.step(
            hypo, pit, adrenal, choroid, cortex, hippo, dt_hours, systemic_sasp,
        );

        // 6. Immune lymphatic cycle
        self.immune_lymphatic.step(marrow, thymus, spleen);

        // 7. Energy & thermoregulatory cycle
        self.energy_thermoregulatory.step(bat, muscle, wat, skin, ambient_temp_c, 0.4);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_all_7_physiological_cycles_synchronized() {
        let atlas = AnatomicalSpatialAtlas::from_somatic_age(25.0);
        let mut circuits = WholeBodySystemicCircuits::new();

        // Run 24 hourly steps (a full circadian day)
        for h in 0..24 {
            circuits.step_synchronized(&atlas, 1.0, 22.0, if h == 8 || h == 13 || h == 19 { 600.0 } else { 0.0 }, 0.05);
        }

        // Verify cardiovascular homeostasis
        assert!(circuits.cardiovascular.mean_arterial_pressure_mmhg > 75.0 && circuits.cardiovascular.mean_arterial_pressure_mmhg < 115.0);
        assert!(circuits.cardiovascular.ejection_fraction_pct > 50.0);

        // Verify cardiopulmonary arterial blood gas
        assert!(circuits.cardiopulmonary.arterial_po2_mmhg > 85.0);
        assert!(circuits.cardiopulmonary.arterial_ph > 7.35 && circuits.cardiopulmonary.arterial_ph < 7.45);

        // Verify enterohepatic glucose control
        assert!(circuits.enterohepatic.blood_glucose_mg_dl > 70.0 && circuits.enterohepatic.blood_glucose_mg_dl < 140.0);

        // Verify renal eGFR and Klotho
        assert!(circuits.renal_raas.egfr_ml_min > 85.0);
        assert!(circuits.renal_raas.circulating_alpha_klotho_pg_ml > 600.0);

        // Verify glymphatic clearance and neurogenesis
        assert!(circuits.neuroendocrine_glymphatic.adult_hippocampal_neurogenesis_cells_day > 400.0);

        // Verify immune patrolling and thermoregulation
        assert!(circuits.immune_lymphatic.peripheral_immune_patrolling_competence > 0.70);
        assert!(circuits.energy_thermoregulatory.core_body_temperature_c > 36.5 && circuits.energy_thermoregulatory.core_body_temperature_c < 37.5);
    }
}
