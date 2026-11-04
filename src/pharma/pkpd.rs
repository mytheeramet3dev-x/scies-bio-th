//! One-compartment Pharmacokinetic (PK) and Hill-equation Pharmacodynamic (PD) modeling.
//!
//! Models drug absorption, elimination half-life ($t_{1/2}$), concentration-dependent biological
//! efficacy ($E_{\max}$ Hill equation), and Therapeutic Index ($\text{TI} = \text{TD}_{50} / \text{EC}_{50}$).

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Pharmacokinetic and pharmacodynamic parameters for a small molecule or biologics intervention.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PkPdParameters {
    /// Compound name (e.g. `"Navitoclax_BH3"`, `"Doxycycline_OSK_Inducer"`).
    pub compound_name: String,
    /// Half-maximal effective concentration $\text{EC}_{50}$ ($\mu\text{M}$).
    pub ec50_um: f64,
    /// Median toxic concentration $\text{TD}_{50}$ ($\mu\text{M}$).
    pub td50_um: f64,
    /// Hill cooperativity coefficient $n_H$ (typically $1.0\text{--}3.0$).
    pub hill_coefficient: f64,
    /// Elimination half-life $t_{1/2}$ in hours.
    pub half_life_hours: f64,
    /// Bioavailability fraction $F \in (0.0, 1.0]$.
    pub bioavailability: f64,
}

impl PkPdParameters {
    /// First-order elimination rate constant $k_e = \frac{\ln 2}{t_{1/2}}$ ($\text{hr}^{-1}$).
    pub fn elimination_rate_ke(&self) -> f64 {
        std::f64::consts::LN_2 / self.half_life_hours.max(0.1)
    }

    /// Classical Therapeutic Index $\text{TI} = \frac{\text{TD}_{50}}{\text{EC}_{50}}$.
    pub fn therapeutic_index(&self) -> f64 {
        self.td50_um / self.ec50_um.max(1e-9)
    }
}

/// Evaluates the Hill equation fractional response $E(C) = \frac{C^n}{\text{EC}_{50}^n + C^n} \in [0.0, 1.0]$.
pub fn evaluate_hill_response(concentration_um: f64, ec50_um: f64, hill_n: f64) -> f64 {
    let c = concentration_um.max(0.0);
    if c <= 0.0 {
        return 0.0;
    }
    let ec50 = ec50_um.max(1e-9);
    let ratio = c / ec50;
    if ratio >= 1e12 {
        return 1.0;
    }
    if ratio <= 1e-12 {
        return 0.0;
    }
    let n = hill_n.clamp(0.2, 10.0);
    let ratio_n = ratio.powf(n);
    if ratio_n.is_infinite() {
        return 1.0;
    }
    (ratio_n / (1.0 + ratio_n)).clamp(0.0, 1.0)
}

/// Summary report of a multi-dose PK/PD simulation over time.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PkPdRegimenReport {
    /// Peak plasma/media concentration $C_{\max}$ ($\mu\text{M}$).
    pub c_max_um: f64,
    /// Trough concentration prior to final dose $C_{\min}$ ($\mu\text{M}$).
    pub c_trough_um: f64,
    /// Time-averaged therapeutic efficacy in $[0.0, 1.0]$.
    pub mean_efficacy: f64,
    /// Peak off-target toxicity probability in $[0.0, 1.0]$.
    pub peak_toxicity: f64,
    /// True if the regimen maintains $\text{mean\_efficacy} \ge 0.50$ while keeping $\text{peak\_toxicity} < 0.20$.
    pub within_therapeutic_window: bool,
}

/// Simulates a repeated bolus dosing regimen (`dose_um` every `interval_hours` for `num_doses`)
/// at 1-hour resolution.
pub fn simulate_pkpd_regimen(
    params: &PkPdParameters,
    dose_um: f64,
    interval_hours: usize,
    num_doses: usize,
) -> PkPdRegimenReport {
    let interval = interval_hours.max(1);
    let doses = num_doses.max(1);
    let total_hours = interval * doses;
    let ke = params.elimination_rate_ke();
    let decay_per_hour = (-ke).exp();
    let effective_bolus = dose_um.max(0.0) * params.bioavailability.clamp(0.0, 1.0);

    let mut conc = 0.0_f64;
    let mut c_max = 0.0_f64;
    let mut c_trough = f64::INFINITY;
    let mut sum_efficacy = 0.0_f64;
    let mut peak_toxicity = 0.0_f64;

    for hr in 0..total_hours {
        if hr % interval == 0 {
            if hr > 0 && conc < c_trough {
                c_trough = conc;
            }
            conc += effective_bolus;
        }
        if conc > c_max {
            c_max = conc;
        }

        let eff = evaluate_hill_response(conc, params.ec50_um, params.hill_coefficient);
        let tox = evaluate_hill_response(conc, params.td50_um, params.hill_coefficient);
        sum_efficacy += eff;
        if tox > peak_toxicity {
            peak_toxicity = tox;
        }

        conc *= decay_per_hour;
    }

    if !c_trough.is_finite() {
        c_trough = conc;
    }
    let mean_efficacy = sum_efficacy / (total_hours as f64);

    PkPdRegimenReport {
        c_max_um: c_max,
        c_trough_um: c_trough,
        mean_efficacy,
        peak_toxicity,
        within_therapeutic_window: mean_efficacy >= 0.50 && peak_toxicity < 0.20,
    }
}

/// Single scheduled drug administration event in hours and dose amount.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct DosingEvent {
    pub time_hours: f64,
    pub dose_amount: f64,
}

/// Flexible or cyclic intermittent dosing schedule.
#[derive(Debug, Clone, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct IntermittentDosingSchedule {
    pub events: Vec<DosingEvent>,
}

impl IntermittentDosingSchedule {
    /// Creates a regular pulsed schedule (e.g., 5.0 mg every 24 hours for N doses).
    pub fn pulsed(dose: f64, interval_hours: f64, num_doses: usize) -> Self {
        let mut events = Vec::with_capacity(num_doses);
        for i in 0..num_doses {
            events.push(DosingEvent {
                time_hours: i as f64 * interval_hours,
                dose_amount: dose,
            });
        }
        Self { events }
    }

    /// Creates an intermittent cyclic regimen (e.g. on_days daily, off_days off, for N cycles).
    pub fn cyclic_on_off(
        daily_dose: f64,
        on_days: usize,
        off_days: usize,
        num_cycles: usize,
    ) -> Self {
        let mut events = Vec::new();
        let cycle_duration_hours = (on_days + off_days) as f64 * 24.0;
        for c in 0..num_cycles {
            let cycle_start = c as f64 * cycle_duration_hours;
            for d in 0..on_days {
                events.push(DosingEvent {
                    time_hours: cycle_start + d as f64 * 24.0,
                    dose_amount: daily_dose,
                });
            }
        }
        Self { events }
    }
}

/// Two-compartment pharmacokinetic model with central and peripheral distribution.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TwoCompartmentPk {
    pub compound_name: String,
    /// Central distribution volume $V_1$ in Liters.
    pub volume_central_l: f64,
    /// Peripheral distribution volume $V_2$ in Liters.
    pub volume_peripheral_l: f64,
    /// Elimination rate constant from central compartment $k_{10}$ ($\text{hr}^{-1}$).
    pub kel: f64,
    /// Distribution rate constant from central to peripheral $k_{12}$ ($\text{hr}^{-1}$).
    pub k12: f64,
    /// Return rate constant from peripheral to central $k_{21}$ ($\text{hr}^{-1}$).
    pub k21: f64,
    /// Oral/subcutaneous bioavailability fraction $F \in (0.0, 1.0]$.
    pub bioavailability: f64,
}

impl TwoCompartmentPk {
    /// Creates a calibrated two-compartment model from clearance and distribution volumes.
    pub fn new(
        compound_name: impl Into<String>,
        volume_central_l: f64,
        volume_peripheral_l: f64,
        kel: f64,
        k12: f64,
        k21: f64,
        bioavailability: f64,
    ) -> Self {
        Self {
            compound_name: compound_name.into(),
            volume_central_l: volume_central_l.max(0.1),
            volume_peripheral_l: volume_peripheral_l.max(0.1),
            kel: kel.max(1e-4),
            k12: k12.max(0.0),
            k21: k21.max(0.0),
            bioavailability: bioavailability.clamp(0.01, 1.0),
        }
    }

    /// Computes distribution half-life $t_{1/2,\alpha}$ and terminal elimination half-life $t_{1/2,\beta}$ in hours.
    pub fn hybrid_half_lives(&self) -> (f64, f64) {
        let b = self.kel + self.k12 + self.k21;
        let c = self.kel * self.k21;
        let disc = (b * b - 4.0 * c).max(0.0).sqrt();
        let alpha = (b + disc) / 2.0;
        let beta = (b - disc) / 2.0;
        let t_half_alpha = std::f64::consts::LN_2 / alpha.max(1e-4);
        let t_half_beta = std::f64::consts::LN_2 / beta.max(1e-4);
        (t_half_alpha, t_half_beta)
    }
}

/// Pharmacodynamic parameters evaluating the Hill equation with explicit $E_{\max}$ and $\gamma$.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HillPharmacodynamics {
    /// Maximum efficacy attainable $E_{\max} \in [0.0, 1.0]$.
    pub emax: f64,
    /// Potency parameter $\text{EC}_{50}$ in $\mu\text{M}$.
    pub ec50_um: f64,
    /// Hill coefficient / slope parameter $\gamma$.
    pub gamma: f64,
    /// Median toxic concentration $\text{TD}_{50}$ in $\mu\text{M}$.
    pub td50_um: f64,
    /// Maximum toxicity probability $T_{\max} \in [0.0, 1.0]$.
    pub toxicity_emax: f64,
    /// Hill slope for toxicity curve $\gamma_{\text{tox}}$.
    pub toxicity_gamma: f64,
}

impl HillPharmacodynamics {
    pub fn new(emax: f64, ec50_um: f64, gamma: f64, td50_um: f64) -> Self {
        Self {
            emax: emax.clamp(0.0, 1.0),
            ec50_um: ec50_um.max(1e-9),
            gamma: gamma.clamp(0.2, 10.0),
            td50_um: td50_um.max(1e-9),
            toxicity_emax: 1.0,
            toxicity_gamma: 2.0,
        }
    }

    pub fn evaluate_efficacy(&self, concentration_um: f64) -> f64 {
        self.emax * evaluate_hill_response(concentration_um, self.ec50_um, self.gamma)
    }

    pub fn evaluate_toxicity(&self, concentration_um: f64) -> f64 {
        self.toxicity_emax
            * evaluate_hill_response(concentration_um, self.td50_um, self.toxicity_gamma)
    }
}

/// Quantified modulation of primary hallmarks / levers of aging resulting from pharmacodynamic response.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct AgingLeverModulation {
    /// Fractional clearance of senescent cell load ($0.0\text{--}1.0$).
    pub senescent_fraction_cleared: f64,
    /// Acceleration or restoration factor for intrinsic DNA/cellular repair ($\ge 1.0$).
    pub repair_rate_enhancement_factor: f64,
    /// Fractional reduction in systemic inflammatory secretome / SASP ($0.0\text{--}1.0$).
    pub systemic_inflammation_reduction: f64,
}

/// Recommended clinical decision rule for the tested dosing regimen.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum RegimenDosingRecommendation {
    /// Maintain current dosing: robust efficacy achieved without approaching toxic threshold.
    MaintainOptimalDose,
    /// Underdosed: mean efficacy is subtherapeutic (< 40%); increase dose or dosing frequency.
    IncreaseDoseUnderdosed,
    /// Overdosed: peak toxicity exceeds safety threshold (> 25%); reduce nominal dose or extend interval.
    ReduceDoseOverdosed,
    /// Severe acute toxicity spike detected (> 45%); pause treatment immediately to prevent organ injury.
    PauseTreatmentToxicitySpike,
}

/// Report produced by advanced multi-compartment PK/PD simulation.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct AdvancedPkPdReport {
    pub peak_central_conc_um: f64,
    pub trough_central_conc_um: f64,
    pub mean_efficacy: f64,
    pub peak_toxicity: f64,
    pub within_therapeutic_window: bool,
    pub aging_lever_effects: AgingLeverModulation,
    pub recommendation: RegimenDosingRecommendation,
    pub decision_rationale: String,
}

/// Simulates a two-compartment PK and Hill-PD model over an arbitrary intermittent dosing schedule.
pub fn simulate_two_compartment_pkpd(
    pk: &TwoCompartmentPk,
    pd: &HillPharmacodynamics,
    schedule: &IntermittentDosingSchedule,
) -> AdvancedPkPdReport {
    if schedule.events.is_empty() {
        return AdvancedPkPdReport {
            peak_central_conc_um: 0.0,
            trough_central_conc_um: 0.0,
            mean_efficacy: 0.0,
            peak_toxicity: 0.0,
            within_therapeutic_window: false,
            aging_lever_effects: AgingLeverModulation {
                senescent_fraction_cleared: 0.0,
                repair_rate_enhancement_factor: 1.0,
                systemic_inflammation_reduction: 0.0,
            },
            recommendation: RegimenDosingRecommendation::IncreaseDoseUnderdosed,
            decision_rationale: "No doses scheduled.".to_string(),
        };
    }

    let last_dose_time = schedule
        .events
        .iter()
        .map(|e| e.time_hours)
        .fold(0.0_f64, f64::max);
    let total_sim_hours = (last_dose_time + 24.0).max(48.0);
    let dt = 0.5_f64;
    let num_steps = (total_sim_hours / dt) as usize;

    let mut c1 = 0.0_f64;
    let mut c2 = 0.0_f64;
    let mut peak_c1 = 0.0_f64;
    let mut min_trough_c1 = f64::INFINITY;
    let mut sum_eff = 0.0_f64;
    let mut peak_tox = 0.0_f64;

    let mut event_idx = 0;
    let mut events_sorted = schedule.events.clone();
    events_sorted
        .sort_by(|a, b| a.time_hours.partial_cmp(&b.time_hours).unwrap());

    for step in 0..num_steps {
        let current_time = step as f64 * dt;

        // Apply any dosing events scheduled at or before current_time
        while event_idx < events_sorted.len()
            && events_sorted[event_idx].time_hours <= current_time
        {
            let dose = events_sorted[event_idx].dose_amount;
            let added_conc = (dose * pk.bioavailability) / pk.volume_central_l;
            if step > 0 && c1 < min_trough_c1 {
                min_trough_c1 = c1;
            }
            c1 += added_conc;
            event_idx += 1;
        }

        if c1 > peak_c1 {
            peak_c1 = c1;
        }

        let eff = pd.evaluate_efficacy(c1);
        let tox = pd.evaluate_toxicity(c1);

        sum_eff += eff;
        if tox > peak_tox {
            peak_tox = tox;
        }

        // Two-compartment linear ODE state evolution:
        // dc1/dt = - (kel + k12) * c1 + k21 * c2 * (V2 / V1)
        // dc2/dt = k12 * c1 * (V1 / V2) - k21 * c2
        let dc1 = -(pk.kel + pk.k12) * c1 + pk.k21 * c2 * (pk.volume_peripheral_l / pk.volume_central_l);
        let dc2 = pk.k12 * c1 * (pk.volume_central_l / pk.volume_peripheral_l) - pk.k21 * c2;

        c1 = (c1 + dc1 * dt).max(0.0);
        c2 = (c2 + dc2 * dt).max(0.0);
    }

    if !min_trough_c1.is_finite() {
        min_trough_c1 = c1;
    }

    let mean_eff = sum_eff / num_steps as f64;
    let within_window = mean_eff >= 0.40 && peak_tox < 0.25;

    let (rec, rationale) = if peak_tox >= 0.45 {
        (
            RegimenDosingRecommendation::PauseTreatmentToxicitySpike,
            format!("Severe toxicity spike ({:.1}%); pause treatment immediately.", peak_tox * 100.0),
        )
    } else if peak_tox >= 0.25 {
        (
            RegimenDosingRecommendation::ReduceDoseOverdosed,
            format!("Excessive off-target toxicity ({:.1}%); reduce dose by 25-50%.", peak_tox * 100.0),
        )
    } else if mean_eff < 0.40 {
        (
            RegimenDosingRecommendation::IncreaseDoseUnderdosed,
            format!("Subtherapeutic mean efficacy ({:.1}%); increase dose or dosing frequency.", mean_eff * 100.0),
        )
    } else {
        (
            RegimenDosingRecommendation::MaintainOptimalDose,
            format!("Optimal therapeutic window achieved (efficacy: {:.1}%, toxicity: {:.1}%).", mean_eff * 100.0, peak_tox * 100.0),
        )
    };

    let aging_levers = AgingLeverModulation {
        senescent_fraction_cleared: (mean_eff * 0.85 * (1.0 - peak_tox * 0.5)).clamp(0.0, 0.95),
        repair_rate_enhancement_factor: 1.0 + (mean_eff * 0.65).clamp(0.0, 1.5),
        systemic_inflammation_reduction: (mean_eff * 0.70).clamp(0.0, 0.90),
    };

    AdvancedPkPdReport {
        peak_central_conc_um: peak_c1,
        trough_central_conc_um: min_trough_c1,
        mean_efficacy: mean_eff,
        peak_toxicity: peak_tox,
        within_therapeutic_window: within_window,
        aging_lever_effects: aging_levers,
        recommendation: rec,
        decision_rationale: rationale,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pkpd_therapeutic_window_simulation() {
        let senolytic = PkPdParameters {
            compound_name: "Fisetin_Senolytic".to_string(),
            ec50_um: 2.0,
            td50_um: 25.0,
            hill_coefficient: 1.8,
            half_life_hours: 8.0,
            bioavailability: 0.85,
        };
        assert!(senolytic.therapeutic_index() > 10.0);

        let report = simulate_pkpd_regimen(&senolytic, 5.0, 12, 4);
        assert!(report.mean_efficacy > 0.50);
        assert!(report.peak_toxicity < 0.20);
        assert!(report.within_therapeutic_window);
    }

    #[test]
    fn test_hill_extreme_concentration_stability() {
        // Astronomical concentration would previously overflow to Infinity and produce NaN
        let eff_huge = evaluate_hill_response(1e50, 1.0, 10.0);
        assert_eq!(eff_huge, 1.0);
        assert!(!eff_huge.is_nan());

        let eff_tiny = evaluate_hill_response(1e-50, 1.0, 10.0);
        assert_eq!(eff_tiny, 0.0);
        assert!(!eff_tiny.is_nan());
    }

    #[test]
    fn test_two_compartment_underdose_optimal_overdose_regimens() {
        let pk = TwoCompartmentPk::new(
            "Dasatinib_Senolytic",
            12.0, // V1 = 12 L central
            28.0, // V2 = 28 L peripheral
            0.15, // kel = 0.15 hr^-1
            0.08, // k12
            0.04, // k21
            0.80, // Bioavailability
        );
        let (t_alpha, t_beta) = pk.hybrid_half_lives();
        assert!(t_alpha < t_beta);
        assert!(t_beta > 4.0);

        let pd = HillPharmacodynamics::new(
            0.95, // Emax = 95%
            1.50, // EC50 = 1.5 uM
            1.80, // gamma = 1.8
            18.0, // TD50 = 18.0 uM
        );

        // Case 1: Underdose (0.5 uM/dose every 48h) -> Subtherapeutic efficacy
        let sched_under = IntermittentDosingSchedule::pulsed(2.0, 48.0, 3);
        let rep_under = simulate_two_compartment_pkpd(&pk, &pd, &sched_under);
        assert!(rep_under.mean_efficacy < 0.40);
        assert!(!rep_under.within_therapeutic_window);
        assert_eq!(
            rep_under.recommendation,
            RegimenDosingRecommendation::IncreaseDoseUnderdosed
        );

        // Case 2: Optimal Regimen (60.0 mg every 6h for 8 doses) -> Efficacy >= 40%, Toxicity < 25%
        let sched_opt = IntermittentDosingSchedule::pulsed(60.0, 6.0, 8);
        let rep_opt = simulate_two_compartment_pkpd(&pk, &pd, &sched_opt);
        assert!(rep_opt.mean_efficacy >= 0.40);
        assert!(rep_opt.peak_toxicity < 0.25);
        assert!(rep_opt.within_therapeutic_window);
        assert_eq!(
            rep_opt.recommendation,
            RegimenDosingRecommendation::MaintainOptimalDose
        );
        assert!(rep_opt.aging_lever_effects.senescent_fraction_cleared > 0.30);
        assert!(rep_opt.aging_lever_effects.repair_rate_enhancement_factor > 1.20);
        assert!(rep_opt.aging_lever_effects.systemic_inflammation_reduction > 0.25);

        // Case 3: Overdose (200.0 uM/dose every 12h) -> High Toxicity Spike (> 45%)
        let sched_over = IntermittentDosingSchedule::pulsed(200.0, 12.0, 6);
        let rep_over = simulate_two_compartment_pkpd(&pk, &pd, &sched_over);
        assert!(rep_over.peak_toxicity > 0.25);
        assert!(!rep_over.within_therapeutic_window);
        assert!(matches!(
            rep_over.recommendation,
            RegimenDosingRecommendation::ReduceDoseOverdosed
                | RegimenDosingRecommendation::PauseTreatmentToxicitySpike
        ));
    }
}
