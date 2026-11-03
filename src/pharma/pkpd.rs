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
}
