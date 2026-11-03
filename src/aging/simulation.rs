//! Multi-generation cell population simulation and long-term replicative trajectory modeling.
//!
//! Simulates the evolutionary trajectory of edited cell populations over multiple generations
//! (20 to 100 cell divisions), tracking:
//! - Chromosomal telomere attrition and telomerase-mediated renewal.
//! - Transition into replicative senescence (M1 arrest).
//! - Neoplastic crisis bypass and malignant transformation risk (M2).
//! - Active clearance of senescent and transformed cells by synthetic circuits.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// State of the cell population at a specific division generation.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CellPopulationState {
    /// Generation index ($0 = \text{initial edited baseline}$).
    pub generation: u32,
    /// Population mean terminal telomere length in base pairs.
    pub mean_telomere_length_bp: f64,
    /// Projected biological age in years.
    pub biological_age_years: f64,
    /// Proportion of healthy, proliferating cells ($0.0$ to $1.0$).
    pub fraction_proliferating: f64,
    /// Proportion of cells arrested in replicative senescence ($0.0$ to $1.0$).
    pub fraction_senescent: f64,
    /// Proportion of cells experiencing neoplastic malignant transformation ($0.0$ to $1.0$).
    pub fraction_transformed: f64,
    /// Cumulative oncogenic transformation risk ($0.0$ to $1.0$).
    pub cumulative_oncogenic_risk: f64,
    /// Shelterin protein complex capping integrity index ($0.0$ to $1.0$).
    pub shelterin_capping: f64,
}

/// Configuration parameters for multi-generation population simulation.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MultiGenerationSimulationConfig {
    /// Number of cell divisions/generations to simulate (typically 20 to 100).
    pub total_generations: u32,
    /// Level of chronic oxidative stress (ROS / 8-oxodG formation, $0.0$ to $3.0$).
    pub oxidative_stress: f64,
    /// Level of lagging strand G-quadruplex replication fork stalling ($0.0$ to $2.0$).
    pub replication_stress_g4: f64,
    /// Regulated telomerase (TERT/TERC) expression level ($0.0$ to $1.0$).
    pub telomerase_expression: f64,
    /// Telomerase repeat addition processivity (nucleotides synthesized per binding event).
    pub telomerase_processivity: f64,
    /// Extra Pachyderm Super-p53 surveillance active.
    pub has_super_p53: bool,
    /// Inducible cell suicide switch (iCasp9) active.
    pub has_suicide_switch: bool,
    /// Rate of senescent cell clearance by senolytic circuits per generation ($0.0$ to $0.5$).
    pub senolytic_clearance_rate: f64,
    /// Critical length threshold triggering replicative arrest in bp (~3,800 bp).
    pub critical_hayflick_limit_bp: f64,
}

impl Default for MultiGenerationSimulationConfig {
    fn default() -> Self {
        Self {
            total_generations: 50,
            oxidative_stress: 1.0,
            replication_stress_g4: 0.8,
            telomerase_expression: 0.80,
            telomerase_processivity: 52.0,
            has_super_p53: true,
            has_suicide_switch: true,
            senolytic_clearance_rate: 0.10,
            critical_hayflick_limit_bp: 3_800.0,
        }
    }
}

/// Complete multi-generation trajectory output.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SimulationTrajectory {
    /// State at each discrete division generation.
    pub history: Vec<CellPopulationState>,
    /// Final cellular population state.
    pub final_state: CellPopulationState,
    /// Total population doublings sustained.
    pub doublings_achieved: u32,
    /// Whether an uncontrolled neoplastic crisis or population extinction occurred.
    pub crisis_alarm_triggered: bool,
    /// Whether the tissue population remains healthy and viable.
    pub is_population_viable: bool,
}

impl SimulationTrajectory {
    /// Formats a concise generation summary for inclusion in longevity reports.
    pub fn summary_overview(&self) -> String {
        let mut out = String::new();
        out.push_str("Gen | Telomere (bp) | Bio Age (yr) | Prolif % | Senescent % | Transformed % | Cumul Risk\n");
        out.push_str("----+---------------+--------------+----------+-------------+---------------+-----------\n");

        // Sample initial, quarter, midpoint, and final generations
        let len = self.history.len();
        let sample_indices = [0, len / 4, len / 2, (3 * len) / 4, len.saturating_sub(1)];

        let mut seen = std::collections::HashSet::new();
        for &idx in &sample_indices {
            if idx < len && seen.insert(idx) {
                let s = &self.history[idx];
                out.push_str(&format!(
                    "{:>3} | {:>13.0} | {:>12.1} | {:>7.1}% | {:>10.1}% | {:>12.4}% | {:>10.4}\n",
                    s.generation,
                    s.mean_telomere_length_bp,
                    s.biological_age_years,
                    s.fraction_proliferating * 100.0,
                    s.fraction_senescent * 100.0,
                    s.fraction_transformed * 100.0,
                    s.cumulative_oncogenic_risk
                ));
            }
        }
        out
    }
}

/// Executes multi-generation cell population simulation over divisions.
pub fn run_multi_generation_simulation(
    initial_telomere_bp: f64,
    initial_bio_age: f64,
    config: &MultiGenerationSimulationConfig,
) -> SimulationTrajectory {
    let mut history = Vec::with_capacity(config.total_generations as usize + 1);

    let mut current_length = initial_telomere_bp;
    let mut current_age = initial_bio_age;
    let mut fraction_proliferating = 0.95_f64;
    let mut fraction_senescent = 0.05_f64;
    let mut fraction_transformed = 0.00_f64;
    let mut cumulative_risk = 0.005_f64;
    let mut shelterin = (current_length / 10_000.0).clamp(0.20, 0.98);

    history.push(CellPopulationState {
        generation: 0,
        mean_telomere_length_bp: current_length,
        biological_age_years: current_age,
        fraction_proliferating,
        fraction_senescent,
        fraction_transformed,
        cumulative_oncogenic_risk: cumulative_risk,
        shelterin_capping: shelterin,
    });

    let basal_attrition = 55.0_f64;
    let ox_factor = 1.0 + config.oxidative_stress.clamp(0.0, 5.0) * 0.35;
    let g4_factor = 1.0 + config.replication_stress_g4.clamp(0.0, 4.0) * 0.25;

    let p53_damping = if config.has_super_p53 { 0.15 } else { 1.00 };
    let suicide_damping = if config.has_suicide_switch { 0.18 } else { 1.00 };

    let mut crisis_alarm = false;
    let mut doublings = 0;

    for gen in 1..=config.total_generations {
        // 1. Telomere attrition dynamics
        let capping_buffer = (1.0 / shelterin.clamp(0.25, 1.0)).min(2.5);
        let loss = basal_attrition * ox_factor * g4_factor * capping_buffer;

        // 2. Telomerase renewal dynamics
        let gain = if config.telomerase_expression > 0.01 {
            let accessibility = (1.1 - shelterin).clamp(0.25, 0.95);
            config.telomerase_expression * config.telomerase_processivity * accessibility
        } else {
            0.0
        };

        let net_change = gain - loss;
        current_length = (current_length + net_change).clamp(600.0, 18_500.0);

        // Update shelterin capping
        shelterin = if current_length > 6_000.0 {
            (0.85 + (current_length - 6_000.0) / 15_000.0).clamp(0.85, 0.98)
        } else {
            ((current_length - 600.0) / 5_400.0).clamp(0.10, 0.85)
        };

        // 3. Compartment transition dynamics
        if current_length <= config.critical_hayflick_limit_bp {
            // Replicative senescence triggers
            let arrest_rate = 0.08_f64;
            let newly_senescent = fraction_proliferating * arrest_rate;
            fraction_proliferating -= newly_senescent;
            fraction_senescent += newly_senescent;
        }

        // Senolytic clearance
        if config.senolytic_clearance_rate > 0.0 {
            let cleared = fraction_senescent * config.senolytic_clearance_rate;
            fraction_senescent -= cleared;
            // Space reclaimed by proliferating healthy cells
            fraction_proliferating += cleared * 0.85;
        }

        // Transformation hazard
        if config.telomerase_expression > 0.01 {
            let gen_hazard = (config.telomerase_expression.powf(1.4) * 0.004) * p53_damping * suicide_damping;
            let newly_transformed = fraction_proliferating * gen_hazard;
            fraction_proliferating -= newly_transformed;
            fraction_transformed += newly_transformed;
            cumulative_risk = (cumulative_risk + gen_hazard).min(0.99);

            // Active suicide clearance eliminates transformed cells
            if config.has_suicide_switch {
                fraction_transformed *= 0.25;
            }
        }

        // Renormalize fractions
        let total = (fraction_proliferating + fraction_senescent + fraction_transformed).max(0.001);
        fraction_proliferating /= total;
        fraction_senescent /= total;
        fraction_transformed /= total;

        // Biological age progression: 1 division ~ 0.2 chronological years drift
        let age_drift = 0.20;
        let telomere_benefit = if net_change > 0.0 { (net_change / 60.0) * 0.15 } else { 0.0 };
        current_age = (current_age + age_drift - telomere_benefit).max(18.0);

        if fraction_transformed > 0.05 {
            crisis_alarm = true;
        }
        if fraction_proliferating > 0.10 {
            doublings += 1;
        }

        history.push(CellPopulationState {
            generation: gen,
            mean_telomere_length_bp: current_length,
            biological_age_years: current_age,
            fraction_proliferating,
            fraction_senescent,
            fraction_transformed,
            cumulative_oncogenic_risk: cumulative_risk,
            shelterin_capping: shelterin,
        });
    }

    let final_state = history.last().unwrap().clone();
    let is_viable = final_state.fraction_proliferating >= 0.50 && !crisis_alarm;

    SimulationTrajectory {
        history,
        final_state,
        doublings_achieved: doublings,
        crisis_alarm_triggered: crisis_alarm,
        is_population_viable: is_viable,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multi_generation_simulation_balanced_viability() {
        let config = MultiGenerationSimulationConfig {
            total_generations: 40,
            telomerase_expression: 0.75,
            telomerase_processivity: 50.0,
            has_super_p53: true,
            has_suicide_switch: true,
            ..Default::default()
        };

        let trajectory = run_multi_generation_simulation(7_000.0, 52.0, &config);

        assert_eq!(trajectory.history.len(), 41);
        assert!(trajectory.is_population_viable);
        assert!(!trajectory.crisis_alarm_triggered);
        assert!(trajectory.final_state.fraction_proliferating > 0.60);
        assert!(trajectory.final_state.fraction_transformed < 0.01);
    }

    #[test]
    fn test_conservative_strategy_senescence_without_telomerase() {
        let config = MultiGenerationSimulationConfig {
            total_generations: 40,
            telomerase_expression: 0.0, // Silenced
            telomerase_processivity: 0.0,
            has_super_p53: true,
            has_suicide_switch: false,
            senolytic_clearance_rate: 0.0,
            critical_hayflick_limit_bp: 3_800.0,
            ..Default::default()
        };

        // Starting with short telomeres (4,200 bp)
        let trajectory = run_multi_generation_simulation(4_200.0, 65.0, &config);

        // Telomeres should decline and senescence should rise
        assert!(trajectory.final_state.mean_telomere_length_bp < 4_200.0);
        assert!(trajectory.final_state.fraction_senescent > 0.10);
        assert!(!trajectory.crisis_alarm_triggered);
    }
}
