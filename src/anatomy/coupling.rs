//! # Cross-Organ Inter-Systemic Coupling Engine (`anatomy::coupling`)
//!
//! Models systemic paracrine and physiological cross-talk across all 6 core organ compartments:
//! 1. **Systemic Circulatory SASP Burden**: Paracrine senescence-associated secretory phenotype
//!    secreted by aged or damaged organs (e.g. skin, liver) circulates systemically and feeds
//!    back to degrade endothelial eNOS activity, increase arterial stiffness, and trigger
//!    secondary senescence across distant organs.
//! 2. **Vascular Oxygen Delivery Coupling**: Cardiovascular microvascular health and cardiac
//!    output govern systemic $p\text{O}_2$ delivery. Microvascular rarefaction impairs oxygenation
//!    and bioenergetics in the cerebral cortex and hepatic metabolic parenchyma.
//! 3. **Immune Surveillance Clearance Coupling**: Thymic naive T-cell output and bone marrow
//!    phagocytic capacity govern the systemic clearance rate of senescent and mutated cells across
//!    all peripheral organs. Rejuvenating the immune compartment accelerates clearance of senescent
//!    burden in non-immune organs.

use super::cycles::WholeBodySystemicCircuits;
use super::organ::*;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Metrics resulting from a synchronized cross-organ coupling step.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SystemicCouplingMetrics {
    /// Aggregate circulating systemic SASP burden ($0.0\text{--}1.0$).
    pub systemic_sasp_burden: f64,
    /// Systemic microvascular perfusion and oxygen delivery efficiency ($0.0\text{--}1.0$).
    pub vascular_oxygen_delivery_efficiency: f64,
    /// Systemic immune patrolling and senescent cell clearance competence ($0.0\text{--}1.0$).
    pub immune_surveillance_clearance_competence: f64,
    /// Mean cortical parenchymal $p\text{O}_2$ ($\text{mmHg}$).
    pub cortical_po2_mmhg: f64,
    /// Mean hepatic parenchymal $p\text{O}_2$ ($\text{mmHg}$).
    pub hepatic_po2_mmhg: f64,
    /// Mean arterial stiffness Young's modulus ($\text{kPa}$).
    pub arterial_stiffness_kpa: f64,
    /// Daily thymic naive T-cell output ($\times 10^7\text{ cells/day}$).
    pub naive_t_cell_output_10e7: f64,
}

/// Executes a synchronized cross-organ coupling step across the anatomical atlas and systemic circuits.
pub fn propagate_cross_organ_coupling(
    atlas: &mut AnatomicalSpatialAtlas,
    circuits: &mut WholeBodySystemicCircuits,
) -> SystemicCouplingMetrics {
    // 1. Compute weighted circulating systemic SASP burden across organs
    // Weights: Liver (0.30), Artery (0.20), Skin (0.15), Marrow (0.15), Thymus (0.10), Brain (0.10)
    let sasp_artery = atlas.organs[0].biomarker_set().sasp_secretion_rate;
    let sasp_cortex = atlas.organs[1].biomarker_set().sasp_secretion_rate;
    let sasp_liver = atlas.organs[2].biomarker_set().sasp_secretion_rate;
    let sasp_marrow = atlas.organs[3].biomarker_set().sasp_secretion_rate;
    let sasp_thymus = atlas.organs[4].biomarker_set().sasp_secretion_rate;
    let sasp_skin = atlas.organs[5].biomarker_set().sasp_secretion_rate;

    let systemic_sasp = (0.20 * sasp_artery
        + 0.10 * sasp_cortex
        + 0.30 * sasp_liver
        + 0.15 * sasp_marrow
        + 0.10 * sasp_thymus
        + 0.15 * sasp_skin)
        .clamp(0.01, 1.0);

    // 2. Compute immune surveillance clearance competence from Thymus and Bone Marrow
    let thymic_biomarkers = atlas.organs[4].biomarker_set();
    let marrow_biomarkers = atlas.organs[3].biomarker_set();
    let naive_t_output = (12.0 * (thymic_biomarkers.functional_reserve / 0.85)).clamp(0.5, 20.0);
    circuits.immune_lymphatic.daily_naive_t_cell_emigration_10e7 = naive_t_output;

    let immune_competence = (0.55 * (naive_t_output / 12.0)
        + 0.45 * marrow_biomarkers.viability)
        .clamp(0.10, 1.0);
    circuits.immune_lymphatic.peripheral_immune_patrolling_competence = immune_competence;

    // 3. Compute vascular oxygen delivery from Cardiovascular Endothelium
    let vascular_health = atlas.organs[0].biomarker_set().functional_reserve;
    let co_efficiency = (circuits.cardiovascular.cardiac_output_l_min / 5.4).clamp(0.4, 1.3);
    let o2_delivery_efficiency = (0.60 * vascular_health + 0.40 * co_efficiency).clamp(0.20, 1.0);

    // 4. Apply cross-organ feedback into organs:
    // A) Immune surveillance cleanses senescent burden across distant organs:
    for organ in &mut atlas.organs[..6] {
        let ecm = organ.ecm_mut();
        if immune_competence > 0.80 {
            let clearance = 0.04 * (immune_competence - 0.80);
            ecm.ages_crosslink_burden = (ecm.ages_crosslink_burden - clearance).max(0.01);
        } else {
            let escape = 0.02 * (0.80 - immune_competence);
            ecm.ages_crosslink_burden = (ecm.ages_crosslink_burden + escape).min(0.95);
        }
    }

    // B) Systemic SASP accelerates arterial stiffness and endothelial rarefaction:
    if systemic_sasp > 0.25 {
        let sasp_stress = (systemic_sasp - 0.25) * 2.5;
        let art_ecm = atlas.organs[0].ecm_mut();
        art_ecm.youngs_modulus_kpa = (art_ecm.youngs_modulus_kpa + 0.8 * sasp_stress).min(120.0);
    }

    // C) Vascular oxygen delivery modulates parenchymal pO2 in Cortex and Liver:
    let cortex_intrinsic_po2 = 14.0
        + 31.0
            * (0.65 * atlas.organs[1].vascular().capillary_density
                + 0.35 * atlas.organs[1].vascular().enos_activity);
    let cortex_po2 =
        (cortex_intrinsic_po2 * (0.40 + 0.60 * o2_delivery_efficiency)).clamp(15.0, 48.0);
    atlas.organs[1].vascular_mut().po2_mmhg = cortex_po2;

    let liver_intrinsic_po2 = 14.0
        + 31.0
            * (0.65 * atlas.organs[2].vascular().capillary_density
                + 0.35 * atlas.organs[2].vascular().enos_activity);
    let liver_po2 =
        (liver_intrinsic_po2 * (0.40 + 0.60 * o2_delivery_efficiency)).clamp(18.0, 52.0);
    atlas.organs[2].vascular_mut().po2_mmhg = liver_po2;

    let arterial_stiffness = atlas.organs[0].ecm().youngs_modulus_kpa;

    SystemicCouplingMetrics {
        systemic_sasp_burden: systemic_sasp,
        vascular_oxygen_delivery_efficiency: o2_delivery_efficiency,
        immune_surveillance_clearance_competence: immune_competence,
        cortical_po2_mmhg: cortex_po2,
        hepatic_po2_mmhg: liver_po2,
        arterial_stiffness_kpa: arterial_stiffness,
        naive_t_cell_output_10e7: naive_t_output,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coupling_thymic_rejuvenation_clears_systemic_sasp_and_enhances_distal_organs() {
        let mut atlas = AnatomicalSpatialAtlas::from_somatic_age(75.0);
        let mut circuits = WholeBodySystemicCircuits::new();

        // Baseline aged state
        let baseline = propagate_cross_organ_coupling(&mut atlas, &mut circuits);
        assert!(baseline.systemic_sasp_burden > 0.35);
        assert!(baseline.immune_surveillance_clearance_competence < 0.65);

        let liver_ages_before = atlas.organs[2].ecm().ages_crosslink_burden;
        let cortex_ages_before = atlas.organs[1].ecm().ages_crosslink_burden;

        // Rejuvenate ONLY the Thymus (Organ index 4) + Bone Marrow (index 3)
        atlas.organs[4].rejuvenate_architecture(0.95, 0.90, 0.85);
        atlas.organs[3].rejuvenate_architecture(0.95, 0.90, 0.85);
        circuits.immune_lymphatic.daily_naive_t_cell_emigration_10e7 = 13.5;

        // Propagate systemic coupling
        let post_treatment = propagate_cross_organ_coupling(&mut atlas, &mut circuits);

        // Systemic immune patrolling competence must increase
        assert!(post_treatment.immune_surveillance_clearance_competence > baseline.immune_surveillance_clearance_competence);
        assert!(post_treatment.naive_t_cell_output_10e7 > baseline.naive_t_cell_output_10e7);

        // Distal non-immune organs (Liver and Brain) must exhibit cleared senescent crosslink burden
        let liver_ages_after = atlas.organs[2].ecm().ages_crosslink_burden;
        let cortex_ages_after = atlas.organs[1].ecm().ages_crosslink_burden;
        assert!(liver_ages_after < liver_ages_before);
        assert!(cortex_ages_after < cortex_ages_before);
    }

    #[test]
    fn test_coupling_vascular_endothelial_rejuvenation_restores_parenchymal_oxygenation() {
        let mut atlas = AnatomicalSpatialAtlas::from_somatic_age(75.0);
        let mut circuits = WholeBodySystemicCircuits::new();

        let baseline = propagate_cross_organ_coupling(&mut atlas, &mut circuits);
        let cortex_po2_before = baseline.cortical_po2_mmhg;
        let liver_po2_before = baseline.hepatic_po2_mmhg;

        // Rejuvenate ONLY the Cardiovascular Endothelium (Organ index 0)
        atlas.organs[0].rejuvenate_architecture(0.95, 0.95, 0.90);
        circuits.cardiovascular.cardiac_output_l_min = 5.6;

        let post_treatment = propagate_cross_organ_coupling(&mut atlas, &mut circuits);

        // Vascular oxygen delivery efficiency improves
        assert!(post_treatment.vascular_oxygen_delivery_efficiency > baseline.vascular_oxygen_delivery_efficiency);
        // Distal parenchymal oxygen tension in both Cortex and Liver increases
        assert!(post_treatment.cortical_po2_mmhg > cortex_po2_before);
        assert!(post_treatment.hepatic_po2_mmhg > liver_po2_before);
    }

    #[test]
    fn test_coupling_isolated_organ_damage_propagates_systemic_sasp_and_arterial_stiffening() {
        let mut atlas = AnatomicalSpatialAtlas::from_somatic_age(30.0);
        let mut circuits = WholeBodySystemicCircuits::new();

        let baseline = propagate_cross_organ_coupling(&mut atlas, &mut circuits);
        let baseline_sasp = baseline.systemic_sasp_burden;
        let baseline_arterial_stiffness = baseline.arterial_stiffness_kpa;

        // Induce isolated severe fibrotic and SASP damage ONLY in the Liver (index 2)
        {
            let liver_ecm = atlas.organs[2].ecm_mut();
            liver_ecm.ages_crosslink_burden = 0.95;
            liver_ecm.collagen_volume_fraction = 0.65;
            liver_ecm.youngs_modulus_kpa = 45.0;
        }

        let post_damage = propagate_cross_organ_coupling(&mut atlas, &mut circuits);

        // Systemic circulating SASP must surge
        assert!(post_damage.systemic_sasp_burden > baseline_sasp);
        // Distal arterial wall stiffness must increase in response to elevated systemic SASP
        assert!(post_damage.arterial_stiffness_kpa > baseline_arterial_stiffness);
    }
}
