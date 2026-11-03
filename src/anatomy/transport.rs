//! Spatial Reaction-Diffusion PDE & Krogh Cylinder Transport Kinetics across anatomical tissue layers.
//!
//! Models radial oxygen diffusion ($p\text{O}_2(r)$) and macromolecular drug/vector concentration
//! gradients ($C(z)/C_0 = \frac{\cosh(\Phi (1 - z/L))}{\cosh(\Phi)}$) across vascular barriers
//! (`Blood-Brain Barrier`, `Continuous`, `Fenestrated`, `Sinusoidal`) and interstitial ECM depth.
//!
//! References:
//! - Krogh (1919) *J Physiol* 52:409–415 (The number and distribution of capillaries in muscles with calculations of the oxygen pressure head).
//! - Truskey, Yuan & Katz (2009) *Transport Phenomena in Biological Systems* (2nd ed., Pearson; Thiele modulus & reaction-diffusion in tissues).
//! - Pardridge (2012) *J Cereb Blood Flow Metab* 32:1959–1972 (Drug transport across the blood-brain barrier).

use super::ecm::ExtracellularMatrixState;
use crate::aging::vascular::{drug_bioavailability, VascularPerfusionState};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Microvascular endothelial barrier phenotype governing trans-vascular extravasation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum EndothelialBarrierType {
    /// Blood-Brain Barrier (`Claudin-5` / `Occludin` tight junctions, pore $< 1\text{ nm}$).
    /// Requires receptor-mediated transcytosis (`TfR` / `LRP1` shuttle) for macromolecular vectors.
    ContinuousTightJunctionBbb,
    /// Continuous non-fenestrated endothelium (`Heart`, `Skeletal Muscle`, `Skin`, `Thymus`, pore $\approx 4\text{--}6\text{ nm}$).
    ContinuousNonFenestrated,
    /// Fenestrated diaphragmed endothelium (`Kidney Cortex`, `Endocrine Glands`, `Choroid Plexus`, pore $\approx 15\text{ nm}$).
    FenestratedDiaphragmed,
    /// Discontinuous sinusoidal endothelium (`Liver`, `Bone Marrow`, `Spleen`, fenestrae $\approx 100\text{--}180\text{ nm}$).
    SinusoidalDiscontinuous,
}

impl EndothelialBarrierType {
    /// Computes trans-endothelial extravasation permeability factor ($0.02$ to $1.0$)
    /// for a molecule/vector of `molecular_weight_da` given optional `receptor_shuttle_efficiency`
    /// ($0.0\text{--}1.0$, e.g., engineered Transferrin Receptor `TfR` / `AAV.PHP.eB` BBB shuttle).
    pub fn trans_endothelial_permeability(
        &self,
        molecular_weight_da: f64,
        receptor_shuttle_efficiency: f64,
    ) -> f64 {
        let mw = molecular_weight_da.max(50.0);
        let shuttle = receptor_shuttle_efficiency.clamp(0.0, 1.0);

        match self {
            EndothelialBarrierType::ContinuousTightJunctionBbb => {
                if mw <= 450.0 {
                    0.82
                } else {
                    // Passive BBB penetration decays steeply above 450 Da; rescued by TfR/LRP1 receptor transcytosis shuttle
                    let passive = (450.0 / mw).powf(0.45).clamp(0.02, 0.65);
                    (passive + (0.88 - passive) * shuttle).clamp(0.02, 0.92)
                }
            }
            EndothelialBarrierType::ContinuousNonFenestrated => {
                let base = if mw <= 1_000.0 {
                    0.92
                } else {
                    (1.0 - (mw.log10() - 3.0).max(0.0) * 0.055).clamp(0.62, 0.92)
                };
                (base + (0.95 - base) * 0.5 * shuttle).clamp(0.50, 0.96)
            }
            EndothelialBarrierType::FenestratedDiaphragmed => {
                (0.95 - (mw.log10() - 3.5).max(0.0) * 0.035).clamp(0.75, 0.97)
            }
            EndothelialBarrierType::SinusoidalDiscontinuous => {
                // Large open fenestrae in liver and bone marrow permit high macromolecular extravasation
                (0.98 - (mw.log10() - 4.0).max(0.0) * 0.020).clamp(0.85, 0.99)
            }
        }
    }
}

/// Anatomical tissue layer within an organ compartment positioned at a radial distance from capillaries.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TissueLayer {
    /// Canonical histological layer name (e.g., `"TunicaIntima"`, `"CorticalNeuropilLayer"`).
    pub layer_name: String,
    /// Mean radial diffusion distance $z$ from nearest capillary lumen ($\mu\text{m}$).
    pub radial_depth_um: f64,
    /// Layer thickness $L$ ($\mu\text{m}$).
    pub thickness_um: f64,
    /// Cellular packing density fraction ($0.0$ to $1.0$).
    pub cellular_packing_density: f64,
    /// Basal tissue oxygen consumption rate parameter ($\text{mmHg}/\mu\text{m}^2$).
    pub o2_consumption_coeff: f64,
}

/// Solved steady-state spatial Reaction-Diffusion PDE profile for a specific `TissueLayer`.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SpatialLayerProfile {
    /// Layer name.
    pub layer_name: String,
    /// Radial distance from capillary ($\mu\text{m}$).
    pub radial_depth_um: f64,
    /// Steady-state interstitial partial pressure of oxygen $p\text{O}_2(r)$ ($\text{mmHg}$).
    pub local_po2_mmhg: f64,
    /// Effective drug/vector bioavailability reaching cells at this radial depth ($0.0$ to $1.0$).
    pub local_drug_bioavailability: f64,
    /// Local `YAP`/`TAZ` mechanotransduction stress ($0.0$ to $1.0$).
    pub mechanotransduction_stress: f64,
    /// Fraction of cells in this layer experiencing micro-hypoxia ($p\text{O}_2 < 20\text{ mmHg}$).
    pub hypoxic_stress_fraction: f64,
}

/// Solves the steady-state 1D/cylindrical Krogh-Erlang Reaction-Diffusion PDE for oxygen ($p\text{O}_2$)
/// and Thiele-modulus interstitial drug diffusion ($C(z)/C_0 = \frac{\cosh(\Phi(1 - z/L))}{\cosh(\Phi)}$)
/// for a given `TissueLayer`.
pub fn solve_layer_reaction_diffusion(
    layer: &TissueLayer,
    vascular: &VascularPerfusionState,
    ecm: &ExtracellularMatrixState,
    barrier: EndothelialBarrierType,
    molecular_weight_da: f64,
    receptor_shuttle_efficiency: f64,
) -> SpatialLayerProfile {
    // 1. Interstitial Krogh cylinder radius R_t increases when capillary density drops (R_t ~ R_0 / sqrt(cap_density))
    let cap_density = vascular.capillary_density.clamp(0.15, 1.0);
    let rarefaction_stretch = 1.0 / cap_density.sqrt();
    let effective_depth_um = layer.radial_depth_um * rarefaction_stretch;
    let max_domain_um = (layer.radial_depth_um + layer.thickness_um) * rarefaction_stretch;

    // 2. Steady-state Krogh-Erlang oxygen reaction-diffusion solution:
    // P(z) = P_cap - (M_0 / (4 * D_eff)) * z^2
    let ecm_o2_conductance =
        (1.0 - 0.30 * ecm.ages_crosslink_burden - 0.25 * (ecm.collagen_volume_fraction - 0.14).max(0.0))
            .clamp(0.45, 1.0);
    let o2_drop_mmhg = (layer.o2_consumption_coeff
        * layer.cellular_packing_density
        * effective_depth_um.powi(2))
        / ecm_o2_conductance;
    let local_po2 = (vascular.po2_mmhg - o2_drop_mmhg).clamp(6.0, 48.0);

    // 3. Steady-state Thiele-modulus drug/vector reaction-diffusion PDE:
    // d^2 C / dz^2 = (k_uptake / D_ecm) * C  =>  C(z)/C_0 = cosh(Phi * (1 - z/L)) / cosh(Phi)
    let base_vascular_bioavail = drug_bioavailability(cap_density, molecular_weight_da);
    let barrier_perm =
        barrier.trans_endothelial_permeability(molecular_weight_da, receptor_shuttle_efficiency);
    let ecm_transmittance = ecm.interstitial_diffusion_transmittance(molecular_weight_da);

    // Dimensionless Thiele modulus Phi = L * sqrt(k / D_eff)
    let z_norm = (effective_depth_um / max_domain_um.max(1.0)).clamp(0.0, 1.0);
    let thiele_phi = (0.35
        * layer.cellular_packing_density
        * rarefaction_stretch
        / ecm_transmittance.sqrt())
    .clamp(0.05, 2.5);
    let pde_depth_attenuation =
        (thiele_phi * (1.0 - z_norm)).cosh() / thiele_phi.cosh().max(1.0);

    let local_bioavail = (base_vascular_bioavail
        * barrier_perm
        * (0.45 + 0.55 * ecm_transmittance)
        * pde_depth_attenuation)
        .clamp(0.02, 0.99);

    let hypoxic_stress = if local_po2 < 24.0 {
        ((24.0 - local_po2) / 18.0).clamp(0.0, 1.0)
    } else {
        0.0
    };

    SpatialLayerProfile {
        layer_name: layer.layer_name.clone(),
        radial_depth_um: effective_depth_um,
        local_po2_mmhg: local_po2,
        local_drug_bioavailability: local_bioavail,
        mechanotransduction_stress: ecm.yap_taz_mechanotransduction_stress(),
        hypoxic_stress_fraction: hypoxic_stress,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_krogh_pde_radial_depth_and_bbb_shuttle_gradient() {
        let shallow = TissueLayer {
            layer_name: "PerivascularLayer".to_string(),
            radial_depth_um: 12.0,
            thickness_um: 20.0,
            cellular_packing_density: 0.55,
            o2_consumption_coeff: 0.0025,
        };
        let deep = TissueLayer {
            layer_name: "DeepNeuropilLayer".to_string(),
            radial_depth_um: 75.0,
            thickness_um: 40.0,
            cellular_packing_density: 0.70,
            o2_consumption_coeff: 0.0025,
        };

        let vasc = VascularPerfusionState::from_somatic_age(72.0);
        let ecm = ExtracellularMatrixState::from_somatic_age(1.5, 72.0);

        // Without BBB receptor shuttle, large AAV (3.7 MDa) is strongly blocked by BBB
        let deep_no_shuttle = solve_layer_reaction_diffusion(
            &deep,
            &vasc,
            &ecm,
            EndothelialBarrierType::ContinuousTightJunctionBbb,
            3_700_000.0,
            0.0,
        );
        // With engineered TfR BBB shuttle (0.90), delivery increases dramatically
        let deep_with_shuttle = solve_layer_reaction_diffusion(
            &deep,
            &vasc,
            &ecm,
            EndothelialBarrierType::ContinuousTightJunctionBbb,
            3_700_000.0,
            0.90,
        );
        let shallow_with_shuttle = solve_layer_reaction_diffusion(
            &shallow,
            &vasc,
            &ecm,
            EndothelialBarrierType::ContinuousTightJunctionBbb,
            3_700_000.0,
            0.90,
        );

        assert!(deep_with_shuttle.local_drug_bioavailability > deep_no_shuttle.local_drug_bioavailability * 4.0);
        // Shallow perivascular layer has higher pO2 and higher bioavailability than deep parenchymal layer
        assert!(shallow_with_shuttle.local_po2_mmhg > deep_with_shuttle.local_po2_mmhg + 4.0);
        assert!(shallow_with_shuttle.local_drug_bioavailability > deep_with_shuttle.local_drug_bioavailability);
    }
}
