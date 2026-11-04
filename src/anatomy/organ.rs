//! Ultra-Detailed *Gray's Anatomy* (42nd Ed.) & *Terminologia Anatomica* (`TA2`)
//! Object-Oriented Organ Hierarchy and Whole-Body 3D Spatial Topological Atlas.
//!
//! Models 12 canonical human organ systems with 6–8 sub-microscopic histological laminae each,
//! complete branching arterial/venous/lymphatic/neural trees, 3D `RAS` stereotaxic coordinates,
//! and inter-organ 3D spatial relationship graphs (`AnatomicalSpatialAtlas`).

use super::ecm::ExtracellularMatrixState;
use super::spatial::{
    classify_directional_relations, AnatomicalPoint3D, DirectionalAxisRelation, FascialBodyCavity,
    SpatialConnectivityConduit, SpatialRelationshipEdge,
};
use super::transport::{
    solve_layer_reaction_diffusion, EndothelialBarrierType, SpatialLayerProfile, TissueLayer,
};
use crate::aging::vascular::VascularPerfusionState;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Biological, biomechanical, and metabolic biomarker state for an organ.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct OrganBiomarkerSet {
    pub viability: f64,
    pub youngs_modulus_kpa: f64,
    pub interstitial_po2_mmhg: f64,
    pub sasp_secretion_rate: f64,
    pub functional_reserve: f64,
    pub senescent_cell_fraction: f64,
}

/// Dynamic kinetic parameters governing chronological degradation and functional decline.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct OrganAgingRateParams {
    pub baseline_senescence_rate_per_year: f64,
    pub ecm_crosslink_accumulation_rate_per_year: f64,
    pub microvascular_rarefaction_rate_per_year: f64,
    pub parenchymal_loss_rate_per_year: f64,
}

/// Pharmacodynamic and epigenetic responsiveness coefficients to rejuvenation therapies.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct OrganInterventionResponseCoefficients {
    pub osk_epigenetic_responsiveness: f64,
    pub senolytic_clearance_selectivity: f64,
    pub metabolic_sirt_responsiveness: f64,
}

/// Object-Oriented interface for an anatomical human organ with *Terminologia Anatomica* (`TA2`)
/// identity, 3D stereotaxic coordinates, neurovascular trees, stratified sub-microscopic layers,
/// and Extracellular Matrix (`ECM`) biomechanics.
pub trait OrganArchitecture {
    /// Canonical anatomical name of this organ.
    fn organ_name(&self) -> &'static str;

    /// Official *Terminologia Anatomica* (2nd Edition, FIPAT) code (e.g., `"TA2:5434"`).
    fn ta2_code(&self) -> &'static str;

    /// 3D somatic stereotaxic centroid in Right-Anterior-Superior (`RAS`, $\text{mm}$) coordinates.
    fn stereotaxic_center_mm(&self) -> AnatomicalPoint3D;

    /// Coelomic or fascial body cavity enclosing this organ.
    fn fascial_cavity(&self) -> FascialBodyCavity;

    /// Complete arterial supply branching tree from the aorta to precapillary arterioles.
    fn arterial_supply_tree(&self) -> &'static [&'static str];

    /// Complete venous drainage tree from postcapillary venules to the venae cavae / portal vein.
    fn venous_drainage_tree(&self) -> &'static [&'static str];

    /// Primary lymphatic / glymphatic drainage conduit and regional lymph node chain.
    fn lymphatic_pathway(&self) -> &'static str;

    /// Autonomic (sympathetic / parasympathetic) and somatic innervation pathway.
    fn innervation_pathway(&self) -> &'static str;

    /// Fraction of resting cardiac output ($0.0\text{--}1.0$) perfusing this organ.
    fn cardiac_output_fraction(&self) -> f64;

    /// Microvascular endothelial barrier phenotype governing trans-vascular permeability.
    fn barrier_type(&self) -> EndothelialBarrierType;

    /// Stratified histological & sub-microscopic layers ordered from capillary lumen outward.
    fn layers(&self) -> &[TissueLayer];

    /// Immutable reference to the organ's Extracellular Matrix (`ECM`) state.
    fn ecm(&self) -> &ExtracellularMatrixState;

    /// Mutable reference to the organ's Extracellular Matrix (`ECM`) state.
    fn ecm_mut(&mut self) -> &mut ExtracellularMatrixState;

    /// Immutable reference to the organ's microvascular perfusion state.
    fn vascular(&self) -> &VascularPerfusionState;

    /// Mutable reference to the organ's microvascular perfusion state.
    fn vascular_mut(&mut self) -> &mut VascularPerfusionState;

    /// Computes exact 3D Euclidean distance ($\text{mm}$) between this organ's centroid and `other`.
    fn euclidean_distance_to_mm(&self, other: &dyn OrganArchitecture) -> f64 {
        self.stereotaxic_center_mm()
            .distance_mm(&other.stereotaxic_center_mm())
    }

    /// Classifies *Terminologia Anatomica* directional relationships (`SuperiorTo`, `AnteriorTo`,
    /// `MedialTo`, `DeepTo`, etc.) of `self` relative to `other`.
    fn directional_relations_to(
        &self,
        other: &dyn OrganArchitecture,
    ) -> Vec<DirectionalAxisRelation> {
        classify_directional_relations(&self.stereotaxic_center_mm(), &other.stereotaxic_center_mm())
    }

    /// Solves the steady-state spatial Reaction-Diffusion PDE across all histological layers
    /// of this organ for a therapeutic molecule/vector of `molecular_weight_da` and
    /// `receptor_shuttle_efficiency` ($0.0\text{--}1.0$).
    fn solve_spatial_pde_profiles(
        &self,
        molecular_weight_da: f64,
        receptor_shuttle_efficiency: f64,
    ) -> Vec<SpatialLayerProfile> {
        self.layers()
            .iter()
            .map(|layer| {
                solve_layer_reaction_diffusion(
                    layer,
                    self.vascular(),
                    self.ecm(),
                    self.barrier_type(),
                    molecular_weight_da,
                    receptor_shuttle_efficiency,
                )
            })
            .collect()
    }

    /// Computes thickness-weighted mean interstitial oxygen tension ($p\text{O}_2$, $\text{mmHg}$) across all layers.
    fn mean_parenchymal_po2_mmhg(&self) -> f64 {
        let profiles = self.solve_spatial_pde_profiles(500.0, 0.0);
        let mut weighted_sum = 0.0;
        let mut total_thickness = 0.0;
        for (layer, prof) in self.layers().iter().zip(profiles.iter()) {
            weighted_sum += prof.local_po2_mmhg * layer.thickness_um;
            total_thickness += layer.thickness_um;
        }
        if total_thickness > 0.0 {
            weighted_sum / total_thickness
        } else {
            self.vascular().po2_mmhg
        }
    }

    /// Computes thickness-weighted mean drug/vector bioavailability ($0.0\text{--}1.0$) across all layers.
    fn mean_parenchymal_bioavailability(
        &self,
        molecular_weight_da: f64,
        receptor_shuttle_efficiency: f64,
    ) -> f64 {
        let profiles =
            self.solve_spatial_pde_profiles(molecular_weight_da, receptor_shuttle_efficiency);
        let mut weighted_sum = 0.0;
        let mut total_thickness = 0.0;
        for (layer, prof) in self.layers().iter().zip(profiles.iter()) {
            weighted_sum += prof.local_drug_bioavailability * layer.thickness_um;
            total_thickness += layer.thickness_um;
        }
        if total_thickness > 0.0 {
            weighted_sum / total_thickness
        } else {
            0.5
        }
    }

    /// Executes a combined microvascular angiogenesis + ECM deglycation/antifibrotic remodeling
    /// protocol on this organ.
    fn rejuvenate_architecture(
        &mut self,
        vegf_sirt1_boost: f64,
        ages_breaker_activity: f64,
        antifibrotic_mmp_remodeling: f64,
    ) {
        self.vascular_mut()
            .restore_microvascular_angiogenesis(vegf_sirt1_boost);
        self.ecm_mut().apply_matrix_rejuvenation(
            ages_breaker_activity,
            antifibrotic_mmp_remodeling,
            vegf_sirt1_boost,
        );
    }

    /// Advances chronological aging by `delta_years` across both the organ's microvascular bed
    /// and its Extracellular Matrix (`ECM`).
    fn apply_chronological_aging(&mut self, delta_years: f64) {
        let dt = delta_years.max(0.0);
        let v = self.vascular_mut();
        v.capillary_density = (v.capillary_density - 0.0055 * dt).clamp(0.22, 0.96);
        v.enos_activity = (v.enos_activity - 0.0060 * dt).clamp(0.20, 0.96);
        v.po2_mmhg =
            (14.0 + 31.0 * (0.65 * v.capillary_density + 0.35 * v.enos_activity)).clamp(15.0, 46.0);
        self.ecm_mut().apply_chronological_aging(dt);
    }

    /// Returns the organ's current biological and biomechanical biomarker set.
    fn biomarker_set(&self) -> OrganBiomarkerSet {
        let viab = self.viability_hook();
        let stiff = self.stiffness_hook_kpa();
        let po2 = self.oxygenation_hook_po2();
        let ecm = self.ecm();
        let sasp = (0.45 * ecm.ages_crosslink_burden
            + 0.35 * (1.0 - viab)
            + 0.20 * ecm.yap_taz_mechanotransduction_stress())
        .clamp(0.01, 1.0);
        let reserve = match self.organ_name() {
            "ThymicLymphoidLobule" => {
                // Thymic involution severely depresses functional thymic reserve in aged adults
                (viab * (1.0 - ecm.ages_crosslink_burden).powi(2)).clamp(0.05, 1.0)
            }
            _ => (0.50 * viab
                + 0.30 * (po2 / 40.0).min(1.0)
                + 0.20 * (1.0 - ecm.collagen_volume_fraction))
            .clamp(0.05, 1.0),
        };
        let senescent = (0.50 * ecm.ages_crosslink_burden + 0.50 * (1.0 - viab)).clamp(0.01, 0.95);
        OrganBiomarkerSet {
            viability: viab,
            youngs_modulus_kpa: stiff,
            interstitial_po2_mmhg: po2,
            sasp_secretion_rate: sasp,
            functional_reserve: reserve,
            senescent_cell_fraction: senescent,
        }
    }

    /// Intrinsic aging rate parameters governing chronological decay.
    fn aging_rate_params(&self) -> OrganAgingRateParams {
        match self.organ_name() {
            "ArterialVascularWall" => OrganAgingRateParams {
                baseline_senescence_rate_per_year: 0.009,
                ecm_crosslink_accumulation_rate_per_year: 0.012,
                microvascular_rarefaction_rate_per_year: 0.010,
                parenchymal_loss_rate_per_year: 0.004,
            },
            "CerebralCortexParenchyma" => OrganAgingRateParams {
                baseline_senescence_rate_per_year: 0.006,
                ecm_crosslink_accumulation_rate_per_year: 0.005,
                microvascular_rarefaction_rate_per_year: 0.008,
                parenchymal_loss_rate_per_year: 0.007,
            },
            "HepaticLobuleAcinus" => OrganAgingRateParams {
                baseline_senescence_rate_per_year: 0.008,
                ecm_crosslink_accumulation_rate_per_year: 0.009,
                microvascular_rarefaction_rate_per_year: 0.007,
                parenchymal_loss_rate_per_year: 0.005,
            },
            "HematopoieticBoneMarrowNiche" => OrganAgingRateParams {
                baseline_senescence_rate_per_year: 0.011,
                ecm_crosslink_accumulation_rate_per_year: 0.008,
                microvascular_rarefaction_rate_per_year: 0.009,
                parenchymal_loss_rate_per_year: 0.008,
            },
            "ThymicLymphoidLobule" => OrganAgingRateParams {
                baseline_senescence_rate_per_year: 0.025,
                ecm_crosslink_accumulation_rate_per_year: 0.015,
                microvascular_rarefaction_rate_per_year: 0.018,
                parenchymal_loss_rate_per_year: 0.022,
            },
            "DermalCutaneousTissue" => OrganAgingRateParams {
                baseline_senescence_rate_per_year: 0.010,
                ecm_crosslink_accumulation_rate_per_year: 0.014,
                microvascular_rarefaction_rate_per_year: 0.008,
                parenchymal_loss_rate_per_year: 0.006,
            },
            _ => OrganAgingRateParams {
                baseline_senescence_rate_per_year: 0.008,
                ecm_crosslink_accumulation_rate_per_year: 0.008,
                microvascular_rarefaction_rate_per_year: 0.007,
                parenchymal_loss_rate_per_year: 0.006,
            },
        }
    }

    /// Responsiveness coefficients to OSK reprogramming, senolytic clearance, and metabolic enhancement.
    fn intervention_response_coeffs(&self) -> OrganInterventionResponseCoefficients {
        match self.organ_name() {
            "ArterialVascularWall" => OrganInterventionResponseCoefficients {
                osk_epigenetic_responsiveness: 0.75,
                senolytic_clearance_selectivity: 0.90,
                metabolic_sirt_responsiveness: 0.88,
            },
            "CerebralCortexParenchyma" => OrganInterventionResponseCoefficients {
                osk_epigenetic_responsiveness: 0.70,
                senolytic_clearance_selectivity: 0.75,
                metabolic_sirt_responsiveness: 0.92,
            },
            "HepaticLobuleAcinus" => OrganInterventionResponseCoefficients {
                osk_epigenetic_responsiveness: 0.95,
                senolytic_clearance_selectivity: 0.92,
                metabolic_sirt_responsiveness: 0.95,
            },
            "HematopoieticBoneMarrowNiche" => OrganInterventionResponseCoefficients {
                osk_epigenetic_responsiveness: 0.90,
                senolytic_clearance_selectivity: 0.85,
                metabolic_sirt_responsiveness: 0.85,
            },
            "ThymicLymphoidLobule" => OrganInterventionResponseCoefficients {
                osk_epigenetic_responsiveness: 0.92,
                senolytic_clearance_selectivity: 0.88,
                metabolic_sirt_responsiveness: 0.90,
            },
            "DermalCutaneousTissue" => OrganInterventionResponseCoefficients {
                osk_epigenetic_responsiveness: 0.88,
                senolytic_clearance_selectivity: 0.95,
                metabolic_sirt_responsiveness: 0.82,
            },
            _ => OrganInterventionResponseCoefficients {
                osk_epigenetic_responsiveness: 0.80,
                senolytic_clearance_selectivity: 0.80,
                metabolic_sirt_responsiveness: 0.80,
            },
        }
    }

    /// Organ parenchymal cellular viability hook ($0.0\text{--}1.0$).
    fn viability_hook(&self) -> f64 {
        let ecm = self.ecm();
        let v = self.vascular();
        (0.40 * (v.capillary_density / 0.90).min(1.0)
            + 0.35 * (1.0 - ecm.ages_crosslink_burden)
            + 0.25 * (v.po2_mmhg / 38.0).min(1.0))
        .clamp(0.10, 1.0)
    }

    /// Elastic stiffness Young's modulus hook ($E$, $\text{kPa}$).
    fn stiffness_hook_kpa(&self) -> f64 {
        self.ecm().youngs_modulus_kpa
    }

    /// Interstitial tissue oxygen tension hook ($p\text{O}_2$, $\text{mmHg}$).
    fn oxygenation_hook_po2(&self) -> f64 {
        self.mean_parenchymal_po2_mmhg()
    }

    /// Applies targeted intervention responses (OSK dosage, senolytics, metabolic activators)
    /// scaled by this organ's specific response coefficients.
    fn apply_targeted_intervention_response(
        &mut self,
        osk_dose: f64,
        senolytic_dose: f64,
        metabolic_dose: f64,
    ) {
        let coeffs = self.intervention_response_coeffs();
        let eff_osk = osk_dose * coeffs.osk_epigenetic_responsiveness;
        let eff_seno = senolytic_dose * coeffs.senolytic_clearance_selectivity;
        let eff_meta = metabolic_dose * coeffs.metabolic_sirt_responsiveness;

        self.rejuvenate_architecture(
            (0.5 * eff_osk + 0.5 * eff_meta).clamp(0.0, 1.0),
            (0.7 * eff_seno + 0.3 * eff_meta).clamp(0.0, 1.0),
            (0.6 * eff_seno + 0.4 * eff_osk).clamp(0.0, 1.0),
        );
    }
}

/// Macro generating ultra-detailed *Gray's Anatomy* / `TA2` organ structs implementing [`OrganArchitecture`].
macro_rules! define_anatomical_organ {
    (
        $(#[$meta:meta])*
        pub struct $struct_name:ident {
            organ_name: $organ_name:expr,
            ta2_code: $ta2:expr,
            center_mm: ($cx:expr, $cy:expr, $cz:expr),
            cavity: $cavity:expr,
            arterial_tree: $art_tree:expr,
            venous_tree: $ven_tree:expr,
            lymphatic: $lymph:expr,
            innervation: $nerve:expr,
            cardiac_output_fraction: $co_frac:expr,
            barrier: $barrier:expr,
            youthful_modulus_kpa: $modulus:expr,
            default_layers: [ $( ($lname:expr, $depth:expr, $thick:expr, $pack:expr, $o2:expr) ),+ $(,)? ]
        }
    ) => {
        $(#[$meta])*
        #[derive(Debug, Clone, PartialEq)]
        #[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
        pub struct $struct_name {
            pub layers: Vec<TissueLayer>,
            pub ecm: ExtracellularMatrixState,
            pub vascular: VascularPerfusionState,
        }

        impl $struct_name {
            /// Creates a youthful organ architecture at baseline physiological compliance.
            pub fn youthful() -> Self {
                Self::from_somatic_age(20.0)
            }

            /// Initializes the organ's stratified histological laminae, microvascular state, and ECM biomechanics
            /// calibrated to chronological `age_years`.
            pub fn from_somatic_age(age_years: f64) -> Self {
                Self {
                    layers: vec![
                        $(
                            TissueLayer {
                                layer_name: $lname.to_string(),
                                radial_depth_um: $depth,
                                thickness_um: $thick,
                                cellular_packing_density: $pack,
                                o2_consumption_coeff: $o2,
                            },
                        )+
                    ],
                    ecm: ExtracellularMatrixState::from_somatic_age($modulus, age_years),
                    vascular: VascularPerfusionState::from_somatic_age(age_years),
                }
            }
        }

        impl OrganArchitecture for $struct_name {
            fn organ_name(&self) -> &'static str {
                $organ_name
            }

            fn ta2_code(&self) -> &'static str {
                $ta2
            }

            fn stereotaxic_center_mm(&self) -> AnatomicalPoint3D {
                AnatomicalPoint3D::new($cx, $cy, $cz)
            }

            fn fascial_cavity(&self) -> FascialBodyCavity {
                $cavity
            }

            fn arterial_supply_tree(&self) -> &'static [&'static str] {
                $art_tree
            }

            fn venous_drainage_tree(&self) -> &'static [&'static str] {
                $ven_tree
            }

            fn lymphatic_pathway(&self) -> &'static str {
                $lymph
            }

            fn innervation_pathway(&self) -> &'static str {
                $nerve
            }

            fn cardiac_output_fraction(&self) -> f64 {
                $co_frac
            }

            fn barrier_type(&self) -> EndothelialBarrierType {
                $barrier
            }

            fn layers(&self) -> &[TissueLayer] {
                &self.layers
            }

            fn ecm(&self) -> &ExtracellularMatrixState {
                &self.ecm
            }

            fn ecm_mut(&mut self) -> &mut ExtracellularMatrixState {
                &mut self.ecm
            }

            fn vascular(&self) -> &VascularPerfusionState {
                &self.vascular
            }

            fn vascular_mut(&mut self) -> &mut VascularPerfusionState {
                &mut self.vascular
            }
        }
    };
}

define_anatomical_organ! {
    /// 1. Arterial & Coronary Vascular Wall (`TA2:3891`, 7 histological laminae from glycocalyx to vasa vasorum).
    pub struct ArterialVascularOrgan {
        organ_name: "ArterialVascularWall",
        ta2_code: "TA2:3891",
        center_mm: (8.0, 22.0, -18.0),
        cavity: FascialBodyCavity::MiddleMediastinumPericardialSac,
        arterial_tree: &[
            "AscendingAorta_SinusOfValsalva",
            "LeftAndRightCoronaryArteries",
            "IntramuralPenetratingArterioles",
            "AdventitialVasaVasorumInterna",
        ],
        venous_tree: &[
            "VenousVasaVasorumPlexus",
            "GreatCardiacVein",
            "CoronarySinus",
            "RightAtrium",
        ],
        lymphatic: "SubepicardialLymphaticPlexus -> TracheobronchialNodes -> BronchomediastinalTrunk",
        innervation: "CardiacPlexus_SympatheticT1T4_And_VagusNerveCNX",
        cardiac_output_fraction: 0.045,
        barrier: EndothelialBarrierType::ContinuousNonFenestrated,
        youthful_modulus_kpa: 28.0,
        default_layers: [
            ("EndothelialHeparanSulfateGlycocalyx", 2.0, 4.0, 0.35, 0.0010),
            ("TunicaIntimaSquamousEndothelium", 6.0, 8.0, 0.72, 0.0018),
            ("SubendothelialLaminaPropriaMyointima", 14.0, 12.0, 0.60, 0.0016),
            ("InternalElasticLaminaFenestratedSheet", 26.0, 8.0, 0.45, 0.0012),
            ("TunicaMediaConcentricElasticSmoothMuscleLamellae", 42.0, 35.0, 0.75, 0.0022),
            ("ExternalElasticLaminaElastinMeshwork", 70.0, 10.0, 0.42, 0.0012),
            ("TunicaAdventitiaVasaVasorumNerviVasorumFibrosa", 85.0, 30.0, 0.52, 0.0015),
        ]
    }
}

define_anatomical_organ! {
    /// 2. Cerebral Cortex & Neurovascular Unit (`TA2:5434`, 8 Brodmann & Virchow-Robin sub-layers).
    pub struct CerebralCortexOrgan {
        organ_name: "CerebralCortexParenchyma",
        ta2_code: "TA2:5434",
        center_mm: (0.0, 15.0, 385.0),
        cavity: FascialBodyCavity::CranialSubarachnoidVault,
        arterial_tree: &[
            "ArchOfAorta",
            "CommonAndInternalCarotidArtery_CerebralPart",
            "CircleOfWillis_MiddleAndAnteriorCerebralArteries",
            "PialLeptomeningealAnastomoticPlexus",
            "PenetratingCorticalArterioles",
        ],
        venous_tree: &[
            "CorticalAscendingVenules",
            "SuperficialMiddleCerebralVein_BridgingVeins",
            "SuperiorSagittalAndTransverseDuralSinuses",
            "InternalJugularVein -> SuperiorVenaCava",
        ],
        lymphatic: "VirchowRobinGlymphaticAquaporin4 -> DuralMeningealLymphatics -> DeepCervicalLymphNodes -> JugularTrunk",
        innervation: "BasalForebrainNucleusBasalisMeynert_LocusCoeruleusNoradrenergic_RapheSerotonergic",
        cardiac_output_fraction: 0.150,
        barrier: EndothelialBarrierType::ContinuousTightJunctionBbb,
        youthful_modulus_kpa: 1.5,
        default_layers: [
            ("LuminalBbbClaudin5OccludinEndothelium", 4.0, 6.0, 0.68, 0.0025),
            ("CapillaryBasalLaminaContractilePericyteSheath", 10.0, 8.0, 0.62, 0.0024),
            ("VirchowRobinPerivascularGlymphaticInterstitium", 18.0, 10.0, 0.30, 0.0012),
            ("AstrocyticGliaLimitansAqp4Endfeet", 28.0, 12.0, 0.70, 0.0028),
            ("LaminaIMolecularCajalRetziusPlexiformLayer", 40.0, 15.0, 0.55, 0.0026),
            ("LaminaII_IIIExternalGranularAndPyramidalCorticoCortical", 55.0, 22.0, 0.74, 0.0031),
            ("LaminaIVInternalGranularThalamocorticalRecipient", 74.0, 18.0, 0.80, 0.0034),
            ("LaminaV_VIDeepBetzPyramidalAndMultiformNeuropil", 92.0, 25.0, 0.76, 0.0033),
        ]
    }
}

define_anatomical_organ! {
    /// 3. Hepatic Lobule & Rappaport Acinus (`TA2:3025`, 6 sinusoidal, Disse, and metabolic zonation laminae).
    pub struct HepaticLobuleOrgan {
        organ_name: "HepaticLobuleAcinus",
        ta2_code: "TA2:3025",
        center_mm: (-68.0, 28.0, -145.0),
        cavity: FascialBodyCavity::IntraperitonealPeritonealCavity,
        arterial_tree: &[
            "AbdominalAorta -> CeliacTrunk",
            "CommonAndProperHepaticArtery",
            "InterlobularHepaticArteriole_Plus_PortalVeinRadicle",
            "PeriportalSinusoidalInlet",
        ],
        venous_tree: &[
            "CentrilobularCentralVein",
            "SublobularCollectingVeins",
            "RightMiddleLeftHepaticVeins",
            "InferiorVenaCava_HepaticSegment",
        ],
        lymphatic: "SpaceOfDisse -> PeriportalSpaceOfMall -> HepaticHilarNodes -> CisternaChyli -> ThoracicDuct",
        innervation: "HepaticPlexus_CeliacSympatheticGanglion_And_PosteriorVagalTrunk",
        cardiac_output_fraction: 0.250,
        barrier: EndothelialBarrierType::SinusoidalDiscontinuous,
        youthful_modulus_kpa: 4.2,
        default_layers: [
            ("FenestratedSinusoidalLsecEndotheliumKupfferLumen", 6.0, 10.0, 0.65, 0.0020),
            ("PerisinusoidalSpaceOfDisseHepaticStellateItoNiche", 16.0, 12.0, 0.52, 0.0018),
            ("PeriportalZone1OxidativeGluconeogenicHepatocytes", 30.0, 22.0, 0.82, 0.0028),
            ("MidzonalZone2TransitionalHepaticPlates", 52.0, 24.0, 0.80, 0.0026),
            ("CentrilobularZone3PericentralCyp450DetoxHepatocytes", 76.0, 24.0, 0.78, 0.0029),
            ("BileCanaliculiCanalsOfHeringCholangiocyteProgenitorNiche", 98.0, 14.0, 0.70, 0.0021),
        ]
    }
}

define_anatomical_organ! {
    /// 4. Hematopoietic Bone Marrow Stem Niche (`TA2:5152`, 6 perivascular CAR, Nestin+ & endosteal laminae).
    pub struct BoneMarrowNicheOrgan {
        organ_name: "HematopoieticBoneMarrowNiche",
        ta2_code: "TA2:5152",
        center_mm: (0.0, 45.0, -35.0),
        cavity: FascialBodyCavity::EndostealMedullaryBoneCavity,
        arterial_tree: &[
            "InternalThoracicAndIntercostalArteries",
            "DiaphysealAndMetaphysealNutrientArteries",
            "CentralLongitudinalMedullaryArtery",
            "EndostealSca1PositiveTransitionalCapillaries",
        ],
        venous_tree: &[
            "FenestratedMedullaryVenousSinusoids",
            "CentralVenousSinus",
            "EmissaryAndNutrientVeins -> AzygosAndVenaCava",
        ],
        lymphatic: "PeriostealCorticalCanalicularDrainage -> InternalThoracicAndIliacLymphNodes",
        innervation: "UnmyelinatedSympatheticTyrosineHydroxylaseAdrenergicBeta3NerveFibers_CircadianCxcl12Control",
        cardiac_output_fraction: 0.040,
        barrier: EndothelialBarrierType::SinusoidalDiscontinuous,
        youthful_modulus_kpa: 8.0,
        default_layers: [
            ("SinusoidalVegfr3FenestratedEndothelialWall", 8.0, 10.0, 0.68, 0.0018),
            ("PerivascularCxcl12AbundantReticularCarCellNiche", 18.0, 15.0, 0.72, 0.0020),
            ("PeriarteriolarNestinPositiveLeprMesenchymalStemNiche", 33.0, 18.0, 0.70, 0.0021),
            ("CentralMedullaryErythroblasticIslandMegakaryocyteCords", 51.0, 22.0, 0.78, 0.0023),
            ("TrabecularEndostealOsteoblasticQuiescentLtHscNiche", 73.0, 18.0, 0.66, 0.0019),
            ("CalcifiedTrabecularSpongiosaOsteocyteLacunarCanalicularMatrix", 92.0, 15.0, 0.55, 0.0014),
        ]
    }
}

define_anatomical_organ! {
    /// 5. Thymic Lymphoid Organ (`TA2:5160`, 6 cortical-medullary AIRE+ mTEC & Hassall laminae).
    pub struct ThymicLymphoidOrgan {
        organ_name: "ThymicLymphoidLobule",
        ta2_code: "TA2:5160",
        center_mm: (0.0, 38.0, 42.0),
        cavity: FascialBodyCavity::SuperiorAnteriorMediastinum,
        arterial_tree: &[
            "SubclavianArtery -> InternalThoracicArtery",
            "AnteriorMediastinalAndPericardiacophrenicBranches",
            "InferiorThyroidArtery_ThymicBranches",
            "CorticomedullaryJunctionArcuateArterioles",
        ],
        venous_tree: &[
            "MedullaryPostcapillaryHighEndothelialVenules",
            "ThymicVeins_KeynesPlexus",
            "LeftBrachiocephalicVein -> SuperiorVenaCava",
        ],
        lymphatic: "EfferentInterlobularLymphatics -> AnteriorMediastinalParasternaNodes -> BronchomediastinalTrunk",
        innervation: "CervicothoracicStellateSympatheticGanglion_And_VagusRecurrentLaryngealBranches",
        cardiac_output_fraction: 0.015,
        barrier: EndothelialBarrierType::ContinuousNonFenestrated,
        youthful_modulus_kpa: 3.5,
        default_layers: [
            ("BloodThymusBarrierPerivascularEpithelialSheath", 8.0, 10.0, 0.68, 0.0018),
            ("SubcapsularCorticalFoxn1NurseEpithelialZone", 18.0, 14.0, 0.72, 0.0020),
            ("OuterCortexCd4Cd8DoubleNegativeNotch1ProliferationZone", 32.0, 18.0, 0.80, 0.0023),
            ("InnerCortexDoublePositiveCtecPositiveSelectionReticulum", 50.0, 20.0, 0.82, 0.0024),
            ("CorticomedullaryJunctionPerivascularTcellEgressGate", 68.0, 12.0, 0.70, 0.0020),
            ("MedullaryAirePositiveMtecNegativeSelectionHassallCorpuscleZone", 82.0, 22.0, 0.68, 0.0019),
        ]
    }
}

define_anatomical_organ! {
    /// 6. Cutaneous Dermal-Epidermal & Pilosebaceous Organ (`TA2:6098`, 7 epidermal, basement & dermal laminae).
    pub struct DermalCutaneousOrgan {
        organ_name: "DermalCutaneousTissue",
        ta2_code: "TA2:6098",
        center_mm: (145.0, 65.0, -20.0),
        cavity: FascialBodyCavity::SuperficialIntegumentaryMantle,
        arterial_tree: &[
            "MusculocutaneousAndFasciocutaneousPerforators",
            "DeepSubdermalCutaneousArterialPlexus",
            "SuperficialSubpapillaryArterialPlexus",
            "DermalPapillaryCapillaryHairpinLoops",
        ],
        venous_tree: &[
            "PostcapillarySuperficialAndDeepVenousPlexuses",
            "SubcutaneousSuperficialVeins",
            "AxillaryAndFemoralDeepVeins -> VenaeCavae",
        ],
        lymphatic: "InitialBlindDermalLymphaticCapillaries -> PrecollectorValvedVessels -> AxillaryAndInguinalNodes",
        innervation: "CutaneousSensoryAfferents_MeissnerMerkelPacinian_And_CholinergicSudomotorSympathetic",
        cardiac_output_fraction: 0.080,
        barrier: EndothelialBarrierType::ContinuousNonFenestrated,
        youthful_modulus_kpa: 12.0,
        default_layers: [
            ("PapillaryDermisCapillaryLoopMeissnerCorpuscleZone", 8.0, 14.0, 0.58, 0.0015),
            ("DermoEpidermalLaminaDensaCollagenVIIAnchoringFibrilJunction", 20.0, 8.0, 0.65, 0.0016),
            ("StratumBasaleGerminativumCol17a1EpidermalStemLayer", 28.0, 12.0, 0.82, 0.0022),
            ("HairFollicleBulgeCd34Lgr5MultipotentStemNiche", 40.0, 16.0, 0.76, 0.0020),
            ("StratumSpinosumAndGranulosumDesmosomalLamellarBodyLayer", 56.0, 18.0, 0.74, 0.0016),
            ("ReticularDermisDenseIrregularCollagenIElastinFibroblastPlexus", 74.0, 28.0, 0.48, 0.0013),
            ("HypodermalSubcutaneousAdiposeVascularFascia", 98.0, 20.0, 0.42, 0.0011),
        ]
    }
}

define_anatomical_organ! {
    /// 7. Left Ventricular Myocardium & His-Purkinje Syncytium (`TA2:3932`, 6 myocardial & conduction laminae).
    pub struct MyocardialVentricularOrgan {
        organ_name: "LeftVentricularMyocardium",
        ta2_code: "TA2:3932",
        center_mm: (32.0, 35.0, -48.0),
        cavity: FascialBodyCavity::MiddleMediastinumPericardialSac,
        arterial_tree: &[
            "AscendingAorta_LeftCoronarySinus",
            "LeftMainCoronaryArtery -> LeftAnteriorDescendingAndCircumflex",
            "TransmuralSubendocardialPerforatingArterioles",
            "HighDensityPerimyocyticCapillaryPlexus",
        ],
        venous_tree: &[
            "ThebesianAndIntramyocardialVenules",
            "GreatAndMiddleCardiacVeins",
            "CoronarySinus -> RightAtrium",
        ],
        lymphatic: "SubendocardialToSubepicardialLymphatics -> LeftCoronaryLymphaticTrunk -> TracheobronchialNodes",
        innervation: "IntrinsicCardiacGanglionicPlexus_SympatheticBeta1Adrenergic_And_VagalM2Cholinergic",
        cardiac_output_fraction: 0.050,
        barrier: EndothelialBarrierType::ContinuousNonFenestrated,
        youthful_modulus_kpa: 18.0,
        default_layers: [
            ("PerimyocyticCoronaryCapillaryEndothelium", 6.0, 8.0, 0.70, 0.0028),
            ("EndocardialEndotheliumSubendocardialConnectiveLayer", 14.0, 12.0, 0.62, 0.0024),
            ("SubendocardialHisPurkinjeFastConductionFiberNetwork", 26.0, 15.0, 0.74, 0.0032),
            ("InnerObliqueHelicalCardiomyocyteIntercalatedDiscSyncytium", 42.0, 24.0, 0.85, 0.0038),
            ("MidwallCircumferentialHighTensionSarcomereMyocardium", 66.0, 28.0, 0.86, 0.0040),
            ("OuterSpiralSubepicardialMyocardiumVisceralPericardium", 92.0, 20.0, 0.72, 0.0030),
        ]
    }
}

define_anatomical_organ! {
    /// 8. Renal Cortex Glomerular-Tubular Nephron (`TA2:3358`, 6 podocyte, tubular & vasa recta laminae).
    pub struct RenalNephronOrgan {
        organ_name: "RenalCorticalNephron",
        ta2_code: "TA2:3358",
        center_mm: (58.0, -42.0, -210.0),
        cavity: FascialBodyCavity::RetroperitonealFascialSpace,
        arterial_tree: &[
            "AbdominalAorta -> RenalArtery",
            "SegmentalAndInterlobarArteries",
            "ArcuateAndCorticalInterlobularArteries",
            "AfferentGlomerularArteriole -> GlomerularTuft -> EfferentArteriole",
        ],
        venous_tree: &[
            "PeritubularCapillaryPlexusAndVasaRectaRecta",
            "StellateAndInterlobularVeins -> ArcuateAndInterlobarVeins",
            "RenalVein -> InferiorVenaCava",
        ],
        lymphatic: "CorticalInterstitialLymphatics -> RenalHilarNodes -> LateralAorticLumbarNodes -> CisternaChyli",
        innervation: "RenalSympatheticPlexus_AorticorenalGanglion_LesserAndLeastSplanchnicNerves_JuxtaglomerularRenin",
        cardiac_output_fraction: 0.200,
        barrier: EndothelialBarrierType::FenestratedDiaphragmed,
        youthful_modulus_kpa: 5.0,
        default_layers: [
            ("FenestratedGlomerularCapillaryEndothelium", 5.0, 8.0, 0.68, 0.0022),
            ("GlomerularBasementMembraneCol4a3LaminaDensa", 13.0, 8.0, 0.50, 0.0015),
            ("VisceralPodocyteNephrinSlitDiaphragmBowmanCapsule", 22.0, 14.0, 0.72, 0.0026),
            ("MaculaDensaJuxtaglomerularGranularReninSensor", 36.0, 14.0, 0.76, 0.0028),
            ("ProximalConvolutedTubuleBrushBorderMitochondrialTransportZone", 54.0, 26.0, 0.84, 0.0036),
            ("OuterMedullaryThickAscendingLimbVasaRectaEpoInterstitium", 82.0, 24.0, 0.80, 0.0035),
        ]
    }
}

define_anatomical_organ! {
    /// 9. Pulmonary Alveolar-Capillary Acinus (`TA2:3195`, 6 blood-gas barrier & AT2 stem laminae).
    pub struct PulmonaryAlveolarOrgan {
        organ_name: "PulmonaryAlveolarAcinus",
        ta2_code: "TA2:3195",
        center_mm: (72.0, 10.0, 12.0),
        cavity: FascialBodyCavity::ThoracicPleuralCavity,
        arterial_tree: &[
            "RightVentricle -> PulmonaryTrunk",
            "LeftAndRightPulmonaryArteries -> LobarAndSegmentalBranches",
            "IntralobularAcinarArterioles",
            "AlveolarSeptalSheetCapillaryBed_Plus_BronchialArteries",
        ],
        venous_tree: &[
            "AlveolarPostcapillaryVenules",
            "InterlobularPulmonaryVenousTributaries",
            "FourSuperiorAndInferiorPulmonaryVeins -> LeftAtrium",
        ],
        lymphatic: "JuxtaAlveolarInterstitialLymphatics -> BronchopulmonaryHilarNodes -> TracheobronchialNodes",
        innervation: "PulmonaryPlexus_VagalBronchomotorParasympathetic_And_T2T5SympatheticBronchodilator",
        cardiac_output_fraction: 1.000,
        barrier: EndothelialBarrierType::ContinuousNonFenestrated,
        youthful_modulus_kpa: 2.2,
        default_layers: [
            ("AlveolarSeptalCapillaryEndothelialSheet", 3.0, 5.0, 0.65, 0.0014),
            ("FusedAlveolarCapillaryBasementMembrane0_5Micron", 8.0, 5.0, 0.45, 0.0010),
            ("SquamousTypeIAlveolarPneumocyteGasExchangeWall", 14.0, 10.0, 0.62, 0.0016),
            ("CuboidalTypeIIAlveolarSurfactantProgenitorStemNiche", 26.0, 14.0, 0.74, 0.0022),
            ("AlveolarInterstitialElastinCollagenSeptalMatrix", 42.0, 18.0, 0.52, 0.0015),
            ("BronchiolarAlveolarDuctJunctionClubCellScgb1a1Niche", 62.0, 20.0, 0.70, 0.0020),
        ]
    }
}

define_anatomical_organ! {
    /// 10. Pancreatic Islets of Langerhans & Acinar Interface (`TA2:3061`, 6 endocrine-exocrine laminae).
    pub struct PancreaticIsletOrgan {
        organ_name: "PancreaticIsletOfLangerhans",
        ta2_code: "TA2:3061",
        center_mm: (18.0, -15.0, -185.0),
        cavity: FascialBodyCavity::RetroperitonealFascialSpace,
        arterial_tree: &[
            "CeliacTrunkAndSuperiorMesentericArtery",
            "SplenicArtery_DorsalAndGreatPancreaticArteries_PancreaticoduodenalArcades",
            "IntralobularInsularAfferentArterioles",
            "GlomerulusLikeIntraIsletFenestratedCapillaryNetwork",
        ],
        venous_tree: &[
            "InsuloAcinarPortalVenules",
            "SplenicAndSuperiorMesentericVeins",
            "HepaticPortalVein -> LiverSinusoids",
        ],
        lymphatic: "PerilobularInterstitialLymphatics -> PancreaticosplenicAndPyloricNodes -> CeliacNodes",
        innervation: "CeliacAndSuperiorMesentericPlexuses_VagalCephalicPhaseInsulin_And_SympatheticSplanchnic",
        cardiac_output_fraction: 0.015,
        barrier: EndothelialBarrierType::FenestratedDiaphragmed,
        youthful_modulus_kpa: 3.0,
        default_layers: [
            ("IntraIsletDoubleBasementMembraneFenestratedCapillary", 6.0, 8.0, 0.70, 0.0024),
            ("CoreInsulinSecretingPdx1MafaBetaCellSyncytium", 16.0, 20.0, 0.84, 0.0034),
            ("MantleGlucagonAlphaAndSomatostatinDeltaParacrineRing", 38.0, 16.0, 0.78, 0.0028),
            ("PancreaticPolypeptideGammaAndEpsilonGhrelinRim", 54.0, 12.0, 0.72, 0.0024),
            ("PeriIsletSchwannGlialAutonomicNeuroinsularSheath", 68.0, 10.0, 0.60, 0.0020),
            ("InsuloAcinarExocrineZymogenCentroacinarDuctalInterface", 82.0, 20.0, 0.80, 0.0029),
        ]
    }
}

define_anatomical_organ! {
    /// 11. Striated Skeletal Muscle Fascicle & Pax7+ Satellite Niche (`TA2:2150`, 6 myofiber & niche laminae).
    pub struct SkeletalMuscleFascicleOrgan {
        organ_name: "SkeletalMuscleFascicle",
        ta2_code: "TA2:2150",
        center_mm: (95.0, 12.0, -320.0),
        cavity: FascialBodyCavity::DeepMuscularFascialCompartment,
        arterial_tree: &[
            "ExternalIliac -> FemoralAndProfundaFemorisArteries",
            "EpimysialAndPerimysialFeedArteries",
            "TransverseTerminalArterioles",
            "LongitudinalEndomysialCapillaryAnastomoses",
        ],
        venous_tree: &[
            "EndomysialPostcapillaryVenules",
            "PerimysialValvedVenousArcades_SkeletalMusclePump",
            "DeepFemoralAndIliacVeins -> InferiorVenaCava",
        ],
        lymphatic: "PerimysialAndEpimysialLymphaticPlexus -> DeepInguinalAndIliacNodes -> LumbarTrunk",
        innervation: "VentralHornAlphaMotorNeuron_CholinergicNeuromuscularJunctionMotorEndplate_And_MuscleSpindleIa",
        cardiac_output_fraction: 0.120,
        barrier: EndothelialBarrierType::ContinuousNonFenestrated,
        youthful_modulus_kpa: 14.0,
        default_layers: [
            ("EndomysialCapillaryMicrovesselWall", 6.0, 8.0, 0.66, 0.0022),
            ("EndomysialReticularLaminaCollagenIIIVIMatrix", 14.0, 10.0, 0.55, 0.0018),
            ("SublaminalPax7PositiveMyf5SatelliteStemCellNiche", 24.0, 12.0, 0.72, 0.0024),
            ("SarcolemmaDystrophinGlycoproteinTtubuleTriadMembrane", 36.0, 12.0, 0.78, 0.0028),
            ("SubsarcolemmalMitochondrialReticulumMyonuclearDomain", 50.0, 20.0, 0.82, 0.0033),
            ("DeepIntermyofibrillarSarcomericActomyosinTitinCore", 74.0, 26.0, 0.86, 0.0035),
        ]
    }
}

define_anatomical_organ! {
    /// 12. Gastrointestinal Crypt of Lieberkühn-Villus Axis (`TA2:2960`, 6 Lgr5+ stem, Paneth & villus laminae).
    pub struct GastrointestinalCryptOrgan {
        organ_name: "IntestinalCryptVillusAxis",
        ta2_code: "TA2:2960",
        center_mm: (12.0, 32.0, -265.0),
        cavity: FascialBodyCavity::IntraperitonealPeritonealCavity,
        arterial_tree: &[
            "AbdominalAorta -> SuperiorMesentericArtery",
            "JejunalAndIlealMesentericArterialArcades -> VasaRecta",
            "SubmucosalArterialPlexusOfHeller",
            "PericryptalAndVillusFountainCapillaryNetwork",
        ],
        venous_tree: &[
            "VillusSubepithelialVenules",
            "SuperiorMesentericVein",
            "HepaticPortalVein -> LiverSinusoids_FirstPass",
        ],
        lymphatic: "CentralVillusChyleLacteal -> SubmucosalAndMesentericLymphNodes -> CisternaChyli -> ThoracicDuct",
        innervation: "EntericNervousSystem_MeissnerSubmucosalAndAuerbachMyentericPlexuses_Plus_VagusNerveCNX",
        cardiac_output_fraction: 0.035,
        barrier: EndothelialBarrierType::FenestratedDiaphragmed,
        youthful_modulus_kpa: 2.8,
        default_layers: [
            ("LaminaPropriaCentralLactealFenestratedCapillaryCore", 6.0, 10.0, 0.65, 0.0020),
            ("PericryptalFoxl1TelocyteMyofibroblastWntRspo3Sheath", 16.0, 12.0, 0.70, 0.0022),
            ("CryptBaseLgr5PositiveColumnarStemAndLysozymePanethNiche", 28.0, 16.0, 0.84, 0.0029),
            ("TransitAmplifyingNotch1ProgenitorMidCryptWall", 44.0, 18.0, 0.82, 0.0028),
            ("DifferentiatedVillusEnterocyteGobletEnteroendocrineEpithelium", 62.0, 22.0, 0.78, 0.0025),
            ("ApicalMicrovilliBrushBorderMuc2GlycocalyxTightJunctionBarrier", 84.0, 12.0, 0.68, 0.0018),
        ]
    }
}

define_anatomical_organ! {
    /// 13. Hypothalamic Neuroendocrine & Circadian Master Pacemaker (`TA2:5562`, 6 SCN, Arcuate, PVN & Median Eminence laminae).
    pub struct HypothalamicNeuroendocrineOrgan {
        organ_name: "HypothalamicNeuroendocrineNuclei",
        ta2_code: "TA2:5562",
        center_mm: (0.0, 8.0, 328.0),
        cavity: FascialBodyCavity::CranialSubarachnoidVault,
        arterial_tree: &[
            "InternalCarotidArtery -> SuperiorHypophysealArtery",
            "CircleOfWillis -> AnteriorAndPosteriorCommunicatingPerforators",
            "MedianEminencePrimaryPortalCapillaryLoops",
        ],
        venous_tree: &[
            "HypothalamoHypophysealLongPortalVenules",
            "CavernousSinus -> InferiorPetrosalSinus -> InternalJugularVein",
        ],
        lymphatic: "ThirdVentricleTanycyteCsfExchange -> BasalMeningealLymphatics -> DeepCervicalNodes",
        innervation: "RetinohypothalamicTract_To_SuprachiasmaticNucleus_And_BrainstemSolitaryTractAfferents",
        cardiac_output_fraction: 0.008,
        barrier: EndothelialBarrierType::FenestratedDiaphragmed,
        youthful_modulus_kpa: 1.4,
        default_layers: [
            ("MedianEminenceFenestratedPrimaryPortalCapillaryPlexus", 5.0, 8.0, 0.65, 0.0024),
            ("ThirdVentricleBetaTanycyteGlialStemBarrier", 13.0, 10.0, 0.72, 0.0025),
            ("ArcuateNucleusAgRP_POMC_KisspeptinGnrhPulseGenerator", 25.0, 16.0, 0.80, 0.0032),
            ("SuprachiasmaticNucleusBmal1ClockMasterPacemakerCore", 42.0, 16.0, 0.82, 0.0033),
            ("ParaventricularNucleusCrhOxytocinVasopressinMagnocellularZone", 60.0, 18.0, 0.78, 0.0031),
            ("DorsomedialAndLateralHypothalamicOrexinSirt1LongevityZone", 80.0, 20.0, 0.74, 0.0029),
        ]
    }
}

define_anatomical_organ! {
    /// 14. Hippocampal Formation & Dentate Gyrus Neurogenic Niche (`TA2:5518`, 6 SGZ stem & Trisynaptic CA3/CA1 laminae).
    pub struct HippocampalDentateGyrusOrgan {
        organ_name: "HippocampalDentateGyrusFormation",
        ta2_code: "TA2:5518",
        center_mm: (26.0, -12.0, 336.0),
        cavity: FascialBodyCavity::CranialSubarachnoidVault,
        arterial_tree: &[
            "PosteriorCerebralArtery_And_AnteriorChoroidalArtery",
            "AnteriorMiddlePosteriorHippocampalArteries",
            "DentateSulcusRakeLikePenetratingArterioles",
            "SubgranularZoneAngiogenicNeurovascularNicheCapillaries",
        ],
        venous_tree: &[
            "IntraventricularSubependymalHippocampalVeins",
            "BasalVeinOfRosenthal -> GreatCerebralVeinOfGalen -> StraightSinus",
        ],
        lymphatic: "PerihippocampalVirchowRobinGlymphatic -> AmbientCistern -> TransverseSinusDuralLymphatics",
        innervation: "MedialSeptalCholinergicFimbriaFornixPathway_And_EntorhinalPerforantPathGlutamatergic",
        cardiac_output_fraction: 0.012,
        barrier: EndothelialBarrierType::ContinuousTightJunctionBbb,
        youthful_modulus_kpa: 1.3,
        default_layers: [
            ("SubgranularZoneBbbPerivascularEndothelialNiche", 5.0, 7.0, 0.68, 0.0026),
            ("SubgranularZoneSox2NestinRadialGliaLikeNeuralStemNiche", 12.0, 12.0, 0.78, 0.0030),
            ("DentateGyrusProx1GranuleCellMossyFiberOriginLayer", 26.0, 18.0, 0.84, 0.0034),
            ("PolymorphicHilarMossyCellGabaergicInterneuronZone", 44.0, 16.0, 0.70, 0.0029),
            ("CornuAmmonisCA3RecurrentAutoassociativePyramidalLayer", 62.0, 20.0, 0.80, 0.0035),
            ("CornuAmmonisCA1SchafferCollateralLtpPyramidalLayer", 84.0, 22.0, 0.82, 0.0036),
        ]
    }
}

define_anatomical_organ! {
    /// 15. Choroid Plexus & Blood-CSF Klotho Secretory Barrier (`TA2:5410`, 6 fenestrated stromal & tight cuboidal ependymal laminae).
    pub struct ChoroidPlexusCsfOrgan {
        organ_name: "ChoroidPlexusBloodCsfBarrier",
        ta2_code: "TA2:5410",
        center_mm: (14.0, -8.0, 355.0),
        cavity: FascialBodyCavity::CranialSubarachnoidVault,
        arterial_tree: &[
            "InternalCarotidArtery -> AnteriorChoroidalArtery",
            "PosteriorCerebralArtery -> PosteriorLateralAndMedialChoroidalArteries",
            "ChoroidalVillousFrondFenestratedCapillaryLoops",
        ],
        venous_tree: &[
            "SuperiorAndInferiorChoroidVeins",
            "InternalCerebralVeins -> GreatCerebralVeinOfGalen",
        ],
        lymphatic: "VentricularCsfSecretion -> ArachnoidGranulationsAndCribriformPlateLymphatics",
        innervation: "SuperiorCervicalSympatheticGanglion_And_GlossopharyngealVagalParasympatheticPlexus",
        cardiac_output_fraction: 0.006,
        barrier: EndothelialBarrierType::FenestratedDiaphragmed,
        youthful_modulus_kpa: 1.8,
        default_layers: [
            ("ChoroidalVillousFenestratedDiaphragmedCapillaryCore", 4.0, 6.0, 0.65, 0.0022),
            ("PerivascularPiaArachnoidLeptomeningealStroma", 10.0, 8.0, 0.55, 0.0018),
            ("StromalHofbauerLikeBorderAssociatedMacrophageSurveillanceZone", 18.0, 10.0, 0.62, 0.0020),
            ("BasolateralInfolding ChoroidalEpithelialBasementMembrane", 28.0, 8.0, 0.70, 0.0024),
            ("CuboidalChoroidEpitheliumClaudin1TightJunctionKlothoTtrFactory", 38.0, 16.0, 0.85, 0.0034),
            ("ApicalMicrovilliCiliaAqp1NaKAtpaseCsfSecretoryBorder", 54.0, 10.0, 0.74, 0.0028),
        ]
    }
}

define_anatomical_organ! {
    /// 16. Pituitary Gland / Hypophysis (`TA2:3851`, 6 adenohypophyseal & neurohypophyseal endocrine laminae).
    pub struct PituitaryHypophysealOrgan {
        organ_name: "PituitaryHypophysealGland",
        ta2_code: "TA2:3851",
        center_mm: (0.0, 12.0, 312.0),
        cavity: FascialBodyCavity::CranialSubarachnoidVault,
        arterial_tree: &[
            "InternalCarotidArtery -> SuperiorAndInferiorHypophysealArteries",
            "HypothalamoHypophysealPortalVenousSinusoids_SecondaryPlexus",
        ],
        venous_tree: &[
            "AdenohypophysealCollectingVenules",
            "CavernousSinus -> PetrosalSinuses -> InternalJugularVein",
        ],
        lymphatic: "SellaTurcicaDiaphragmaSellaeDuralLymphatics -> RetropharyngealAndDeepCervicalNodes",
        innervation: "HypothalamoHypophysealTract_SupraopticAndParaventricularAxons",
        cardiac_output_fraction: 0.004,
        barrier: EndothelialBarrierType::SinusoidalDiscontinuous,
        youthful_modulus_kpa: 2.4,
        default_layers: [
            ("SecondaryPortalFenestratedSinusoidalEndothelium", 5.0, 8.0, 0.68, 0.0022),
            ("FolliculostellateSox2PositivePituitaryStemCellReticulum", 13.0, 10.0, 0.72, 0.0024),
            ("ParsDistalisAcidophilSomatotrophGhAndLactotrophCords", 25.0, 18.0, 0.84, 0.0031),
            ("ParsDistalisBasophilCorticotrophActhGonadotrophThyrotrophZone", 43.0, 18.0, 0.82, 0.0030),
            ("ParsIntermediaMelanotrophPomcRathkeCleftRemnant", 61.0, 10.0, 0.70, 0.0022),
            ("ParsNervosaNeurohypophysisHerringBodyPituicyteAxonalPlexus", 73.0, 18.0, 0.75, 0.0026),
        ]
    }
}

define_anatomical_organ! {
    /// 17. Adrenal / Suprarenal Gland (`TA2:3865`, 6 cortical steroidogenic & medullary chromaffin laminae).
    pub struct AdrenalSuprarenalOrgan {
        organ_name: "AdrenalSuprarenalGland",
        ta2_code: "TA2:3865",
        center_mm: (45.0, -35.0, -172.0),
        cavity: FascialBodyCavity::RetroperitonealFascialSpace,
        arterial_tree: &[
            "InferiorPhrenicArtery -> SuperiorSuprarenalArteries",
            "AbdominalAorta -> MiddleSuprarenalArtery",
            "RenalArtery -> InferiorSuprarenalArtery -> SubcapsularArterialPlexus",
        ],
        venous_tree: &[
            "CorticalFenestratedSinusoids -> MedullaryVenousSinusoids",
            "CentralAdrenalVein -> InferiorVenaCava_Right_Or_RenalVein_Left",
        ],
        lymphatic: "CapsularAndMedullaryLymphatics -> LateralAorticLumbarNodes -> CisternaChyli",
        innervation: "PreganglionicSympatheticGreaterSplanchnicNerveT5T9_DirectToMedullaryChromaffinCells",
        cardiac_output_fraction: 0.005,
        barrier: EndothelialBarrierType::SinusoidalDiscontinuous,
        youthful_modulus_kpa: 3.8,
        default_layers: [
            ("SubcapsularSonicHedgehogGli1AdrenocorticalStemCapsule", 6.0, 10.0, 0.70, 0.0020),
            ("ZonaGlomerulosaCyp11b2AldosteroneMineralocorticoidBalls", 16.0, 14.0, 0.78, 0.0026),
            ("ZonaFasciculataSpongiosaCyp17CortisolGlucocorticoidCords", 32.0, 24.0, 0.82, 0.0031),
            ("ZonaReticularisSult2a1DheaAndrogenAgingSensitiveZone", 56.0, 18.0, 0.80, 0.0030),
            ("CorticomedullaryPortalVenousSinusoidalInterface", 74.0, 10.0, 0.68, 0.0022),
            ("AdrenalMedullaPnmtChromaffinEpinephrineNorepinephrineGanglion", 86.0, 20.0, 0.84, 0.0033),
        ]
    }
}

define_anatomical_organ! {
    /// 18. Splenic White & Red Pulp Lymphoid Organ (`TA2:5172`, 6 marginal sinus, PALS & Billroth cord laminae).
    pub struct SplenicLymphoidPulpOrgan {
        organ_name: "SplenicWhiteAndRedPulp",
        ta2_code: "TA2:5172",
        center_mm: (88.0, -18.0, -152.0),
        cavity: FascialBodyCavity::IntraperitonealPeritonealCavity,
        arterial_tree: &[
            "CeliacTrunk -> TortuousSplenicArtery",
            "TrabecularArteries -> CentralArteriolesWithPalsSheath",
            "PenicillarArterioles -> SheathedEllipsoids_OpenAndClosedCirculation",
        ],
        venous_tree: &[
            "RedPulpVenousSinusoidsWithStaveEndothelium",
            "PulpAndTrabecularVeins -> SplenicVein",
            "HepaticPortalVein -> LiverSinusoids",
        ],
        lymphatic: "WhitePulpPeriarteriolarEfferentLymphatics -> SplenicHilarNodes -> CeliacNodes",
        innervation: "CeliacPlexus_SplenicSympatheticNerve_To_ChATPositiveTcells_CholinergicAntiInflammatoryReflex",
        cardiac_output_fraction: 0.030,
        barrier: EndothelialBarrierType::SinusoidalDiscontinuous,
        youthful_modulus_kpa: 3.2,
        default_layers: [
            ("MarginalSinusMadcam1MetallophilicMacrophageGate", 6.0, 10.0, 0.72, 0.0021),
            ("MarginalZoneIgM_MemoryBcellAndMarcoMacrophageTrap", 16.0, 14.0, 0.78, 0.0023),
            ("PeriarteriolarLymphoidSheathPalsCcl19Ccl21TcellZone", 30.0, 18.0, 0.84, 0.0026),
            ("SplenicFollicularGerminalCenterCxcl13BcellZone", 48.0, 18.0, 0.82, 0.0025),
            ("RedPulpCordsOfBillrothSenescentErythrocyteMacrophageFilter", 66.0, 22.0, 0.76, 0.0022),
            ("RedPulpStaveEndothelialInterendothelialSlitVenousSinusoids", 88.0, 14.0, 0.70, 0.0019),
        ]
    }
}

define_anatomical_organ! {
    /// 19. Visceral Omental/Mesenteric White Adipose Tissue (`TA2:6115`, 6 SVF stem, crown macrophage & adipocyte laminae).
    pub struct VisceralWhiteAdiposeOrgan {
        organ_name: "VisceralWhiteAdiposeDepot",
        ta2_code: "TA2:6115",
        center_mm: (24.0, 48.0, -240.0),
        cavity: FascialBodyCavity::IntraperitonealPeritonealCavity,
        arterial_tree: &[
            "SuperiorAndInferiorMesentericArteries_Plus_GastroepiploicArteries",
            "OmentalAndMesentericArterialArcades",
            "PeriadipocyteContinuousCapillaryMeshwork",
        ],
        venous_tree: &[
            "OmentalAndMesentericVenules",
            "SuperiorMesentericVein -> HepaticPortalVein_FreeFattyAcidAndSaspDirectToLiver",
        ],
        lymphatic: "OmentalMilkySpotsAndMesentericLymphatics -> MesentericNodes -> CisternaChyli",
        innervation: "SympatheticNoradrenergicVaricosities_Beta3LipolyticInnervation",
        cardiac_output_fraction: 0.025,
        barrier: EndothelialBarrierType::ContinuousNonFenestrated,
        youthful_modulus_kpa: 2.0,
        default_layers: [
            ("PeriadipocyteCapillaryEndotheliumCd36FattyAcidTransporter", 6.0, 8.0, 0.64, 0.0016),
            ("StromalVascularFractionPdgfraCd34AdiposeStemProgenitorNiche", 14.0, 12.0, 0.70, 0.0018),
            ("OmentalMilkySpotB1LymphocyteAndM2MacrophageSurveillanceZone", 26.0, 14.0, 0.68, 0.0019),
            ("PericellularCollagenVIMatrixAndCrownLikeCd38MacrophageRing", 40.0, 16.0, 0.58, 0.0017),
            ("UnilocularLipidDropletLeptinAdiponectinEndocrineAdipocytes", 58.0, 26.0, 0.52, 0.0015),
            ("MesothelialVisceralPeritoneumSlipperyLubricinSurface", 84.0, 10.0, 0.65, 0.0014),
        ]
    }
}

define_anatomical_organ! {
    /// 20. Supraclavicular & Paravertebral Thermogenic Brown Adipose Tissue (`TA2:6118`, 6 UCP1+ mitochondrial laminae).
    pub struct ThermogenicBrownAdiposeOrgan {
        organ_name: "ThermogenicBrownAdiposeDepot",
        ta2_code: "TA2:6118",
        center_mm: (42.0, 14.0, 78.0),
        cavity: FascialBodyCavity::SuperiorAnteriorMediastinum,
        arterial_tree: &[
            "SubclavianArtery -> ThyrocervicalTrunk -> SuprascapularAndTransverseCervicalArteries",
            "InternalThoracicAndIntercostalArterioles",
            "HyperdenseThermogenicCapillaryBed",
        ],
        venous_tree: &[
            "SulzerVeinAndSupraclavicularThermogenicVenules",
            "BrachiocephalicAndAzygosVeins -> SuperiorVenaCava_WarmBloodDistribution",
        ],
        lymphatic: "SupraclavicularLymphatics -> InferiorDeepCervicalNodes -> ThoracicOrRightLymphaticDuct",
        innervation: "DenseUnmyelinatedSympatheticTyrosineHydroxylaseFibers_From_StellateGanglion",
        cardiac_output_fraction: 0.010,
        barrier: EndothelialBarrierType::ContinuousNonFenestrated,
        youthful_modulus_kpa: 2.5,
        default_layers: [
            ("HyperdenseThermogenicCapillaryEndothelium", 4.0, 6.0, 0.72, 0.0026),
            ("PerivascularSympatheticNoradrenergicBeta3VaricosityPlexus", 10.0, 8.0, 0.68, 0.0024),
            ("Myf5PositivePrdm16BrownAdipocyteProgenitorStemNiche", 18.0, 12.0, 0.75, 0.0028),
            ("MultilocularLipidDropletUcp1UncouplingMitochondrialCore", 32.0, 22.0, 0.86, 0.0042),
            ("BeigeTmem26_Cd137InducibleThermogenicAdipocyteTransitionZone", 54.0, 18.0, 0.78, 0.0034),
            ("BatokineFgf21Nrg4Il6EndocrineSecretoryInterstitium", 72.0, 14.0, 0.66, 0.0025),
        ]
    }
}

/// Polymorphic enum wrapper dispatching [`OrganArchitecture`] across all 20 canonical *Gray's Anatomy* / `TA2` organs.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum AnatomicalOrgan {
    ArterialVascular(ArterialVascularOrgan),
    CerebralCortex(CerebralCortexOrgan),
    HepaticLobule(HepaticLobuleOrgan),
    BoneMarrowNiche(BoneMarrowNicheOrgan),
    ThymicLymphoid(ThymicLymphoidOrgan),
    DermalCutaneous(DermalCutaneousOrgan),
    MyocardialVentricular(MyocardialVentricularOrgan),
    RenalNephron(RenalNephronOrgan),
    PulmonaryAlveolar(PulmonaryAlveolarOrgan),
    PancreaticIslet(PancreaticIsletOrgan),
    SkeletalMuscleFascicle(SkeletalMuscleFascicleOrgan),
    GastrointestinalCrypt(GastrointestinalCryptOrgan),
    HypothalamicNeuroendocrine(HypothalamicNeuroendocrineOrgan),
    HippocampalDentateGyrus(HippocampalDentateGyrusOrgan),
    ChoroidPlexusCsf(ChoroidPlexusCsfOrgan),
    PituitaryHypophyseal(PituitaryHypophysealOrgan),
    AdrenalSuprarenal(AdrenalSuprarenalOrgan),
    SplenicLymphoidPulp(SplenicLymphoidPulpOrgan),
    VisceralWhiteAdipose(VisceralWhiteAdiposeOrgan),
    ThermogenicBrownAdipose(ThermogenicBrownAdiposeOrgan),
}

macro_rules! delegate_anatomical_organ {
    ($self:ident, $method:ident $(, $arg:expr)*) => {
        match $self {
            AnatomicalOrgan::ArterialVascular(o) => o.$method($($arg),*),
            AnatomicalOrgan::CerebralCortex(o) => o.$method($($arg),*),
            AnatomicalOrgan::HepaticLobule(o) => o.$method($($arg),*),
            AnatomicalOrgan::BoneMarrowNiche(o) => o.$method($($arg),*),
            AnatomicalOrgan::ThymicLymphoid(o) => o.$method($($arg),*),
            AnatomicalOrgan::DermalCutaneous(o) => o.$method($($arg),*),
            AnatomicalOrgan::MyocardialVentricular(o) => o.$method($($arg),*),
            AnatomicalOrgan::RenalNephron(o) => o.$method($($arg),*),
            AnatomicalOrgan::PulmonaryAlveolar(o) => o.$method($($arg),*),
            AnatomicalOrgan::PancreaticIslet(o) => o.$method($($arg),*),
            AnatomicalOrgan::SkeletalMuscleFascicle(o) => o.$method($($arg),*),
            AnatomicalOrgan::GastrointestinalCrypt(o) => o.$method($($arg),*),
            AnatomicalOrgan::HypothalamicNeuroendocrine(o) => o.$method($($arg),*),
            AnatomicalOrgan::HippocampalDentateGyrus(o) => o.$method($($arg),*),
            AnatomicalOrgan::ChoroidPlexusCsf(o) => o.$method($($arg),*),
            AnatomicalOrgan::PituitaryHypophyseal(o) => o.$method($($arg),*),
            AnatomicalOrgan::AdrenalSuprarenal(o) => o.$method($($arg),*),
            AnatomicalOrgan::SplenicLymphoidPulp(o) => o.$method($($arg),*),
            AnatomicalOrgan::VisceralWhiteAdipose(o) => o.$method($($arg),*),
            AnatomicalOrgan::ThermogenicBrownAdipose(o) => o.$method($($arg),*),
        }
    };
}

impl OrganArchitecture for AnatomicalOrgan {
    fn organ_name(&self) -> &'static str {
        delegate_anatomical_organ!(self, organ_name)
    }

    fn ta2_code(&self) -> &'static str {
        delegate_anatomical_organ!(self, ta2_code)
    }

    fn stereotaxic_center_mm(&self) -> AnatomicalPoint3D {
        delegate_anatomical_organ!(self, stereotaxic_center_mm)
    }

    fn fascial_cavity(&self) -> FascialBodyCavity {
        delegate_anatomical_organ!(self, fascial_cavity)
    }

    fn arterial_supply_tree(&self) -> &'static [&'static str] {
        delegate_anatomical_organ!(self, arterial_supply_tree)
    }

    fn venous_drainage_tree(&self) -> &'static [&'static str] {
        delegate_anatomical_organ!(self, venous_drainage_tree)
    }

    fn lymphatic_pathway(&self) -> &'static str {
        delegate_anatomical_organ!(self, lymphatic_pathway)
    }

    fn innervation_pathway(&self) -> &'static str {
        delegate_anatomical_organ!(self, innervation_pathway)
    }

    fn cardiac_output_fraction(&self) -> f64 {
        delegate_anatomical_organ!(self, cardiac_output_fraction)
    }

    fn barrier_type(&self) -> EndothelialBarrierType {
        delegate_anatomical_organ!(self, barrier_type)
    }

    fn layers(&self) -> &[TissueLayer] {
        delegate_anatomical_organ!(self, layers)
    }

    fn ecm(&self) -> &ExtracellularMatrixState {
        delegate_anatomical_organ!(self, ecm)
    }

    fn ecm_mut(&mut self) -> &mut ExtracellularMatrixState {
        delegate_anatomical_organ!(self, ecm_mut)
    }

    fn vascular(&self) -> &VascularPerfusionState {
        delegate_anatomical_organ!(self, vascular)
    }

    fn vascular_mut(&mut self) -> &mut VascularPerfusionState {
        delegate_anatomical_organ!(self, vascular_mut)
    }

    fn biomarker_set(&self) -> OrganBiomarkerSet {
        delegate_anatomical_organ!(self, biomarker_set)
    }

    fn aging_rate_params(&self) -> OrganAgingRateParams {
        delegate_anatomical_organ!(self, aging_rate_params)
    }

    fn intervention_response_coeffs(&self) -> OrganInterventionResponseCoefficients {
        delegate_anatomical_organ!(self, intervention_response_coeffs)
    }

    fn viability_hook(&self) -> f64 {
        delegate_anatomical_organ!(self, viability_hook)
    }

    fn stiffness_hook_kpa(&self) -> f64 {
        delegate_anatomical_organ!(self, stiffness_hook_kpa)
    }

    fn oxygenation_hook_po2(&self) -> f64 {
        delegate_anatomical_organ!(self, oxygenation_hook_po2)
    }

    fn apply_targeted_intervention_response(
        &mut self,
        osk_dose: f64,
        senolytic_dose: f64,
        metabolic_dose: f64,
    ) {
        delegate_anatomical_organ!(
            self,
            apply_targeted_intervention_response,
            osk_dose,
            senolytic_dose,
            metabolic_dose
        )
    }
}

/// Complete 3D Whole-Body Spatial Atlas containing all 20 *Gray's Anatomy* / `TA2` organs and their
/// directed hemodynamic, portal venous, lymphatic, fascial, and neuroendocrine spatial edges.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct AnatomicalSpatialAtlas {
    pub organs: Vec<AnatomicalOrgan>,
    pub spatial_edges: Vec<SpatialRelationshipEdge>,
}

impl AnatomicalSpatialAtlas {
    /// Builds the complete 20-organ human anatomical atlas at chronological `age_years`
    /// with exact 3D Euclidean distances, `TA2` directional vectors, and vascular/lymphatic/neural conduits.
    pub fn from_somatic_age(age_years: f64) -> Self {
        let organs = vec![
            AnatomicalOrgan::ArterialVascular(ArterialVascularOrgan::from_somatic_age(age_years)),
            AnatomicalOrgan::CerebralCortex(CerebralCortexOrgan::from_somatic_age(age_years)),
            AnatomicalOrgan::HepaticLobule(HepaticLobuleOrgan::from_somatic_age(age_years)),
            AnatomicalOrgan::BoneMarrowNiche(BoneMarrowNicheOrgan::from_somatic_age(age_years)),
            AnatomicalOrgan::ThymicLymphoid(ThymicLymphoidOrgan::from_somatic_age(age_years)),
            AnatomicalOrgan::DermalCutaneous(DermalCutaneousOrgan::from_somatic_age(age_years)),
            AnatomicalOrgan::MyocardialVentricular(MyocardialVentricularOrgan::from_somatic_age(
                age_years,
            )),
            AnatomicalOrgan::RenalNephron(RenalNephronOrgan::from_somatic_age(age_years)),
            AnatomicalOrgan::PulmonaryAlveolar(PulmonaryAlveolarOrgan::from_somatic_age(age_years)),
            AnatomicalOrgan::PancreaticIslet(PancreaticIsletOrgan::from_somatic_age(age_years)),
            AnatomicalOrgan::SkeletalMuscleFascicle(SkeletalMuscleFascicleOrgan::from_somatic_age(
                age_years,
            )),
            AnatomicalOrgan::GastrointestinalCrypt(GastrointestinalCryptOrgan::from_somatic_age(
                age_years,
            )),
            AnatomicalOrgan::HypothalamicNeuroendocrine(
                HypothalamicNeuroendocrineOrgan::from_somatic_age(age_years),
            ),
            AnatomicalOrgan::HippocampalDentateGyrus(
                HippocampalDentateGyrusOrgan::from_somatic_age(age_years),
            ),
            AnatomicalOrgan::ChoroidPlexusCsf(ChoroidPlexusCsfOrgan::from_somatic_age(age_years)),
            AnatomicalOrgan::PituitaryHypophyseal(PituitaryHypophysealOrgan::from_somatic_age(
                age_years,
            )),
            AnatomicalOrgan::AdrenalSuprarenal(AdrenalSuprarenalOrgan::from_somatic_age(
                age_years,
            )),
            AnatomicalOrgan::SplenicLymphoidPulp(SplenicLymphoidPulpOrgan::from_somatic_age(
                age_years,
            )),
            AnatomicalOrgan::VisceralWhiteAdipose(VisceralWhiteAdiposeOrgan::from_somatic_age(
                age_years,
            )),
            AnatomicalOrgan::ThermogenicBrownAdipose(
                ThermogenicBrownAdiposeOrgan::from_somatic_age(age_years),
            ),
        ];

        let make_edge = |src_idx: usize,
                         dst_idx: usize,
                         conduit: SpatialConnectivityConduit|
         -> SpatialRelationshipEdge {
            let src = &organs[src_idx];
            let dst = &organs[dst_idx];
            SpatialRelationshipEdge {
                source_organ_ta2: src.ta2_code().to_string(),
                source_organ_name: src.organ_name().to_string(),
                target_organ_ta2: dst.ta2_code().to_string(),
                target_organ_name: dst.organ_name().to_string(),
                euclidean_distance_mm: src.euclidean_distance_to_mm(dst),
                directional_relations: src.directional_relations_to(dst),
                conduit,
            }
        };

        // Canonical Whole-Body 20-Organ Spatial & Closed-Loop Hemodynamic/Endocrine Connectivity Graph:
        // Indices: 0=Artery, 1=BrainCortex, 2=Liver, 3=BoneMarrow, 4=Thymus, 5=Skin,
        //          6=HeartMyocardium, 7=Kidney, 8=Lung, 9=Pancreas, 10=Muscle, 11=Gut,
        //          12=Hypothalamus, 13=HippocampusDG, 14=ChoroidPlexus, 15=Pituitary,
        //          16=Adrenal, 17=Spleen, 18=VisceralWAT, 19=SupraclavicularBAT
        let spatial_edges = vec![
            // 1. Pulmonary-Cardiac Oxygenation & Closed Venous Return Loop (Lung 8 <-> Heart 6 -> Artery 0)
            make_edge(
                8,
                6,
                SpatialConnectivityConduit::ArterialHemodynamicConduit {
                    vessel_pathway: "FourPulmonaryVeins -> LeftAtrium -> MitralValve -> LeftVentricle".to_string(),
                    cardiac_output_fraction: 1.0,
                    transit_time_sec: 4.5,
                },
            ),
            make_edge(
                6,
                0,
                SpatialConnectivityConduit::ContiguousFascialPlane {
                    shared_cavity: FascialBodyCavity::MiddleMediastinumPericardialSac,
                    contact_surface_cm2: 185.0,
                },
            ),
            make_edge(
                2,
                8,
                SpatialConnectivityConduit::VenousFirstPassConduit {
                    vein_pathway: "HepaticVeins -> InferiorVenaCava -> RightHeart -> PulmonaryArteryToAlveolarCapillaries".to_string(),
                    flow_fraction: 0.28,
                    transit_time_sec: 6.8,
                },
            ),
            make_edge(
                7,
                8,
                SpatialConnectivityConduit::VenousFirstPassConduit {
                    vein_pathway: "RenalVeins -> InferiorVenaCava -> RightHeart -> PulmonaryCirculation".to_string(),
                    flow_fraction: 0.20,
                    transit_time_sec: 7.2,
                },
            ),
            // 2. Systemic Arterial Distribution from Ascending/Descending Aorta (0 -> All Target Organs)
            make_edge(
                0,
                1,
                SpatialConnectivityConduit::ArterialHemodynamicConduit {
                    vessel_pathway: "ArchOfAorta -> InternalCarotidAndVertebralArteries -> CircleOfWillis".to_string(),
                    cardiac_output_fraction: 0.14,
                    transit_time_sec: 6.2,
                },
            ),
            make_edge(
                0,
                2,
                SpatialConnectivityConduit::ArterialHemodynamicConduit {
                    vessel_pathway: "AbdominalAorta -> CeliacTrunk -> ProperHepaticArtery".to_string(),
                    cardiac_output_fraction: 0.065,
                    transit_time_sec: 8.5,
                },
            ),
            make_edge(
                0,
                7,
                SpatialConnectivityConduit::ArterialHemodynamicConduit {
                    vessel_pathway: "AbdominalAorta -> BilateralRenalArteries -> GlomerularAfferents".to_string(),
                    cardiac_output_fraction: 0.20,
                    transit_time_sec: 7.0,
                },
            ),
            make_edge(
                0,
                3,
                SpatialConnectivityConduit::ArterialHemodynamicConduit {
                    vessel_pathway: "InternalThoracicAndIntercostalArteries -> MedullaryNutrientArteries".to_string(),
                    cardiac_output_fraction: 0.04,
                    transit_time_sec: 9.0,
                },
            ),
            make_edge(
                0,
                4,
                SpatialConnectivityConduit::ArterialHemodynamicConduit {
                    vessel_pathway: "InternalThoracicArtery -> AnteriorMediastinalThymicBranches".to_string(),
                    cardiac_output_fraction: 0.015,
                    transit_time_sec: 4.8,
                },
            ),
            make_edge(
                0,
                5,
                SpatialConnectivityConduit::ArterialHemodynamicConduit {
                    vessel_pathway: "MusculocutaneousPerforators -> SubpapillaryDermalPlexus".to_string(),
                    cardiac_output_fraction: 0.075,
                    transit_time_sec: 14.0,
                },
            ),
            make_edge(
                0,
                10,
                SpatialConnectivityConduit::ArterialHemodynamicConduit {
                    vessel_pathway: "FemoralAndProfundaFemorisArteries -> EndomysialCapillaryPlexus".to_string(),
                    cardiac_output_fraction: 0.15,
                    transit_time_sec: 12.5,
                },
            ),
            make_edge(
                0,
                17,
                SpatialConnectivityConduit::ArterialHemodynamicConduit {
                    vessel_pathway: "AbdominalAorta -> CeliacTrunk -> TortuousSplenicArtery -> MarginalZoneSinuses".to_string(),
                    cardiac_output_fraction: 0.04,
                    transit_time_sec: 7.5,
                },
            ),
            // 3. Hepatic Portal Venous First-Pass Axis (Gut 11, Pancreas 9, Spleen 17, Visceral WAT 18 -> Liver 2)
            make_edge(
                11,
                2,
                SpatialConnectivityConduit::VenousFirstPassConduit {
                    vein_pathway: "SuperiorMesentericVein -> HepaticPortalVein -> PeriportalZone1Sinusoids".to_string(),
                    flow_fraction: 0.75,
                    transit_time_sec: 5.5,
                },
            ),
            make_edge(
                9,
                2,
                SpatialConnectivityConduit::VenousFirstPassConduit {
                    vein_pathway: "PancreaticoduodenalAndSplenicVeins -> HepaticPortalVein (DirectFirstPassInsulinGlucagon)".to_string(),
                    flow_fraction: 0.95,
                    transit_time_sec: 3.8,
                },
            ),
            make_edge(
                17,
                2,
                SpatialConnectivityConduit::VenousFirstPassConduit {
                    vein_pathway: "SplenicVein -> HepaticPortalVein (HemeBilirubinAndCytokinePortalClearance)".to_string(),
                    flow_fraction: 0.90,
                    transit_time_sec: 4.4,
                },
            ),
            make_edge(
                18,
                2,
                SpatialConnectivityConduit::VenousFirstPassConduit {
                    vein_pathway: "RightAndLeftGastroepiploicOmentalVeins -> HepaticPortalVein (DirectPortalFFAsAndSaspIL6)".to_string(),
                    flow_fraction: 0.85,
                    transit_time_sec: 5.1,
                },
            ),
            // 4. Choroid Plexus CSF-Glymphatic-Hippocampus-Meningeal Axis
            make_edge(
                14,
                13,
                SpatialConnectivityConduit::LymphaticDrainageConduit {
                    lymphatic_chain: "ChoroidPlexusVentricularCSF -> SubependymalSubgranularZoneKlothoAndBDNFBath".to_string(),
                    clearance_rate_ml_hr: 21.0,
                },
            ),
            make_edge(
                14,
                1,
                SpatialConnectivityConduit::LymphaticDrainageConduit {
                    lymphatic_chain: "SubarachnoidCSF -> VirchowRobinPeriarterialInflux -> AQP4GlymphaticParenchymalSweep".to_string(),
                    clearance_rate_ml_hr: 19.5,
                },
            ),
            make_edge(
                1,
                4,
                SpatialConnectivityConduit::LymphaticDrainageConduit {
                    lymphatic_chain: "VirchowRobinGlymphatic -> MeningealLymphatics -> DeepCervicalNodes -> SuperiorMediastinum".to_string(),
                    clearance_rate_ml_hr: 18.5,
                },
            ),
            make_edge(
                3,
                4,
                SpatialConnectivityConduit::VenousFirstPassConduit {
                    vein_pathway: "BoneMarrowSinusoids -> ETP_CommonLymphoidProgenitorHoming -> ThymicCorticomedullaryHEV".to_string(),
                    flow_fraction: 0.65,
                    transit_time_sec: 11.0,
                },
            ),
            // 5. Hypothalamic-Pituitary-Adrenal (HPA), Somatotropic (GH/IGF-1), Klotho & Sympathetic BAT Axes
            make_edge(
                12,
                15,
                SpatialConnectivityConduit::VenousFirstPassConduit {
                    vein_pathway: "MedianEminencePrimaryCapillaryPlexus -> LongHypophysealPortalVessels -> Adenohypophysis (CRH_GHRH_GnRH)".to_string(),
                    flow_fraction: 0.98,
                    transit_time_sec: 1.2,
                },
            ),
            make_edge(
                15,
                16,
                SpatialConnectivityConduit::NeuroendocrineReflexArc {
                    nerve_or_hormone_tract: "CorticotrophACTH_SystemicEndocrineCascade -> AdrenalZonaFasciculataMC2R_CortisolDHEAS".to_string(),
                    coupling_gain: 0.91,
                },
            ),
            make_edge(
                15,
                2,
                SpatialConnectivityConduit::NeuroendocrineReflexArc {
                    nerve_or_hormone_tract: "SomatotrophPulsatileGH -> HepatocyteGHR_JAK2_STAT5b_SystemicIGF1Synthesis".to_string(),
                    coupling_gain: 0.87,
                },
            ),
            make_edge(
                7,
                13,
                SpatialConnectivityConduit::NeuroendocrineReflexArc {
                    nerve_or_hormone_tract: "RenalDistalTubuleSecretedAlphaKlotho -> HippocampalDentateGyrusGluN2BAndFGF23Protection".to_string(),
                    coupling_gain: 0.84,
                },
            ),
            make_edge(
                12,
                19,
                SpatialConnectivityConduit::NeuroendocrineReflexArc {
                    nerve_or_hormone_tract: "PreopticHypothalamus -> SympatheticStellateGanglionBeta3Adrenergic -> SupraclavicularBAT_UCP1Thermogenesis".to_string(),
                    coupling_gain: 0.88,
                },
            ),
            make_edge(
                1,
                11,
                SpatialConnectivityConduit::NeuroendocrineReflexArc {
                    nerve_or_hormone_tract: "DorsalMotorNucleusOfVagusCNX -> EntericMyentericPlexus_MicrobiomeBrainGutAxis".to_string(),
                    coupling_gain: 0.82,
                },
            ),
            make_edge(
                12,
                3,
                SpatialConnectivityConduit::NeuroendocrineReflexArc {
                    nerve_or_hormone_tract: "SuprachiasmaticNucleus -> SympatheticBeta3Adrenergic_CircadianHscCxcl12NicheControl".to_string(),
                    coupling_gain: 0.76,
                },
            ),
        ];

        Self {
            organs,
            spatial_edges,
        }
    }

    /// Looks up all outgoing spatial/vascular/lymphatic/neural edges originating from `source_ta2`.
    pub fn edges_from_ta2(&self, source_ta2: &str) -> Vec<&SpatialRelationshipEdge> {
        self.spatial_edges
            .iter()
            .filter(|e| e.source_organ_ta2 == source_ta2)
            .collect()
    }

    /// Advances chronological aging by `delta_years` across all 20 organs in the spatial atlas.
    pub fn apply_chronological_aging(&mut self, delta_years: f64) {
        for organ in &mut self.organs {
            organ.apply_chronological_aging(delta_years);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_oop_organ_architecture_biodistribution_and_spatial_layers() {
        let brain = CerebralCortexOrgan::from_somatic_age(72.0);
        let liver = HepaticLobuleOrgan::from_somatic_age(72.0);
        let artery = ArterialVascularOrgan::from_somatic_age(72.0);

        // Verify Gray's Anatomy / TA2 metadata and ultra-fine histological layer counts (6-8 layers)
        assert_eq!(brain.ta2_code(), "TA2:5434");
        assert_eq!(brain.layers().len(), 8);
        assert_eq!(artery.ta2_code(), "TA2:3891");
        assert_eq!(artery.layers().len(), 7);
        assert_eq!(liver.ta2_code(), "TA2:3025");
        assert_eq!(liver.layers().len(), 6);

        // 1. Organ-specific Young's modulus hierarchy: Artery (>> 28 kPa) > Liver (~4.2 kPa) > Brain (~1.5 kPa)
        assert!(artery.ecm().youngs_modulus_kpa > liver.ecm().youngs_modulus_kpa * 4.0);
        assert!(liver.ecm().youngs_modulus_kpa > brain.ecm().youngs_modulus_kpa * 2.0);

        // 2. Without BBB receptor shuttle, sinusoidal liver has vastly higher LNP (800 kDa) uptake than BBB-gated brain
        let lnp_mw = 800_000.0;
        let liver_bioavail = liver.mean_parenchymal_bioavailability(lnp_mw, 0.0);
        let brain_unshuttled = brain.mean_parenchymal_bioavailability(lnp_mw, 0.0);
        let brain_tfr_shuttled = brain.mean_parenchymal_bioavailability(lnp_mw, 0.92);

        assert!(liver_bioavail > brain_unshuttled * 3.5);
        assert!(brain_tfr_shuttled > brain_unshuttled * 3.0);

        // 3. Rejuvenating organ architecture restores pO2 and lowers Young's modulus & YAP/TAZ stress
        let mut aged_artery = AnatomicalOrgan::ArterialVascular(artery);
        let before_po2 = aged_artery.mean_parenchymal_po2_mmhg();
        let before_e = aged_artery.ecm().youngs_modulus_kpa;
        let before_yap = aged_artery.ecm().yap_taz_mechanotransduction_stress();

        aged_artery.rejuvenate_architecture(1.0, 1.0, 1.0);
        assert!(aged_artery.mean_parenchymal_po2_mmhg() > before_po2 + 6.0);
        assert!(aged_artery.ecm().youngs_modulus_kpa < before_e * 0.70);
        assert!(aged_artery.ecm().yap_taz_mechanotransduction_stress() < before_yap * 0.35);
    }

    #[test]
    fn test_whole_body_20_organ_spatial_atlas_and_ta2_directional_relations() {
        let atlas = AnatomicalSpatialAtlas::from_somatic_age(72.0);
        assert_eq!(atlas.organs.len(), 20);
        assert_eq!(atlas.spatial_edges.len(), 27);

        let brain = &atlas.organs[1];
        let liver = &atlas.organs[2];
        let skin = &atlas.organs[5];
        let hypothalamus = &atlas.organs[12];
        let pituitary = &atlas.organs[15];

        // Brain (Z = +385 mm) is SuperiorTo Liver (Z = -145 mm)
        let brain_to_liver = brain.directional_relations_to(liver);
        assert!(brain_to_liver.contains(&DirectionalAxisRelation::SuperiorTo));
        assert!(brain_to_liver.contains(&DirectionalAxisRelation::MedialTo));
        assert!(brain.euclidean_distance_to_mm(liver) > 500.0);

        // Hypothalamus (Z = +322 mm) is SuperiorTo Pituitary (Z = +308 mm)
        let hypo_to_pit = hypothalamus.directional_relations_to(pituitary);
        assert!(hypo_to_pit.contains(&DirectionalAxisRelation::SuperiorTo));
        assert!(hypothalamus.euclidean_distance_to_mm(pituitary) < 25.0);

        // Skin (radial distance > 150 mm) is SuperficialTo Thymus (radial distance 38 mm)
        let thymus = &atlas.organs[4];
        let skin_to_thymus = skin.directional_relations_to(thymus);
        assert!(skin_to_thymus.contains(&DirectionalAxisRelation::SuperficialTo));
        assert!(skin_to_thymus.contains(&DirectionalAxisRelation::LateralTo));
    }

    #[test]
    fn test_six_humanz_organs_completeness_and_intervention_response() {
        let mut atlas = AnatomicalSpatialAtlas::from_somatic_age(65.0);
        // Indices: 0: Artery, 1: Cortex, 2: Liver, 3: Marrow, 4: Thymus, 5: Skin
        for i in 0..6 {
            let organ = &atlas.organs[i];
            let bm = organ.biomarker_set();
            assert!(bm.viability > 0.0 && bm.viability <= 1.0);
            assert!(bm.youngs_modulus_kpa > 0.5);
            assert!(bm.interstitial_po2_mmhg > 10.0 && bm.interstitial_po2_mmhg < 60.0);
            assert!(bm.sasp_secretion_rate >= 0.01 && bm.sasp_secretion_rate <= 1.0);
            assert!(bm.functional_reserve >= 0.05 && bm.functional_reserve <= 1.0);
            assert!(bm.senescent_cell_fraction >= 0.01 && bm.senescent_cell_fraction <= 1.0);

            let aging_params = organ.aging_rate_params();
            assert!(aging_params.baseline_senescence_rate_per_year > 0.0);
            assert!(aging_params.ecm_crosslink_accumulation_rate_per_year > 0.0);
            assert!(aging_params.microvascular_rarefaction_rate_per_year > 0.0);
            assert!(aging_params.parenchymal_loss_rate_per_year > 0.0);

            let resp = organ.intervention_response_coeffs();
            assert!(resp.osk_epigenetic_responsiveness > 0.5);
            assert!(resp.senolytic_clearance_selectivity > 0.5);
            assert!(resp.metabolic_sirt_responsiveness > 0.5);
        }

        // Test applying targeted intervention response to liver
        let liver = &mut atlas.organs[2];
        let stiff_before = liver.stiffness_hook_kpa();
        let viab_before = liver.viability_hook();
        liver.apply_targeted_intervention_response(0.9, 0.9, 0.9);
        assert!(liver.viability_hook() > viab_before);
        assert!(liver.stiffness_hook_kpa() < stiff_before);
    }
}
