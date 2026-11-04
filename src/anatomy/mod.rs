//! # Anatomical Architecture, Sub-Microscopic Histology & 3D Spatial Atlas (`anatomy`)
//!
//! Provides an Object-Oriented (`OrganArchitecture` trait + 20 concrete *Gray's Anatomy* / *Terminologia Anatomica* `TA2` organ classes)
//! framework modeling:
//! - **Extracellular Matrix (`ECM`) Biomechanics** ([`ecm`]): Collagen I/III fibrosis, elastin
//!   fragmentation, glucosepane Advanced Glycation End-product (`AGE`) cross-links, Young's elastic
//!   modulus ($E$, $\text{kPa}$), `YAP`/`TAZ` mechanotransduction stress, and Ogston porous-gel
//!   macromolecular steric hindrance.
//! - **Spatial Reaction-Diffusion PDE & Barrier Transport** ([`transport`]): Krogh-Erlang radial
//!   oxygenation ($p\text{O}_2(r)$), Thiele-modulus interstitial drug/vector gradients, and
//!   organ-specific endothelial permeability (`Blood-Brain Barrier`, `Continuous`, `Fenestrated`,
//!   `Sinusoidal`).
//! - **3D Stereotaxic Coordinate Frame & Spatial Conduits** ([`spatial`]): Right-Anterior-Superior (`RAS`, $\text{mm}$)
//!   coordinates, `TA2` directional relation classification (`SuperiorTo`, `AnteriorTo`, `MedialTo`, `SuperficialTo`),
//!   coelomic/fascial cavities, and arterial/portal-venous/glymphatic/neuroendocrine conduits.
//! - **Whole-Body PBPK Biodistribution, Moens-Korteweg Hemodynamics & Endocrine Axes** ([`physiology`]):
//!   Closed-loop Pulse Wave Velocity (`PWV`), hepatic/splenic `RES` macrophage first-pass trapping,
//!   and `HPA`, `GH`/`IGF-1`, `RAAS`, and `α-Klotho` systemic feedback loops.
//! - **20 Ultra-Stratified Organ Classes (6–8 Histological Laminae Each)** ([`organ`]):
//!   [`ArterialVascularOrgan`], [`CerebralCortexOrgan`], [`HepaticLobuleOrgan`], [`BoneMarrowNicheOrgan`],
//!   [`ThymicLymphoidOrgan`], [`DermalCutaneousOrgan`], [`MyocardialVentricularOrgan`],
//!   [`RenalNephronOrgan`], [`PulmonaryAlveolarOrgan`], [`PancreaticIsletOrgan`],
//!   [`SkeletalMuscleFascicleOrgan`], [`GastrointestinalCryptOrgan`], [`HypothalamicNeuroendocrineOrgan`],
//!   [`HippocampalDentateGyrusOrgan`], [`ChoroidPlexusCsfOrgan`], [`PituitaryHypophysealOrgan`],
//!   [`AdrenalSuprarenalOrgan`], [`SplenicLymphoidPulpOrgan`], [`VisceralWhiteAdiposeOrgan`],
//!   [`ThermogenicBrownAdiposeOrgan`], and [`AnatomicalSpatialAtlas`].
//! - **Organ-Specific Physiological Behaviors** ([`behavior`]): Frank-Starling cardiac cycle,
//!   alveolar gas exchange, zonated hepatic metabolism, glomerular eGFR, and SCN circadian regulation.
//! - **Closed-Loop Whole-Body Systemic Cycles** ([`cycles`]): 7 major physiological feedback circuits
//!   integrating cardiovascular, cardiopulmonary, enterohepatic, renal-RAAS, neuroendocrine-glymphatic,
//!   immune-lymphatic, and energy-thermoregulatory systems ([`WholeBodySystemicCircuits`]).

pub mod anchors;
pub mod behavior;
pub mod coupling;
pub mod cycles;
pub mod ecm;
pub mod organ;
pub mod physiology;
pub mod spatial;
pub mod transport;

pub use anchors::{
    audit_clinical_reference_anchors, ClinicalAnchorCheck, ClinicalBiomarkerValues,
    ClinicalReferenceAuditReport, ReferenceRange, REF_EGFR_YOUTHFUL, REF_FASTING_GLUCOSE_NORMAL,
    REF_OGTT_2H_GLUCOSE_NORMAL, REF_SASP_YOUTHFUL, REF_SKIN_STIFFNESS_YOUTH,
    REF_SYSTOLIC_BP_NORMAL, REF_TELOMERE_AGED, REF_TELOMERE_MIDDLE_AGE, REF_TELOMERE_YOUNG_ADULT,
    REF_THYMIC_NAIVE_T_OUTPUT_YOUTH,
};
pub use behavior::{
    AdrenalSteroidMetrics, ArterialHemodynamicsMetrics, BrownAdiposeMetrics, CardiacCycleMetrics,
    ChoroidCsfMetrics, CorticalCognitiveMetrics, CutaneousBarrierMetrics,
    GastrointestinalFunctionMetrics, HematopoieticNicheMetrics, HepaticMetabolicMetrics,
    HippocampalFunctionMetrics, HypothalamicRegulationMetrics, MuscleFunctionMetrics,
    OrganBehavioralSnapshot, PancreaticIsletMetrics, PituitaryEndocrineMetrics,
    PulmonaryGasExchangeMetrics, RenalFiltrationMetrics, SplenicImmunityMetrics,
    ThymicSelectionMetrics, VisceralAdiposeMetrics,
};
pub use coupling::{propagate_cross_organ_coupling, SystemicCouplingMetrics};
pub use cycles::{
    CardiopulmonaryCycle, CardiovascularCycle, EnergyThermoregulatoryCycle,
    EnterohepaticMetabolicCycle, ImmuneLymphaticCycle, NeuroendocrineGlymphaticCycle,
    RenalCardiovascularRaasCycle, WholeBodySystemicCircuits,
};
pub use ecm::ExtracellularMatrixState;
pub use organ::{
    AdrenalSuprarenalOrgan, AnatomicalOrgan, AnatomicalSpatialAtlas, ArterialVascularOrgan,
    BoneMarrowNicheOrgan, CerebralCortexOrgan, ChoroidPlexusCsfOrgan, DermalCutaneousOrgan,
    GastrointestinalCryptOrgan, HepaticLobuleOrgan, HippocampalDentateGyrusOrgan,
    HypothalamicNeuroendocrineOrgan, MyocardialVentricularOrgan, OrganAgingRateParams,
    OrganArchitecture, OrganBiomarkerSet, OrganInterventionResponseCoefficients,
    PancreaticIsletOrgan, PituitaryHypophysealOrgan, PulmonaryAlveolarOrgan, RenalNephronOrgan,
    SkeletalMuscleFascicleOrgan, SplenicLymphoidPulpOrgan, ThermogenicBrownAdiposeOrgan,
    ThymicLymphoidOrgan, VisceralWhiteAdiposeOrgan,
};
pub use physiology::{
    compute_whole_body_pbpk_biodistribution, HemodynamicBiomechanicsState,
    OrganBiodistributionShare, SystemicEndocrineAxes,
};
pub use spatial::{
    classify_directional_relations, AnatomicalPoint3D, DirectionalAxisRelation, FascialBodyCavity,
    SpatialConnectivityConduit, SpatialRelationshipEdge,
};
pub use transport::{
    solve_layer_reaction_diffusion, EndothelialBarrierType, SpatialLayerProfile, TissueLayer,
};
