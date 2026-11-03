//! 3D Stereotaxic Coordinate Frame (`RAS` in $\text{mm}$), Topological Anatomical Relations,
//! Fascial Compartments, and Whole-Body Vascular / Lymphatic / Autonomic Neural Conduits.
//!
//! Provides spatial relationship modeling beyond *Gray's Anatomy* (42nd Edition) and *Terminologia Anatomica* (`TA2`):
//! - **3D Stereotaxic Coordinate System**: Right-Anterior-Superior (`RAS`) somatic coordinates in $\text{mm}$
//!   with exact Euclidean distance, bounding ellipsoids, and directional classification
//!   (`Superior`/`Inferior`, `Anterior`/`Posterior`, `Medial`/`Lateral`, `Superficial`/`Deep`).
//! - **Anatomical Conduits**: Arterial supply trees (`Aorta` $\to$ `Arterioles`), venous drainage
//!   (`Portal` vs `Caval` first-pass routing), lymphatic drainage (`Glymphatic`, `Thoracic Duct`,
//!   `Right Lymphatic Duct`), autonomic innervation (`Vagal Cholinergic` vs `Sympathetic`), and
//!   contiguous serous/fascial cavities.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// 3D anatomical coordinates in standard Right-Anterior-Superior (`RAS`) somatic stereotaxic space ($\text{mm}$).
/// Origin $(0, 0, 0)$ is calibrated at the mid-sagittal sternal angle / T4 vertebral body plane.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct AnatomicalPoint3D {
    /// Medio-lateral coordinate ($+X = \text{Left lateral}$, $-X = \text{Right lateral}$, $0 = \text{Mid-sagittal plane}$, in $\text{mm}$).
    pub x_mediolateral_mm: f64,
    /// Ventro-dorsal coordinate ($+Y = \text{Anterior/Ventral}$, $-Y = \text{Posterior/Dorsal}$, in $\text{mm}$).
    pub y_anteroposterior_mm: f64,
    /// Cranio-caudal coordinate ($+Z = \text{Superior/Cranial}$, $-Z = \text{Inferior/Caudal}$, in $\text{mm}$).
    pub z_superoinferior_mm: f64,
}

impl AnatomicalPoint3D {
    pub const fn new(x_mm: f64, y_mm: f64, z_mm: f64) -> Self {
        Self {
            x_mediolateral_mm: x_mm,
            y_anteroposterior_mm: y_mm,
            z_superoinferior_mm: z_mm,
        }
    }

    /// Computes exact 3D Euclidean distance ($\text{mm}$) between two anatomical loci.
    pub fn distance_mm(&self, other: &Self) -> f64 {
        let dx = self.x_mediolateral_mm - other.x_mediolateral_mm;
        let dy = self.y_anteroposterior_mm - other.y_anteroposterior_mm;
        let dz = self.z_superoinferior_mm - other.z_superoinferior_mm;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Radial distance ($\text{mm}$) from the central longitudinal body axis ($X=0, Y=0$).
    pub fn radial_distance_from_midline_mm(&self) -> f64 {
        (self.x_mediolateral_mm.powi(2) + self.y_anteroposterior_mm.powi(2)).sqrt()
    }
}

/// Canonical *Terminologia Anatomica* (`TA2`) directional relationship between two structures.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum DirectionalAxisRelation {
    /// Structure A is cranial / superior to Structure B ($\Delta Z > +15\text{ mm}$).
    SuperiorTo,
    /// Structure A is caudal / inferior to Structure B ($\Delta Z < -15\text{ mm}$).
    InferiorTo,
    /// Structure A is ventral / anterior to Structure B ($\Delta Y > +10\text{ mm}$).
    AnteriorTo,
    /// Structure A is dorsal / posterior to Structure B ($\Delta Y < -10\text{ mm}$).
    PosteriorTo,
    /// Structure A is closer to the median sagittal plane than Structure B.
    MedialTo,
    /// Structure A is farther from the median sagittal plane than Structure B.
    LateralTo,
    /// Structure A is superficial (closer to skin surface) relative to Structure B.
    SuperficialTo,
    /// Structure A is deep (closer to visceral core) relative to Structure B.
    DeepTo,
}

/// Classifies all applicable *Terminologia Anatomica* directional relations of `source` relative to `target`.
pub fn classify_directional_relations(
    source: &AnatomicalPoint3D,
    target: &AnatomicalPoint3D,
) -> Vec<DirectionalAxisRelation> {
    let mut rels = Vec::with_capacity(4);
    let dz = source.z_superoinferior_mm - target.z_superoinferior_mm;
    if dz > 15.0 {
        rels.push(DirectionalAxisRelation::SuperiorTo);
    } else if dz < -15.0 {
        rels.push(DirectionalAxisRelation::InferiorTo);
    }

    let dy = source.y_anteroposterior_mm - target.y_anteroposterior_mm;
    if dy > 10.0 {
        rels.push(DirectionalAxisRelation::AnteriorTo);
    } else if dy < -10.0 {
        rels.push(DirectionalAxisRelation::PosteriorTo);
    }

    let dm = source.x_mediolateral_mm.abs() - target.x_mediolateral_mm.abs();
    if dm < -10.0 {
        rels.push(DirectionalAxisRelation::MedialTo);
    } else if dm > 10.0 {
        rels.push(DirectionalAxisRelation::LateralTo);
    }

    let dr = source.radial_distance_from_midline_mm() - target.radial_distance_from_midline_mm();
    if dr > 25.0 {
        rels.push(DirectionalAxisRelation::SuperficialTo);
    } else if dr < -25.0 {
        rels.push(DirectionalAxisRelation::DeepTo);
    }

    rels
}

/// Coelomic / fascial body cavity enclosure defined in *Gray's Anatomy*.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum FascialBodyCavity {
    /// Supratentorial & infratentorial cranial vault enclosed by dura/arachnoid/pia mater.
    CranialSubarachnoidVault,
    /// Superior & anterior mediastinum retrosternal compartment.
    SuperiorAnteriorMediastinum,
    /// Middle mediastinum enclosed by fibrous and serous pericardium.
    MiddleMediastinumPericardialSac,
    /// Bilateral pleural cavities enclosed by visceral and parietal pleura.
    ThoracicPleuralCavity,
    /// Intraperitoneal subphrenic and mesenteric cavity.
    IntraperitonealPeritonealCavity,
    /// Retroperitoneal pararenal and perirenal fascial space (Gerota's fascia).
    RetroperitonealFascialSpace,
    /// Endosteal and medullary trabecular cavity inside axial/appendicular bones.
    EndostealMedullaryBoneCavity,
    /// Deep muscular fascial compartment (epimysium / deepinvesting fascia).
    DeepMuscularFascialCompartment,
    /// Superficial cutaneous & subcutaneous integumentary mantle.
    SuperficialIntegumentaryMantle,
}

/// Physical, vascular, lymphatic, fascial, or neural spatial connection between two organs.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum SpatialConnectivityConduit {
    /// Direct arterial hemodynamic perfusion conduit (`source_organ` $\to$ `target_organ`).
    ArterialHemodynamicConduit {
        vessel_pathway: String,
        cardiac_output_fraction: f64,
        transit_time_sec: f64,
    },
    /// Hepatic portal or systemic caval venous drainage conduit (`source_organ` $\to$ `target_organ`).
    VenousFirstPassConduit {
        vein_pathway: String,
        flow_fraction: f64,
        transit_time_sec: f64,
    },
    /// Lymphatic & glymphatic interstitial drainage (`source_organ` $\to$ `target_organ`).
    LymphaticDrainageConduit {
        lymphatic_chain: String,
        clearance_rate_ml_hr: f64,
    },
    /// Contiguous fascial plane or serous cavity allowing direct paracrine SASP diffusion.
    ContiguousFascialPlane {
        shared_cavity: FascialBodyCavity,
        contact_surface_cm2: f64,
    },
    /// Autonomic / neuroendocrine efferent-afferent reflex arc (`source_organ` $\to$ `target_organ`).
    NeuroendocrineReflexArc {
        nerve_or_hormone_tract: String,
        coupling_gain: f64,
    },
}

/// Directed edge in the whole-body 3D anatomical spatial graph.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SpatialRelationshipEdge {
    pub source_organ_ta2: String,
    pub source_organ_name: String,
    pub target_organ_ta2: String,
    pub target_organ_name: String,
    pub euclidean_distance_mm: f64,
    pub directional_relations: Vec<DirectionalAxisRelation>,
    pub conduit: SpatialConnectivityConduit,
}

impl SpatialRelationshipEdge {
    /// Computes the paracrine / hemodynamic signal transmission weight ($0.0$ to $1.0$)
    /// across this anatomical relationship edge.
    pub fn signal_transmission_weight(&self) -> f64 {
        let distance_decay = (-self.euclidean_distance_mm / 450.0).exp().clamp(0.10, 1.0);
        match &self.conduit {
            SpatialConnectivityConduit::ArterialHemodynamicConduit {
                cardiac_output_fraction,
                transit_time_sec,
                ..
            } => {
                let transit_factor = (-transit_time_sec / 30.0).exp().clamp(0.25, 1.0);
                (0.35 + 1.8 * cardiac_output_fraction) * transit_factor
            }
            SpatialConnectivityConduit::VenousFirstPassConduit {
                flow_fraction,
                transit_time_sec,
                ..
            } => {
                let transit_factor = (-transit_time_sec / 25.0).exp().clamp(0.30, 1.0);
                (0.45 * flow_fraction * transit_factor).clamp(0.05, 0.95)
            }
            SpatialConnectivityConduit::LymphaticDrainageConduit {
                clearance_rate_ml_hr,
                ..
            } => (0.15 + 0.008 * clearance_rate_ml_hr).clamp(0.05, 0.75) * distance_decay,
            SpatialConnectivityConduit::ContiguousFascialPlane {
                contact_surface_cm2,
                ..
            } => ((contact_surface_cm2 / 120.0).clamp(0.08, 0.85)) * distance_decay,
            SpatialConnectivityConduit::NeuroendocrineReflexArc { coupling_gain, .. } => {
                coupling_gain.clamp(0.05, 0.95)
            }
        }
        .clamp(0.02, 0.99)
    }
}
