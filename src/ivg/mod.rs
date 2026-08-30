//! In Vitro Gametogenesis (IVG), Somatic Cell Nuclear Transfer (SCNT), and Synthetic Embryology.
//!
//! This module provides models and simulation tools for:
//! - Reprogramming stem cells into synthetic gametes (PGCLCs -> synthetic eggs & sperm).
//! - Epigenetic imprint erasure and DNA demethylation modeling.
//! - Somatic Cell Nuclear Transfer (SCNT) and nuclear reprogramming barriers (H3K9me3 removal via KDM4).
//! - Synthetic embryogenesis (Zygote -> Blastocyst -> Gastrula) and artificial womb compatibility.

pub mod embryo;
pub mod epigenetics;
pub mod gamete;
pub mod scnt;

pub use embryo::{
    ArtificialWombChamber, DevelopmentalStage, EmbryoCultureEnvironment, EmbryoOrigin,
    SyntheticEmbryo,
};
pub use epigenetics::{EpigeneticProfile, HistoneMark, ImprintPattern};
pub use gamete::{
    GameteQuality, GameteType, InductionProtocol, PgclcCulture, StemCellState, SyntheticGamete,
    ZonaPellucida,
};
pub use scnt::{ActivationMethod, EnucleatedOocyte, ScntConstruct, SomaticDonorCell};
