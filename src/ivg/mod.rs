//! In Vitro Gametogenesis (IVG), Somatic Cell Nuclear Transfer (SCNT), and Synthetic Embryology.
//!
//! This module provides models and simulation tools for:
//! - Reprogramming stem cells into synthetic gametes (PGCLCs -> synthetic eggs & sperm).
//! - Reconstituted Ovarian Organoid (`rOvary`) 3-stage folliculogenesis (`IVD -> IVG -> IVM`).
//! - Epigenetic imprint erasure and DNA demethylation modeling.
//! - Somatic Cell Nuclear Transfer (SCNT) and nuclear reprogramming barriers (H3K9me3 removal via KDM4).
//! - Synthetic embryogenesis (Zygote -> Blastocyst -> Gastrula) and artificial womb compatibility.

pub mod embryo;
pub mod epigenetics;
pub mod gamete;
pub mod mitochondria;
pub mod organoid;
pub mod scnt;

pub use embryo::{
    ArtificialWombChamber, DevelopmentalStage, EmbryoCultureEnvironment, EmbryoOrigin,
    SyntheticEmbryo,
};
pub use epigenetics::{EpigeneticProfile, HistoneMark, ImprintPattern};
pub use gamete::{
    derive_homologous_partner, recombine_homologous_chromosomes,
    simulate_karyotype_meiotic_crossovers, species_bivalent_pairs, GameteQuality, GameteType,
    InductionProtocol, PgclcCulture, StemCellState, SyntheticGamete, ZonaPellucida,
};
pub use mitochondria::{
    simulate_germline_mitochondrial_bottleneck, simulate_pink1_parkin_mitophagy,
    MitochondrialBottleneckReport, MitochondrialPopulation, MitochondrialPopulationState,
};
pub use organoid::{
    simulate_rovary_folliculogenesis, FollicleStage, OogenesisOrganoidReport,
    OvarianOrganoidConfig,
};
pub use scnt::{ActivationMethod, EnucleatedOocyte, ScntConstruct, SomaticDonorCell};
