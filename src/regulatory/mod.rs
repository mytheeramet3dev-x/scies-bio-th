//! Gene regulatory logic: JASPAR TFBS scanning, core promoter elements, H3K27ac enhancer scoring,
//! and Gene Regulatory Networks (GRN DAGs).

pub mod enhancer;
pub mod grn;
pub mod promoter;
pub mod tfbs;

pub use enhancer::{score_enhancer_activity, EnhancerAssessment, EnhancerClass};
pub use grn::{GeneRegulatoryNetwork, RegulationType, RegulatoryEdge};
pub use promoter::{
    detect_core_promoter_elements, CorePromoterElement, CorePromoterMotifType,
};
pub use tfbs::{
    jaspar_ma0003_1_tfap2a, jaspar_ma0079_3_sp1, scan_jaspar_tfbs, JasparMatrix, TfbsHit,
};
