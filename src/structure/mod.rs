//! Macromolecular 3D structure analysis, PDB/mmCIF parsing, RMSD, TM-score, and contact maps.

pub mod contact;
pub mod pdb;
pub mod rmsd;
pub mod secondary;

pub use contact::{calculate_ca_contact_map, calculate_ca_distance_matrix, ContactMap3D, DistanceMatrix3D};
pub use pdb::{parse_mmcif, parse_pdb, Atom3D, Chain3D, ProteinStructure3D, Residue3D};
pub use rmsd::{apply_rotation, calculate_rmsd_and_tm_score, SuperpositionResult};
pub use secondary::{assign_secondary_structure, SecondaryStructure};
