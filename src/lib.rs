//! scies-bio-th — A high-performance bio toolkit for Rust.

#[cfg(feature = "aging")]
pub mod aging;
#[cfg(feature = "aging")]
pub mod anatomy;
pub mod align;
pub mod annotation;
pub mod codon;
pub mod edit;
pub mod error;
pub mod genome;
pub mod genomics;
pub mod hmm;
pub mod io;
pub mod metagenomics;
pub mod motif;
pub mod phylo;
pub mod popgen;
pub mod proteomics;
pub mod seq;
pub mod transcriptomics;
pub mod vaccine;
pub mod variant;

#[cfg(feature = "structure")]
pub mod structure;

#[cfg(feature = "network")]
pub mod network;

#[cfg(feature = "stats")]
pub mod stats;

#[cfg(feature = "docking")]
pub mod docking;

#[cfg(feature = "epigenome")]
pub mod epigenome;

#[cfg(feature = "regulatory")]
pub mod regulatory;

#[cfg(feature = "multiomics")]
pub mod multiomics;

#[cfg(feature = "safety")]
pub mod safety;

#[cfg(feature = "imaging")]
pub mod imaging;

#[cfg(feature = "singlecell")]
pub mod singlecell;

#[cfg(feature = "pharma")]
pub mod pharma;

#[cfg(feature = "experimental")]
pub mod ivg;

pub use error::{BioError, Result};
pub use seq::{Dna, Kmer, Protein, Rna};
