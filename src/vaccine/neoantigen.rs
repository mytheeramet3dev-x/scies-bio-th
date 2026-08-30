//! Tumor-specific neoantigen identification and mutant peptide processing.
//!
//! Neoantigens arise from somatic mutations in cancer cells and represent foreign,
//! tumor-exclusive epitopes that can be recognized by cytotoxic T cells without inducing autoimmunity.

use crate::error::{BioError, Result};
use crate::seq::Protein;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A somatic cancer mutation producing a novel peptide sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Neoantigen {
    /// Associated gene symbol (e.g. "TP53", "KRAS", "BRAF").
    pub gene_symbol: String,
    /// Wild-type flanking peptide window (e.g. 25-27 amino acids).
    pub wt_peptide: Protein,
    /// Mutated neoantigen peptide window containing the altered amino acid.
    pub mut_peptide: Protein,
    /// 0-based position of the mutated residue within the window peptide.
    pub mutation_offset_in_window: usize,
    /// Original wild-type amino acid single-letter code.
    pub wt_aa: char,
    /// Novel mutant amino acid single-letter code.
    pub mut_aa: char,
    /// 1-based original position within the full-length protein.
    pub original_protein_pos: usize,
}

impl Neoantigen {
    /// Constructs a `Neoantigen` by extracting a window of amino acids around the mutation.
    ///
    /// `flank_size`: number of residues on each side (e.g. 13 for a 27-mer peptide).
    pub fn from_mutation(
        gene_symbol: &str,
        full_wt_protein: &Protein,
        aa_pos_1based: usize,
        mut_aa: char,
        flank_size: usize,
    ) -> Result<Self> {
        let seq = &full_wt_protein.sequence;
        let total_len = seq.len();

        if aa_pos_1based == 0 || aa_pos_1based > total_len {
            return Err(BioError::ParseError(format!(
                "Mutation position {} is out of bounds for protein of length {}",
                aa_pos_1based, total_len
            )));
        }

        let idx_0based = aa_pos_1based - 1;
        let wt_aa = seq[idx_0based] as char;

        let start = idx_0based.saturating_sub(flank_size);
        let end = (idx_0based + flank_size + 1).min(total_len);

        let wt_window_slice = &seq[start..end];
        let offset = idx_0based - start;

        let mut mut_window_slice = wt_window_slice.to_vec();
        mut_window_slice[offset] = mut_aa as u8;

        let wt_peptide = Protein::new(wt_window_slice)?;
        let mut_peptide = Protein::new(&mut_window_slice)?;

        Ok(Self {
            gene_symbol: gene_symbol.to_string(),
            wt_peptide,
            mut_peptide,
            mutation_offset_in_window: offset,
            wt_aa,
            mut_aa,
            original_protein_pos: aa_pos_1based,
        })
    }

    /// Generates all overlapping 9-mer peptides that encompass the mutated amino acid.
    ///
    /// Returns a vector of tuples: `(offset_in_window, wt_9mer, mut_9mer)`.
    pub fn generate_9mers(&self) -> Vec<(usize, Protein, Protein)> {
        let mut results = Vec::new();
        let mut_bytes = &self.mut_peptide.sequence;
        let wt_bytes = &self.wt_peptide.sequence;
        let len = mut_bytes.len();
        let k = 9;

        if len < k {
            return results;
        }

        for i in 0..=(len - k) {
            let kmer_end = i + k;
            // The 9-mer must contain the mutated residue
            if i <= self.mutation_offset_in_window && self.mutation_offset_in_window < kmer_end {
                if let (Ok(wt_9mer), Ok(mut_9mer)) = (
                    Protein::new(&wt_bytes[i..kmer_end]),
                    Protein::new(&mut_bytes[i..kmer_end]),
                ) {
                    results.push((i, wt_9mer, mut_9mer));
                }
            }
        }

        results
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kras_g12d_neoantigen_extraction() {
        // KRAS first 20 aa: MTEYKLVVVG AGGVGKSALT
        let kras_wt = Protein::new(b"MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSY").unwrap();
        // Gly at pos 12 mutated to Asp (G12D)
        let neo = Neoantigen::from_mutation("KRAS", &kras_wt, 12, 'D', 6).unwrap();

        assert_eq!(neo.wt_aa, 'G');
        assert_eq!(neo.mut_aa, 'D');
        assert_eq!(neo.gene_symbol, "KRAS");

        let nine_mers = neo.generate_9mers();
        assert!(!nine_mers.is_empty());

        // Check that all 9-mers contain 'D' in the mutant peptide
        for (_, wt, mut_kmer) in &nine_mers {
            assert!(mut_kmer.sequence.contains(&(b'D')));
            assert!(!wt.sequence.contains(&(b'D')));
        }
    }
}
