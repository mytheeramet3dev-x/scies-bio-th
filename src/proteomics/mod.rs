//! Proteomics, enzymatic digestion simulation, and peptide mass spectrometry tools.
//!
//! Provides in silico proteolytic cleavage (Trypsin, Chymotrypsin, LysC),
//! monoisotopic and average molecular mass calculation, and isoelectric point (pI) estimation.

use crate::seq::Protein;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Proteolytic enzymes used in mass spectrometry sample preparation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Enzyme {
    /// Trypsin: cleaves C-terminal to Lys (K) or Arg (R), except when followed by Pro (P).
    Trypsin,
    /// Chymotrypsin: cleaves C-terminal to Phe (F), Trp (W), or Tyr (Y), except before Pro (P).
    Chymotrypsin,
    /// Endoproteinase Lys-C: cleaves strictly C-terminal to Lys (K).
    LysC,
}

/// A digested peptide fragment with sequence, coordinates, and physical properties.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Peptide {
    pub sequence: Protein,
    pub start: usize,
    pub end: usize,
    pub missed_cleavages: usize,
    pub mass_da: f64,
}

/// Monoisotopic residue mass table (Da) for standard 20 amino acids.
fn amino_acid_monoisotopic_mass(aa: u8) -> f64 {
    match aa.to_ascii_uppercase() {
        b'A' => 71.03711,
        b'R' => 156.10111,
        b'N' => 114.04293,
        b'D' => 115.02694,
        b'C' => 103.00919,
        b'E' => 129.04259,
        b'Q' => 128.05858,
        b'G' => 57.02146,
        b'H' => 137.05891,
        b'I' | b'L' => 113.08406,
        b'K' => 128.09496,
        b'M' => 131.04049,
        b'F' => 147.06841,
        b'P' => 97.05276,
        b'S' => 87.03203,
        b'T' => 101.04768,
        b'W' => 186.07931,
        b'Y' => 163.06333,
        b'V' => 99.06841,
        _ => 110.0, // Average generic placeholder
    }
}

/// Calculates the monoisotopic mass of a peptide or protein sequence in Daltons (Da).
///
/// Sums residue masses plus terminal $H_2O$ (+18.01056 Da).
pub fn monoisotopic_mass(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let residue_sum: f64 = seq.iter().map(|&aa| amino_acid_monoisotopic_mass(aa)).sum();
    residue_sum + 18.010565 // Terminal H- and -OH
}

/// Estimates the isoelectric point (pI) of a protein sequence using the bisection method.
pub fn isoelectric_point(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 7.0;
    }

    let mut low = 2.0;
    let mut high = 14.0;

    for _ in 0..50 {
        let mid = (low + high) / 2.0;
        let charge = calculate_net_charge(seq, mid);
        if charge > 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }

    (low + high) / 2.0
}

/// Calculates net charge of a peptide at a given pH.
fn calculate_net_charge(seq: &[u8], ph: f64) -> f64 {
    // pKa values: Lehninger biochemistry standards
    let n_term_pka = 9.69;
    let c_term_pka = 2.34;

    let mut charge = 0.0;
    // N-terminus
    charge += 1.0 / (1.0 + 10.0_f64.powf(ph - n_term_pka));
    // C-terminus
    charge -= 1.0 / (1.0 + 10.0_f64.powf(c_term_pka - ph));

    for &aa in seq {
        match aa.to_ascii_uppercase() {
            b'D' => charge -= 1.0 / (1.0 + 10.0_f64.powf(3.86 - ph)), // Asp
            b'E' => charge -= 1.0 / (1.0 + 10.0_f64.powf(4.25 - ph)), // Glu
            b'C' => charge -= 1.0 / (1.0 + 10.0_f64.powf(8.33 - ph)), // Cys
            b'Y' => charge -= 1.0 / (1.0 + 10.0_f64.powf(10.07 - ph)), // Tyr
            b'H' => charge += 1.0 / (1.0 + 10.0_f64.powf(ph - 6.00)), // His
            b'K' => charge += 1.0 / (1.0 + 10.0_f64.powf(ph - 10.53)), // Lys
            b'R' => charge += 1.0 / (1.0 + 10.0_f64.powf(ph - 12.48)), // Arg
            _ => {}
        }
    }

    charge
}

/// Simulates in silico proteolytic cleavage with support for missed cleavages.
pub fn digest(protein: &Protein, enzyme: Enzyme, max_missed_cleavages: usize) -> Vec<Peptide> {
    let bytes = &protein.sequence;
    let len = bytes.len();
    if len == 0 {
        return Vec::new();
    }

    // 1. Identify all cleavage sites (end index of cleaved fragment)
    let mut cut_indices = Vec::new();
    cut_indices.push(0); // Start of protein

    for i in 0..len {
        let curr = bytes[i];
        let next = if i + 1 < len { Some(bytes[i + 1]) } else { None };

        let is_cut = match enzyme {
            Enzyme::Trypsin => {
                (curr == b'K' || curr == b'R') && next != Some(b'P')
            }
            Enzyme::Chymotrypsin => {
                (curr == b'F' || curr == b'W' || curr == b'Y') && next != Some(b'P')
            }
            Enzyme::LysC => curr == b'K',
        };

        if is_cut && (i + 1 < len) {
            cut_indices.push(i + 1);
        }
    }
    cut_indices.push(len); // End of protein

    let num_fragments = cut_indices.len() - 1;
    let mut peptides = Vec::new();

    for missed in 0..=max_missed_cleavages {
        for i in 0..num_fragments {
            if i + missed < num_fragments {
                let start = cut_indices[i];
                let end = cut_indices[i + missed + 1];
                let slice = &bytes[start..end];
                if let Ok(pep_prot) = Protein::new(slice) {
                    let mass = monoisotopic_mass(slice);
                    peptides.push(Peptide {
                        sequence: pep_prot,
                        start,
                        end,
                        missed_cleavages: missed,
                        mass_da: mass,
                    });
                }
            }
        }
    }

    peptides
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trypsin_cleavage_rule() {
        // M-K-A-R-P-A-R-G:
        // Cleavage after K (pos 2).
        // Arg at pos 4 is followed by Pro (P) -> blocked!
        // Cleavage after Arg at pos 7.
        // Fragments: MK, ARPAR, G
        let protein = Protein::new(b"MKARPARG").unwrap();
        let peps = digest(&protein, Enzyme::Trypsin, 0);

        assert_eq!(peps.len(), 3);
        assert_eq!(peps[0].sequence.sequence, b"MK");
        assert_eq!(peps[1].sequence.sequence, b"ARPAR");
        assert_eq!(peps[2].sequence.sequence, b"G");
    }

    #[test]
    fn test_missed_cleavages() {
        let protein = Protein::new(b"AKARAG").unwrap();
        // Fully cleaved (missed 0): AK, AR, AG (3)
        // Missed 1: AKAR, ARAG (2)
        let peps = digest(&protein, Enzyme::Trypsin, 1);
        assert_eq!(peps.len(), 5);
    }

    #[test]
    fn test_monoisotopic_mass() {
        // Glycine (57.02146) + Alanine (71.03711) + H2O (18.010565) = 146.069135 Da
        let mass = monoisotopic_mass(b"GA");
        assert!((mass - 146.069135).abs() < 1e-4);
    }

    #[test]
    fn test_isoelectric_point_basic() {
        // Poly-Lysine should have very basic pI (> 9.0)
        let poly_k = b"KKKKKKKK";
        let pi = isoelectric_point(poly_k);
        assert!(pi > 9.5);

        // Poly-Aspartic acid should have very acidic pI (< 4.5)
        let poly_d = b"DDDDDDDD";
        let pi_d = isoelectric_point(poly_d);
        assert!(pi_d < 4.5);
    }
}
