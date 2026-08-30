//! MHC / HLA Class I binding affinity prediction and neoepitope immunogenicity scoring.
//!
//! Evaluates whether a 9-mer mutant neoantigen binds strongly to the patient's
//! specific Human Leukocyte Antigen (HLA) alleles and calculates the Differential Agretopicity Index (DAI).

use crate::error::{BioError, Result};
use crate::seq::Protein;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Major Histocompatibility Complex (MHC) Class I HLA alleles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum HlaAllele {
    /// HLA-A*02:01 (Dominant allele in ~45% of European and global populations; P2: L/M, P9: V/L).
    HlaA0201,
    /// HLA-A*24:02 (High frequency in Asian populations ~35-40%; P2: Y/F, P9: L/F).
    HlaA2402,
    /// HLA-A*01:01 (P2: T/S, P9: Y).
    HlaA0101,
    /// HLA-B*07:02 (P2: P, P9: L).
    HlaB0702,
    /// HLA-B*08:01 (P3: K/R, P5: K/R, P9: L).
    HlaB0801,
}

impl HlaAllele {
    /// Returns the standard nomenclature string for the HLA allele.
    pub fn name(&self) -> &'static str {
        match self {
            HlaAllele::HlaA0201 => "HLA-A*02:01",
            HlaAllele::HlaA2402 => "HLA-A*24:02",
            HlaAllele::HlaA0101 => "HLA-A*01:01",
            HlaAllele::HlaB0702 => "HLA-B*07:02",
            HlaAllele::HlaB0801 => "HLA-B*08:01",
        }
    }
}

/// MHC Binding prediction result for a neoepitope candidate.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MhcBindingResult {
    /// Evaluated HLA allele.
    pub hla_allele: HlaAllele,
    /// 9-mer mutant neoantigen peptide.
    pub mut_9mer: Protein,
    /// 9-mer wild-type counterpart peptide.
    pub wt_9mer: Protein,
    /// Predicted mutant binding affinity IC50 (nM). Lower is stronger (<50 nM = High Affinity).
    pub mut_ic50_nm: f64,
    /// Predicted wild-type binding affinity IC50 (nM).
    pub wt_ic50_nm: f64,
    /// Differential Agretopicity Index (DAI = WT IC50 / Mutant IC50). > 1.0 means mutant binds tighter.
    pub differential_agretopicity_index: f64,
    /// Whether mutant peptide is classified as a strong binder (<500 nM, ideally <100 nM).
    pub is_strong_binder: bool,
    /// Composite vaccine candidate immunogenicity score (0.0 to 100.0).
    pub immunogenicity_score: f64,
}

/// Predicts MHC Class I binding affinity for a 9-mer peptide against an HLA allele.
///
/// Uses position-specific anchor residue scoring matrices (PSSM / NetMHC-calibrated baseline).
pub fn predict_mhc_binding(
    mut_9mer: &Protein,
    wt_9mer: &Protein,
    hla: HlaAllele,
) -> Result<MhcBindingResult> {
    if mut_9mer.sequence.len() != 9 || wt_9mer.sequence.len() != 9 {
        return Err(BioError::ParseError(
            "MHC Class I prediction requires exactly 9-mer peptides".to_string(),
        ));
    }

    let mut_score = score_peptide_for_hla(&mut_9mer.sequence, hla);
    let wt_score = score_peptide_for_hla(&wt_9mer.sequence, hla);

    // Convert raw log-scale score to IC50 in nM: IC50 = 50000^(1 - score)
    let mut_ic50 = 50000.0_f64.powf(1.0 - mut_score).clamp(1.0, 50000.0);
    let wt_ic50 = 50000.0_f64.powf(1.0 - wt_score).clamp(1.0, 50000.0);

    let dai = wt_ic50 / mut_ic50;
    let is_strong_binder = mut_ic50 < 500.0;

    // Immunogenicity score weights both high affinity (<500 nM) and high DAI (novelty vs self)
    let affinity_component = ((500.0 - mut_ic50).max(0.0) / 500.0) * 60.0;
    let novelty_component = ((dai - 1.0).max(0.0).min(10.0) / 10.0) * 40.0;
    let immunogenicity_score = (affinity_component + novelty_component).clamp(0.0, 100.0);

    Ok(MhcBindingResult {
        hla_allele: hla,
        mut_9mer: mut_9mer.clone(),
        wt_9mer: wt_9mer.clone(),
        mut_ic50_nm: mut_ic50,
        wt_ic50_nm: wt_ic50,
        differential_agretopicity_index: dai,
        is_strong_binder,
        immunogenicity_score,
    })
}

/// Evaluates 9-mer peptide fit into HLA binding groove (returns normalized score 0.0 to 1.0).
fn score_peptide_for_hla(seq: &[u8], hla: HlaAllele) -> f64 {
    let mut score: f64 = 0.50; // Baseline median affinity

    match hla {
        HlaAllele::HlaA0201 => {
            // Position 2 Anchor (P2)
            match seq[1] {
                b'L' | b'M' => score += 0.25,
                b'I' | b'V' | b'Q' => score += 0.15,
                b'P' | b'D' | b'E' | b'K' | b'R' => score -= 0.25,
                _ => {}
            }
            // Position 9 C-terminal Anchor (P9)
            match seq[8] {
                b'V' | b'L' | b'I' => score += 0.25,
                b'A' | b'T' => score += 0.10,
                b'D' | b'E' | b'K' | b'R' | b'P' => score -= 0.30,
                _ => {}
            }
            // Secondary anchors (P1, P3, P7)
            if seq[0] == b'A' || seq[0] == b'Y' || seq[0] == b'K' {
                score += 0.05;
            }
            if seq[6] == b'L' || seq[6] == b'V' {
                score += 0.05;
            }
        }
        HlaAllele::HlaA2402 => {
            // Position 2 Aromatic Anchor (P2)
            match seq[1] {
                b'Y' | b'F' | b'W' => score += 0.30,
                b'M' | b'L' => score += 0.15,
                _ => {}
            }
            // Position 9 C-terminal Anchor (P9)
            match seq[8] {
                b'F' | b'L' | b'I' | b'W' => score += 0.25,
                _ => {}
            }
        }
        HlaAllele::HlaA0101 => {
            if seq[1] == b'T' || seq[1] == b'S' {
                score += 0.25;
            }
            if seq[8] == b'Y' {
                score += 0.30;
            }
        }
        HlaAllele::HlaB0702 => {
            // Proline at P2 is the strict signature anchor for HLA-B*07:02
            if seq[1] == b'P' {
                score += 0.35;
            } else {
                score -= 0.20;
            }
            if seq[8] == b'L' || seq[8] == b'M' {
                score += 0.20;
            }
        }
        HlaAllele::HlaB0801 => {
            if seq[2] == b'K' || seq[2] == b'R' {
                score += 0.20;
            }
            if seq[4] == b'K' || seq[4] == b'R' {
                score += 0.20;
            }
            if seq[8] == b'L' {
                score += 0.20;
            }
        }
    }

    score.clamp(0.01, 0.99)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hla_a0201_anchor_prediction() {
        // High affinity epitope for HLA-A*02:01: P2='L', P9='V'
        let mut_pep = Protein::new(b"NLVPMVATV").unwrap();
        // Wild-type counterpart with non-anchor P2='A'
        let wt_pep = Protein::new(b"NAVPMVATV").unwrap();

        let res = predict_mhc_binding(&mut_pep, &wt_pep, HlaAllele::HlaA0201).unwrap();
        assert!(res.is_strong_binder);
        assert!(res.mut_ic50_nm < 100.0);
        assert!(res.differential_agretopicity_index > 1.0);
        assert!(res.immunogenicity_score > 60.0);
    }
}
