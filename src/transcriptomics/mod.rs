//! Transcriptomics quantification and expression normalization tools.
//!
//! Provides functions for calculating Transcripts Per Million (TPM),
//! Reads Per Kilobase Million (RPKM/FPKM), and Counts Per Million (CPM).

use crate::error::{BioError, Result};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Quantified expression level of a transcript.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TranscriptExpression {
    pub id: String,
    pub raw_count: u64,
    pub length_bp: u64,
    pub tpm: f64,
    pub rpkm: f64,
    pub cpm: f64,
}

/// Calculates Transcripts Per Million (TPM) from raw read counts and transcript lengths (bp).
///
/// Formula:
/// 1. $r_i = \text{count}_i / (\text{length}_i / 1000.0) = \text{rate (RPK)}$
/// 2. $S = \sum r_i$
/// 3. $\text{TPM}_i = (r_i / S) \times 10^6$
pub fn calculate_tpm(counts: &[u64], lengths_bp: &[u64]) -> Result<Vec<f64>> {
    if counts.len() != lengths_bp.len() {
        return Err(BioError::InvalidInput(format!(
            "Counts length ({}) does not match lengths ({})",
            counts.len(),
            lengths_bp.len()
        )));
    }
    if counts.is_empty() {
        return Ok(Vec::new());
    }

    let mut rates = Vec::with_capacity(counts.len());
    let mut sum_rates = 0.0;

    for (&count, &len) in counts.iter().zip(lengths_bp.iter()) {
        if len == 0 {
            return Err(BioError::InvalidInput(
                "Transcript length cannot be zero".to_string(),
            ));
        }
        let rate = count as f64 / (len as f64 / 1000.0);
        rates.push(rate);
        sum_rates += rate;
    }

    if sum_rates == 0.0 {
        return Ok(vec![0.0; counts.len()]);
    }

    let scaling_factor = 1_000_000.0 / sum_rates;
    let tpm = rates.into_iter().map(|r| r * scaling_factor).collect();
    Ok(tpm)
}

/// Calculates Reads Per Kilobase Million (RPKM / FPKM).
///
/// Formula:
/// $\text{RPKM}_i = \frac{\text{count}_i}{(\text{length}_i / 1000.0) \times (N / 10^6)}$
/// where $N = \sum \text{count}_j$.
pub fn calculate_rpkm(counts: &[u64], lengths_bp: &[u64]) -> Result<Vec<f64>> {
    if counts.len() != lengths_bp.len() {
        return Err(BioError::InvalidInput(format!(
            "Counts length ({}) does not match lengths ({})",
            counts.len(),
            lengths_bp.len()
        )));
    }
    if counts.is_empty() {
        return Ok(Vec::new());
    }

    let total_mapped_reads: u64 = counts.iter().sum();
    if total_mapped_reads == 0 {
        return Ok(vec![0.0; counts.len()]);
    }

    let per_million_factor = total_mapped_reads as f64 / 1_000_000.0;

    let mut rpkm = Vec::with_capacity(counts.len());
    for (&count, &len) in counts.iter().zip(lengths_bp.iter()) {
        if len == 0 {
            return Err(BioError::InvalidInput(
                "Transcript length cannot be zero".to_string(),
            ));
        }
        let kb = len as f64 / 1000.0;
        let value = (count as f64) / (kb * per_million_factor);
        rpkm.push(value);
    }

    Ok(rpkm)
}

/// Calculates Counts Per Million (CPM).
///
/// Formula:
/// $\text{CPM}_i = (\text{count}_i / N) \times 10^6$
pub fn calculate_cpm(counts: &[u64]) -> Result<Vec<f64>> {
    let total: u64 = counts.iter().sum();
    if total == 0 {
        return Ok(vec![0.0; counts.len()]);
    }
    let factor = 1_000_000.0 / total as f64;
    Ok(counts.iter().map(|&c| c as f64 * factor).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_tpm_sum_to_million() {
        let counts = [10, 20, 30];
        let lengths = [1000, 2000, 3000];
        let tpm = calculate_tpm(&counts, &lengths).unwrap();

        assert_eq!(tpm.len(), 3);
        // Each RPK: 10/1 = 10, 20/2 = 10, 30/3 = 10 -> total 30.
        // Each TPM: (10/30) * 1,000,000 = 333333.33...
        for &val in &tpm {
            assert!((val - 333_333.333_3).abs() < 1.0);
        }
        let total_tpm: f64 = tpm.iter().sum();
        assert!((total_tpm - 1_000_000.0).abs() < 1e-6);
    }

    #[test]
    fn test_calculate_rpkm() {
        let counts = [50, 150];
        let lengths = [1000, 3000];
        let rpkm = calculate_rpkm(&counts, &lengths).unwrap();

        // total reads = 200 -> per million = 200 / 1e6 = 0.0002
        // RPKM1 = 50 / (1 * 0.0002) = 250,000
        // RPKM2 = 150 / (3 * 0.0002) = 250,000
        assert_eq!(rpkm.len(), 2);
        assert!((rpkm[0] - 250_000.0).abs() < 1e-6);
        assert!((rpkm[1] - 250_000.0).abs() < 1e-6);
    }

    #[test]
    fn test_calculate_cpm() {
        let counts = [100, 300];
        let cpm = calculate_cpm(&counts).unwrap();
        assert_eq!(cpm.len(), 2);
        assert_eq!(cpm[0], 250_000.0);
        assert_eq!(cpm[1], 750_000.0);
    }

    #[test]
    fn test_invalid_lengths() {
        let counts = [10, 20];
        let lengths = [1000];
        assert!(calculate_tpm(&counts, &lengths).is_err());
        assert!(calculate_rpkm(&counts, &lengths).is_err());
    }
}
