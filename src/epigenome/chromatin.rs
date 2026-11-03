//! Chromatin accessibility modeling and ATAC-seq peak analysis.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// An open chromatin peak (e.g. from ATAC-seq or DNase-seq).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ChromatinPeak {
    /// Chromosome or contig name.
    pub chromosome: String,
    /// 0-based start coordinate.
    pub start: usize,
    /// 0-based end coordinate.
    pub end: usize,
    /// Peak enrichment signal value (e.g. fold enrichment).
    pub signal_value: f64,
    /// Statistical significance (-log10 p-value).
    pub p_value_neg_log10: f64,
    /// False discovery rate (-log10 q-value).
    pub q_value_neg_log10: f64,
}

impl ChromatinPeak {
    /// Peak width in base pairs.
    pub fn length(&self) -> usize {
        self.end.saturating_sub(self.start)
    }

    /// Summit coordinate approximation (center of peak).
    pub fn summit(&self) -> usize {
        (self.start + self.end) / 2
    }
}

/// Chromatin accessibility profile of a locus or genome region.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ChromatinAccessibilityProfile {
    /// Total region length analyzed in base pairs.
    pub total_region_length: usize,
    /// Total accessible base pairs covered by peaks.
    pub open_chromatin_length: usize,
    /// Fraction of open chromatin (0.0 to 1.0).
    pub accessibility_fraction: f64,
    /// Mean signal value across all peaks.
    pub mean_signal: f64,
    /// Peak records.
    pub peaks: Vec<ChromatinPeak>,
}

/// Evaluates chromatin accessibility profile across a target genomic region.
pub fn evaluate_chromatin_accessibility(
    total_region_length: usize,
    peaks: Vec<ChromatinPeak>,
    q_value_threshold: f64,
) -> ChromatinAccessibilityProfile {
    let significant_peaks: Vec<ChromatinPeak> = peaks
        .into_iter()
        .filter(|p| p.q_value_neg_log10 >= q_value_threshold)
        .collect();

    if significant_peaks.is_empty() || total_region_length == 0 {
        return ChromatinAccessibilityProfile {
            total_region_length,
            open_chromatin_length: 0,
            accessibility_fraction: 0.0,
            mean_signal: 0.0,
            peaks: Vec::new(),
        };
    }

    let mut open_len = 0usize;
    let mut sum_signal = 0.0;

    for peak in &significant_peaks {
        open_len += peak.length();
        sum_signal += peak.signal_value;
    }

    let mean_signal = sum_signal / significant_peaks.len() as f64;
    let fraction = (open_len as f64 / total_region_length as f64).min(1.0);

    ChromatinAccessibilityProfile {
        total_region_length,
        open_chromatin_length: open_len,
        accessibility_fraction: fraction,
        mean_signal,
        peaks: significant_peaks,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chromatin_accessibility_profile() {
        let peaks = vec![
            ChromatinPeak {
                chromosome: "chr5".to_string(),
                start: 100,
                end: 500,
                signal_value: 12.5,
                p_value_neg_log10: 8.0,
                q_value_neg_log10: 5.5,
            },
            ChromatinPeak {
                chromosome: "chr5".to_string(),
                start: 700,
                end: 800,
                signal_value: 1.2,
                p_value_neg_log10: 0.8,
                q_value_neg_log10: 0.3, // Below threshold
            },
        ];

        let profile = evaluate_chromatin_accessibility(1000, peaks, 2.0);
        assert_eq!(profile.peaks.len(), 1);
        assert_eq!(profile.open_chromatin_length, 400);
        assert!((profile.accessibility_fraction - 0.4).abs() < 1e-6);
        assert!((profile.mean_signal - 12.5).abs() < 1e-6);
        assert_eq!(profile.peaks[0].summit(), 300);
    }
}
