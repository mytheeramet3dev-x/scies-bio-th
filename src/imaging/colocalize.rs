//! Dual-channel fluorescence co-localization analysis (Pearson's $r$ and Manders' $M_1, M_2$).

use crate::BioError;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Co-localization metrics between two microscopy fluorescence channels.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ColocalizationMetrics {
    /// Pearson Correlation Coefficient (PCC) $r \in [-1.0, 1.0]$.
    pub pearson_r: f64,
    /// Manders' Overlap Coefficient $R \in [0.0, 1.0]$.
    pub manders_overlap: f64,
    /// Manders' split colocalization coefficient $M_1$ (fraction of channel 1 overlapping channel 2).
    pub manders_m1: f64,
    /// Manders' split colocalization coefficient $M_2$ (fraction of channel 2 overlapping channel 1).
    pub manders_m2: f64,
}

/// Calculates Pearson's correlation coefficient and Manders' split coefficients ($M_1, M_2$)
/// between two equal-sized fluorescence intensity arrays.
pub fn calculate_colocalization(
    channel_1: &[f64],
    channel_2: &[f64],
    threshold_1: f64,
    threshold_2: f64,
) -> crate::Result<ColocalizationMetrics> {
    if channel_1.len() != channel_2.len() || channel_1.is_empty() {
        return Err(BioError::InvalidInput(
            "Fluorescence channels must be non-empty and have identical pixel count".to_string(),
        ));
    }

    let n = channel_1.len() as f64;
    let mean_1 = channel_1.iter().sum::<f64>() / n;
    let mean_2 = channel_2.iter().sum::<f64>() / n;

    let mut cov = 0.0;
    let mut var_1 = 0.0;
    let mut var_2 = 0.0;

    let mut sum_prod = 0.0;
    let mut sum_sq_1 = 0.0;
    let mut sum_sq_2 = 0.0;

    let mut total_c1 = 0.0;
    let mut total_c2 = 0.0;
    let mut coloc_c1 = 0.0;
    let mut coloc_c2 = 0.0;

    for (&c1, &c2) in channel_1.iter().zip(channel_2.iter()) {
        let d1 = c1 - mean_1;
        let d2 = c2 - mean_2;
        cov += d1 * d2;
        var_1 += d1 * d1;
        var_2 += d2 * d2;

        let p1 = c1.max(0.0);
        let p2 = c2.max(0.0);
        sum_prod += p1 * p2;
        sum_sq_1 += p1 * p1;
        sum_sq_2 += p2 * p2;

        total_c1 += p1;
        total_c2 += p2;
        if p2 > threshold_2 {
            coloc_c1 += p1;
        }
        if p1 > threshold_1 {
            coloc_c2 += p2;
        }
    }

    let denom_pcc = (var_1 * var_2).sqrt();
    let pearson_r = if denom_pcc > 1e-12 {
        (cov / denom_pcc).clamp(-1.0, 1.0)
    } else {
        0.0
    };

    let denom_moc = (sum_sq_1 * sum_sq_2).sqrt();
    let manders_overlap = if denom_moc > 1e-12 {
        (sum_prod / denom_moc).clamp(0.0, 1.0)
    } else {
        0.0
    };

    let manders_m1 = if total_c1 > 1e-12 {
        (coloc_c1 / total_c1).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let manders_m2 = if total_c2 > 1e-12 {
        (coloc_c2 / total_c2).clamp(0.0, 1.0)
    } else {
        0.0
    };

    Ok(ColocalizationMetrics {
        pearson_r,
        manders_overlap,
        manders_m1,
        manders_m2,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_colocalization_pearson_above_095_validation() {
        // Validation target: Pearson coefficient on synthetic co-localized channels >= 0.95
        let ch1 = [10.0, 15.0, 180.0, 240.0, 195.0, 12.0, 210.0, 14.0];
        let ch2 = [12.0, 18.0, 175.0, 250.0, 190.0, 10.0, 205.0, 16.0];

        let metrics = calculate_colocalization(&ch1, &ch2, 50.0, 50.0).unwrap();
        assert!(metrics.pearson_r >= 0.95);
        assert!(metrics.manders_m1 > 0.90);
        assert!(metrics.manders_m2 > 0.90);
    }
}
