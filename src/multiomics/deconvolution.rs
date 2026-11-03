//! Reference-based bulk tissue cell-type deconvolution via multiplicative-update
//! Non-Negative Least Squares / NMF-lite (`multiomics::deconvolution`).

use crate::BioError;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Estimated cell-type composition of a bulk multi-omics sample.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct DeconvolutionResult {
    /// Normalized cell-type fractions summing to $1.0$.
    pub proportions: Vec<f64>,
    /// Residual L2 reconstruction error $\|v - W \cdot h\|_2$.
    pub reconstruction_residual: f64,
}

/// Estimates non-negative cell-type proportions $h$ ($k$ cell types, $\sum h_c = 1$)
/// from a bulk expression vector `bulk_profile` ($m$ genes) and a reference signature matrix
/// `reference_signatures` ($k$ cell types $\times m$ genes) using Lee-Seung multiplicative updates.
pub fn deconvolve_cell_types_nmf(
    bulk_profile: &[f64],
    reference_signatures: &[Vec<f64>],
    iterations: usize,
) -> crate::Result<DeconvolutionResult> {
    let m = bulk_profile.len();
    let k = reference_signatures.len();
    if m == 0 || k == 0 {
        return Err(BioError::InvalidInput(
            "Bulk profile and reference signature list must be non-empty".to_string(),
        ));
    }
    if reference_signatures.iter().any(|r| r.len() != m) {
        return Err(BioError::InvalidInput(
            "All reference cell-type signatures must match bulk profile gene count".to_string(),
        ));
    }

    let mut h = vec![1.0 / (k as f64); k];
    let eps = 1e-12;

    // Multiplicative update: h_c <- h_c * (W^T v)_c / (W^T W h)_c
    for _ in 0..iterations.max(25) {
        // Compute reconstructed bulk hat_v = W * h (length m)
        let mut hat_v = vec![0.0; m];
        for (c, ref_sig) in reference_signatures.iter().enumerate() {
            let hc = h[c];
            for (g, &w_gc) in ref_sig.iter().enumerate() {
                hat_v[g] += w_gc.max(0.0) * hc;
            }
        }

        for (c, ref_sig) in reference_signatures.iter().enumerate() {
            let mut numer = 0.0;
            let mut denom = 0.0;
            for g in 0..m {
                let w_gc = ref_sig[g].max(0.0);
                numer += w_gc * bulk_profile[g].max(0.0);
                denom += w_gc * hat_v[g];
            }
            h[c] *= (numer + eps) / (denom + eps);
        }
    }

    // Normalize proportions to sum to 1.0
    let sum_h: f64 = h.iter().sum::<f64>().max(eps);
    for hc in &mut h {
        *hc /= sum_h;
    }

    // Compute residual
    let mut sq_err = 0.0;
    for (g, &obs) in bulk_profile.iter().enumerate() {
        let pred: f64 = reference_signatures
            .iter()
            .enumerate()
            .map(|(c, sig)| sig[g] * h[c])
            .sum();
        let diff = obs - pred;
        sq_err += diff * diff;
    }

    Ok(DeconvolutionResult {
        proportions: h,
        reconstruction_residual: sq_err.sqrt(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_deconvolve_cell_types_nmf_recovery() {
        // Reference cell types: C1 = [10, 1, 0, 8], C2 = [0, 9, 12, 1]
        let ref_c1 = vec![10.0, 1.0, 0.0, 8.0];
        let ref_c2 = vec![0.0, 9.0, 12.0, 1.0];
        // Synthetic bulk mixture: 70% C1 + 30% C2
        let bulk: Vec<f64> = (0..4).map(|g| 0.70 * ref_c1[g] + 0.30 * ref_c2[g]).collect();

        let res = deconvolve_cell_types_nmf(&bulk, &[ref_c1, ref_c2], 100).unwrap();
        assert!((res.proportions[0] - 0.70).abs() < 0.02);
        assert!((res.proportions[1] - 0.30).abs() < 0.02);
    }
}
