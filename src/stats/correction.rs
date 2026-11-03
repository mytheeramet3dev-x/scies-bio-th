//! Multiple hypothesis testing corrections for high-throughput genomics and omics.

/// Bonferroni family-wise error rate (FWER) correction:
/// $p'_i = \min(1.0, p_i \times m)$
pub fn bonferroni_correction(p_values: &[f64]) -> Vec<f64> {
    let m = p_values.len() as f64;
    p_values
        .iter()
        .map(|&p| (p * m).clamp(0.0, 1.0))
        .collect()
}

/// Benjamini-Hochberg (BH) False Discovery Rate (FDR) controlling procedure:
/// Converts raw p-values into adjusted FDR q-values:
/// $q_{(i)} = \min_{j \ge i} \left( \frac{p_{(j)} \cdot m}{j} \right)$
pub fn benjamini_hochberg(p_values: &[f64]) -> Vec<f64> {
    let m = p_values.len();
    if m == 0 {
        return Vec::new();
    }

    // Pair p-value with its original index
    let mut indexed: Vec<(usize, f64)> = p_values.iter().copied().enumerate().collect();
    // Sort ascending by p-value
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut q_values = vec![1.0; m];
    let mut min_so_far = 1.0_f64;

    // Traverse from largest p-value down to smallest (1-based rank j from m down to 1)
    for rank in (1..=m).rev() {
        let i = rank - 1;
        let p = indexed[i].1;
        let q = (p * (m as f64) / (rank as f64)).min(1.0);
        min_so_far = min_so_far.min(q);
        let orig_idx = indexed[i].0;
        q_values[orig_idx] = min_so_far;
    }

    q_values
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bonferroni() {
        let p_vals = [0.01, 0.04, 0.20];
        let corrected = bonferroni_correction(&p_vals);
        // m = 3
        assert!((corrected[0] - 0.03).abs() < 1e-9);
        assert!((corrected[1] - 0.12).abs() < 1e-9);
        assert!((corrected[2] - 0.60).abs() < 1e-9);
    }

    #[test]
    fn test_benjamini_hochberg() {
        let p_vals = [0.001, 0.01, 0.02, 0.05];
        let q_vals = benjamini_hochberg(&p_vals);
        // m = 4
        // p[0] = 0.001 -> 0.001 * 4 / 1 = 0.004
        // p[1] = 0.01  -> 0.01 * 4 / 2 = 0.02
        // p[2] = 0.02  -> 0.02 * 4 / 3 = 0.026666...
        // p[3] = 0.05  -> 0.05 * 4 / 4 = 0.05
        assert!((q_vals[0] - 0.004).abs() < 1e-5);
        assert!((q_vals[1] - 0.02).abs() < 1e-5);
        assert!((q_vals[2] - 0.026666).abs() < 1e-4);
        assert!((q_vals[3] - 0.05).abs() < 1e-5);
    }
}
