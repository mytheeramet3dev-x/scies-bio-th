//! Metagenomics and ecological biodiversity metrics.
//!
//! Provides functions for computing Alpha Diversity indices:
//! Shannon-Wiener entropy index, Simpson diversity, Pielou's evenness, and Chao1 richness estimator.

/// Computes the Shannon-Wiener diversity index ($H$).
///
/// Formula:
/// $H = -\sum_{i=1}^S p_i \ln(p_i)$
/// where $p_i = \text{abundance}_i / \sum \text{abundance}_j$.
pub fn shannon_index(abundances: &[f64]) -> f64 {
    let total: f64 = abundances.iter().filter(|&&x| x > 0.0).sum();
    if total <= 0.0 {
        return 0.0;
    }

    let mut h = 0.0;
    for &val in abundances {
        if val > 0.0 {
            let p = val / total;
            h -= p * p.ln();
        }
    }
    h
}

/// Computes the Simpson concentration index ($\lambda = \sum p_i^2$).
///
/// Ranges from $0.0$ to $1.0$. Represents the probability that two randomly
/// sampled individuals belong to the same species.
pub fn simpson_index(abundances: &[f64]) -> f64 {
    let total: f64 = abundances.iter().filter(|&&x| x > 0.0).sum();
    if total <= 0.0 {
        return 0.0;
    }

    abundances
        .iter()
        .filter(|&&x| x > 0.0)
        .map(|&val| {
            let p = val / total;
            p * p
        })
        .sum()
}

/// Computes the Gini-Simpson diversity index ($1 - \lambda$).
///
/// Ranges from $0.0$ (no diversity) to $1.0$ (infinite diversity).
pub fn gini_simpson_index(abundances: &[f64]) -> f64 {
    1.0 - simpson_index(abundances)
}

/// Computes Pielou's species evenness ($J = H / \ln(S)$).
///
/// Ranges from $0.0$ (highly skewed) to $1.0$ (perfectly uniform abundance).
pub fn pielou_evenness(abundances: &[f64]) -> f64 {
    let s = abundances.iter().filter(|&&x| x > 0.0).count();
    if s <= 1 {
        return 1.0;
    }
    let h = shannon_index(abundances);
    h / (s as f64).ln()
}

/// Computes the Chao1 species richness estimator.
///
/// Uses the bias-corrected formula:
/// $S_{\text{chao1}} = S_{\text{obs}} + \frac{f_1(f_1 - 1)}{2(f_2 + 1)}$
/// where $S_{\text{obs}}$ is the number of observed species,
/// $f_1$ is the number of singletons (count = 1),
/// and $f_2$ is the number of doubletons (count = 2).
pub fn chao1(counts: &[u64]) -> f64 {
    let s_obs = counts.iter().filter(|&&c| c > 0).count() as f64;
    let f1 = counts.iter().filter(|&&c| c == 1).count() as f64;
    let f2 = counts.iter().filter(|&&c| c == 2).count() as f64;

    if f2 > 0.0 {
        s_obs + (f1 * (f1 - 1.0)) / (2.0 * (f2 + 1.0))
    } else {
        s_obs + (f1 * (f1 - 1.0)) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shannon_index_uniform() {
        // Equal abundance of 4 species: H = ln(4) = 1.38629...
        let abundances = [25.0, 25.0, 25.0, 25.0];
        let h = shannon_index(&abundances);
        let expected = 4.0_f64.ln();
        assert!((h - expected).abs() < 1e-6);

        // Evenness should be exactly 1.0 for uniform distribution
        assert!((pielou_evenness(&abundances) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_simpson_and_gini_simpson() {
        let abundances = [50.0, 50.0];
        let d = simpson_index(&abundances);
        // (0.5)^2 + (0.5)^2 = 0.5
        assert!((d - 0.5).abs() < 1e-6);
        assert!((gini_simpson_index(&abundances) - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_chao1_richness() {
        // 5 observed species: three singletons (count=1), one doubleton (count=2), one high (count=10)
        let counts = [1, 1, 1, 2, 10];
        let c1 = chao1(&counts);
        // S_obs = 5, f1 = 3, f2 = 1
        // Chao1 = 5 + (3 * 2) / (2 * (1 + 1)) = 5 + 6 / 4 = 6.5
        assert!((c1 - 6.5).abs() < 1e-6);
    }
}
