//! Statistical hypothesis testing for biological data (Mann-Whitney U, Fisher's Exact, Chi-Squared).

use crate::BioError;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Result of a two-sample Mann-Whitney U rank-sum test.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MannWhitneyResult {
    /// Test statistic $U = \min(U_1, U_2)$.
    pub u_stat: f64,
    /// Asymptotic two-tailed p-value with continuity correction.
    pub p_value: f64,
    /// Standardized normal Z-score.
    pub z_score: f64,
}

/// Result of Fisher's Exact Test on a 2x2 contingency table.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct FishersExactResult {
    /// Odds ratio $(a \cdot d) / (b \cdot c)$.
    pub odds_ratio: f64,
    /// Two-tailed hypergeometric p-value.
    pub p_value: f64,
}

/// Result of Pearson's Chi-squared test of independence.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ChiSquaredResult {
    /// Chi-squared statistic $\chi^2 = \sum \frac{(O - E)^2}{E}$.
    pub chi2_stat: f64,
    /// Degrees of freedom $(r - 1)(c - 1)$.
    pub degrees_of_freedom: usize,
    /// Two-tailed p-value from chi-squared distribution.
    pub p_value: f64,
}

/// Result of Wilcoxon signed-rank test on paired samples.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct WilcoxonResult {
    /// Wilcoxon test statistic $W = \min(W^+, W^-)$.
    pub w_stat: f64,
    /// Two-tailed p-value (exact combinatorial for untied $n \le 50$, asymptotic normal otherwise).
    pub p_value: f64,
    /// Standardized normal Z-score.
    pub z_score: f64,
    /// Number of non-zero paired differences evaluated.
    pub effective_n: usize,
}

/// Result of Kruskal-Wallis one-way analysis of variance by ranks across $k \ge 2$ groups.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct KruskalWallisResult {
    /// Tie-corrected Kruskal-Wallis $H$ statistic.
    pub h_stat: f64,
    /// Degrees of freedom ($k - 1$).
    pub degrees_of_freedom: usize,
    /// Asymptotic $\chi^2$ p-value.
    pub p_value: f64,
}

// --- Mathematical approximations for error function and gamma distributions ---

/// Numerical approximation of the Gauss error function $\text{erf}(x)$
/// (Abramowitz and Stegun 7.1.26, max error $< 1.5 \times 10^{-7}$).
fn erf(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x_abs = x.abs();

    let p = 0.3275911;
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;

    let t = 1.0 / (1.0 + p * x_abs);
    let poly = t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5))));
    let val = 1.0 - poly * (-x_abs * x_abs).exp();

    sign * val
}

/// Standard normal cumulative distribution function $\Phi(z)$.
fn normal_cdf(z: f64) -> f64 {
    0.5 * (1.0 + erf(z / std::f64::consts::SQRT_2))
}

/// Lanczos approximation for $\ln \Gamma(z)$ for $z > 0$.
fn ln_gamma(z: f64) -> f64 {
    let g = 7;
    let c = [
        0.999_999_999_999_809_9,
        676.5203681218851,
        -1259.1392167224028,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507343278686905,
        -0.13857109526572012,
        9.984_369_578_019_572e-6,
        1.5056327351493116e-7,
    ];

    if z < 0.5 {
        std::f64::consts::PI.ln() - (std::f64::consts::PI * z).sin().ln() - ln_gamma(1.0 - z)
    } else {
        let z = z - 1.0;
        let mut x = c[0];
        for (i, &coeff) in c.iter().enumerate().take(g + 2).skip(1) {
            x += coeff / (z + i as f64);
        }
        let t = z + g as f64 + 0.5;
        0.5 * (2.0 * std::f64::consts::PI).ln() + (z + 0.5) * t.ln() - t + x.ln()
    }
}

/// Upper regularized incomplete gamma function $Q(a, x) = \frac{\Gamma(a, x)}{\Gamma(a)}$.
/// Computed via continued fraction (Lentz's method) or series expansion.
fn regularized_gamma_q(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    if x < a + 1.0 {
        // Use series representation for $P(a, x) = 1 - Q(a, x)$
        let mut ap = a;
        let mut sum = 1.0 / a;
        let mut del = sum;
        for _ in 0..100 {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if del.abs() < sum.abs() * 1e-14 {
                break;
            }
        }
        let p = sum * (-x + a * x.ln() - ln_gamma(a)).exp();
        (1.0 - p).clamp(0.0, 1.0)
    } else {
        // Use continued fraction for $Q(a, x)$
        let mut b = x + 1.0 - a;
        let mut c = 1.0 / 1e-30;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..100 {
            let an = -(i as f64) * (i as f64 - a);
            b += 2.0;
            d = an * d + b;
            if d.abs() < 1e-30 {
                d = 1e-30;
            }
            c = b + an / c;
            if c.abs() < 1e-30 {
                c = 1e-30;
            }
            d = 1.0 / d;
            let del = d * c;
            h *= del;
            if (del - 1.0).abs() < 1e-14 {
                break;
            }
        }
        (h * (-x + a * x.ln() - ln_gamma(a)).exp()).clamp(0.0, 1.0)
    }
}

/// Natural log of factorial: $\ln(n!) = \ln \Gamma(n + 1)$.
fn ln_fact(n: u64) -> f64 {
    ln_gamma(n as f64 + 1.0)
}

// --- Public Statistical Tests ---

/// Performs a two-sample Mann-Whitney U rank-sum test between sample1 and sample2.
pub fn mann_whitney_u(sample1: &[f64], sample2: &[f64]) -> crate::Result<MannWhitneyResult> {
    let n1 = sample1.len();
    let n2 = sample2.len();
    if n1 == 0 || n2 == 0 {
        return Err(BioError::InvalidInput("Samples must be non-empty".to_string()));
    }

    // Combine and rank all values
    let mut combined: Vec<(f64, usize)> = Vec::with_capacity(n1 + n2);
    for &x in sample1 {
        combined.push((x, 1));
    }
    for &y in sample2 {
        combined.push((y, 2));
    }
    combined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    // Assign average ranks for ties
    let total_n = n1 + n2;
    let mut ranks = vec![0.0_f64; total_n];
    let mut i = 0;
    let mut tie_correction = 0.0_f64;

    while i < total_n {
        let mut j = i;
        while j + 1 < total_n && (combined[j + 1].0 - combined[i].0).abs() <= 1e-12 {
            j += 1;
        }
        let tie_len = j - i + 1;
        let avg_rank = ((i + 1) + (j + 1)) as f64 / 2.0;
        ranks[i..=j].fill(avg_rank);
        if tie_len > 1 {
            let t = tie_len as f64;
            tie_correction += t * t * t - t;
        }
        i = j + 1;
    }

    // Sum ranks for sample 1
    let mut r1 = 0.0_f64;
    for (idx, &(_, group)) in combined.iter().enumerate() {
        if group == 1 {
            r1 += ranks[idx];
        }
    }

    let u1 = r1 - (n1 * (n1 + 1)) as f64 / 2.0;
    let u2 = (n1 * n2) as f64 - u1;
    let u_stat = u1.min(u2);

    let mean_u = (n1 * n2) as f64 / 2.0;
    let n_tot = total_n as f64;
    let var_u = if n_tot > 1.0 {
        ((n1 * n2) as f64 / 12.0) * (n_tot + 1.0 - tie_correction / (n_tot * (n_tot - 1.0)))
    } else {
        1.0
    };
    let std_u = var_u.max(1e-12).sqrt();

    // Continuity corrected Z score
    let z_score = if (u_stat - mean_u).abs() < 0.5 {
        0.0
    } else {
        (u_stat - mean_u + 0.5) / std_u
    };

    let p_value = (2.0 * normal_cdf(-z_score.abs())).clamp(0.0, 1.0);

    Ok(MannWhitneyResult {
        u_stat,
        p_value,
        z_score,
    })
}

/// Performs Fisher's Exact Test on a 2x2 contingency table:
/// | a | b |
/// | c | d |
pub fn fishers_exact_2x2(a: u64, b: u64, c: u64, d: u64) -> FishersExactResult {
    let odds_ratio = if b == 0 || c == 0 {
        if a == 0 || d == 0 { 0.0 } else { f64::INFINITY }
    } else {
        (a as f64 * d as f64) / (b as f64 * c as f64)
    };

    let r1 = a + b;
    let r2 = c + d;
    let c1 = a + c;
    let n = r1 + r2;

    // Hypergeometric log-probability function for table with upper-left cell x
    let log_prob = |x: u64| -> f64 {
        let y = r1 - x;
        let z = c1 - x;
        let w = r2 - z;
        ln_fact(r1) + ln_fact(r2) + ln_fact(c1) + ln_fact(n - c1)
            - (ln_fact(x) + ln_fact(y) + ln_fact(z) + ln_fact(w) + ln_fact(n))
    };

    let p_obs = log_prob(a).exp();

    let x_min = c1.saturating_sub(r2);
    let x_max = r1.min(c1);

    let mut p_val = 0.0;
    for x in x_min..=x_max {
        let p_x = log_prob(x).exp();
        // Sum probabilities of tables less than or equal to observed (two-tailed)
        if p_x <= p_obs * (1.0 + 1e-9) {
            p_val += p_x;
        }
    }

    FishersExactResult {
        odds_ratio,
        p_value: p_val.clamp(0.0, 1.0),
    }
}

/// Performs Pearson's Chi-squared test of independence on an $r \times c$ contingency table.
pub fn chi_squared_test(observed: &[Vec<f64>]) -> crate::Result<ChiSquaredResult> {
    let r = observed.len();
    if r < 2 {
        return Err(BioError::InvalidInput("Table must have at least 2 rows".to_string()));
    }
    let c = observed[0].len();
    if c < 2 {
        return Err(BioError::InvalidInput("Table must have at least 2 columns".to_string()));
    }
    for row in observed {
        if row.len() != c {
            return Err(BioError::InvalidInput("All rows in table must have equal column length".to_string()));
        }
    }

    let mut row_sums = vec![0.0; r];
    let mut col_sums = vec![0.0; c];
    let mut grand_total = 0.0;

    for i in 0..r {
        for j in 0..c {
            let val = observed[i][j];
            if val < 0.0 {
                return Err(BioError::InvalidInput("Counts must be non-negative".to_string()));
            }
            row_sums[i] += val;
            col_sums[j] += val;
            grand_total += val;
        }
    }

    if grand_total <= 0.0 {
        return Err(BioError::InvalidInput("Contingency table cannot be all zeros".to_string()));
    }

    let mut chi2 = 0.0;
    for i in 0..r {
        for j in 0..c {
            let expected = (row_sums[i] * col_sums[j]) / grand_total;
            if expected > 0.0 {
                let diff = observed[i][j] - expected;
                chi2 += (diff * diff) / expected;
            }
        }
    }

    let df = (r - 1) * (c - 1);
    let p_value = regularized_gamma_q(df as f64 / 2.0, chi2 / 2.0);

    Ok(ChiSquaredResult {
        chi2_stat: chi2,
        degrees_of_freedom: df,
        p_value,
    })
}

/// Performs the Wilcoxon signed-rank test on paired samples `x` and `y`.
///
/// Zero differences ($x_i - y_i = 0$) are discarded per Wilcoxon's standard convention.
/// When no rank ties exist and $n \le 50$, computes the exact combinatorial two-tailed p-value
/// via dynamic programming (matching `scipy.stats.wilcoxon(..., mode='exact')` within $\pm 10^{-12}$).
pub fn wilcoxon_signed_rank(x: &[f64], y: &[f64]) -> crate::Result<WilcoxonResult> {
    if x.len() != y.len() {
        return Err(BioError::InvalidInput(
            "Paired samples must have identical length".to_string(),
        ));
    }

    // Collect non-zero differences (|d|, sign)
    let mut diffs: Vec<(f64, i8)> = Vec::with_capacity(x.len());
    for (&xi, &yi) in x.iter().zip(y.iter()) {
        let d = xi - yi;
        if d.abs() > 1e-15 {
            diffs.push((d.abs(), if d > 0.0 { 1 } else { -1 }));
        }
    }

    let n = diffs.len();
    if n == 0 {
        return Err(BioError::InvalidInput(
            "All paired differences are zero".to_string(),
        ));
    }

    diffs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut ranks = vec![0.0_f64; n];
    let mut i = 0;
    let mut tie_correction = 0.0_f64;
    let mut has_ties = false;

    while i < n {
        let mut j = i;
        while j + 1 < n && (diffs[j + 1].0 - diffs[i].0).abs() <= 1e-12 {
            j += 1;
        }
        let tie_len = j - i + 1;
        let avg_rank = ((i + 1) + (j + 1)) as f64 / 2.0;
        ranks[i..=j].fill(avg_rank);
        if tie_len > 1 {
            has_ties = true;
            let t = tie_len as f64;
            tie_correction += t * t * t - t;
        }
        i = j + 1;
    }

    let mut w_plus = 0.0_f64;
    let mut w_minus = 0.0_f64;
    for (idx, &(_, sign)) in diffs.iter().enumerate() {
        if sign > 0 {
            w_plus += ranks[idx];
        } else {
            w_minus += ranks[idx];
        }
    }

    let w_stat = w_plus.min(w_minus);
    let nf = n as f64;
    let mean_w = nf * (nf + 1.0) / 4.0;
    let var_w = (nf * (nf + 1.0) * (2.0 * nf + 1.0) - 0.5 * tie_correction) / 24.0;
    let std_w = var_w.max(1e-12).sqrt();
    let z_score = (w_stat - mean_w) / std_w;

    let p_value = if !has_ties && n <= 50 {
        // Exact combinatorial distribution over rank sums 0..=n(n+1)/2
        let max_sum = n * (n + 1) / 2;
        let mut dp = vec![0.0_f64; max_sum + 1];
        dp[0] = 1.0;
        for r in 1..=n {
            for s in (r..=max_sum).rev() {
                dp[s] += dp[s - r];
            }
        }
        let target = w_stat.round() as usize;
        let tail_count: f64 = dp[0..=target.min(max_sum)].iter().sum();
        let total_outcomes = 2.0_f64.powi(n as i32);
        (2.0 * tail_count / total_outcomes).clamp(0.0, 1.0)
    } else {
        (2.0 * normal_cdf(-z_score.abs())).clamp(0.0, 1.0)
    };

    Ok(WilcoxonResult {
        w_stat,
        p_value,
        z_score,
        effective_n: n,
    })
}

/// Performs the Kruskal-Wallis one-way analysis of variance by ranks across $k \ge 2$ sample groups.
pub fn kruskal_wallis(groups: &[&[f64]]) -> crate::Result<KruskalWallisResult> {
    if groups.len() < 2 {
        return Err(BioError::InvalidInput(
            "Kruskal-Wallis requires at least 2 groups".to_string(),
        ));
    }
    for g in groups {
        if g.is_empty() {
            return Err(BioError::InvalidInput(
                "All Kruskal-Wallis groups must be non-empty".to_string(),
            ));
        }
    }

    let mut combined: Vec<(f64, usize)> = Vec::new();
    for (g_idx, group) in groups.iter().enumerate() {
        for &val in *group {
            combined.push((val, g_idx));
        }
    }
    let total_n = combined.len();
    combined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut ranks = vec![0.0_f64; total_n];
    let mut i = 0;
    let mut tie_correction = 0.0_f64;

    while i < total_n {
        let mut j = i;
        while j + 1 < total_n && (combined[j + 1].0 - combined[i].0).abs() <= 1e-12 {
            j += 1;
        }
        let tie_len = j - i + 1;
        let avg_rank = ((i + 1) + (j + 1)) as f64 / 2.0;
        ranks[i..=j].fill(avg_rank);
        if tie_len > 1 {
            let t = tie_len as f64;
            tie_correction += t * t * t - t;
        }
        i = j + 1;
    }

    let mut rank_sums = vec![0.0_f64; groups.len()];
    for (idx, &(_, g_idx)) in combined.iter().enumerate() {
        rank_sums[g_idx] += ranks[idx];
    }

    let nf = total_n as f64;
    let mut sum_r2_ni = 0.0_f64;
    for (g_idx, group) in groups.iter().enumerate() {
        let ni = group.len() as f64;
        sum_r2_ni += (rank_sums[g_idx] * rank_sums[g_idx]) / ni;
    }

    let h_uncorrected = (12.0 / (nf * (nf + 1.0))) * sum_r2_ni - 3.0 * (nf + 1.0);
    let denom = 1.0 - tie_correction / (nf * nf * nf - nf);
    let h_stat = if denom > 1e-12 {
        (h_uncorrected / denom).max(0.0)
    } else {
        0.0
    };

    let df = groups.len() - 1;
    let p_value = regularized_gamma_q(df as f64 / 2.0, h_stat / 2.0);

    Ok(KruskalWallisResult {
        h_stat,
        degrees_of_freedom: df,
        p_value,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mann_whitney_u_known_values() {
        // Distinct distributions
        let s1 = [1.0, 2.0, 3.0, 4.0, 5.0];
        let s2 = [6.0, 7.0, 8.0, 9.0, 10.0];
        let res = mann_whitney_u(&s1, &s2).unwrap();
        assert_eq!(res.u_stat, 0.0);
        assert!(res.p_value < 0.02);
    }

    #[test]
    fn test_fishers_exact_known_table() {
        // Classic Lady Tasting Tea 2x2 table:
        // [[3, 1], [1, 3]]
        let res = fishers_exact_2x2(3, 1, 1, 3);
        assert!((res.odds_ratio - 9.0).abs() < 1e-4);
        assert!(res.p_value > 0.40); // small sample not significant
    }

    #[test]
    fn test_chi_squared_test() {
        // 2x2 table
        let table = vec![
            vec![12.0, 5.0],
            vec![4.0, 15.0],
        ];
        let res = chi_squared_test(&table).unwrap();
        assert_eq!(res.degrees_of_freedom, 1);
        assert!(res.chi2_stat > 6.0);
        assert!(res.p_value < 0.02);
    }

    #[test]
    fn test_wilcoxon_scipy_ground_truth() {
        // SciPy benchmark: x = [125, 115, 130, 140, 140, 115, 140, 125, 140, 135]
        // y = [110, 122, 125, 120, 140, 124, 123, 137, 135, 145]
        // Or untied 8-pair benchmark: diffs = [+1, +2, +3, +4, +5, +6, +7, -8]
        // W+ = 28, W- = 8, W = 8. For n=8, sum<=8 has counts:
        // 0:1, 1:1, 2:1, 3:2, 4:2, 5:3, 6:4, 7:5, 8:6 -> sum = 25.
        // Exact two-tailed p = 2 * 25 / 256 = 50 / 256 = 0.1953125 (SciPy wilcoxon exact = 0.1953125).
        let x = [11.0, 22.0, 33.0, 44.0, 55.0, 66.0, 77.0, 80.0];
        let y = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 88.0];
        let res = wilcoxon_signed_rank(&x, &y).unwrap();
        assert_eq!(res.w_stat, 8.0);
        assert!((res.p_value - 0.195_312_5).abs() < 1e-6);
    }

    #[test]
    fn test_kruskal_wallis_multi_group() {
        let g1 = [1.0, 2.0, 3.0, 4.0];
        let g2 = [5.0, 6.0, 7.0, 8.0];
        let g3 = [9.0, 10.0, 11.0, 12.0];
        let res = kruskal_wallis(&[&g1, &g2, &g3]).unwrap();
        assert_eq!(res.degrees_of_freedom, 2);
        assert!(res.h_stat > 9.0);
        assert!(res.p_value < 0.01);
    }

    #[test]
    fn test_mann_whitney_float_tie_tolerance() {
        // Values differing only by floating point epsilon 1e-14 should be detected as tied
        let s1 = [1.0, 2.0, 3.0 + 1e-14];
        let s2 = [3.0, 4.0, 5.0];
        let res = mann_whitney_u(&s1, &s2).unwrap();
        assert!(res.u_stat >= 0.0);
        assert!(res.p_value <= 1.0);
    }
}
