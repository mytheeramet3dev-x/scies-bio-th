//! Population genetics, allele frequencies, and Hardy-Weinberg equilibrium modeling.
//!
//! Provides calculations for allele frequencies, observed and expected heterozygosity,
//! Wright's inbreeding coefficient ($F_{\text{IS}}$), and Chi-square Hardy-Weinberg Equilibrium (HWE) testing.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Diploid genotype counts in a sample population for a biallelic locus (Alleles A and B).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct GenotypeCounts {
    /// Number of Homozygous Major (AA) individuals.
    pub aa: u64,
    /// Number of Heterozygous (AB) individuals.
    pub ab: u64,
    /// Number of Homozygous Minor (BB) individuals.
    pub bb: u64,
}

/// Result of a Chi-square Hardy-Weinberg Equilibrium test.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HweTestResult {
    /// Chi-square statistic ($\chi^2$).
    pub chi_square: f64,
    /// Approximate p-value (degrees of freedom = 1).
    pub p_value: f64,
    /// Whether the locus is in Hardy-Weinberg equilibrium at the tested alpha significance level.
    pub in_equilibrium: bool,
    /// Observed heterozygosity ($H_o$).
    pub observed_heterozygosity: f64,
    /// Expected heterozygosity ($H_e = 2pq$).
    pub expected_heterozygosity: f64,
    /// Wright's inbreeding coefficient ($F = 1 - H_o / H_e$).
    pub inbreeding_coefficient: f64,
}

impl GenotypeCounts {
    /// Creates a new `GenotypeCounts` struct.
    pub fn new(aa: u64, ab: u64, bb: u64) -> Self {
        Self { aa, ab, bb }
    }

    /// Total number of individuals in the sample.
    pub fn total_individuals(&self) -> u64 {
        self.aa + self.ab + self.bb
    }

    /// Total number of alleles in the sample ($2 \times N$).
    pub fn total_alleles(&self) -> u64 {
        self.total_individuals() * 2
    }

    /// Calculates allele frequencies $(p, q)$ for Alleles A and B.
    pub fn allele_frequencies(&self) -> (f64, f64) {
        let total_alleles = self.total_alleles();
        if total_alleles == 0 {
            return (0.0, 0.0);
        }

        let p = (2 * self.aa + self.ab) as f64 / total_alleles as f64;
        let q = 1.0 - p;
        (p, q)
    }

    /// Observed heterozygosity ($H_o = N_{AB} / N$).
    pub fn observed_heterozygosity(&self) -> f64 {
        let n = self.total_individuals();
        if n == 0 {
            return 0.0;
        }
        self.ab as f64 / n as f64
    }

    /// Expected heterozygosity under Hardy-Weinberg assumptions ($H_e = 2pq$).
    pub fn expected_heterozygosity(&self) -> f64 {
        let (p, q) = self.allele_frequencies();
        2.0 * p * q
    }

    /// Wright's Inbreeding Coefficient ($F = 1 - H_o / H_e$).
    pub fn inbreeding_coefficient(&self) -> f64 {
        let he = self.expected_heterozygosity();
        if he <= 0.0 {
            return 0.0;
        }
        1.0 - (self.observed_heterozygosity() / he)
    }

    /// Expected genotype counts $(E_{AA}, E_{AB}, E_{BB})$ under HWE.
    pub fn expected_counts(&self) -> (f64, f64, f64) {
        let n = self.total_individuals() as f64;
        let (p, q) = self.allele_frequencies();
        (p * p * n, 2.0 * p * q * n, q * q * n)
    }

    /// Performs a Chi-square goodness-of-fit test for Hardy-Weinberg Equilibrium.
    ///
    /// Tests the null hypothesis that genotype frequencies conform to $p^2 : 2pq : q^2$ ($df = 1$).
    pub fn test_hardy_weinberg(&self, alpha: f64) -> HweTestResult {
        let n = self.total_individuals();
        if n == 0 {
            return HweTestResult {
                chi_square: 0.0,
                p_value: 1.0,
                in_equilibrium: true,
                observed_heterozygosity: 0.0,
                expected_heterozygosity: 0.0,
                inbreeding_coefficient: 0.0,
            };
        }

        let (exp_aa, exp_ab, exp_bb) = self.expected_counts();
        let obs_aa = self.aa as f64;
        let obs_ab = self.ab as f64;
        let obs_bb = self.bb as f64;

        let chi_square = ((obs_aa - exp_aa).powi(2) / exp_aa.max(1e-9))
            + ((obs_ab - exp_ab).powi(2) / exp_ab.max(1e-9))
            + ((obs_bb - exp_bb).powi(2) / exp_bb.max(1e-9));

        // For df=1, p-value = erfc(sqrt(chi_square / 2))
        let p_value = chi_square_df1_pvalue(chi_square);
        let in_equilibrium = p_value >= alpha;

        HweTestResult {
            chi_square,
            p_value,
            in_equilibrium,
            observed_heterozygosity: self.observed_heterozygosity(),
            expected_heterozygosity: self.expected_heterozygosity(),
            inbreeding_coefficient: self.inbreeding_coefficient(),
        }
    }
}

/// Approximate p-value for Chi-Square distribution with df = 1.
///
/// Uses complementary error function: $P(X \ge \chi^2) = \text{erfc}(\sqrt{\chi^2 / 2})$.
fn chi_square_df1_pvalue(chi2: f64) -> f64 {
    if chi2 <= 0.0 {
        return 1.0;
    }
    let z = (chi2 / 2.0).sqrt();
    erfc_approx(z).clamp(0.0, 1.0)
}

/// Numerical approximation of complementary error function erfc(x) with high precision (Abramowitz & Stegun).
fn erfc_approx(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc_approx(-x);
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736
                + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    poly * (-x * x).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hwe_perfect_equilibrium() {
        // p = 0.6, q = 0.4, N = 100
        // AA: 36, AB: 48, BB: 16
        let counts = GenotypeCounts::new(36, 48, 16);
        let (p, q) = counts.allele_frequencies();
        assert!((p - 0.6).abs() < 1e-6);
        assert!((q - 0.4).abs() < 1e-6);

        let hwe = counts.test_hardy_weinberg(0.05);
        assert!(hwe.chi_square < 1e-6);
        assert!((hwe.p_value - 1.0).abs() < 1e-4);
        assert!(hwe.in_equilibrium);
        assert!((hwe.inbreeding_coefficient).abs() < 1e-4);
    }

    #[test]
    fn test_hwe_severe_inbreeding_departure() {
        // Complete deficiency of heterozygotes: AA: 50, AB: 0, BB: 50
        let counts = GenotypeCounts::new(50, 0, 50);
        let (p, q) = counts.allele_frequencies();
        assert_eq!(p, 0.5);
        assert_eq!(q, 0.5);

        let hwe = counts.test_hardy_weinberg(0.05);
        assert!(hwe.chi_square > 50.0);
        assert!(hwe.p_value < 1e-5);
        assert!(!hwe.in_equilibrium);
        // F should be 1.0 (complete inbreeding)
        assert!((hwe.inbreeding_coefficient - 1.0).abs() < 1e-4);
    }
}
