//! Vector distance and dissimilarity metrics for expression and phenotypic profiles.

use crate::BioError;

/// Calculates Euclidean distance between two numerical feature vectors:
/// $d(a, b) = \sqrt{\sum (a_i - b_i)^2}$
pub fn euclidean_distance(a: &[f64], b: &[f64]) -> crate::Result<f64> {
    if a.len() != b.len() {
        return Err(BioError::InvalidInput(format!(
            "Vectors must have equal length ({} != {})",
            a.len(),
            b.len()
        )));
    }
    if a.is_empty() {
        return Ok(0.0);
    }
    let sum_sq: f64 = a.iter().zip(b.iter()).map(|(x, y)| (x - y).powi(2)).sum();
    Ok(sum_sq.sqrt())
}

/// Calculates Manhattan (taxicab) distance between two vectors:
/// $d(a, b) = \sum |a_i - b_i|$
pub fn manhattan_distance(a: &[f64], b: &[f64]) -> crate::Result<f64> {
    if a.len() != b.len() {
        return Err(BioError::InvalidInput(format!(
            "Vectors must have equal length ({} != {})",
            a.len(),
            b.len()
        )));
    }
    let sum_abs: f64 = a.iter().zip(b.iter()).map(|(x, y)| (x - y).abs()).sum();
    Ok(sum_abs)
}

/// Calculates Cosine distance ($1.0 - \text{cosine similarity}$):
/// $d_{\cos}(a, b) = 1.0 - \frac{a \cdot b}{\|a\| \|b\|}$
pub fn cosine_distance(a: &[f64], b: &[f64]) -> crate::Result<f64> {
    if a.len() != b.len() {
        return Err(BioError::InvalidInput(format!(
            "Vectors must have equal length ({} != {})",
            a.len(),
            b.len()
        )));
    }
    if a.is_empty() {
        return Ok(0.0);
    }

    let mut dot = 0.0;
    let mut norm_a = 0.0;
    let mut norm_b = 0.0;

    for (&x, &y) in a.iter().zip(b.iter()) {
        dot += x * y;
        norm_a += x * x;
        norm_b += y * y;
    }

    if norm_a <= 0.0 || norm_b <= 0.0 {
        return Ok(1.0); // Orthogonal/zero vector convention
    }

    let sim = dot / (norm_a.sqrt() * norm_b.sqrt());
    let dist = (1.0 - sim).clamp(0.0, 2.0);
    Ok(dist)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_distances() {
        let v1 = [1.0, 2.0, 3.0];
        let v2 = [4.0, 6.0, 8.0];

        let euc = euclidean_distance(&v1, &v2).unwrap();
        // sqrt((3^2) + (4^2) + (5^2)) = sqrt(9 + 16 + 25) = sqrt(50) = 7.0710678
        assert!((euc - 7.0710678).abs() < 1e-5);

        let man = manhattan_distance(&v1, &v2).unwrap();
        assert_eq!(man, 12.0); // 3 + 4 + 5 = 12

        let cos = cosine_distance(&v1, &v1).unwrap();
        assert!(cos < 1e-6); // Identity is 0 distance
    }
}
