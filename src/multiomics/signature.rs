//! Multi-omics molecular signature comparison (Cosine similarity & 2D PCA projection).

use crate::stats::dist::cosine_distance;
use crate::BioError;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A multi-omics feature vector representing a cell type or intervention state.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MolecularSignature {
    /// Sample or cell-type label (e.g. `"Young_Fibroblast"`, `"Rejuvenated_iPSC"`).
    pub label: String,
    /// Standardized multi-omics feature values.
    pub features: Vec<f64>,
}

/// Result of comparing two [`MolecularSignature`] profiles.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SignatureComparison {
    /// Cosine similarity in $[-1.0, 1.0]$ ($> 0.85$ indicates strong phenotypic identity).
    pub cosine_similarity: f64,
    /// Root mean square feature deviation.
    pub rms_deviation: f64,
    /// Whether the two signatures belong to the same phenotypic class ($\text{similarity} > 0.85$).
    pub is_same_class: bool,
}

/// Compares two [`MolecularSignature`] vectors using cosine similarity and RMSD.
pub fn compare_molecular_signatures(
    sig_a: &MolecularSignature,
    sig_b: &MolecularSignature,
) -> crate::Result<SignatureComparison> {
    if sig_a.features.len() != sig_b.features.len() || sig_a.features.is_empty() {
        return Err(BioError::InvalidInput(
            "Molecular signatures must be non-empty and have equal feature length".to_string(),
        ));
    }

    let cos_dist = cosine_distance(&sig_a.features, &sig_b.features)?;
    let cosine_similarity = (1.0 - cos_dist).clamp(-1.0, 1.0);

    let n = sig_a.features.len() as f64;
    let sum_sq: f64 = sig_a
        .features
        .iter()
        .zip(sig_b.features.iter())
        .map(|(a, b)| (a - b) * (a - b))
        .sum();
    let rms_deviation = (sum_sq / n).sqrt();

    Ok(SignatureComparison {
        cosine_similarity,
        rms_deviation,
        is_same_class: cosine_similarity > 0.85,
    })
}

/// Projects a collection of [`MolecularSignature`] profiles onto the top 2 principal components (PC1, PC2)
/// using mean-centering and deflation power iteration.
pub fn project_pca_2d(signatures: &[MolecularSignature]) -> crate::Result<Vec<[f64; 2]>> {
    if signatures.len() < 2 {
        return Err(BioError::InvalidInput(
            "At least 2 signatures are required for PCA projection".to_string(),
        ));
    }
    let d = signatures[0].features.len();
    if d < 2 || signatures.iter().any(|s| s.features.len() != d) {
        return Err(BioError::InvalidInput(
            "All signatures must share at least 2 features".to_string(),
        ));
    }

    let n = signatures.len();
    let mut mean = vec![0.0; d];
    for s in signatures {
        for (j, &v) in s.features.iter().enumerate() {
            mean[j] += v;
        }
    }
    for m in &mut mean {
        *m /= n as f64;
    }

    let centered: Vec<Vec<f64>> = signatures
        .iter()
        .map(|s| s.features.iter().enumerate().map(|(j, &v)| v - mean[j]).collect())
        .collect();

    // Extract top 2 principal axes via power iteration on X^T X
    let extract_axis = |data: &[Vec<f64>]| -> Vec<f64> {
        let mut v: Vec<f64> = (0..d).map(|j| 1.0 / ((j + 1) as f64)).collect();
        for _ in 0..40 {
            let mut next_v = vec![0.0; d];
            for row in data {
                let dot: f64 = row.iter().zip(v.iter()).map(|(a, b)| a * b).sum();
                for j in 0..d {
                    next_v[j] += row[j] * dot;
                }
            }
            let norm = next_v.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-12);
            for x in &mut next_v {
                *x /= norm;
            }
            v = next_v;
        }
        v
    };

    let pc1 = extract_axis(&centered);
    let deflated: Vec<Vec<f64>> = centered
        .iter()
        .map(|row| {
            let proj: f64 = row.iter().zip(pc1.iter()).map(|(a, b)| a * b).sum();
            row.iter().zip(pc1.iter()).map(|(a, p)| a - proj * p).collect()
        })
        .collect();
    let pc2 = extract_axis(&deflated);

    let coords = centered
        .iter()
        .map(|row| {
            let c1 = row.iter().zip(pc1.iter()).map(|(a, b)| a * b).sum();
            let c2 = row.iter().zip(pc2.iter()).map(|(a, b)| a * b).sum();
            [c1, c2]
        })
        .collect();

    Ok(coords)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multiomics_signature_within_class_similarity() {
        // Validation target: cosine similarity between known cell type signatures > 0.85 within class
        let young_fib_1 = MolecularSignature {
            label: "Young_Fibroblast_Rep1".to_string(),
            features: vec![0.92, 0.88, 0.15, 0.12, 0.78, 0.84],
        };
        let young_fib_2 = MolecularSignature {
            label: "Young_Fibroblast_Rep2".to_string(),
            features: vec![0.89, 0.85, 0.18, 0.14, 0.81, 0.82],
        };
        let senescent_fib = MolecularSignature {
            label: "Senescent_Fibroblast".to_string(),
            features: vec![0.12, 0.15, 0.94, 0.91, 0.10, 0.18],
        };

        let within = compare_molecular_signatures(&young_fib_1, &young_fib_2).unwrap();
        assert!(within.cosine_similarity > 0.85);
        assert!(within.is_same_class);

        let across = compare_molecular_signatures(&young_fib_1, &senescent_fib).unwrap();
        assert!(across.cosine_similarity < 0.50);
        assert!(!across.is_same_class);

        let pca = project_pca_2d(&[young_fib_1, young_fib_2, senescent_fib]).unwrap();
        assert_eq!(pca.len(), 3);
    }
}
