//! Root Mean Square Deviation (RMSD) and TM-score with optimal Kabsch / Quaternion superposition.

use crate::BioError;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Result of 3D macromolecular optimal superposition and alignment.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SuperpositionResult {
    /// Coordinate Root Mean Square Deviation in Ångströms.
    pub rmsd: f64,
    /// Template Modeling score ($0.0$ to $1.0$, $> 0.5$ indicates same global fold).
    pub tm_score: f64,
    /// Centroid of mobile structure P.
    pub centroid_p: [f64; 3],
    /// Centroid of target structure Q.
    pub centroid_q: [f64; 3],
    /// $3 \times 3$ optimal rotation matrix.
    pub rotation_matrix: [[f64; 3]; 3],
}

/// Solves the largest eigenvalue and corresponding eigenvector of a $4 \times 4$ real symmetric matrix
/// using the classical cyclic Jacobi rotation algorithm.
#[allow(clippy::needless_range_loop)]
fn largest_eigen_symmetric_4x4(m: &[[f64; 4]; 4]) -> [f64; 4] {
    let mut a = *m;
    let mut v = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];

    for _ in 0..50 {
        let mut max_off = 0.0;
        let mut p = 0;
        let mut q = 1;

        for i in 0..4 {
            for j in (i + 1)..4 {
                let off = a[i][j].abs();
                if off > max_off {
                    max_off = off;
                    p = i;
                    q = j;
                }
            }
        }

        if max_off < 1e-12 {
            break;
        }

        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];

        let theta = (aqq - app) / (2.0 * apq);
        let t = if theta >= 0.0 {
            1.0 / (theta + (1.0 + theta * theta).sqrt())
        } else {
            -1.0 / (-theta + (1.0 + theta * theta).sqrt())
        };

        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;
        let tau = s / (1.0 + c);

        a[p][p] -= t * apq;
        a[q][q] += t * apq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;

        for k in 0..4 {
            if k != p && k != q {
                let akp = a[k][p];
                let akq = a[k][q];
                a[k][p] = akp - s * (akq + tau * akp);
                a[p][k] = a[k][p];
                a[k][q] = akq + s * (akp - tau * akq);
                a[q][k] = a[k][q];
            }
        }

        for k in 0..4 {
            let vkp = v[k][p];
            let vkq = v[k][q];
            v[k][p] = vkp - s * (vkq + tau * vkp);
            v[k][q] = vkq + s * (vkp - tau * vkq);
        }
    }

    let mut max_eig = a[0][0];
    let mut best_col = 0;
    for i in 1..4 {
        if a[i][i] > max_eig {
            max_eig = a[i][i];
            best_col = i;
        }
    }

    let mut q_vec = [v[0][best_col], v[1][best_col], v[2][best_col], v[3][best_col]];
    if q_vec[0] < 0.0 {
        q_vec[0] = -q_vec[0];
        q_vec[1] = -q_vec[1];
        q_vec[2] = -q_vec[2];
        q_vec[3] = -q_vec[3];
    }
    q_vec
}

/// Converts a unit quaternion $q = [w, x, y, z]$ into a $3 \times 3$ orthogonal rotation matrix.
fn quaternion_to_rotation_matrix(q: &[f64; 4]) -> [[f64; 3]; 3] {
    let (w, x, y, z) = (q[0], q[1], q[2], q[3]);

    [
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y - z * w),
            2.0 * (x * z + y * w),
        ],
        [
            2.0 * (x * y + z * w),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z - x * w),
        ],
        [
            2.0 * (x * z - y * w),
            2.0 * (y * z + x * w),
            1.0 - 2.0 * (x * x + y * y),
        ],
    ]
}

/// Applies a $3 \times 3$ rotation matrix to a 3D vector.
pub fn apply_rotation(r: &[[f64; 3]; 3], v: &[f64; 3]) -> [f64; 3] {
    [
        r[0][0] * v[0] + r[0][1] * v[1] + r[0][2] * v[2],
        r[1][0] * v[0] + r[1][1] * v[1] + r[1][2] * v[2],
        r[2][0] * v[0] + r[2][1] * v[1] + r[2][2] * v[2],
    ]
}

/// Computes the optimal superposition RMSD and TM-score between coordinate sets P and Q.
pub fn calculate_rmsd_and_tm_score(
    coords_p: &[[f64; 3]],
    coords_q: &[[f64; 3]],
    target_length: usize,
) -> crate::Result<SuperpositionResult> {
    let n = coords_p.len();
    if n != coords_q.len() {
        return Err(BioError::InvalidInput(format!(
            "Coordinate sets must have identical size ({} != {})",
            n,
            coords_q.len()
        )));
    }
    if n < 3 {
        return Err(BioError::InvalidInput("At least 3 coordinate pairs required for 3D superposition".to_string()));
    }

    // 1. Calculate centroids
    let mut cp = [0.0; 3];
    let mut cq = [0.0; 3];
    for i in 0..n {
        cp[0] += coords_p[i][0];
        cp[1] += coords_p[i][1];
        cp[2] += coords_p[i][2];

        cq[0] += coords_q[i][0];
        cq[1] += coords_q[i][1];
        cq[2] += coords_q[i][2];
    }
    let nf = n as f64;
    cp[0] /= nf;
    cp[1] /= nf;
    cp[2] /= nf;

    cq[0] /= nf;
    cq[1] /= nf;
    cq[2] /= nf;

    // 2. Compute cross-covariance matrix H = sum( x_i * y_i^T )
    let mut h = [[0.0; 3]; 3];
    for i in 0..n {
        let x = [coords_p[i][0] - cp[0], coords_p[i][1] - cp[1], coords_p[i][2] - cp[2]];
        let y = [coords_q[i][0] - cq[0], coords_q[i][1] - cq[1], coords_q[i][2] - cq[2]];

        for r in 0..3 {
            for c in 0..3 {
                h[r][c] += x[r] * y[c];
            }
        }
    }

    // 3. Build Horn's 4x4 matrix for optimal quaternion extraction
    let m = [
        [
            h[0][0] + h[1][1] + h[2][2],
            h[1][2] - h[2][1],
            h[2][0] - h[0][2],
            h[0][1] - h[1][0],
        ],
        [
            h[1][2] - h[2][1],
            h[0][0] - h[1][1] - h[2][2],
            h[0][1] + h[1][0],
            h[2][0] + h[0][2],
        ],
        [
            h[2][0] - h[0][2],
            h[0][1] + h[1][0],
            -h[0][0] + h[1][1] - h[2][2],
            h[1][2] + h[2][1],
        ],
        [
            h[0][1] - h[1][0],
            h[2][0] + h[0][2],
            h[1][2] + h[2][1],
            -h[0][0] - h[1][1] + h[2][2],
        ],
    ];

    let q = largest_eigen_symmetric_4x4(&m);
    let rot = quaternion_to_rotation_matrix(&q);

    // 4. Calculate RMSD and distances for TM-score
    let mut sum_sq_diff = 0.0;
    let mut tm_sum = 0.0;

    let l_target = target_length.max(n);
    let d0 = if l_target > 15 {
        (1.24 * ((l_target as f64 - 15.0).cbrt()) - 1.8).max(0.5)
    } else {
        0.5
    };
    let d0_sq = d0 * d0;

    for i in 0..n {
        let x = [coords_p[i][0] - cp[0], coords_p[i][1] - cp[1], coords_p[i][2] - cp[2]];
        let y = [coords_q[i][0] - cq[0], coords_q[i][1] - cq[1], coords_q[i][2] - cq[2]];

        let rx = apply_rotation(&rot, &x);
        let dx = rx[0] - y[0];
        let dy = rx[1] - y[1];
        let dz = rx[2] - y[2];
        let dist_sq = dx * dx + dy * dy + dz * dz;

        sum_sq_diff += dist_sq;
        tm_sum += 1.0 / (1.0 + dist_sq / d0_sq);
    }

    let rmsd = (sum_sq_diff / nf).sqrt();
    let tm_score = (tm_sum / l_target as f64).clamp(0.0, 1.0);

    Ok(SuperpositionResult {
        rmsd,
        tm_score,
        centroid_p: cp,
        centroid_q: cq,
        rotation_matrix: rot,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical_superposition() {
        let coords = [
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 10.0],
            [2.0, 3.0, 4.0],
        ];

        let res = calculate_rmsd_and_tm_score(&coords, &coords, 4).unwrap();
        assert!(res.rmsd < 1e-6);
        assert!((res.tm_score - 1.0).abs() < 1e-4);
    }

    #[test]
    fn test_rigid_translation_and_rotation_rmsd() {
        let c1 = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];

        // Translate by (10, 10, 10)
        let c2: Vec<[f64; 3]> = c1.iter().map(|p| [p[0] + 10.0, p[1] + 10.0, p[2] + 10.0]).collect();

        let res = calculate_rmsd_and_tm_score(&c1, &c2, 4).unwrap();
        assert!(res.rmsd < 1e-6); // Pure translation has 0 RMSD
        assert!((res.tm_score - 1.0).abs() < 1e-4);
    }

    #[test]
    fn test_benchmark_1crn_vs_1ubq_crystallographic_segments() {
        // Real crystallographic Cα coordinates from PDB 1CRN (residues 7-14 alpha-helix)
        let crambin_1crn_helix = [
            [13.462, 10.916,  9.520],
            [13.197,  7.185,  8.951],
            [ 9.671,  6.774, 10.370],
            [ 8.953,  9.711, 12.662],
            [12.316, 11.464, 13.093],
            [13.662,  8.243, 11.584],
            [10.806,  6.357, 13.214],
            [ 8.587,  8.627, 15.275],
        ];

        // Real crystallographic Cα coordinates from PDB 1UBQ (residues 23-30 alpha-helix)
        let ubiquitin_1ubq_helix = [
            [25.632, 32.893, 14.375],
            [24.918, 29.381, 13.080],
            [28.504, 28.187, 12.569],
            [29.839, 30.892, 14.885],
            [27.274, 33.428, 16.059],
            [25.568, 30.506, 17.752],
            [28.513, 28.431, 16.597],
            [30.923, 30.461, 18.749],
        ];

        // Comparing two distinct alpha-helices from 1CRN and 1UBQ yields finite local RMSD (< 3.5 Å)
        let res = calculate_rmsd_and_tm_score(&crambin_1crn_helix, &ubiquitin_1ubq_helix, 46).unwrap();
        assert!(res.rmsd > 0.5 && res.rmsd < 3.5);
        assert!(res.tm_score > 0.05 && res.tm_score < 1.0);
    }
}
