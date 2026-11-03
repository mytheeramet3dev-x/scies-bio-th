//! Fluorescent microscopy spot detection and sub-pixel 2D Gaussian centroid fitting (`imaging::spot`).

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A 2D single-channel fluorescence intensity image.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ImageChannel2D {
    /// Image width in pixels.
    pub width: usize,
    /// Image height in pixels.
    pub height: usize,
    /// Row-major pixel intensities (`height * width`).
    pub pixels: Vec<f64>,
}

impl ImageChannel2D {
    /// Gets pixel intensity at `(x, y)`.
    #[inline]
    pub fn get(&self, x: usize, y: usize) -> f64 {
        self.pixels[y * self.width + x]
    }
}

/// Detected fluorescent puncta / focus (e.g. $\gamma$-H2AX DNA damage focus or smFISH spot).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct FluorescentSpot {
    /// Sub-pixel fitted X centroid coordinate.
    pub x_centroid: f64,
    /// Sub-pixel fitted Y centroid coordinate.
    pub y_centroid: f64,
    /// Peak pixel intensity above local background.
    pub peak_intensity: f64,
    /// Estimated Gaussian width $\sigma$ in pixels.
    pub sigma_xy: f64,
    /// Integrated spot intensity in a $3 \times 3$ neighborhood.
    pub integrated_intensity: f64,
}

/// Detects local $3 \times 3$ intensity maxima above `threshold` and refines sub-pixel centroids
/// via 2D log-parabolic (Gaussian) peak interpolation.
pub fn detect_fluorescent_spots(image: &ImageChannel2D, threshold: f64) -> Vec<FluorescentSpot> {
    if image.width < 3 || image.height < 3 || image.pixels.len() < image.width * image.height {
        return Vec::new();
    }

    let mut spots = Vec::new();

    for y in 1..(image.height - 1) {
        for x in 1..(image.width - 1) {
            let center = image.get(x, y);
            if center < threshold {
                continue;
            }

            let mut is_local_max = true;
            let mut sum_3x3 = 0.0;
            for dy in -1isize..=1 {
                for dx in -1isize..=1 {
                    let v = image.get((x as isize + dx) as usize, (y as isize + dy) as usize);
                    sum_3x3 += v;
                    if (dx != 0 || dy != 0) && v >= center {
                        is_local_max = false;
                    }
                }
            }

            if !is_local_max {
                continue;
            }

            // Sub-pixel Gaussian log-parabola refinement along X and Y
            let left = image.get(x - 1, y).max(1.0).ln();
            let mid = center.max(1.0).ln();
            let right = image.get(x + 1, y).max(1.0).ln();
            let denom_x = 2.0 * (left - 2.0 * mid + right);
            let dx_sub = if denom_x.abs() > 1e-6 {
                ((left - right) / denom_x).clamp(-0.5, 0.5)
            } else {
                0.0
            };

            let up = image.get(x, y - 1).max(1.0).ln();
            let down = image.get(x, y + 1).max(1.0).ln();
            let denom_y = 2.0 * (up - 2.0 * mid + down);
            let dy_sub = if denom_y.abs() > 1e-6 {
                ((up - down) / denom_y).clamp(-0.5, 0.5)
            } else {
                0.0
            };

            let curv = (-(denom_x + denom_y) / 4.0).max(0.15);
            let sigma_xy = (1.0 / curv.sqrt()).clamp(0.5, 4.0);

            spots.push(FluorescentSpot {
                x_centroid: (x as f64) + dx_sub,
                y_centroid: (y as f64) + dy_sub,
                peak_intensity: center,
                sigma_xy,
                integrated_intensity: sum_3x3,
            });
        }
    }

    spots
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_fluorescent_spots_gaussian() {
        let mut pixels = vec![10.0; 25];
        // Place a synthetic 2D Gaussian spot centered at (2, 2) on a 5x5 grid
        pixels[2 * 5 + 2] = 250.0;
        pixels[2 * 5 + 1] = 120.0;
        pixels[2 * 5 + 3] = 120.0;
        pixels[5 + 2] = 120.0;
        pixels[3 * 5 + 2] = 120.0;

        let img = ImageChannel2D {
            width: 5,
            height: 5,
            pixels,
        };
        let spots = detect_fluorescent_spots(&img, 100.0);
        assert_eq!(spots.len(), 1);
        assert!((spots[0].x_centroid - 2.0).abs() < 1e-4);
        assert!((spots[0].y_centroid - 2.0).abs() < 1e-4);
    }
}
