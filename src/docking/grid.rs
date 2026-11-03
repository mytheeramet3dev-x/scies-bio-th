//! 3D spatial grid, bounding boxes, and binding pocket cavity detection.

use crate::structure::Atom3D;
use crate::BioError;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// 3D Axis-Aligned Bounding Box (AABB).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct BoundingBox3D {
    /// Minimum coordinates [min_x, min_y, min_z].
    pub min: [f64; 3],
    /// Maximum coordinates [max_x, max_y, max_z].
    pub max: [f64; 3],
}

impl BoundingBox3D {
    /// Creates a bounding box around a set of atoms with a specified padding in Ångströms.
    pub fn from_atoms(atoms: &[Atom3D], padding: f64) -> crate::Result<Self> {
        if atoms.is_empty() {
            return Err(BioError::InvalidInput("Cannot create bounding box from empty atom list".to_string()));
        }

        let mut min_x = f64::INFINITY;
        let mut min_y = f64::INFINITY;
        let mut min_z = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        let mut max_z = f64::NEG_INFINITY;

        for a in atoms {
            if a.x < min_x { min_x = a.x; }
            if a.x > max_x { max_x = a.x; }
            if a.y < min_y { min_y = a.y; }
            if a.y > max_y { max_y = a.y; }
            if a.z < min_z { min_z = a.z; }
            if a.z > max_z { max_z = a.z; }
        }

        Ok(Self {
            min: [min_x - padding, min_y - padding, min_z - padding],
            max: [max_x + padding, max_y + padding, max_z + padding],
        })
    }

    /// Extents / dimensions along [X, Y, Z] axes.
    pub fn dimensions(&self) -> [f64; 3] {
        [
            self.max[0] - self.min[0],
            self.max[1] - self.min[1],
            self.max[2] - self.min[2],
        ]
    }
}

/// A detected binding pocket cavity in a macromolecule.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct BindingPocket3D {
    /// Identifier index of the pocket.
    pub id: usize,
    /// Center of mass of the cavity [X, Y, Z].
    pub center: [f64; 3],
    /// Estimated cavity volume in cubic Ångströms ($Å^3$).
    pub volume_angstrom3: f64,
    /// Number of grid points inside the cavity pocket.
    pub grid_point_count: usize,
}

/// Regular 3D grid for cavity detection and docking grid potentials.
#[derive(Debug, Clone)]
pub struct CavityGrid3D {
    /// Origin coordinate of grid index (0,0,0).
    pub origin: [f64; 3],
    /// Grid spacing / voxel resolution in Ångströms (e.g. 1.0 Å).
    pub spacing: f64,
    /// Grid dimensions (nx, ny, nz).
    pub dims: [usize; 3],
    /// Voxel values (true if occupied by protein atom).
    pub occupied: Vec<bool>,
}

impl CavityGrid3D {
    /// Initializes a 3D grid encompassing the specified bounding box.
    pub fn new(bbox: &BoundingBox3D, spacing: f64) -> Self {
        let dims_f = bbox.dimensions();
        let nx = (dims_f[0] / spacing).ceil() as usize + 1;
        let ny = (dims_f[1] / spacing).ceil() as usize + 1;
        let nz = (dims_f[2] / spacing).ceil() as usize + 1;

        let total_voxels = nx * ny * nz;
        Self {
            origin: bbox.min,
            spacing,
            dims: [nx, ny, nz],
            occupied: vec![false; total_voxels],
        }
    }

    /// Converts 3D grid coordinate (i, j, k) to 1D index.
    #[inline]
    pub fn index_of(&self, i: usize, j: usize, k: usize) -> usize {
        i + j * self.dims[0] + k * (self.dims[0] * self.dims[1])
    }

    /// Converts 1D index to Cartesian coordinate [x, y, z].
    pub fn coord_of(&self, idx: usize) -> [f64; 3] {
        let area = self.dims[0] * self.dims[1];
        let k = idx / area;
        let rem = idx % area;
        let j = rem / self.dims[0];
        let i = rem % self.dims[0];

        [
            self.origin[0] + (i as f64) * self.spacing,
            self.origin[1] + (j as f64) * self.spacing,
            self.origin[2] + (k as f64) * self.spacing,
        ]
    }

    /// Marks grid points within Van der Waals radius of atoms as occupied.
    pub fn populate_atoms(&mut self, atoms: &[Atom3D], probe_radius: f64) {
        for atom in atoms {
            let r_atom = match atom.element.as_str() {
                "H" => 1.20,
                "C" => 1.70,
                "N" => 1.55,
                "O" => 1.52,
                "P" => 1.80,
                "S" => 1.80,
                _ => 1.65,
            };
            let cutoff = r_atom + probe_radius;
            let cutoff_sq = cutoff * cutoff;

            let min_i = (((atom.x - cutoff - self.origin[0]) / self.spacing).floor() as isize).max(0) as usize;
            let max_i = (((atom.x + cutoff - self.origin[0]) / self.spacing).ceil() as usize).min(self.dims[0] - 1);
            let min_j = (((atom.y - cutoff - self.origin[1]) / self.spacing).floor() as isize).max(0) as usize;
            let max_j = (((atom.y + cutoff - self.origin[1]) / self.spacing).ceil() as usize).min(self.dims[1] - 1);
            let min_k = (((atom.z - cutoff - self.origin[2]) / self.spacing).floor() as isize).max(0) as usize;
            let max_k = (((atom.z + cutoff - self.origin[2]) / self.spacing).ceil() as usize).min(self.dims[2] - 1);

            for i in min_i..=max_i {
                let gx = self.origin[0] + (i as f64) * self.spacing;
                let dx = gx - atom.x;
                for j in min_j..=max_j {
                    let gy = self.origin[1] + (j as f64) * self.spacing;
                    let dy = gy - atom.y;
                    for k in min_k..=max_k {
                        let gz = self.origin[2] + (k as f64) * self.spacing;
                        let dz = gz - atom.z;
                        if dx * dx + dy * dy + dz * dz <= cutoff_sq {
                            let idx = self.index_of(i, j, k);
                            self.occupied[idx] = true;
                        }
                    }
                }
            }
        }
    }

    /// Detects candidate binding cavities / pockets using geometric ray casting.
    ///
    /// Finds unoccupied points that are buried in at least 4 out of 6 cardinal directions
    /// within `max_buried_dist` Ångströms.
    pub fn detect_pockets(&self, min_cluster_points: usize, max_buried_dist: f64) -> Vec<BindingPocket3D> {
        let max_steps = (max_buried_dist / self.spacing).ceil() as usize;
        let mut cavity_mask = vec![false; self.occupied.len()];

        let (nx, ny, nz) = (self.dims[0], self.dims[1], self.dims[2]);

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                for k in 1..nz - 1 {
                    let idx = self.index_of(i, j, k);
                    if self.occupied[idx] {
                        continue;
                    }

                    // Check 6 cardinal directions
                    let mut buried_directions = 0;

                    // -X
                    for step in 1..=max_steps {
                        if i < step { break; }
                        if self.occupied[self.index_of(i - step, j, k)] {
                            buried_directions += 1;
                            break;
                        }
                    }
                    // +X
                    for step in 1..=max_steps {
                        if i + step >= nx { break; }
                        if self.occupied[self.index_of(i + step, j, k)] {
                            buried_directions += 1;
                            break;
                        }
                    }
                    // -Y
                    for step in 1..=max_steps {
                        if j < step { break; }
                        if self.occupied[self.index_of(i, j - step, k)] {
                            buried_directions += 1;
                            break;
                        }
                    }
                    // +Y
                    for step in 1..=max_steps {
                        if j + step >= ny { break; }
                        if self.occupied[self.index_of(i, j + step, k)] {
                            buried_directions += 1;
                            break;
                        }
                    }
                    // -Z
                    for step in 1..=max_steps {
                        if k < step { break; }
                        if self.occupied[self.index_of(i, j, k - step)] {
                            buried_directions += 1;
                            break;
                        }
                    }
                    // +Z
                    for step in 1..=max_steps {
                        if k + step >= nz { break; }
                        if self.occupied[self.index_of(i, j, k + step)] {
                            buried_directions += 1;
                            break;
                        }
                    }

                    if buried_directions >= 4 {
                        cavity_mask[idx] = true;
                    }
                }
            }
        }

        // Connected component labeling (BFS)
        let mut visited = vec![false; cavity_mask.len()];
        let mut pockets = Vec::new();
        let voxel_volume = self.spacing.powi(3);

        for idx in 0..cavity_mask.len() {
            if !cavity_mask[idx] || visited[idx] {
                continue;
            }

            let mut cluster = Vec::new();
            let mut queue = std::collections::VecDeque::new();
            visited[idx] = true;
            queue.push_back(idx);

            while let Some(curr) = queue.pop_front() {
                cluster.push(curr);

                let area = nx * ny;
                let k = curr / area;
                let rem = curr % area;
                let j = rem / nx;
                let i = rem % nx;

                for (di, dj, dk) in [
                    (-1isize, 0, 0), (1, 0, 0),
                    (0, -1, 0), (0, 1, 0),
                    (0, 0, -1), (0, 0, 1),
                ] {
                    let ni = i as isize + di;
                    let nj = j as isize + dj;
                    let nk = k as isize + dk;

                    if ni >= 0 && ni < nx as isize &&
                       nj >= 0 && nj < ny as isize &&
                       nk >= 0 && nk < nz as isize {
                        let n_idx = self.index_of(ni as usize, nj as usize, nk as usize);
                        if cavity_mask[n_idx] && !visited[n_idx] {
                            visited[n_idx] = true;
                            queue.push_back(n_idx);
                        }
                    }
                }
            }

            if cluster.len() >= min_cluster_points {
                let mut sum_x = 0.0;
                let mut sum_y = 0.0;
                let mut sum_z = 0.0;

                for &p_idx in &cluster {
                    let c = self.coord_of(p_idx);
                    sum_x += c[0];
                    sum_y += c[1];
                    sum_z += c[2];
                }

                let count = cluster.len() as f64;
                let pocket = BindingPocket3D {
                    id: pockets.len() + 1,
                    center: [sum_x / count, sum_y / count, sum_z / count],
                    volume_angstrom3: count * voxel_volume,
                    grid_point_count: cluster.len(),
                };
                pockets.push(pocket);
            }
        }

        pockets.sort_by(|a, b| b.volume_angstrom3.partial_cmp(&a.volume_angstrom3).unwrap());
        pockets
    }
}
