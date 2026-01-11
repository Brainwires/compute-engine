//! Bohm Quantum Potential Calculations
//!
//! Implements the Bohm quantum potential (quantum stress tensor contribution)
//! arising from the Madelung transformation of the Schrödinger equation:
//!
//! Q = (ℏ²/2m) × ∇²√ρ / √ρ
//!
//! This term appears in the quantum Navier-Stokes equations as a regularizing
//! force that acts at small scales. The potential can be rewritten in a
//! numerically stable form:
//!
//! ∇²√ρ / √ρ = (1/2) × [∇²ρ/ρ - (1/2) × |∇ρ|²/ρ²]
//!
//! which avoids division by √ρ when ρ → 0.

use super::H_BAR;
use ndarray::{Array2, Array3};
use serde::{Deserialize, Serialize};

/// Configuration for Bohm potential calculations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BohmPotentialConfig {
    /// Reduced Planck constant (J·s), defaults to physical ℏ
    pub hbar: f64,

    /// Particle mass (kg)
    pub mass: f64,

    /// Minimum density threshold to avoid numerical instability
    /// Default: 1e-15
    pub min_density: f64,
}

impl Default for BohmPotentialConfig {
    fn default() -> Self {
        Self {
            hbar: H_BAR,
            mass: 1.0,
            min_density: 1e-15,
        }
    }
}

/// Bohm quantum potential calculator
///
/// Provides methods to calculate the Bohm potential Q on 1D, 2D, and 3D grids.
/// The potential acts as a quantum correction term in the Navier-Stokes equations.
#[derive(Debug, Clone)]
pub struct BohmPotential {
    config: BohmPotentialConfig,
    /// Prefactor: ℏ²/(2m)
    prefactor: f64,
}

impl BohmPotential {
    /// Create a new Bohm potential calculator
    pub fn new(config: BohmPotentialConfig) -> Self {
        let prefactor = config.hbar * config.hbar / (2.0 * config.mass);
        Self { config, prefactor }
    }

    /// Create with physical constants for a given particle mass
    pub fn with_mass(mass: f64) -> Self {
        Self::new(BohmPotentialConfig {
            hbar: H_BAR,
            mass,
            min_density: 1e-15,
        })
    }

    /// Get the ℏ²/(2m) prefactor
    pub fn prefactor(&self) -> f64 {
        self.prefactor
    }

    /// Calculate Bohm potential on a 1D grid
    ///
    /// Uses the numerically stable form:
    /// Q = (ℏ²/2m) × (1/2) × [∇²ρ/ρ - (1/2) × (∇ρ)²/ρ²]
    ///
    /// # Arguments
    /// * `density` - Density field ρ(x)
    /// * `dx` - Grid spacing
    ///
    /// # Returns
    /// Bohm potential Q(x) on the grid
    pub fn calculate_1d(&self, density: &[f64], dx: f64) -> Vec<f64> {
        let n = density.len();
        if n < 3 {
            return vec![0.0; n];
        }

        let mut q = vec![0.0; n];
        let dx2 = dx * dx;

        for i in 1..n - 1 {
            let rho = density[i].max(self.config.min_density);
            let rho_l = density[i - 1].max(self.config.min_density);
            let rho_r = density[i + 1].max(self.config.min_density);

            // ∇²ρ using central differences
            let laplacian_rho = (rho_l - 2.0 * rho + rho_r) / dx2;

            // (∇ρ)² using central differences
            let grad_rho = (rho_r - rho_l) / (2.0 * dx);
            let grad_rho_sq = grad_rho * grad_rho;

            // Q = (ℏ²/2m) × (1/2) × [∇²ρ/ρ - (1/2) × (∇ρ)²/ρ²]
            q[i] = self.prefactor * 0.5 * (laplacian_rho / rho - 0.5 * grad_rho_sq / (rho * rho));
        }

        // Boundary conditions: zero gradient (Neumann)
        q[0] = q[1];
        q[n - 1] = q[n - 2];

        q
    }

    /// Calculate the quantum force F_Q = -∇Q on a 1D grid
    ///
    /// This is the force term that appears in the QNS momentum equation:
    /// F_Q = -ρ ∂Q/∂x
    ///
    /// # Arguments
    /// * `density` - Density field ρ(x)
    /// * `dx` - Grid spacing
    ///
    /// # Returns
    /// Quantum force per unit volume on the grid
    pub fn quantum_force_1d(&self, density: &[f64], dx: f64) -> Vec<f64> {
        let q = self.calculate_1d(density, dx);
        let n = q.len();
        let mut force = vec![0.0; n];

        for i in 1..n - 1 {
            let rho = density[i].max(self.config.min_density);
            // F_Q = -ρ ∂Q/∂x
            let grad_q = (q[i + 1] - q[i - 1]) / (2.0 * dx);
            force[i] = -rho * grad_q;
        }

        force[0] = force[1];
        force[n - 1] = force[n - 2];

        force
    }

    /// Calculate Bohm potential on a 2D grid
    ///
    /// # Arguments
    /// * `density` - 2D density field ρ(x, y)
    /// * `dx` - Grid spacing in x
    /// * `dy` - Grid spacing in y
    ///
    /// # Returns
    /// Bohm potential Q(x, y) on the grid
    pub fn calculate_2d(&self, density: &Array2<f64>, dx: f64, dy: f64) -> Array2<f64> {
        let (nx, ny) = density.dim();
        let mut q = Array2::<f64>::zeros((nx, ny));

        if nx < 3 || ny < 3 {
            return q;
        }

        let dx2 = dx * dx;
        let dy2 = dy * dy;

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let rho = density[[i, j]].max(self.config.min_density);
                let rho_xm = density[[i - 1, j]].max(self.config.min_density);
                let rho_xp = density[[i + 1, j]].max(self.config.min_density);
                let rho_ym = density[[i, j - 1]].max(self.config.min_density);
                let rho_yp = density[[i, j + 1]].max(self.config.min_density);

                // ∇²ρ = ∂²ρ/∂x² + ∂²ρ/∂y²
                let laplacian_rho = (rho_xm - 2.0 * rho + rho_xp) / dx2
                    + (rho_ym - 2.0 * rho + rho_yp) / dy2;

                // |∇ρ|² = (∂ρ/∂x)² + (∂ρ/∂y)²
                let grad_rho_x = (rho_xp - rho_xm) / (2.0 * dx);
                let grad_rho_y = (rho_yp - rho_ym) / (2.0 * dy);
                let grad_rho_sq = grad_rho_x * grad_rho_x + grad_rho_y * grad_rho_y;

                // Q = (ℏ²/2m) × (1/2) × [∇²ρ/ρ - (1/2) × |∇ρ|²/ρ²]
                q[[i, j]] = self.prefactor * 0.5 * (laplacian_rho / rho - 0.5 * grad_rho_sq / (rho * rho));
            }
        }

        // Apply Neumann boundary conditions
        for i in 0..nx {
            q[[i, 0]] = q[[i, 1.min(ny - 1)]];
            q[[i, ny - 1]] = q[[i, (ny - 2).max(0)]];
        }
        for j in 0..ny {
            q[[0, j]] = q[[1.min(nx - 1), j]];
            q[[nx - 1, j]] = q[[(nx - 2).max(0), j]];
        }

        q
    }

    /// Calculate quantum force components on a 2D grid
    ///
    /// Returns (F_Q_x, F_Q_y) where F_Q = -ρ∇Q
    pub fn quantum_force_2d(
        &self,
        density: &Array2<f64>,
        dx: f64,
        dy: f64,
    ) -> (Array2<f64>, Array2<f64>) {
        let q = self.calculate_2d(density, dx, dy);
        let (nx, ny) = q.dim();

        let mut fx = Array2::<f64>::zeros((nx, ny));
        let mut fy = Array2::<f64>::zeros((nx, ny));

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let rho = density[[i, j]].max(self.config.min_density);

                // F_Q = -ρ ∇Q
                let grad_q_x = (q[[i + 1, j]] - q[[i - 1, j]]) / (2.0 * dx);
                let grad_q_y = (q[[i, j + 1]] - q[[i, j - 1]]) / (2.0 * dy);

                fx[[i, j]] = -rho * grad_q_x;
                fy[[i, j]] = -rho * grad_q_y;
            }
        }

        (fx, fy)
    }

    /// Calculate Bohm potential on a 3D grid
    ///
    /// # Arguments
    /// * `density` - 3D density field ρ(x, y, z)
    /// * `dx`, `dy`, `dz` - Grid spacings
    ///
    /// # Returns
    /// Bohm potential Q(x, y, z) on the grid
    pub fn calculate_3d(
        &self,
        density: &Array3<f64>,
        dx: f64,
        dy: f64,
        dz: f64,
    ) -> Array3<f64> {
        let (nx, ny, nz) = density.dim();
        let mut q = Array3::<f64>::zeros((nx, ny, nz));

        if nx < 3 || ny < 3 || nz < 3 {
            return q;
        }

        let dx2 = dx * dx;
        let dy2 = dy * dy;
        let dz2 = dz * dz;

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                for k in 1..nz - 1 {
                    let rho = density[[i, j, k]].max(self.config.min_density);
                    let rho_xm = density[[i - 1, j, k]].max(self.config.min_density);
                    let rho_xp = density[[i + 1, j, k]].max(self.config.min_density);
                    let rho_ym = density[[i, j - 1, k]].max(self.config.min_density);
                    let rho_yp = density[[i, j + 1, k]].max(self.config.min_density);
                    let rho_zm = density[[i, j, k - 1]].max(self.config.min_density);
                    let rho_zp = density[[i, j, k + 1]].max(self.config.min_density);

                    // ∇²ρ = ∂²ρ/∂x² + ∂²ρ/∂y² + ∂²ρ/∂z²
                    let laplacian_rho = (rho_xm - 2.0 * rho + rho_xp) / dx2
                        + (rho_ym - 2.0 * rho + rho_yp) / dy2
                        + (rho_zm - 2.0 * rho + rho_zp) / dz2;

                    // |∇ρ|²
                    let grad_rho_x = (rho_xp - rho_xm) / (2.0 * dx);
                    let grad_rho_y = (rho_yp - rho_ym) / (2.0 * dy);
                    let grad_rho_z = (rho_zp - rho_zm) / (2.0 * dz);
                    let grad_rho_sq = grad_rho_x * grad_rho_x
                        + grad_rho_y * grad_rho_y
                        + grad_rho_z * grad_rho_z;

                    // Q = (ℏ²/2m) × (1/2) × [∇²ρ/ρ - (1/2) × |∇ρ|²/ρ²]
                    q[[i, j, k]] = self.prefactor
                        * 0.5
                        * (laplacian_rho / rho - 0.5 * grad_rho_sq / (rho * rho));
                }
            }
        }

        // Apply Neumann boundary conditions on all faces
        for i in 0..nx {
            for j in 0..ny {
                q[[i, j, 0]] = q[[i, j, 1.min(nz - 1)]];
                q[[i, j, nz - 1]] = q[[i, j, (nz - 2).max(0)]];
            }
        }
        for i in 0..nx {
            for k in 0..nz {
                q[[i, 0, k]] = q[[i, 1.min(ny - 1), k]];
                q[[i, ny - 1, k]] = q[[i, (ny - 2).max(0), k]];
            }
        }
        for j in 0..ny {
            for k in 0..nz {
                q[[0, j, k]] = q[[1.min(nx - 1), j, k]];
                q[[nx - 1, j, k]] = q[[(nx - 2).max(0), j, k]];
            }
        }

        q
    }

    /// Calculate quantum force components on a 3D grid
    ///
    /// Returns (F_Q_x, F_Q_y, F_Q_z) where F_Q = -ρ∇Q
    pub fn quantum_force_3d(
        &self,
        density: &Array3<f64>,
        dx: f64,
        dy: f64,
        dz: f64,
    ) -> (Array3<f64>, Array3<f64>, Array3<f64>) {
        let q = self.calculate_3d(density, dx, dy, dz);
        let (nx, ny, nz) = q.dim();

        let mut fx = Array3::<f64>::zeros((nx, ny, nz));
        let mut fy = Array3::<f64>::zeros((nx, ny, nz));
        let mut fz = Array3::<f64>::zeros((nx, ny, nz));

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                for k in 1..nz - 1 {
                    let rho = density[[i, j, k]].max(self.config.min_density);

                    let grad_q_x = (q[[i + 1, j, k]] - q[[i - 1, j, k]]) / (2.0 * dx);
                    let grad_q_y = (q[[i, j + 1, k]] - q[[i, j - 1, k]]) / (2.0 * dy);
                    let grad_q_z = (q[[i, j, k + 1]] - q[[i, j, k - 1]]) / (2.0 * dz);

                    fx[[i, j, k]] = -rho * grad_q_x;
                    fy[[i, j, k]] = -rho * grad_q_y;
                    fz[[i, j, k]] = -rho * grad_q_z;
                }
            }
        }

        (fx, fy, fz)
    }
}

/// Calculate Bohm potential for a Gaussian wavepacket (analytical validation)
///
/// For ψ(x) = A exp(-x²/2σ²), we have ρ = |ψ|² = A² exp(-x²/σ²)
///
/// Derivation:
/// √ρ = A exp(-x²/2σ²)
/// ∇√ρ = -x/σ² × √ρ
/// ∇²√ρ = √ρ × (-1/σ² + x²/σ⁴)
/// ∇²√ρ / √ρ = -1/σ² + x²/σ⁴
///
/// The analytical Bohm potential is:
/// Q(x) = (ℏ²/2m) × (-1/σ² + x²/σ⁴)
///
/// At x = 0: Q(0) = -(ℏ²/2m)/σ² < 0
///
/// # Arguments
/// * `x` - Position
/// * `sigma` - Width of Gaussian
/// * `hbar` - Reduced Planck constant
/// * `mass` - Particle mass
pub fn bohm_potential_gaussian(x: f64, sigma: f64, hbar: f64, mass: f64) -> f64 {
    let prefactor = hbar * hbar / (2.0 * mass);
    let sigma2 = sigma * sigma;
    let sigma4 = sigma2 * sigma2;

    // ∇²√ρ / √ρ = -1/σ² + x²/σ⁴
    prefactor * (-1.0 / sigma2 + x * x / sigma4)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    const ELECTRON_MASS: f64 = 9.109e-31; // kg

    #[test]
    fn test_bohm_potential_1d_uniform() {
        // Uniform density should give zero Bohm potential
        let density = vec![1.0; 100];
        let bohm = BohmPotential::with_mass(ELECTRON_MASS);
        let q = bohm.calculate_1d(&density, 1e-9);

        for i in 1..99 {
            assert!(
                q[i].abs() < 1e-20,
                "Uniform density should give Q ≈ 0, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn test_bohm_potential_1d_gaussian() {
        // Compare numerical to analytical for Gaussian
        // Use high resolution (small dx) for accurate numerical derivatives
        let n = 1001;
        let dx = 1e-10; // 0.1 nm spacing (finer grid)
        let sigma = 10e-9; // 10 nm width (100 points per sigma)
        let mass = ELECTRON_MASS;

        let mut density = vec![0.0; n];
        let x0 = (n / 2) as f64 * dx;

        for i in 0..n {
            let x = i as f64 * dx - x0;
            // ρ = exp(-x²/σ²)
            density[i] = (-x * x / (sigma * sigma)).exp();
        }

        let bohm = BohmPotential::with_mass(mass);
        let q_numerical = bohm.calculate_1d(&density, dx);

        // Check center point against analytical
        let x_center = 0.0;
        let q_analytical = bohm_potential_gaussian(x_center, sigma, H_BAR, mass);
        let q_numerical_center = q_numerical[n / 2];

        // Allow some numerical error due to discretization (~1% with fine grid)
        let rel_error = ((q_numerical_center - q_analytical) / q_analytical).abs();
        assert!(
            rel_error < 0.05,
            "Relative error at center: {:.2}%, numerical={:.3e}, analytical={:.3e}",
            rel_error * 100.0,
            q_numerical_center,
            q_analytical
        );
    }

    #[test]
    fn test_bohm_potential_2d_uniform() {
        let density = Array2::<f64>::from_elem((50, 50), 1.0);
        let bohm = BohmPotential::with_mass(ELECTRON_MASS);
        let q = bohm.calculate_2d(&density, 1e-9, 1e-9);

        for i in 2..48 {
            for j in 2..48 {
                assert!(
                    q[[i, j]].abs() < 1e-20,
                    "Uniform density should give Q ≈ 0"
                );
            }
        }
    }

    #[test]
    fn test_quantum_force_1d() {
        // Gaussian density should produce antisymmetric force around center
        let n = 101;
        let dx = 1e-9;
        let sigma = 10e-9;

        let mut density = vec![0.0; n];
        let x0 = (n / 2) as f64 * dx;

        for i in 0..n {
            let x = i as f64 * dx - x0;
            density[i] = 1.0 + 0.5 * (-x * x / (sigma * sigma)).exp();
        }

        let bohm = BohmPotential::with_mass(ELECTRON_MASS);
        let force = bohm.quantum_force_1d(&density, dx);

        // Force should be approximately antisymmetric around center
        // Due to numerical discretization, we check that center force is much smaller
        // than off-center forces
        let center = n / 2;
        let off_center = n / 4;

        // Force at center should be small relative to off-center forces
        let force_scale = force[off_center].abs().max(force[center + off_center].abs()).max(1e-50);
        let relative_center_force = force[center].abs() / force_scale;

        assert!(
            relative_center_force < 0.1,
            "Force at center should be << off-center forces, got ratio {}",
            relative_center_force
        );

        // Force should be approximately antisymmetric: F(x) ≈ -F(-x)
        let f_left = force[center - 10];
        let f_right = force[center + 10];
        let antisymmetry_error = (f_left + f_right).abs() / (f_left.abs() + f_right.abs() + 1e-50);

        assert!(
            antisymmetry_error < 0.1,
            "Force should be antisymmetric, got error {}",
            antisymmetry_error
        );
    }

    #[test]
    fn test_bohm_prefactor() {
        let bohm = BohmPotential::with_mass(ELECTRON_MASS);
        let expected = H_BAR * H_BAR / (2.0 * ELECTRON_MASS);
        assert!((bohm.prefactor() - expected).abs() < 1e-50);
    }

    #[test]
    fn test_min_density_protection() {
        // Test that very small densities don't cause numerical issues
        let mut density = vec![1e-20; 100];
        density[50] = 1.0; // spike

        let bohm = BohmPotential::with_mass(ELECTRON_MASS);
        let q = bohm.calculate_1d(&density, 1e-9);

        // Should not contain NaN or Inf
        for q_val in &q {
            assert!(q_val.is_finite(), "Q should be finite, got {}", q_val);
        }
    }
}
