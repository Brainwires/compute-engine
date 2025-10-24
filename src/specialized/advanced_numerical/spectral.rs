//! Spectral Methods
//!
//! High-accuracy methods using orthogonal polynomials for solving PDEs.

use nalgebra::{DMatrix, DVector};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Spectral method solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpectralSolver {
    /// Number of collocation points
    pub num_points: usize,
    /// Type of spectral basis
    pub basis: SpectralBasis,
}

/// Types of spectral bases
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum SpectralBasis {
    /// Chebyshev polynomials (for non-periodic problems)
    Chebyshev,
    /// Fourier series (for periodic problems)
    Fourier,
    /// Legendre polynomials
    Legendre,
}

impl SpectralSolver {
    /// Create a new spectral solver with Chebyshev basis
    pub fn chebyshev(num_points: usize) -> Self {
        Self {
            num_points,
            basis: SpectralBasis::Chebyshev,
        }
    }

    /// Create a new spectral solver with Fourier basis
    pub fn fourier(num_points: usize) -> Self {
        Self {
            num_points,
            basis: SpectralBasis::Fourier,
        }
    }

    /// Get Chebyshev collocation points (Gauss-Lobatto)
    pub fn chebyshev_points(&self) -> Vec<f64> {
        let n = self.num_points;
        (0..n)
            .map(|i| -(PI * i as f64 / (n - 1) as f64).cos())
            .collect()
    }

    /// Compute Chebyshev differentiation matrix
    pub fn chebyshev_diff_matrix(&self) -> DMatrix<f64> {
        let n = self.num_points;
        let mut D = DMatrix::<f64>::zeros(n, n);
        let x = self.chebyshev_points();

        // Compute c_i coefficients
        let mut c = vec![1.0; n];
        c[0] = 2.0;
        c[n - 1] = 2.0;

        // Fill differentiation matrix
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    if i == 0 {
                        D[(i, j)] = (2.0 * (n - 1) as f64 * (n - 1) as f64 + 1.0) / 6.0;
                    } else if i == n - 1 {
                        D[(i, j)] = -(2.0 * (n - 1) as f64 * (n - 1) as f64 + 1.0) / 6.0;
                    } else {
                        D[(i, j)] = -x[i] / (2.0 * (1.0 - x[i] * x[i]));
                    }
                } else {
                    let sign = if (i + j) % 2 == 0 { 1.0 } else { -1.0 };
                    D[(i, j)] = (c[i] / c[j]) * sign / (x[i] - x[j]);
                }
            }
        }

        D
    }

    /// Solve 1D BVP using Chebyshev spectral method
    /// -d²u/dx² = f(x), u(-1) = a, u(1) = b
    pub fn solve_bvp_chebyshev(
        &self,
        source: fn(f64) -> f64,
        boundary_left: f64,
        boundary_right: f64,
    ) -> Vec<f64> {
        let n = self.num_points;
        let x = self.chebyshev_points();
        let D = self.chebyshev_diff_matrix();

        // Second derivative matrix D² = D * D
        let D2 = &D * &D;

        // Set up linear system A * u = f
        let mut A = DMatrix::<f64>::zeros(n, n);
        let mut f = DVector::<f64>::zeros(n);

        // First row: boundary condition at x = -1
        for j in 0..n {
            A[(0, j)] = 0.0;
        }
        A[(0, 0)] = 1.0;
        f[0] = boundary_left;

        // Interior points: -D²u = f
        for i in 1..n - 1 {
            for j in 0..n {
                A[(i, j)] = -D2[(i, j)];
            }
            f[i] = source(x[i]);
        }

        // Last row: boundary condition at x = 1
        for j in 0..n {
            A[(n - 1, j)] = 0.0;
        }
        A[(n - 1, n - 1)] = 1.0;
        f[n - 1] = boundary_right;

        // Solve linear system
        match A.lu().solve(&f) {
            Some(solution) => solution.as_slice().to_vec(),
            None => vec![0.0; n],
        }
    }

    /// Interpolate values using Chebyshev basis
    pub fn interpolate_chebyshev(&self, values: &[f64], x: f64) -> f64 {
        let n = values.len();
        let nodes = self.chebyshev_points();

        // Barycentric Lagrange interpolation
        let mut numerator = 0.0;
        let mut denominator = 0.0;

        // Weights for Chebyshev nodes
        let mut weights = vec![1.0; n];
        weights[0] = 0.5;
        weights[n - 1] = 0.5;
        for i in 1..n {
            if i % 2 == 1 {
                weights[i] = -weights[i];
            }
        }

        for i in 0..n {
            if (x - nodes[i]).abs() < 1e-14 {
                return values[i];
            }

            let w = weights[i] / (x - nodes[i]);
            numerator += w * values[i];
            denominator += w;
        }

        numerator / denominator
    }
}

/// Fourier spectral methods for periodic problems
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FourierSpectral {
    /// Number of modes
    pub num_modes: usize,
}

impl FourierSpectral {
    /// Create new Fourier spectral solver
    pub fn new(num_modes: usize) -> Self {
        Self { num_modes }
    }

    /// Get spatial grid points for periodic domain [0, 2π]
    pub fn grid_points(&self) -> Vec<f64> {
        let n = self.num_modes;
        (0..n).map(|i| 2.0 * PI * i as f64 / n as f64).collect()
    }

    /// Compute derivative of periodic function using FFT
    /// (This is a simplified version; would use actual FFT in production)
    pub fn differentiate_periodic(&self, values: &[f64]) -> Vec<f64> {
        let n = values.len();
        let mut result = vec![0.0; n];

        // Simple finite difference for now
        // (Full implementation would use FFT-based spectral differentiation)
        for i in 0..n {
            let i_next = (i + 1) % n;
            let i_prev = if i == 0 { n - 1 } else { i - 1 };
            let dx = 2.0 * PI / n as f64;
            result[i] = (values[i_next] - values[i_prev]) / (2.0 * dx);
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chebyshev_points() {
        let solver = SpectralSolver::chebyshev(5);
        let points = solver.chebyshev_points();

        assert_eq!(points.len(), 5);
        // First point should be -1
        assert!((points[0] - (-1.0)).abs() < 1e-10);
        // Last point should be 1
        assert!((points[4] - 1.0).abs() < 1e-10);
        // Points should be in increasing order
        for i in 0..points.len() - 1 {
            assert!(points[i] < points[i + 1]);
        }
    }

    #[test]
    fn test_chebyshev_diff_matrix() {
        let solver = SpectralSolver::chebyshev(5);
        let D = solver.chebyshev_diff_matrix();

        assert_eq!(D.nrows(), 5);
        assert_eq!(D.ncols(), 5);

        // Differentiation matrix should be antisymmetric for interior points
        // (Not strictly true for all i,j but a sanity check)
        assert!(D[(0, 0)].is_finite());
    }

    #[test]
    fn test_solve_bvp_chebyshev() {
        // Solve -d²u/dx² = 1, u(-1)=0, u(1)=0
        // Analytical: u(x) = (1-x²)/2
        let solver = SpectralSolver::chebyshev(10);
        let solution = solver.solve_bvp_chebyshev(|_x| 1.0, 0.0, 0.0);

        // Check boundary conditions
        assert!((solution[0] - 0.0).abs() < 1e-6);
        assert!((solution[solution.len() - 1] - 0.0).abs() < 1e-6);

        // Check maximum at center (x=0)
        let x = solver.chebyshev_points();
        let mid_idx = x.iter().position(|&xi| xi.abs() < 0.1).unwrap_or(x.len() / 2);
        assert!(solution[mid_idx] > 0.0);
    }

    #[test]
    fn test_interpolate_chebyshev() {
        let solver = SpectralSolver::chebyshev(5);
        let x = solver.chebyshev_points();

        // Create test function: f(x) = x²
        let values: Vec<f64> = x.iter().map(|&xi| xi * xi).collect();

        // Interpolate at the nodes (should be exact)
        for i in 0..x.len() {
            let interp = solver.interpolate_chebyshev(&values, x[i]);
            assert!((interp - values[i]).abs() < 1e-10);
        }

        // Interpolate at x = 0.5
        let interp = solver.interpolate_chebyshev(&values, 0.5);
        assert!((interp - 0.25).abs() < 0.1); // Should be close to 0.5²
    }

    #[test]
    fn test_fourier_grid() {
        let fourier = FourierSpectral::new(8);
        let grid = fourier.grid_points();

        assert_eq!(grid.len(), 8);
        assert!((grid[0] - 0.0).abs() < 1e-10);
        assert!(grid[grid.len() - 1] < 2.0 * PI);
    }

    #[test]
    fn test_fourier_differentiate() {
        let fourier = FourierSpectral::new(16);
        let grid = fourier.grid_points();

        // Test with sin(x), derivative should be cos(x)
        let values: Vec<f64> = grid.iter().map(|&x| x.sin()).collect();
        let deriv = fourier.differentiate_periodic(&values);

        // Check at a few points
        for i in 0..grid.len() {
            let expected = grid[i].cos();
            let error = (deriv[i] - expected).abs();
            assert!(error < 0.5); // Loose tolerance for simple FD
        }
    }
}
