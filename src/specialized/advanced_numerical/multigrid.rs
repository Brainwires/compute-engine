//! Multigrid Methods
//!
//! Efficient iterative solvers for large linear systems arising from PDEs.

use serde::{Deserialize, Serialize};

/// Multigrid solver for 2D Poisson equation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultigridSolver {
    /// Number of V-cycles
    pub num_cycles: usize,
    /// Number of pre-smoothing iterations
    pub pre_smooth: usize,
    /// Number of post-smoothing iterations
    pub post_smooth: usize,
    /// Coarsest grid size
    pub coarsest_size: usize,
}

impl MultigridSolver {
    /// Create a new multigrid solver with default parameters
    pub fn new() -> Self {
        Self {
            num_cycles: 5,
            pre_smooth: 2,
            post_smooth: 2,
            coarsest_size: 4,
        }
    }

    /// Create with custom parameters
    pub fn with_params(
        num_cycles: usize,
        pre_smooth: usize,
        post_smooth: usize,
        coarsest_size: usize,
    ) -> Self {
        Self {
            num_cycles,
            pre_smooth,
            post_smooth,
            coarsest_size,
        }
    }

    /// Solve 2D Poisson equation: -∇²u = f
    /// Using V-cycle multigrid with Gauss-Seidel smoothing
    pub fn solve_poisson_2d(&self, f: &[Vec<f64>], initial_guess: Option<Vec<Vec<f64>>>) -> Vec<Vec<f64>> {
        let n = f.len();
        let m = f[0].len();

        let mut u = initial_guess.unwrap_or_else(|| vec![vec![0.0; m]; n]);

        // Perform V-cycles
        for _ in 0..self.num_cycles {
            u = self.v_cycle(&u, f, 0);
        }

        u
    }

    /// Single V-cycle
    fn v_cycle(&self, u: &[Vec<f64>], f: &[Vec<f64>], level: usize) -> Vec<Vec<f64>> {
        let n = u.len();
        let m = u[0].len();

        // Base case: solve directly on coarsest grid
        if n <= self.coarsest_size || m <= self.coarsest_size {
            return self.solve_direct(u, f);
        }

        // Pre-smoothing
        let mut u_smooth = u.to_vec();
        for _ in 0..self.pre_smooth {
            u_smooth = self.gauss_seidel_step(&u_smooth, f);
        }

        // Compute residual: r = f - A*u
        let residual = self.compute_residual(&u_smooth, f);

        // Restrict residual to coarser grid
        let f_coarse = self.restrict(&residual);

        // Solve on coarser grid
        let n_coarse = f_coarse.len();
        let m_coarse = f_coarse[0].len();
        let u_coarse_init = vec![vec![0.0; m_coarse]; n_coarse];
        let correction_coarse = self.v_cycle(&u_coarse_init, &f_coarse, level + 1);

        // Prolongate correction to fine grid
        let correction = self.prolongate(&correction_coarse, n, m);

        // Apply correction
        let mut u_corrected = vec![vec![0.0; m]; n];
        for i in 0..n {
            for j in 0..m {
                u_corrected[i][j] = u_smooth[i][j] + correction[i][j];
            }
        }

        // Post-smoothing
        for _ in 0..self.post_smooth {
            u_corrected = self.gauss_seidel_step(&u_corrected, f);
        }

        u_corrected
    }

    /// Gauss-Seidel smoothing iteration
    fn gauss_seidel_step(&self, u: &[Vec<f64>], f: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n = u.len();
        let m = u[0].len();
        let h = 1.0 / (n - 1) as f64;
        let h2 = h * h;

        let mut u_new = u.to_vec();

        // Gauss-Seidel iteration
        for i in 1..n - 1 {
            for j in 1..m - 1 {
                // For -∇²u = f, the discrete equation is:
                // (4*u[i,j] - u[i-1,j] - u[i+1,j] - u[i,j-1] - u[i,j+1])/h² = -f[i,j]
                // Solving for u[i,j]:
                u_new[i][j] = 0.25
                    * (u_new[i - 1][j] + u_new[i + 1][j] + u_new[i][j - 1] + u_new[i][j + 1]
                        + h2 * f[i][j]);
            }
        }

        u_new
    }

    /// Compute residual: r = f - A*u
    fn compute_residual(&self, u: &[Vec<f64>], f: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n = u.len();
        let m = u[0].len();
        let h = 1.0 / (n - 1) as f64;
        let h2 = h * h;

        let mut residual = vec![vec![0.0; m]; n];

        for i in 1..n - 1 {
            for j in 1..m - 1 {
                // For -∇²u = f, the discrete form is:
                // -(u[i-1,j] + u[i+1,j] + u[i,j-1] + u[i,j+1] - 4*u[i,j])/h² = f
                // So A*u = (u[i-1,j] + u[i+1,j] + u[i,j-1] + u[i,j+1] - 4*u[i,j])/h²
                let laplacian = (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - 4.0 * u[i][j]) / h2;
                residual[i][j] = f[i][j] + laplacian;
            }
        }

        residual
    }

    /// Restriction operator (fine to coarse): full weighting
    fn restrict(&self, fine: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n_fine = fine.len();
        let m_fine = fine[0].len();
        let n_coarse = (n_fine + 1) / 2;
        let m_coarse = (m_fine + 1) / 2;

        let mut coarse = vec![vec![0.0; m_coarse]; n_coarse];

        for i in 0..n_coarse {
            for j in 0..m_coarse {
                let i_fine = 2 * i.min(n_fine / 2);
                let j_fine = 2 * j.min(m_fine / 2);

                // Full weighting: weighted average of 9 points
                if i_fine > 0 && i_fine < n_fine - 1 && j_fine > 0 && j_fine < m_fine - 1 {
                    coarse[i][j] = 0.25 * fine[i_fine][j_fine]
                        + 0.125
                            * (fine[i_fine - 1][j_fine]
                                + fine[i_fine + 1][j_fine]
                                + fine[i_fine][j_fine - 1]
                                + fine[i_fine][j_fine + 1])
                        + 0.0625
                            * (fine[i_fine - 1][j_fine - 1]
                                + fine[i_fine - 1][j_fine + 1]
                                + fine[i_fine + 1][j_fine - 1]
                                + fine[i_fine + 1][j_fine + 1]);
                } else {
                    // Boundary: simple injection
                    coarse[i][j] = fine[i_fine][j_fine];
                }
            }
        }

        coarse
    }

    /// Prolongation operator (coarse to fine): bilinear interpolation
    fn prolongate(&self, coarse: &[Vec<f64>], n_fine: usize, m_fine: usize) -> Vec<Vec<f64>> {
        let n_coarse = coarse.len();
        let m_coarse = coarse[0].len();
        let mut fine = vec![vec![0.0; m_fine]; n_fine];

        for i in 0..n_fine {
            for j in 0..m_fine {
                let i_coarse = i / 2;
                let j_coarse = j / 2;

                if i % 2 == 0 && j % 2 == 0 {
                    // Coincident points
                    fine[i][j] = coarse[i_coarse.min(n_coarse - 1)][j_coarse.min(m_coarse - 1)];
                } else if i % 2 == 0 {
                    // Horizontal edge
                    let j1 = j_coarse.min(m_coarse - 1);
                    let j2 = (j_coarse + 1).min(m_coarse - 1);
                    fine[i][j] = 0.5
                        * (coarse[i_coarse.min(n_coarse - 1)][j1]
                            + coarse[i_coarse.min(n_coarse - 1)][j2]);
                } else if j % 2 == 0 {
                    // Vertical edge
                    let i1 = i_coarse.min(n_coarse - 1);
                    let i2 = (i_coarse + 1).min(n_coarse - 1);
                    fine[i][j] = 0.5
                        * (coarse[i1][j_coarse.min(m_coarse - 1)]
                            + coarse[i2][j_coarse.min(m_coarse - 1)]);
                } else {
                    // Interior point: bilinear
                    let i1 = i_coarse.min(n_coarse - 1);
                    let i2 = (i_coarse + 1).min(n_coarse - 1);
                    let j1 = j_coarse.min(m_coarse - 1);
                    let j2 = (j_coarse + 1).min(m_coarse - 1);
                    fine[i][j] =
                        0.25 * (coarse[i1][j1] + coarse[i1][j2] + coarse[i2][j1] + coarse[i2][j2]);
                }
            }
        }

        fine
    }

    /// Direct solver for small grids
    fn solve_direct(&self, _u: &[Vec<f64>], f: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n = f.len();
        let m = f[0].len();

        // Use many Gauss-Seidel iterations for small grids
        let mut u = vec![vec![0.0; m]; n];
        for _ in 0..100 {
            u = self.gauss_seidel_step(&u, f);
        }

        u
    }
}

impl Default for MultigridSolver {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multigrid_creation() {
        let solver = MultigridSolver::new();
        assert_eq!(solver.num_cycles, 5);
        assert_eq!(solver.pre_smooth, 2);
        assert_eq!(solver.post_smooth, 2);
    }

    #[test]
    fn test_multigrid_custom_params() {
        let solver = MultigridSolver::with_params(10, 3, 3, 8);
        assert_eq!(solver.num_cycles, 10);
        assert_eq!(solver.pre_smooth, 3);
    }

    #[test]
    fn test_multigrid_solve_simple() {
        // Solve -∇²u = 1 with u=0 on boundary
        let n = 17; // 2^4 + 1 for good multigrid performance
        let f = vec![vec![1.0; n]; n];

        let mut solver = MultigridSolver::new();
        solver.num_cycles = 20; // Increased cycles for better convergence

        let solution = solver.solve_poisson_2d(&f, None);

        // Check boundary conditions
        for i in 0..n {
            assert!((solution[0][i] - 0.0).abs() < 1e-2);
            assert!((solution[n - 1][i] - 0.0).abs() < 1e-2);
            assert!((solution[i][0] - 0.0).abs() < 1e-2);
            assert!((solution[i][n - 1] - 0.0).abs() < 1e-2);
        }

        // Check that interior is positive (with small tolerance near boundary)
        for i in 2..n - 2 {
            for j in 2..n - 2 {
                assert!(solution[i][j] > 0.0,
                    "Solution at ({}, {}) should be positive, got {}", i, j, solution[i][j]);
            }
        }
    }

    #[test]
    fn test_gauss_seidel_step() {
        let solver = MultigridSolver::new();
        let n = 5;
        let u = vec![vec![0.0; n]; n];
        let f = vec![vec![1.0; n]; n];

        let u_new = solver.gauss_seidel_step(&u, &f);

        // Interior values should change
        assert!(u_new[2][2].abs() > 1e-10);
    }

    #[test]
    fn test_restrict_prolongate() {
        let solver = MultigridSolver::new();

        // Create fine grid
        let n_fine = 9;
        let fine = vec![vec![1.0; n_fine]; n_fine];

        // Restrict
        let coarse = solver.restrict(&fine);
        assert_eq!(coarse.len(), 5);

        // Prolongate back
        let fine_again = solver.prolongate(&coarse, n_fine, n_fine);
        assert_eq!(fine_again.len(), n_fine);
        assert_eq!(fine_again[0].len(), n_fine);
    }

    #[test]
    fn test_compute_residual() {
        let solver = MultigridSolver::new();
        let n = 5;

        // If u is the exact solution to -∇²u = 0, residual should be small
        let u = vec![vec![0.0; n]; n];
        let f = vec![vec![0.0; n]; n];

        let residual = solver.compute_residual(&u, &f);

        // All residuals should be zero (or very small)
        for i in 1..n - 1 {
            for j in 1..n - 1 {
                assert!(residual[i][j].abs() < 1e-10);
            }
        }
    }
}
