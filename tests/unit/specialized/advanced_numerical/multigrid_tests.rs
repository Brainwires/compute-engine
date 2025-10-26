// Unit tests for specialized::advanced_numerical::multigrid
use computational_engine::specialized::advanced_numerical::multigrid::*;

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
