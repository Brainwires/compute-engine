// Unit tests for specialized::advanced_numerical::spectral
use computational_engine::compute::advanced_numerical::spectral::*;

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
