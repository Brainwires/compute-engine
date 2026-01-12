// Unit tests for specialized::advanced_numerical::finite_element
use computational_engine::compute::advanced_numerical::finite_element::*;

use super::*;

    #[test]
    fn test_fem_1d_constant_source() {
        // Solve -d²u/dx² = 1, u(0)=0, u(1)=0
        // Analytical solution: u(x) = x(1-x)/2
        let solver = FiniteElementSolver1D::new(
            10,
            (0.0, 1.0),
            (
                BoundaryCondition::Dirichlet(0.0),
                BoundaryCondition::Dirichlet(0.0),
            ),
        );

        let solution = solver.solve(|_x| 1.0);

        // Check boundary conditions
        assert!((solution[0] - 0.0).abs() < 1e-10);
        assert!((solution[solution.len() - 1] - 0.0).abs() < 1e-10);

        // Check maximum at midpoint
        let mid_idx = solution.len() / 2;
        let max_value = solution[mid_idx];
        assert!(max_value > 0.0);
        assert!(max_value < 0.2); // Should be close to 0.125
    }

    #[test]
    fn test_fem_1d_linear_source() {
        // Solve -d²u/dx² = x, u(0)=0, u(1)=0
        let solver = FiniteElementSolver1D::new(
            20,
            (0.0, 1.0),
            (
                BoundaryCondition::Dirichlet(0.0),
                BoundaryCondition::Dirichlet(0.0),
            ),
        );

        let solution = solver.solve(|x| x);

        // Check boundary conditions
        assert!((solution[0] - 0.0).abs() < 1e-10);
        assert!((solution[solution.len() - 1] - 0.0).abs() < 1e-10);

        // Solution should be positive in interior
        for i in 1..solution.len() - 1 {
            assert!(solution[i] > 0.0);
        }
    }

    #[test]
    fn test_fem_2d_poisson() {
        // Solve -∇²u = 1 with u=0 on boundary
        let solver = FiniteElementSolver2D::new(10, 10, ((0.0, 1.0), (0.0, 1.0)));

        let solution = solver.solve_poisson(|_x, _y| 1.0);

        // Check boundary conditions
        let n = solution.len();
        let m = solution[0].len();

        // Boundaries should be zero
        for i in 0..m {
            assert!((solution[0][i] - 0.0).abs() < 1e-10);
            assert!((solution[n - 1][i] - 0.0).abs() < 1e-10);
        }
        for j in 0..n {
            assert!((solution[j][0] - 0.0).abs() < 1e-10);
            assert!((solution[j][m - 1] - 0.0).abs() < 1e-10);
        }

        // Interior should be positive
        for j in 1..n - 1 {
            for i in 1..m - 1 {
                assert!(solution[j][i] > 0.0);
            }
        }
    }

    #[test]
    fn test_fem_get_nodes() {
        let solver = FiniteElementSolver1D::new(
            4,
            (0.0, 1.0),
            (
                BoundaryCondition::Dirichlet(0.0),
                BoundaryCondition::Dirichlet(0.0),
            ),
        );

        let nodes = solver.get_nodes();
        assert_eq!(nodes.len(), 5);
        assert!((nodes[0] - 0.0).abs() < 1e-10);
        assert!((nodes[4] - 1.0).abs() < 1e-10);
        assert!((nodes[2] - 0.5).abs() < 1e-10);
    }
