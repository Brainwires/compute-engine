// Unit tests for specialized::linear_programming::dual
use computational_engine::specialized::linear_programming::dual::*;

use super::*;
    use crate::specialized::linear_programming::simplex;

    #[test]
    fn test_dual_construction() {
        // Primal: max 3x + 2y subject to x + y ≤ 4, x + 3y ≤ 6
        let primal = LinearProgram::new(vec![3.0, 2.0], true).with_inequality_constraints(
            vec![vec![1.0, 1.0], vec![1.0, 3.0]],
            vec![4.0, 6.0],
        );

        let dual = compute_dual(&primal);

        // Dual should have 2 variables (one for each primal constraint)
        assert_eq!(dual.num_variables(), 2);

        // Dual objective should be [-4, -6] (negated primal RHS for maximization form)
        assert_eq!(dual.objective, vec![-4.0, -6.0]);

        // Dual should be maximization (we converted from minimization)
        assert!(dual.maximize);
    }

    #[test]
    #[ignore] // TODO: Requires two-phase simplex or Big-M method for >= constraints
    fn test_strong_duality() {
        // The dual problem involves >= constraints which our basic Simplex
        // implementation doesn't handle yet. This requires either:
        // 1. Two-phase Simplex method
        // 2. Big-M method
        // 3. Revised Simplex with artificial variables
        //
        // For now, we skip this test and mark it as a future enhancement.

        let primal = LinearProgram::new(vec![3.0, 2.0], true).with_inequality_constraints(
            vec![vec![1.0, 1.0], vec![1.0, 3.0]],
            vec![4.0, 6.0],
        );

        let primal_solution = simplex::simplex(&primal).unwrap();
        assert_eq!(primal_solution.status, SolutionStatus::Optimal);

        let dual = compute_dual(&primal);
        let dual_solution = simplex::simplex(&dual).unwrap();

        // Strong duality: primal and dual optimal values should be equal (within tolerance)
        assert!(
            (primal_solution.objective_value - dual_solution.objective_value).abs() < 0.1,
            "Primal: {}, Dual: {}",
            primal_solution.objective_value,
            dual_solution.objective_value
        );
    }
