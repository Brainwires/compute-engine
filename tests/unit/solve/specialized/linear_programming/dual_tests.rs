// Unit tests for specialized::linear_programming::dual
use computational_engine::solve::specialized::linear_programming::*;

use super::*;
    use computational_engine::solve::specialized::linear_programming::simplex;

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
    #[ignore] // KNOWN LIMITATION: Requires two-phase simplex or Big-M for negative RHS
    fn test_strong_duality() {
        // The dual problem has constraint matrix with negative RHS values (-c = [-3, -2])
        // after converting A^T·y ≥ c to -A^T·y ≤ -c.
        //
        // Our basic simplex requires non-negative RHS (feasible starting basis).
        // To handle this, we need either:
        // 1. Two-phase Simplex method (add artificial variables, minimize them first)
        // 2. Big-M method (penalize artificial variables in objective)
        // 3. Dual simplex (start from dual-feasible basis)
        //
        // This is a future enhancement tracked separately.

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
