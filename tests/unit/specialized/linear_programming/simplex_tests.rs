// Unit tests for specialized::linear_programming::simplex
use computational_engine::specialized::linear_programming::*;

use super::*;

    #[test]
    fn test_simple_maximization() {
        // Maximize: 3x + 2y
        // Subject to: x + y ≤ 4, x + 3y ≤ 6, x ≥ 0, y ≥ 0
        let lp = LinearProgram::new(vec![3.0, 2.0], true).with_inequality_constraints(
            vec![vec![1.0, 1.0], vec![1.0, 3.0]],
            vec![4.0, 6.0],
        );

        let solution = simplex(&lp).unwrap();

        assert_eq!(solution.status, SolutionStatus::Optimal);
        // Optimal at (4, 0): 3*4 + 2*0 = 12
        assert!(solution.objective_value > 11.9 && solution.objective_value < 12.1,
                "Expected ~12, got {}", solution.objective_value);
        assert!(is_feasible(&lp, &solution.solution));
    }

    #[test]
    fn test_minimization() {
        // Minimize: -x - y  (equivalent to maximizing x + y but using minimize flag)
        // Subject to: x + y ≤ 4, x + 3y ≤ 6, x ≥ 0, y ≥ 0
        // This should find the minimum of -x - y, which is at (0, 0) with value 0
        // But actually we want to test minimization finds the right answer
        // Better test: Minimize x + y subject to x + y ≥ 4, which in ≤ form is infeasible
        // Or: minimize with constraints that force us away from origin

        // Simpler: Minimize x subject to x ≥ 2, y ≥ 0
        // In standard form: minimize x subject to -x ≤ -2, which needs two-phase
        // Even simpler: minimize -x subject to x ≤ 4, should give x=4, value=-4
        let lp = LinearProgram::new(vec![-1.0, 0.0], false).with_inequality_constraints(
            vec![vec![1.0, 0.0]],
            vec![4.0],
        );

        let solution = simplex(&lp).unwrap();

        assert_eq!(solution.status, SolutionStatus::Optimal);
        // Minimize -x with x ≤ 4: optimal at x=4, objective = -4
        assert!(solution.objective_value > -4.1 && solution.objective_value < -3.9,
                "Expected ~-4, got {}", solution.objective_value);
        assert!(is_feasible(&lp, &solution.solution));
    }

    #[test]
    fn test_unbounded() {
        // Maximize: x + y
        // Subject to: -x + y ≤ 1, x ≥ 0, y ≥ 0
        let lp = LinearProgram::new(vec![1.0, 1.0], true)
            .with_inequality_constraints(vec![vec![-1.0, 1.0]], vec![1.0]);

        let solution = simplex(&lp).unwrap();

        assert_eq!(solution.status, SolutionStatus::Unbounded);
    }

    #[test]
    fn test_feasibility() {
        let lp = LinearProgram::new(vec![3.0, 2.0], true).with_inequality_constraints(
            vec![vec![1.0, 1.0], vec![1.0, 3.0]],
            vec![4.0, 6.0],
        );

        assert!(is_feasible(&lp, &vec![3.0, 1.0]));
        assert!(is_feasible(&lp, &vec![0.0, 0.0]));
        assert!(!is_feasible(&lp, &vec![5.0, 5.0])); // Violates constraints
        assert!(!is_feasible(&lp, &vec![-1.0, 1.0])); // Negative variable
    }

    #[test]
    fn test_objective_value() {
        let lp = LinearProgram::new(vec![3.0, 2.0], true);

        assert_eq!(objective_value(&lp, &vec![3.0, 1.0]), 11.0);
        assert_eq!(objective_value(&lp, &vec![0.0, 0.0]), 0.0);
        assert_eq!(objective_value(&lp, &vec![1.0, 2.0]), 7.0);
    }
