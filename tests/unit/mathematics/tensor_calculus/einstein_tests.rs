// Unit tests for mathematics::tensor_calculus::einstein
use computational_engine::mathematics::tensor_calculus::einstein::*;

use super::*;

    #[test]
    fn test_schwarzschild_solution() {
        let coords = vec![
            "t".to_string(),
            "r".to_string(),
            "theta".to_string(),
            "phi".to_string(),
        ];
        let boundary_conditions = vec![];

        let solutions = solve_spherically_symmetric_vacuum(&coords, &boundary_conditions).unwrap();

        assert!(!solutions.is_empty());
        assert_eq!(solutions[0].solution_type, "exact");
        assert_eq!(solutions[0].coordinates, coords);
    }

    #[test]
    fn test_flrw_solution() {
        let coords = vec![
            "t".to_string(),
            "r".to_string(),
            "theta".to_string(),
            "phi".to_string(),
        ];
        let boundary_conditions = vec![];

        let solutions = solve_flrw_universe(&coords, &boundary_conditions).unwrap();

        assert!(!solutions.is_empty());
        assert_eq!(solutions[0].solution_type, "exact");
        assert!(solutions[0].physical_parameters.contains_key("H"));
    }
