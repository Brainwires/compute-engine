// Unit tests for specialized::control_theory::analysis
use computational_engine::compute::physics::control_systems::analysis::*;

use super::*;

    #[test]
    fn test_controllability() {
        let a = vec![vec![0.0, 1.0], vec![-1.0, 0.0]];
        let b = vec![vec![0.0], vec![1.0]];

        let result = controllability(&a, &b);
        assert!(result.rank > 0);
    }

    #[test]
    fn test_observability() {
        let a = vec![vec![0.0, 1.0], vec![-1.0, 0.0]];
        let c = vec![vec![1.0, 0.0]];

        let result = observability(&a, &c);
        assert!(result.rank > 0);
    }
