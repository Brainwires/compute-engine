// Unit tests for specialized::control_theory::lqr
use computational_engine::specialized::control_theory::lqr::*;

use super::*;

    #[test]
    fn test_lqr_dimensions() {
        let a = vec![vec![0.0, 1.0], vec![-1.0, 0.0]];
        let b = vec![vec![0.0], vec![1.0]];
        let q = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let r = vec![vec![1.0]];

        let result = lqr(&a, &b, &q, &r);
        assert!(result.is_ok());
    }
