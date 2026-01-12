// Unit tests for physics::control_systems::mod
use computational_engine::compute::physics::control_systems::mod::*;

use super::*;

    #[test]
    fn test_transfer_function_dc_gain() {
        let result = transfer_function(TransferFunctionRequest {
            numerator: vec![1.0],
            denominator: vec![1.0, 2.0],
            operation: "evaluate".to_string(),
            frequency: None,
            second_tf: None,
        })
        .unwrap();

        assert!((result.gain - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_routh_hurwitz_stable() {
        let result = routh_hurwitz(RouthHurwitzRequest {
            characteristic_polynomial: vec![1.0, 3.0, 3.0, 1.0], // (s+1)³
        })
        .unwrap();

        assert!(result.stable);
        assert_eq!(result.sign_changes, 0);
    }
