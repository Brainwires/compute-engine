// Unit tests for mathematics::special_functions::error
use computational_engine::compute::special_functions::error::*;

use super::*;

    // ==================== Error Function Tests ====================

    #[test]
    fn test_erf_zero() {
        let result = error_function(ErrorFunctionRequest {
            x: 0.0,
            function: "erf".to_string(),
        })
        .unwrap();
        // erf(0) = 0
        assert!(result.value.abs() < 1e-6);
    }

    #[test]
    fn test_erf_one() {
        let result = error_function(ErrorFunctionRequest {
            x: 1.0,
            function: "erf".to_string(),
        })
        .unwrap();
        // erf(1) ≈ 0.8427
        assert!((result.value - 0.8427).abs() < 0.001);
    }

    #[test]
    fn test_erf_negative() {
        let result = error_function(ErrorFunctionRequest {
            x: -1.0,
            function: "erf".to_string(),
        })
        .unwrap();
        // erf(-1) ≈ -0.8427 (odd function)
        assert!((result.value - (-0.8427)).abs() < 0.001);
    }

    #[test]
    fn test_erf_two() {
        let result = error_function(ErrorFunctionRequest {
            x: 2.0,
            function: "erf".to_string(),
        })
        .unwrap();
        // erf(2) ≈ 0.9953
        assert!((result.value - 0.9953).abs() < 0.001);
    }

    #[test]
    fn test_erf_large() {
        let result = error_function(ErrorFunctionRequest {
            x: 3.0,
            function: "erf".to_string(),
        })
        .unwrap();
        // erf(3) ≈ 0.99998 (approaches 1)
        assert!((result.value - 0.99998).abs() < 0.0001);
    }

    #[test]
    fn test_erfc_zero() {
        let result = error_function(ErrorFunctionRequest {
            x: 0.0,
            function: "erfc".to_string(),
        })
        .unwrap();
        // erfc(0) = 1 - erf(0) = 1
        assert!((result.value - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_erfc_one() {
        let result = error_function(ErrorFunctionRequest {
            x: 1.0,
            function: "erfc".to_string(),
        })
        .unwrap();
        // erfc(1) ≈ 0.1573
        assert!((result.value - 0.1573).abs() < 0.001);
    }

    #[test]
    fn test_erf_erfc_complement() {
        let x = 1.5;
        let erf_result = error_function(ErrorFunctionRequest {
            x,
            function: "erf".to_string(),
        })
        .unwrap();
        let erfc_result = error_function(ErrorFunctionRequest {
            x,
            function: "erfc".to_string(),
        })
        .unwrap();
        // erf(x) + erfc(x) = 1
        assert!((erf_result.value + erfc_result.value - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_erfcx_one() {
        let result = error_function(ErrorFunctionRequest {
            x: 1.0,
            function: "erfcx".to_string(),
        })
        .unwrap();
        // erfcx(1) = exp(1) * erfc(1) ≈ 0.4276
        assert!((result.value - 0.4276).abs() < 0.01);
    }

    #[test]
    fn test_erfi_zero() {
        let result = error_function(ErrorFunctionRequest {
            x: 0.0,
            function: "erfi".to_string(),
        })
        .unwrap();
        // erfi(0) = 0
        assert!(result.value.abs() < 1e-6);
    }

    // ==================== Additional Comprehensive Tests ====================

    #[test]
    fn test_erf_symmetry() {
        // erf(-x) = -erf(x)
        let x = 1.5;
        let erf_pos = error_function(ErrorFunctionRequest {
            x,
            function: "erf".to_string(),
        })
        .unwrap();
        let erf_neg = error_function(ErrorFunctionRequest {
            x: -x,
            function: "erf".to_string(),
        })
        .unwrap();
        assert!((erf_pos.value + erf_neg.value).abs() < 1e-6);
    }

    #[test]
    fn test_error_function_limit() {
        // lim_{x→∞} erf(x) = 1
        let result = error_function(ErrorFunctionRequest {
            x: 5.0,
            function: "erf".to_string(),
        })
        .unwrap();
        assert!((result.value - 1.0).abs() < 0.001);
    }
