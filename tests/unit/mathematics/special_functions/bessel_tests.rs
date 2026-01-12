// Unit tests for mathematics::special_functions::bessel
use computational_engine::compute::special_functions::bessel::*;

use super::*;

    // ==================== Bessel Functions Tests ====================

    #[test]
    fn test_bessel_j0_zero() {
        let result = bessel_function(BesselRequest {
            function_type: "J".to_string(),
            order: 0.0,
            x: 0.0,
        })
        .unwrap();
        // J_0(0) = 1
        assert!((result.value - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_bessel_j0_first_zero() {
        let result = bessel_function(BesselRequest {
            function_type: "J".to_string(),
            order: 0.0,
            x: 2.4048, // First zero of J_0
        })
        .unwrap();
        assert!(result.value.abs() < 1e-2);
    }

    #[test]
    fn test_bessel_j0_pi() {
        let result = bessel_function(BesselRequest {
            function_type: "J".to_string(),
            order: 0.0,
            x: PI,
        })
        .unwrap();
        // J_0(π) ≈ -0.304242
        assert!((result.value - (-0.304242)).abs() < 1e-3);
    }

    #[test]
    fn test_bessel_j1_zero() {
        let result = bessel_function(BesselRequest {
            function_type: "J".to_string(),
            order: 1.0,
            x: 0.0,
        })
        .unwrap();
        // J_1(0) = 0
        assert!(result.value.abs() < 1e-6);
    }

    #[test]
    fn test_bessel_j1_one() {
        let result = bessel_function(BesselRequest {
            function_type: "J".to_string(),
            order: 1.0,
            x: 1.0,
        })
        .unwrap();
        // J_1(1) ≈ 0.4400506
        assert!((result.value - 0.4400506).abs() < 1e-3);
    }

    #[test]
    fn test_bessel_j2_two() {
        let result = bessel_function(BesselRequest {
            function_type: "J".to_string(),
            order: 2.0,
            x: 2.0,
        })
        .unwrap();
        // J_2(2) ≈ 0.3528340
        assert!((result.value - 0.3528340).abs() < 1e-3);
    }

    #[test]
    fn test_bessel_i0_zero() {
        let result = bessel_function(BesselRequest {
            function_type: "I".to_string(),
            order: 0.0,
            x: 0.0,
        })
        .unwrap();
        // I_0(0) = 1
        assert!((result.value - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_bessel_i0_one() {
        let result = bessel_function(BesselRequest {
            function_type: "I".to_string(),
            order: 0.0,
            x: 1.0,
        })
        .unwrap();
        // I_0(1) ≈ 1.266066
        assert!((result.value - 1.266066).abs() < 1e-3);
    }

    #[test]
    fn test_bessel_i1_zero() {
        let result = bessel_function(BesselRequest {
            function_type: "I".to_string(),
            order: 1.0,
            x: 0.0,
        })
        .unwrap();
        // I_1(0) = 0
        assert!(result.value.abs() < 1e-6);
    }

    #[test]
    fn test_bessel_i1_one() {
        let result = bessel_function(BesselRequest {
            function_type: "I".to_string(),
            order: 1.0,
            x: 1.0,
        })
        .unwrap();
        // I_1(1) ≈ 0.565159
        assert!((result.value - 0.565159).abs() < 1e-3);
    }

    #[test]
    fn test_bessel_y0_positive() {
        let result = bessel_function(BesselRequest {
            function_type: "Y".to_string(),
            order: 0.0,
            x: 5.0,
        })
        .unwrap();
        // Y_0(5) ≈ -0.308518
        assert!((result.value - (-0.308518)).abs() < 0.1);
    }

    #[test]
    fn test_bessel_y_negative_x_error() {
        let result = bessel_function(BesselRequest {
            function_type: "Y".to_string(),
            order: 0.0,
            x: -1.0,
        });
        assert!(result.is_err());
    }

    #[test]
    fn test_bessel_j_recurrence() {
        // J_{n-1}(x) + J_{n+1}(x) = (2n/x)J_n(x)
        let x = 2.0;
        let n = 2.0;
        let j1 = bessel_function(BesselRequest {
            function_type: "J".to_string(),
            order: n - 1.0,
            x,
        })
        .unwrap();
        let j2 = bessel_function(BesselRequest {
            function_type: "J".to_string(),
            order: n,
            x,
        })
        .unwrap();
        let j3 = bessel_function(BesselRequest {
            function_type: "J".to_string(),
            order: n + 1.0,
            x,
        })
        .unwrap();
        let left = j1.value + j3.value;
        let right = (2.0 * n / x) * j2.value;
        assert!((left - right).abs() < 0.01);
    }

    #[test]
    fn test_bessel_i_positivity() {
        // I_n(x) > 0 for x > 0
        for n in 0..4 {
            let result = bessel_function(BesselRequest {
                function_type: "I".to_string(),
                order: n as f64,
                x: 1.0,
            })
            .unwrap();
            assert!(result.value > 0.0);
        }
    }
