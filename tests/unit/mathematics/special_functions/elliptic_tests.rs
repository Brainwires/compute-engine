// Unit tests for mathematics::special_functions::elliptic
use computational_engine::mathematics::special_functions::elliptic::*;

use super::*;

    #[test]
    fn test_elliptic_k_zero() {
        let result = elliptic_integral(EllipticIntegralRequest {
            integral_type: "K".to_string(),
            k: 0.0,
            phi: None,
            n: None,
        }).unwrap();
        // K(0) = π/2
        assert!((result.value - PI / 2.0).abs() < 0.01);
    }

    #[test]
    fn test_elliptic_k_half() {
        let result = elliptic_integral(EllipticIntegralRequest {
            integral_type: "K".to_string(),
            k: 0.5,
            phi: None,
            n: None,
        }).unwrap();
        // K(0.5) ≈ 1.685750
        assert!((result.value - 1.685750).abs() < 0.01);
    }

    #[test]
    fn test_elliptic_e_zero() {
        let result = elliptic_integral(EllipticIntegralRequest {
            integral_type: "E".to_string(),
            k: 0.0,
            phi: None,
            n: None,
        }).unwrap();
        // E(0) = π/2
        assert!((result.value - PI / 2.0).abs() < 0.01);
    }

    #[test]
    fn test_elliptic_k_invalid() {
        let result = elliptic_integral(EllipticIntegralRequest {
            integral_type: "K".to_string(),
            k: 1.5,
            phi: None,
            n: None,
        });
        assert!(result.is_err());
    }

    #[test]
    fn test_elliptic_f_incomplete() {
        let result = elliptic_integral(EllipticIntegralRequest {
            integral_type: "F".to_string(),
            k: 0.5,
            phi: Some(PI / 4.0),
            n: None,
        }).unwrap();
        // F(π/4, 0.5) should be a finite value
        assert!(result.value > 0.0 && result.value < 10.0);
    }

    #[test]
    fn test_elliptic_k_e_inequality() {
        // E(k) <= K(k) for 0 <= k < 1
        let k = 0.7;
        let k_result = elliptic_integral(EllipticIntegralRequest {
            integral_type: "K".to_string(),
            k,
            phi: None,
            n: None,
        }).unwrap();
        let e_result = elliptic_integral(EllipticIntegralRequest {
            integral_type: "E".to_string(),
            k,
            phi: None,
            n: None,
        }).unwrap();
        assert!(e_result.value <= k_result.value);
    }
