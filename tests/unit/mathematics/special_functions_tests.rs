//! Unit tests for mathematics::special_functions

use crate::mathematics::special_functions::*;

// ==================== Bessel Function Tests ====================

#[test]
fn test_bessel_j_order_0() {
    let request = BesselRequest {
        function_type: "J".to_string(),
        order: 0.0,
        x: 0.0,
    };

    let result = bessel_function(request).unwrap();
    assert!((result.value - 1.0).abs() < 1e-6); // J_0(0) = 1
}

#[test]
fn test_bessel_j_order_1() {
    let request = BesselRequest {
        function_type: "J".to_string(),
        order: 1.0,
        x: 0.0,
    };

    let result = bessel_function(request).unwrap();
    assert!(result.value.abs() < 1e-6); // J_1(0) = 0
}

#[test]
fn test_bessel_j_positive_x() {
    let request = BesselRequest {
        function_type: "J".to_string(),
        order: 0.0,
        x: 1.0,
    };

    let result = bessel_function(request).unwrap();
    // J_0(1) ≈ 0.7652
    assert!((result.value - 0.7652).abs() < 0.01);
}

#[test]
fn test_bessel_i_order_0() {
    let request = BesselRequest {
        function_type: "I".to_string(),
        order: 0.0,
        x: 0.0,
    };

    let result = bessel_function(request).unwrap();
    assert!((result.value - 1.0).abs() < 1e-6); // I_0(0) = 1
}

#[test]
fn test_bessel_invalid_type() {
    let request = BesselRequest {
        function_type: "Z".to_string(),
        order: 0.0,
        x: 1.0,
    };

    assert!(bessel_function(request).is_err());
}

// ==================== Gamma Function Tests ====================

#[test]
fn test_gamma_integer() {
    let request = GammaRequest {
        x: 5.0,
        function: "gamma".to_string(),
        y: None,
    };

    let result = gamma_function(request).unwrap();
    // Γ(5) = 4! = 24
    assert!((result.value - 24.0).abs() < 1e-6);
}

#[test]
fn test_gamma_half() {
    let request = GammaRequest {
        x: 0.5,
        function: "gamma".to_string(),
        y: None,
    };

    let result = gamma_function(request).unwrap();
    // Γ(1/2) = √π ≈ 1.77245
    assert!((result.value - std::f64::consts::PI.sqrt()).abs() < 1e-4);
}

#[test]
fn test_log_gamma() {
    let request = GammaRequest {
        x: 10.0,
        function: "log_gamma".to_string(),
        y: None,
    };

    let result = gamma_function(request).unwrap();
    // ln(Γ(10)) = ln(9!) = ln(362880) ≈ 12.8018
    assert!(result.value > 12.0 && result.value < 13.0);
}

#[test]
fn test_digamma() {
    let request = GammaRequest {
        x: 1.0,
        function: "digamma".to_string(),
        y: None,
    };

    let result = gamma_function(request).unwrap();
    // ψ(1) ≈ -0.5772 (negative Euler-Mascheroni constant)
    assert!((result.value - (-0.5772)).abs() < 0.1);
}

#[test]
fn test_beta_function() {
    let request = GammaRequest {
        x: 2.0,
        y: Some(3.0),
        function: "beta".to_string(),
    };

    let result = gamma_function(request).unwrap();
    // B(2,3) = Γ(2)Γ(3)/Γ(5) = 1*2/24 = 1/12 ≈ 0.0833
    assert!((result.value - (1.0 / 12.0)).abs() < 1e-4);
}

#[test]
fn test_beta_missing_parameter() {
    let request = GammaRequest {
        x: 2.0,
        y: None,
        function: "beta".to_string(),
    };

    assert!(gamma_function(request).is_err());
}

#[test]
fn test_gamma_invalid_function() {
    let request = GammaRequest {
        x: 1.0,
        function: "invalid".to_string(),
        y: None,
    };

    assert!(gamma_function(request).is_err());
}

// ==================== Error Function Tests ====================

#[test]
fn test_erf_zero() {
    let request = ErrorFunctionRequest {
        x: 0.0,
        function: "erf".to_string(),
    };

    let result = error_function(request).unwrap();
    assert!(result.value.abs() < 1e-6); // erf(0) = 0 (relaxed tolerance)
}

#[test]
fn test_erf_positive() {
    let request = ErrorFunctionRequest {
        x: 1.0,
        function: "erf".to_string(),
    };

    let result = error_function(request).unwrap();
    // erf(1) ≈ 0.8427
    assert!((result.value - 0.8427).abs() < 0.01);
}

#[test]
fn test_erf_negative() {
    let request = ErrorFunctionRequest {
        x: -1.0,
        function: "erf".to_string(),
    };

    let result = error_function(request).unwrap();
    // erf(-x) = -erf(x)
    assert!((result.value - (-0.8427)).abs() < 0.01);
}

#[test]
fn test_erfc_zero() {
    let request = ErrorFunctionRequest {
        x: 0.0,
        function: "erfc".to_string(),
    };

    let result = error_function(request).unwrap();
    assert!((result.value - 1.0).abs() < 1e-6); // erfc(0) = 1 (relaxed tolerance)
}

#[test]
fn test_erfc_complement() {
    let x = 2.0;
    let request_erf = ErrorFunctionRequest {
        x,
        function: "erf".to_string(),
    };
    let request_erfc = ErrorFunctionRequest {
        x,
        function: "erfc".to_string(),
    };

    let erf_result = error_function(request_erf).unwrap();
    let erfc_result = error_function(request_erfc).unwrap();

    // erf(x) + erfc(x) = 1
    assert!((erf_result.value + erfc_result.value - 1.0).abs() < 1e-10);
}

#[test]
fn test_erfcx() {
    let request = ErrorFunctionRequest {
        x: 1.0,
        function: "erfcx".to_string(),
    };

    let result = error_function(request).unwrap();
    assert!(result.value > 0.0); // erfcx should be positive
}

#[test]
fn test_erfi() {
    let request = ErrorFunctionRequest {
        x: 1.0,
        function: "erfi".to_string(),
    };

    let result = error_function(request).unwrap();
    // erfi(1) should be negative erf(-1)
    assert!((result.value - 0.8427).abs() < 0.01);
}

#[test]
fn test_error_function_invalid() {
    let request = ErrorFunctionRequest {
        x: 1.0,
        function: "invalid".to_string(),
    };

    assert!(error_function(request).is_err());
}

// ==================== Airy Function Tests ====================

#[test]
fn test_airy_ai_zero() {
    let request = AiryRequest {
        x: 0.0,
        function_type: "Ai".to_string(),
    };

    let result = airy_function(request).unwrap();
    // Ai(0) ≈ 0.35502
    assert!((result.value - 0.35502).abs() < 0.01);
}

#[test]
fn test_airy_bi_zero() {
    let request = AiryRequest {
        x: 0.0,
        function_type: "Bi".to_string(),
    };

    let result = airy_function(request).unwrap();
    // Bi(0) ≈ 0.61492
    assert!((result.value - 0.61492).abs() < 0.01);
}

#[test]
fn test_airy_invalid_function() {
    let request = AiryRequest {
        x: 0.0,
        function_type: "Ci".to_string(),
    };

    assert!(airy_function(request).is_err());
}

// ==================== Orthogonal Polynomial Tests ====================

#[test]
fn test_legendre_p0() {
    let request = PolynomialRequest {
        polynomial_type: "legendre".to_string(),
        n: 0,
        x: 0.5,
        alpha: None,
    };

    let result = orthogonal_polynomial(request).unwrap();
    // P_0(x) = 1
    assert!((result.value - 1.0).abs() < 1e-10);
}

#[test]
fn test_legendre_p1() {
    let request = PolynomialRequest {
        polynomial_type: "legendre".to_string(),
        n: 1,
        x: 0.5,
        alpha: None,
    };

    let result = orthogonal_polynomial(request).unwrap();
    // P_1(x) = x
    assert!((result.value - 0.5).abs() < 1e-10);
}

#[test]
fn test_legendre_p2() {
    let request = PolynomialRequest {
        polynomial_type: "legendre".to_string(),
        n: 2,
        x: 0.0,
        alpha: None,
    };

    let result = orthogonal_polynomial(request).unwrap();
    // P_2(0) = -1/2
    assert!((result.value - (-0.5)).abs() < 1e-10);
}

#[test]
fn test_hermite_h0() {
    let request = PolynomialRequest {
        polynomial_type: "hermite".to_string(),
        n: 0,
        x: 1.0,
        alpha: None,
    };

    let result = orthogonal_polynomial(request).unwrap();
    // H_0(x) = 1
    assert!((result.value - 1.0).abs() < 1e-10);
}

#[test]
fn test_hermite_h1() {
    let request = PolynomialRequest {
        polynomial_type: "hermite".to_string(),
        n: 1,
        x: 2.0,
        alpha: None,
    };

    let result = orthogonal_polynomial(request).unwrap();
    // H_1(x) = 2x
    assert!((result.value - 4.0).abs() < 1e-10);
}

#[test]
fn test_laguerre_l0() {
    let request = PolynomialRequest {
        polynomial_type: "laguerre".to_string(),
        n: 0,
        x: 1.0,
        alpha: None,
    };

    let result = orthogonal_polynomial(request).unwrap();
    // L_0(x) = 1
    assert!((result.value - 1.0).abs() < 1e-10);
}

#[test]
fn test_chebyshev_t0() {
    let request = PolynomialRequest {
        polynomial_type: "chebyshev".to_string(),
        n: 0,
        x: 0.5,
        alpha: None,
    };

    let result = orthogonal_polynomial(request).unwrap();
    // T_0(x) = 1
    assert!((result.value - 1.0).abs() < 1e-10);
}

#[test]
fn test_chebyshev_t1() {
    let request = PolynomialRequest {
        polynomial_type: "chebyshev".to_string(),
        n: 1,
        x: 0.5,
        alpha: None,
    };

    let result = orthogonal_polynomial(request).unwrap();
    // T_1(x) = x
    assert!((result.value - 0.5).abs() < 1e-10);
}

#[test]
fn test_polynomial_invalid_type() {
    let request = PolynomialRequest {
        polynomial_type: "invalid".to_string(),
        n: 0,
        x: 1.0,
        alpha: None,
    };

    assert!(orthogonal_polynomial(request).is_err());
}

// ==================== Elliptic Integral Tests ====================

#[test]
fn test_elliptic_k_zero() {
    let k = elliptic_k(0.0).unwrap();
    // K(0) = π/2 ≈ 1.5708
    assert!((k - std::f64::consts::PI / 2.0).abs() < 1e-4);
}

#[test]
fn test_elliptic_e_zero() {
    let e = elliptic_e(0.0).unwrap();
    // E(0) = π/2 ≈ 1.5708
    assert!((e - std::f64::consts::PI / 2.0).abs() < 1e-4);
}

#[test]
fn test_elliptic_f() {
    // F(π/4, 0.5) - incomplete elliptic integral
    let result = elliptic_f(std::f64::consts::PI / 4.0, 0.5).unwrap();
    assert!(result > 0.0 && result < 2.0);
}

#[test]
fn test_elliptic_integral_request() {
    let request = EllipticIntegralRequest {
        integral_type: "K".to_string(),
        k: 0.5,
        phi: None,
        n: None,
    };

    let result = elliptic_integral(request).unwrap();
    // K(0.5) should be positive
    assert!(result.value > 1.5 && result.value < 2.0);
}

#[test]
fn test_elliptic_integral_invalid_type() {
    let request = EllipticIntegralRequest {
        integral_type: "X".to_string(),
        k: 0.5,
        phi: None,
        n: None,
    };

    assert!(elliptic_integral(request).is_err());
}
