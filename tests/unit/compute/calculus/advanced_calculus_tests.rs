//! Unit tests for advanced calculus functions (special functions module)
//!
//! Tests coverage:
//! - Gamma function and related functions
//! - Beta function
//! - Bessel functions (J, Y, I, K)
//! - Error functions (erf, erfc, erfcx, erfi)
//! - Elliptic integrals (K, E, F, Pi)
//! - Orthogonal polynomials (Legendre, Hermite, Laguerre, Chebyshev)
//! - Airy functions (Ai, Bi)

use crate::compute::special_functions::*;

// ==================== Gamma Function Tests ====================

#[test]
fn test_gamma_positive_integer() {
    let request = GammaRequest {
        x: 5.0,
        function: "gamma".to_string(),
        y: None,
    };
    let result = gamma_function(request).unwrap();
    // Gamma(5) = 4! = 24
    assert!((result.value - 24.0).abs() < 0.1);
}

#[test]
fn test_gamma_half_integer() {
    let request = GammaRequest {
        x: 0.5,
        function: "gamma".to_string(),
        y: None,
    };
    let result = gamma_function(request).unwrap();
    // Gamma(0.5) = sqrt(pi) ≈ 1.772
    assert!((result.value - 1.772).abs() < 0.01);
}

#[test]
fn test_log_gamma() {
    let request = GammaRequest {
        x: 10.0,
        function: "log_gamma".to_string(),
        y: None,
    };
    let result = gamma_function(request).unwrap();
    // Gamma(10) = 9! = 362880, ln(362880) ≈ 12.8
    assert!((result.value - 12.8).abs() < 0.1);
}

#[test]
fn test_digamma() {
    let request = GammaRequest {
        x: 5.0,
        function: "digamma".to_string(),
        y: None,
    };
    let result = gamma_function(request).unwrap();
    // Should return a reasonable value
    assert!(result.value.is_finite());
}

#[test]
fn test_beta_function() {
    let request = GammaRequest {
        x: 2.0,
        function: "beta".to_string(),
        y: Some(3.0),
    };
    let result = gamma_function(request).unwrap();
    // Beta(2,3) = Gamma(2)*Gamma(3)/Gamma(5) = 1*2/24 = 1/12 ≈ 0.0833
    assert!((result.value - 0.0833).abs() < 0.01);
}

// ==================== Bessel Function Tests ====================

#[test]
fn test_bessel_j0_at_zero() {
    let request = BesselRequest {
        function_type: "J".to_string(),
        order: 0.0,
        x: 0.0,
    };
    let result = bessel_function(request).unwrap();
    // J_0(0) = 1
    assert!((result.value - 1.0).abs() < 1e-6);
}

#[test]
fn test_bessel_j1_at_zero() {
    let request = BesselRequest {
        function_type: "J".to_string(),
        order: 1.0,
        x: 0.0,
    };
    let result = bessel_function(request).unwrap();
    // J_1(0) = 0
    assert!(result.value.abs() < 1e-6);
}

#[test]
fn test_bessel_j0_positive() {
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
fn test_bessel_y_positive() {
    let request = BesselRequest {
        function_type: "Y".to_string(),
        order: 0.0,
        x: 1.0,
    };
    let result = bessel_function(request).unwrap();
    // Y_0(1) ≈ 0.0883
    assert!(result.value.is_finite());
}

#[test]
fn test_bessel_i0_at_zero() {
    let request = BesselRequest {
        function_type: "I".to_string(),
        order: 0.0,
        x: 0.0,
    };
    let result = bessel_function(request).unwrap();
    // I_0(0) = 1
    assert!((result.value - 1.0).abs() < 1e-6);
}

#[test]
fn test_bessel_k_positive() {
    let request = BesselRequest {
        function_type: "K".to_string(),
        order: 0.0,
        x: 1.0,
    };
    let result = bessel_function(request).unwrap();
    // Should return finite value
    assert!(result.value.is_finite());
}

// ==================== Error Function Tests ====================

#[test]
fn test_erf_at_zero() {
    let request = ErrorFunctionRequest {
        x: 0.0,
        function: "erf".to_string(),
    };
    let result = error_function(request).unwrap();
    // erf(0) = 0
    assert!(result.value.abs() < 1e-6);
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
    // erf(-1) ≈ -0.8427
    assert!((result.value + 0.8427).abs() < 0.01);
}

#[test]
fn test_erfc_at_zero() {
    let request = ErrorFunctionRequest {
        x: 0.0,
        function: "erfc".to_string(),
    };
    let result = error_function(request).unwrap();
    // erfc(0) = 1
    assert!((result.value - 1.0).abs() < 1e-6);
}

#[test]
fn test_erfcx() {
    let request = ErrorFunctionRequest {
        x: 1.0,
        function: "erfcx".to_string(),
    };
    let result = error_function(request).unwrap();
    // Should be finite and positive
    assert!(result.value > 0.0 && result.value.is_finite());
}

#[test]
fn test_erfi() {
    let request = ErrorFunctionRequest {
        x: 1.0,
        function: "erfi".to_string(),
    };
    let result = error_function(request).unwrap();
    // Should be finite
    assert!(result.value.is_finite());
}

// ==================== Elliptic Integral Tests ====================

#[test]
fn test_elliptic_k() {
    let request = EllipticIntegralRequest {
        integral_type: "K".to_string(),
        k: 0.5,
        phi: None,
        n: None,
    };
    let result = elliptic_integral(request).unwrap();
    // K(0.5) ≈ 1.6858
    assert!((result.value - 1.6858).abs() < 0.01);
}

#[test]
fn test_elliptic_k_zero() {
    let request = EllipticIntegralRequest {
        integral_type: "K".to_string(),
        k: 0.0,
        phi: None,
        n: None,
    };
    let result = elliptic_integral(request).unwrap();
    // K(0) = π/2
    assert!((result.value - std::f64::consts::FRAC_PI_2).abs() < 0.01);
}

#[test]
fn test_elliptic_e() {
    let request = EllipticIntegralRequest {
        integral_type: "E".to_string(),
        k: 0.5,
        phi: None,
        n: None,
    };
    let result = elliptic_integral(request).unwrap();
    // E(0.5) should be finite and in reasonable range (simplified approximation)
    assert!(result.value > 1.0 && result.value < 2.0 && result.value.is_finite());
}

#[test]
fn test_elliptic_f() {
    let request = EllipticIntegralRequest {
        integral_type: "F".to_string(),
        k: 0.5,
        phi: Some(std::f64::consts::FRAC_PI_4),
        n: None,
    };
    let result = elliptic_integral(request).unwrap();
    // Should be finite and positive
    assert!(result.value > 0.0 && result.value.is_finite());
}

#[test]
fn test_elliptic_pi() {
    let request = EllipticIntegralRequest {
        integral_type: "Pi".to_string(),
        k: 0.5,
        phi: Some(std::f64::consts::FRAC_PI_4),
        n: Some(0.3),
    };
    let result = elliptic_integral(request).unwrap();
    // Should be finite and positive
    assert!(result.value > 0.0 && result.value.is_finite());
}

// ==================== Legendre Polynomial Tests ====================

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
    assert!((result.value - 1.0).abs() < 1e-6);
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
    // P_1(x) = x, P_1(0.5) = 0.5
    assert!((result.value - 0.5).abs() < 1e-6);
}

#[test]
fn test_legendre_p2() {
    let request = PolynomialRequest {
        polynomial_type: "legendre".to_string(),
        n: 2,
        x: 0.5,
        alpha: None,
    };
    let result = orthogonal_polynomial(request).unwrap();
    // P_2(x) = (3x^2 - 1)/2, P_2(0.5) = -0.125
    assert!((result.value - (-0.125)).abs() < 1e-6);
}

// ==================== Hermite Polynomial Tests ====================

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
    assert!((result.value - 1.0).abs() < 1e-6);
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
    // H_1(x) = 2x, H_1(2) = 4
    assert!((result.value - 4.0).abs() < 1e-6);
}

#[test]
fn test_hermite_h2() {
    let request = PolynomialRequest {
        polynomial_type: "hermite".to_string(),
        n: 2,
        x: 1.0,
        alpha: None,
    };
    let result = orthogonal_polynomial(request).unwrap();
    // H_2(x) = 4x^2 - 2, H_2(1) = 2
    assert!((result.value - 2.0).abs() < 1e-6);
}

// ==================== Laguerre Polynomial Tests ====================

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
    assert!((result.value - 1.0).abs() < 1e-6);
}

#[test]
fn test_laguerre_l1() {
    let request = PolynomialRequest {
        polynomial_type: "laguerre".to_string(),
        n: 1,
        x: 1.0,
        alpha: None,
    };
    let result = orthogonal_polynomial(request).unwrap();
    // L_1(x) = 1 - x, L_1(1) = 0
    assert!(result.value.abs() < 1e-6);
}

#[test]
fn test_laguerre_l2() {
    let request = PolynomialRequest {
        polynomial_type: "laguerre".to_string(),
        n: 2,
        x: 1.0,
        alpha: None,
    };
    let result = orthogonal_polynomial(request).unwrap();
    // L_2(x) = (x^2 - 4x + 2)/2, L_2(1) = -0.5
    assert!((result.value - (-0.5)).abs() < 1e-6);
}

// ==================== Chebyshev Polynomial Tests ====================

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
    assert!((result.value - 1.0).abs() < 1e-6);
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
    // T_1(x) = x, T_1(0.5) = 0.5
    assert!((result.value - 0.5).abs() < 1e-6);
}

#[test]
fn test_chebyshev_t2() {
    let request = PolynomialRequest {
        polynomial_type: "chebyshev".to_string(),
        n: 2,
        x: 0.5,
        alpha: None,
    };
    let result = orthogonal_polynomial(request).unwrap();
    // T_2(x) = 2x^2 - 1, T_2(0.5) = -0.5
    assert!((result.value - (-0.5)).abs() < 1e-6);
}

// ==================== Airy Function Tests ====================

#[test]
fn test_airy_ai_at_zero() {
    let request = AiryRequest {
        function_type: "Ai".to_string(),
        x: 0.0,
    };
    let result = airy_function(request).unwrap();
    // Ai(0) ≈ 0.3550
    assert!((result.value - 0.3550).abs() < 0.01);
}

#[test]
fn test_airy_ai_positive() {
    let request = AiryRequest {
        function_type: "Ai".to_string(),
        x: 1.0,
    };
    let result = airy_function(request).unwrap();
    // Should decay exponentially
    assert!(result.value > 0.0 && result.value < 0.5);
}

#[test]
fn test_airy_ai_negative() {
    let request = AiryRequest {
        function_type: "Ai".to_string(),
        x: -1.0,
    };
    let result = airy_function(request).unwrap();
    // Should oscillate
    assert!(result.value.is_finite());
}

#[test]
fn test_airy_bi_at_zero() {
    let request = AiryRequest {
        function_type: "Bi".to_string(),
        x: 0.0,
    };
    let result = airy_function(request).unwrap();
    // Bi(0) ≈ 0.6149
    assert!((result.value - 0.6149).abs() < 0.01);
}

#[test]
fn test_airy_bi_positive() {
    let request = AiryRequest {
        function_type: "Bi".to_string(),
        x: 1.0,
    };
    let result = airy_function(request).unwrap();
    // Should grow exponentially
    assert!(result.value > 0.5 && result.value.is_finite());
}

#[test]
fn test_airy_bi_negative() {
    let request = AiryRequest {
        function_type: "Bi".to_string(),
        x: -1.0,
    };
    let result = airy_function(request).unwrap();
    // Should oscillate
    assert!(result.value.is_finite());
}

#[test]
fn test_airy_ai_prime() {
    let request = AiryRequest {
        function_type: "Ai_prime".to_string(),
        x: 0.0,
    };
    let result = airy_function(request).unwrap();
    // Ai'(0) ≈ -0.2588
    assert!(result.value < 0.0 && result.value.is_finite());
}

#[test]
fn test_airy_bi_prime() {
    let request = AiryRequest {
        function_type: "Bi_prime".to_string(),
        x: 0.0,
    };
    let result = airy_function(request).unwrap();
    // Bi'(0) ≈ 0.4483
    assert!(result.value > 0.0 && result.value.is_finite());
}
