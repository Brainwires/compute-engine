//! Unit tests for special_functions API handler
//!
//! Tests all handler operations including:
//! - Bessel functions (J, Y, I, K types)
//! - Gamma function and variants (gamma, log_gamma, digamma, beta)
//! - Error functions (erf, erfc, erfcx, erfi)
//! - Elliptic integrals (K, E, F, Pi types)
//! - Orthogonal polynomials (Legendre, Hermite, Laguerre, Chebyshev)
//! - Airy functions (Ai, Bi, derivatives)
//! - Invalid parameters and error cases

use super::*; // Import from parent module (special_functions handler)
use crate::api::types::ComputationRequest;
use serde_json::json;
use std::collections::HashMap;

fn create_request(operation: &str, parameters: serde_json::Value) -> ComputationRequest {
    let params_map: HashMap<String, serde_json::Value> = serde_json::from_value(parameters)
        .expect("Failed to convert parameters to HashMap");

    ComputationRequest {
        module: "special_functions".to_string(),
        operation: operation.to_string(),
        parameters: params_map,
    }
}

// ============================================================================
// Bessel Functions Tests
// ============================================================================

#[test]
fn test_bessel_j0() {
    // Test Bessel J_0(1.0)
    let params = json!({
        "function_type": "J",
        "order": 0.0,
        "x": 1.0
    });
    let request = create_request("bessel", params);
    let response = handle(&request);

    assert!(response.success, "Bessel J_0 computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let value = result.get("value").unwrap().as_f64().unwrap();

    // J_0(1) ≈ 0.7652
    assert!((value - 0.7652).abs() < 0.01, "J_0(1) should be approximately 0.7652");
}

#[test]
fn test_bessel_j1() {
    // Test Bessel J_1(2.0)
    let params = json!({
        "function_type": "J",
        "order": 1.0,
        "x": 2.0
    });
    let request = create_request("bessel", params);
    let response = handle(&request);

    assert!(response.success, "Bessel J_1 computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let value = result.get("value").unwrap().as_f64().unwrap();
    let function_type = result.get("function_type").unwrap().as_str().unwrap();

    assert_eq!(function_type, "J");
    assert!(value.is_finite(), "J_1(2) should be finite");
}

#[test]
fn test_bessel_y() {
    // Test Bessel Y_0(1.0)
    let params = json!({
        "function_type": "Y",
        "order": 0.0,
        "x": 1.0
    });
    let request = create_request("bessel", params);
    let response = handle(&request);

    assert!(response.success, "Bessel Y_0 computation should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_bessel_invalid_parameters() {
    // Test missing parameters
    let params = json!({
        "function_type": "J"
        // Missing order and x
    });
    let request = create_request("bessel", params);
    let response = handle(&request);

    assert!(!response.success, "Missing parameters should fail");
    assert!(response.error.is_some());
}

// ============================================================================
// Gamma Function Tests
// ============================================================================

#[test]
fn test_gamma_function() {
    // Test Γ(5) = 24 (4!)
    let params = json!({
        "x": 5.0,
        "function": "gamma"
    });
    let request = create_request("gamma", params);
    let response = handle(&request);

    assert!(response.success, "Gamma function computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let value = result.get("value").unwrap().as_f64().unwrap();

    // Γ(5) = 4! = 24
    assert!((value - 24.0).abs() < 0.01, "Γ(5) should be 24");
}

#[test]
fn test_gamma_half() {
    // Test Γ(0.5) = √π
    let params = json!({
        "x": 0.5,
        "function": "gamma"
    });
    let request = create_request("gamma", params);
    let response = handle(&request);

    assert!(response.success, "Gamma(0.5) computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let value = result.get("value").unwrap().as_f64().unwrap();

    // Γ(0.5) = √π ≈ 1.7725
    assert!((value - std::f64::consts::PI.sqrt()).abs() < 0.01);
}

#[test]
fn test_log_gamma() {
    // Test ln(Γ(10))
    let params = json!({
        "x": 10.0,
        "function": "log_gamma"
    });
    let request = create_request("gamma", params);
    let response = handle(&request);

    assert!(response.success, "Log gamma computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let value = result.get("value").unwrap().as_f64().unwrap();

    // ln(Γ(10)) = ln(9!) = ln(362880)
    let expected = 362880.0_f64.ln();
    assert!((value - expected).abs() < 0.1);
}

#[test]
fn test_digamma() {
    // Test digamma function ψ(x)
    let params = json!({
        "x": 1.0,
        "function": "digamma"
    });
    let request = create_request("gamma", params);
    let response = handle(&request);

    assert!(response.success, "Digamma computation should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_beta_function() {
    // Test Beta(2, 3) = Γ(2)Γ(3)/Γ(5) = 1*2/24 = 1/12
    let params = json!({
        "x": 2.0,
        "y": 3.0,
        "function": "beta"
    });
    let request = create_request("gamma", params);
    let response = handle(&request);

    assert!(response.success, "Beta function computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let value = result.get("value").unwrap().as_f64().unwrap();

    // Beta(2,3) = 1/12 ≈ 0.0833
    assert!((value - 1.0/12.0).abs() < 0.01);
}

// ============================================================================
// Error Function Tests
// ============================================================================

#[test]
fn test_erf() {
    // Test erf(0) = 0
    let params = json!({
        "x": 0.0,
        "function": "erf"
    });
    let request = create_request("erf", params);
    let response = handle(&request);

    assert!(response.success, "Error function computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let value = result.get("value").unwrap().as_f64().unwrap();

    assert!((value - 0.0).abs() < 1e-6, "erf(0) should be 0");
}

#[test]
fn test_erf_one() {
    // Test erf(1.0) ≈ 0.8427
    let params = json!({
        "x": 1.0,
        "function": "erf"
    });
    let request = create_request("erf", params);
    let response = handle(&request);

    assert!(response.success, "Error function computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let value = result.get("value").unwrap().as_f64().unwrap();

    assert!((value - 0.8427).abs() < 0.01, "erf(1) should be approximately 0.8427");
}

#[test]
fn test_erfc() {
    // Test erfc(0) = 1 - erf(0) = 1
    let params = json!({
        "x": 0.0,
        "function": "erfc"
    });
    let request = create_request("erf", params);
    let response = handle(&request);

    assert!(response.success, "Complementary error function should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let value = result.get("value").unwrap().as_f64().unwrap();

    assert!((value - 1.0).abs() < 1e-6, "erfc(0) should be 1");
}

#[test]
fn test_error_function_alias() {
    // Test "error_function" operation alias
    let params = json!({
        "x": 0.5,
        "function": "erf"
    });
    let request = create_request("error_function", params);
    let response = handle(&request);

    assert!(response.success, "error_function alias should work");
    assert!(response.result.is_some());
}

// ============================================================================
// Elliptic Integral Tests
// ============================================================================

#[test]
fn test_elliptic_k() {
    // Test complete elliptic integral of first kind K(0.5)
    let params = json!({
        "integral_type": "K",
        "k": 0.5
    });
    let request = create_request("elliptic", params);
    let response = handle(&request);

    assert!(response.success, "Elliptic K computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let value = result.get("value").unwrap().as_f64().unwrap();

    // K(0.5) ≈ 1.6858
    assert!((value - 1.6858).abs() < 0.01, "K(0.5) should be approximately 1.6858");
}

#[test]
fn test_elliptic_e() {
    // Test complete elliptic integral of second kind E(0.5)
    let params = json!({
        "integral_type": "E",
        "k": 0.5
    });
    let request = create_request("elliptic", params);
    let response = handle(&request);

    assert!(response.success, "Elliptic E computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let value = result.get("value").unwrap().as_f64().unwrap();

    // E(0.5) ≈ 1.5751 (complete elliptic integral of second kind)
    assert!((value - 1.5751).abs() < 0.01, "E(0.5) should be approximately 1.5751");
}

#[test]
fn test_elliptic_f() {
    // Test incomplete elliptic integral of first kind F(π/4, 0.5)
    let params = json!({
        "integral_type": "F",
        "k": 0.5,
        "phi": std::f64::consts::PI / 4.0
    });
    let request = create_request("elliptic", params);
    let response = handle(&request);

    assert!(response.success, "Elliptic F computation should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_elliptic_invalid_k() {
    // Test with invalid modulus k >= 1
    let params = json!({
        "integral_type": "K",
        "k": 1.5
    });
    let request = create_request("elliptic", params);
    let response = handle(&request);

    // Should fail or return NaN for invalid k
    if response.success {
        let result = response.result.unwrap();
        let value = result.get("value").unwrap().as_f64().unwrap();
        assert!(value.is_nan() || value.is_infinite(), "Invalid k should produce NaN or infinity");
    }
}

// ============================================================================
// Orthogonal Polynomial Tests
// ============================================================================

#[test]
fn test_legendre_polynomial() {
    // Test P_2(0.5) = (3*0.5^2 - 1)/2 = -0.125
    let params = json!({
        "polynomial_type": "legendre",
        "n": 2,
        "x": 0.5
    });
    let request = create_request("polynomial", params);
    let response = handle(&request);

    assert!(response.success, "Legendre polynomial computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let value = result.get("value").unwrap().as_f64().unwrap();

    // P_2(0.5) = -0.125
    assert!((value - (-0.125)).abs() < 1e-6);
}

#[test]
fn test_legendre_p0() {
    // Test P_0(x) = 1 for any x
    let params = json!({
        "polynomial_type": "legendre",
        "n": 0,
        "x": 0.7
    });
    let request = create_request("polynomial", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();
    let value = result.get("value").unwrap().as_f64().unwrap();

    assert!((value - 1.0).abs() < 1e-6, "P_0(x) should be 1");
}

#[test]
fn test_hermite_polynomial() {
    // Test Hermite polynomial H_n(x)
    let params = json!({
        "polynomial_type": "hermite",
        "n": 2,
        "x": 1.0
    });
    let request = create_request("polynomial", params);
    let response = handle(&request);

    assert!(response.success, "Hermite polynomial computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let value = result.get("value").unwrap().as_f64().unwrap();

    // H_2(1) = 4*1^2 - 2 = 2
    assert!((value - 2.0).abs() < 0.1);
}

#[test]
fn test_laguerre_polynomial() {
    // Test Laguerre polynomial L_n(x)
    let params = json!({
        "polynomial_type": "laguerre",
        "n": 2,
        "x": 1.0
    });
    let request = create_request("polynomial", params);
    let response = handle(&request);

    assert!(response.success, "Laguerre polynomial computation should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_chebyshev_polynomial() {
    // Test Chebyshev T_2(0.5) = 2*0.5^2 - 1 = -0.5
    let params = json!({
        "polynomial_type": "chebyshev",
        "n": 2,
        "x": 0.5
    });
    let request = create_request("polynomial", params);
    let response = handle(&request);

    assert!(response.success, "Chebyshev polynomial computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let value = result.get("value").unwrap().as_f64().unwrap();

    // T_2(0.5) = -0.5
    assert!((value - (-0.5)).abs() < 0.01);
}

#[test]
fn test_orthogonal_polynomial_alias() {
    // Test "orthogonal_polynomial" operation alias
    let params = json!({
        "polynomial_type": "legendre",
        "n": 1,
        "x": 0.5
    });
    let request = create_request("orthogonal_polynomial", params);
    let response = handle(&request);

    assert!(response.success, "orthogonal_polynomial alias should work");
    assert!(response.result.is_some());
}

// ============================================================================
// Airy Function Tests
// ============================================================================

#[test]
fn test_airy_ai() {
    // Test Airy Ai function
    let params = json!({
        "function_type": "Ai",
        "x": 0.0
    });
    let request = create_request("airy", params);
    let response = handle(&request);

    assert!(response.success, "Airy Ai computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let value = result.get("value").unwrap().as_f64().unwrap();

    // Ai(0) ≈ 0.3550
    assert!((value - 0.3550).abs() < 0.01, "Ai(0) should be approximately 0.3550");
}

#[test]
fn test_airy_bi() {
    // Test Airy Bi function
    let params = json!({
        "function_type": "Bi",
        "x": 0.0
    });
    let request = create_request("airy", params);
    let response = handle(&request);

    assert!(response.success, "Airy Bi computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let value = result.get("value").unwrap().as_f64().unwrap();

    // Bi(0) ≈ 0.6149
    assert!((value - 0.6149).abs() < 0.01, "Bi(0) should be approximately 0.6149");
}

#[test]
fn test_airy_ai_prime() {
    // Test Airy Ai derivative
    let params = json!({
        "function_type": "Ai_prime",
        "x": 0.0
    });
    let request = create_request("airy", params);
    let response = handle(&request);

    assert!(response.success, "Airy Ai' computation should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_airy_bi_prime() {
    // Test Airy Bi derivative
    let params = json!({
        "function_type": "Bi_prime",
        "x": 1.0
    });
    let request = create_request("airy", params);
    let response = handle(&request);

    assert!(response.success, "Airy Bi' computation should succeed");
    assert!(response.result.is_some());
}

// ============================================================================
// Error Handling Tests
// ============================================================================

#[test]
fn test_unknown_operation() {
    // Test unknown operation
    let params = json!({
        "x": 1.0
    });
    let request = create_request("unknown_function", params);
    let response = handle(&request);

    assert!(!response.success, "Unknown operation should fail");
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Unknown"));
}

#[test]
fn test_missing_required_parameters() {
    // Test gamma function without required parameter
    let params = json!({
        "function": "gamma"
        // Missing x parameter
    });
    let request = create_request("gamma", params);
    let response = handle(&request);

    assert!(!response.success, "Missing required parameters should fail");
    assert!(response.error.is_some());
}

#[test]
fn test_invalid_function_type() {
    // Test with invalid function type
    let params = json!({
        "function_type": "invalid_type",
        "order": 0.0,
        "x": 1.0
    });
    let request = create_request("bessel", params);
    let response = handle(&request);

    assert!(!response.success, "Invalid function type should fail");
    assert!(response.error.is_some());
}

#[test]
fn test_beta_missing_second_parameter() {
    // Test beta function without y parameter
    let params = json!({
        "x": 2.0,
        "function": "beta"
        // Missing y parameter
    });
    let request = create_request("gamma", params);
    let response = handle(&request);

    assert!(!response.success, "Beta function without y should fail");
    assert!(response.error.is_some());
}

#[test]
fn test_malformed_json_parameters() {
    // Test with parameters that can't be deserialized
    let params_map: HashMap<String, serde_json::Value> =
        [("x".to_string(), json!("not_a_number"))]
        .iter()
        .cloned()
        .collect();

    let request = ComputationRequest {
        module: "special_functions".to_string(),
        operation: "gamma".to_string(),
        parameters: params_map,
    };

    let response = handle(&request);
    assert!(!response.success, "Malformed parameters should fail");
}
