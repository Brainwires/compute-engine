//! Unit tests for advanced_calculus API handler
//!
//! Tests all handler operations including:
//! - Derivative (polynomial derivatives with coefficients)
//! - Integral (polynomial integration, definite and indefinite)
//! - Limit (polynomial and rational function limits with L'Hôpital's rule)

use super::*; // Import from parent module (advanced_calculus handler)
use crate::api::types::ComputationRequest;
use serde_json::json;
use std::collections::HashMap;

fn create_request(operation: &str, parameters: serde_json::Value) -> ComputationRequest {
    let params_map: HashMap<String, serde_json::Value> = serde_json::from_value(parameters)
        .expect("Failed to convert parameters to HashMap");

    ComputationRequest {
        module: "advanced_calculus".to_string(),
        operation: operation.to_string(),
        parameters: params_map,
    }
}

// ============================================================================
// Derivative Tests
// ============================================================================

#[test]
fn test_derivative_coefficients_only() {
    // d/dx (2 + 3x + 4x²) = 3 + 8x
    let params = json!({
        "coefficients": [2.0, 3.0, 4.0]
    });
    let request = create_request("derivative", params);
    let response = handle(&request);

    assert!(response.success, "Derivative computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let deriv_coeffs = result.get("derivative_coefficients").unwrap().as_array().unwrap();

    assert_eq!(deriv_coeffs.len(), 2);
    assert_eq!(deriv_coeffs[0].as_f64().unwrap(), 3.0);
    assert_eq!(deriv_coeffs[1].as_f64().unwrap(), 8.0);
}

#[test]
fn test_derivative_constant() {
    // d/dx (5) = 0
    let params = json!({
        "coefficients": [5.0]
    });
    let request = create_request("derivative", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();
    let deriv_coeffs = result.get("derivative_coefficients").unwrap().as_array().unwrap();
    assert_eq!(deriv_coeffs.len(), 0); // Derivative of constant is 0
}

#[test]
fn test_derivative_linear() {
    // d/dx (7 + 2x) = 2
    let params = json!({
        "coefficients": [7.0, 2.0]
    });
    let request = create_request("derivative", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();
    let deriv_coeffs = result.get("derivative_coefficients").unwrap().as_array().unwrap();

    assert_eq!(deriv_coeffs.len(), 1);
    assert_eq!(deriv_coeffs[0].as_f64().unwrap(), 2.0);
}

#[test]
fn test_derivative_with_evaluation() {
    // d/dx (1 + 2x + 3x²) = 2 + 6x, evaluated at x=2: 2 + 6(2) = 14
    let params = json!({
        "coefficients": [1.0, 2.0, 3.0],
        "x": 2.0
    });
    let request = create_request("derivative", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();

    let derivative = result.get("derivative").unwrap().as_f64().unwrap();
    assert!((derivative - 14.0).abs() < 1e-10);

    let x = result.get("x").unwrap().as_f64().unwrap();
    assert_eq!(x, 2.0);

    assert!(result.get("derivative_coefficients").is_some());
}

#[test]
fn test_derivative_cubic_polynomial() {
    // d/dx (x³ + 2x² + 3x + 4) = 3x² + 4x + 3
    let params = json!({
        "coefficients": [4.0, 3.0, 2.0, 1.0]
    });
    let request = create_request("derivative", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();
    let deriv_coeffs = result.get("derivative_coefficients").unwrap().as_array().unwrap();

    assert_eq!(deriv_coeffs.len(), 3);
    assert_eq!(deriv_coeffs[0].as_f64().unwrap(), 3.0);
    assert_eq!(deriv_coeffs[1].as_f64().unwrap(), 4.0);
    assert_eq!(deriv_coeffs[2].as_f64().unwrap(), 3.0);
}

#[test]
fn test_derivative_invalid_params() {
    let params = json!({
        "invalid_field": "test"
    });
    let request = create_request("derivative", params);
    let response = handle(&request);

    assert!(!response.success, "Should fail with invalid parameters");
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Invalid parameters"));
}

// ============================================================================
// Integral Tests
// ============================================================================

#[test]
fn test_integral_indefinite() {
    // ∫(2 + 3x)dx = 2x + 3x²/2 + C
    let params = json!({
        "coefficients": [2.0, 3.0]
    });
    let request = create_request("integral", params);
    let response = handle(&request);

    assert!(response.success, "Integral computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let integral_coeffs = result.get("antiderivative_coefficients").unwrap().as_array().unwrap();

    assert_eq!(integral_coeffs.len(), 3);
    assert_eq!(integral_coeffs[0].as_f64().unwrap(), 0.0); // C = 0
    assert_eq!(integral_coeffs[1].as_f64().unwrap(), 2.0);
    assert_eq!(integral_coeffs[2].as_f64().unwrap(), 1.5); // 3/2
}

#[test]
fn test_integral_constant() {
    // ∫5 dx = 5x + C
    let params = json!({
        "coefficients": [5.0]
    });
    let request = create_request("integral", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();
    let integral_coeffs = result.get("antiderivative_coefficients").unwrap().as_array().unwrap();

    assert_eq!(integral_coeffs.len(), 2);
    assert_eq!(integral_coeffs[0].as_f64().unwrap(), 0.0); // C
    assert_eq!(integral_coeffs[1].as_f64().unwrap(), 5.0);
}

#[test]
fn test_integral_definite() {
    // ∫₁³ (x²) dx = [x³/3]₁³ = 27/3 - 1/3 = 26/3 ≈ 8.667
    let params = json!({
        "coefficients": [0.0, 0.0, 1.0], // x²
        "lower": 1.0,
        "upper": 3.0
    });
    let request = create_request("integral", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();

    let integral_value = result.get("integral").unwrap().as_f64().unwrap();
    let expected = 26.0 / 3.0;
    assert!((integral_value - expected).abs() < 1e-10,
            "Expected {}, got {}", expected, integral_value);

    assert_eq!(result.get("lower").unwrap().as_f64().unwrap(), 1.0);
    assert_eq!(result.get("upper").unwrap().as_f64().unwrap(), 3.0);
    assert!(result.get("antiderivative_coefficients").is_some());
}

#[test]
fn test_integral_definite_linear() {
    // ∫₀² (2x + 1) dx = [x² + x]₀² = 4 + 2 = 6
    let params = json!({
        "coefficients": [1.0, 2.0],
        "lower": 0.0,
        "upper": 2.0
    });
    let request = create_request("integral", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();

    let integral_value = result.get("integral").unwrap().as_f64().unwrap();
    assert!((integral_value - 6.0).abs() < 1e-10);
}

#[test]
fn test_integral_cubic_polynomial() {
    // ∫(1 + x + x² + x³)dx = x + x²/2 + x³/3 + x⁴/4 + C
    let params = json!({
        "coefficients": [1.0, 1.0, 1.0, 1.0]
    });
    let request = create_request("integral", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();
    let integral_coeffs = result.get("antiderivative_coefficients").unwrap().as_array().unwrap();

    assert_eq!(integral_coeffs.len(), 5);
    assert_eq!(integral_coeffs[0].as_f64().unwrap(), 0.0); // C
    assert_eq!(integral_coeffs[1].as_f64().unwrap(), 1.0);
    assert_eq!(integral_coeffs[2].as_f64().unwrap(), 0.5);
    assert!((integral_coeffs[3].as_f64().unwrap() - 1.0/3.0).abs() < 1e-10);
    assert_eq!(integral_coeffs[4].as_f64().unwrap(), 0.25);
}

#[test]
fn test_integral_invalid_params() {
    let params = json!({
        "invalid_field": "test"
    });
    let request = create_request("integral", params);
    let response = handle(&request);

    assert!(!response.success, "Should fail with invalid parameters");
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Invalid parameters"));
}

// ============================================================================
// Limit Tests
// ============================================================================

#[test]
fn test_limit_polynomial_simple() {
    // lim(x→2) (1 + 2x + x²) = 1 + 4 + 4 = 9
    let params = json!({
        "coefficients": [1.0, 2.0, 1.0],
        "x": 2.0
    });
    let request = create_request("limit", params);
    let response = handle(&request);

    assert!(response.success, "Limit computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    let limit = result.get("limit").unwrap().as_f64().unwrap();
    assert!((limit - 9.0).abs() < 1e-10);
    assert_eq!(result.get("x").unwrap().as_f64().unwrap(), 2.0);
}

#[test]
fn test_limit_rational_function() {
    // lim(x→2) (x² - 4)/(x - 2) = lim(x→2) (x² - 4)/(x - 2)
    // Numerator: x² - 4 = [−4, 0, 1]
    // Denominator: x - 2 = [−2, 1]
    // At x=2: (4-4)/(2-2) = 0/0 → L'Hôpital's rule
    // d/dx(x² - 4) = 2x, at x=2: 4
    // d/dx(x - 2) = 1
    // Result: 4/1 = 4
    let params = json!({
        "function_type": "rational",
        "numerator": [-4.0, 0.0, 1.0],
        "denominator": [-2.0, 1.0],
        "point": 2.0
    });
    let request = create_request("limit", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();
    let limit = result.get("limit").unwrap().as_f64().unwrap();
    assert!((limit - 4.0).abs() < 1e-10);
    assert_eq!(result.get("method").unwrap().as_str().unwrap(), "lhopital");
}

#[test]
fn test_limit_rational_no_indeterminate() {
    // lim(x→1) (x + 1)/(x - 2) = 2/(-1) = -2
    let params = json!({
        "function_type": "rational",
        "numerator": [1.0, 1.0],
        "denominator": [-2.0, 1.0],
        "point": 1.0
    });
    let request = create_request("limit", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();
    let limit = result.get("limit").unwrap().as_f64().unwrap();
    assert!((limit + 2.0).abs() < 1e-10);
    assert!(result.get("method").is_none()); // No L'Hôpital needed
}

#[test]
fn test_limit_rational_infinite() {
    // lim(x→2) 1/(x - 2) as x→2⁺ = +∞
    // Note: f64::INFINITY serializes as null in JSON, so we check the "infinite" flag
    let params = json!({
        "function_type": "rational",
        "numerator": [1.0],
        "denominator": [-2.0, 1.0],
        "point": 2.0
    });
    let request = create_request("limit", params);
    let response = handle(&request);

    assert!(response.success, "Response should be successful");
    let result = response.result.as_ref().expect("Should have result");

    // The limit field will be null because f64::INFINITY doesn't serialize to JSON properly
    assert!(result.get("limit").is_some(), "Should have 'limit' field");
    assert!(result.get("infinite").is_some(), "Should have 'infinite' field");
    assert_eq!(result.get("infinite").unwrap().as_bool().unwrap(), true);
    assert_eq!(result.get("point").unwrap().as_f64().unwrap(), 2.0);
}

#[test]
fn test_limit_rational_lhopital_x_squared() {
    // lim(x→0) (x²)/(2x) = 0/0 → L'Hôpital
    // d/dx(x²) = 2x, at x=0: 0
    // d/dx(2x) = 2
    // Result: 0/2 = 0
    let params = json!({
        "function_type": "rational",
        "numerator": [0.0, 0.0, 1.0], // x²
        "denominator": [0.0, 2.0],     // 2x
        "point": 0.0
    });
    let request = create_request("limit", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();
    let limit = result.get("limit").unwrap().as_f64().unwrap();
    assert!(limit.abs() < 1e-10); // Should be 0
    assert_eq!(result.get("method").unwrap().as_str().unwrap(), "lhopital");
}

#[test]
fn test_limit_polynomial_at_zero() {
    // lim(x→0) (3 + 2x²) = 3
    let params = json!({
        "coefficients": [3.0, 0.0, 2.0],
        "x": 0.0
    });
    let request = create_request("limit", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();
    let limit = result.get("limit").unwrap().as_f64().unwrap();
    assert!((limit - 3.0).abs() < 1e-10);
}

#[test]
fn test_limit_invalid_params() {
    let params = json!({
        "invalid_field": "test"
    });
    let request = create_request("limit", params);
    let response = handle(&request);

    assert!(!response.success, "Should fail with invalid parameters");
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Invalid parameters"));
}

// ============================================================================
// Operation Not Found Tests
// ============================================================================

#[test]
fn test_unknown_operation() {
    let params = json!({
        "coefficients": [1.0, 2.0]
    });
    let request = create_request("nonexistent_operation", params);
    let response = handle(&request);

    assert!(!response.success, "Should fail for unknown operation");
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Unknown operation"));
}

#[test]
fn test_operation_case_sensitivity() {
    // Operations are case-sensitive, "Derivative" should fail
    let params = json!({
        "coefficients": [1.0, 2.0]
    });
    let request = create_request("Derivative", params);
    let response = handle(&request);

    assert!(!response.success, "Operations should be case-sensitive");
    assert!(response.error.is_some());
}
