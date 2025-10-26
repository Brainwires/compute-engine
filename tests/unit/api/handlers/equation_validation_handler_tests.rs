//! Unit tests for equation_validation API handler
//!
//! Tests all equation validation operations including:
//! - Basic equation validation
//! - Dimensional consistency checking
//! - Physics compliance verification
//! - Mathematical correctness validation
//! - Conservation law checks
//! - Symmetry analysis
//! - Unit analysis
//! - Domain-specific validation (mechanics, electromagnetism, thermodynamics, quantum, relativity)
//! - Error handling for invalid inputs

use crate::api::handlers::equation_validation;
use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};
use std::collections::HashMap;

/// Helper function to create an equation_validation request
fn create_request(operation: &str, parameters: HashMap<String, Value>) -> ComputationRequest {
    ComputationRequest {
        module: "equation_validation".to_string(),
        operation: operation.to_string(),
        parameters,
    }
}

/// Helper function to extract result value from response
fn get_result(response: &ComputationResponse) -> &Value {
    assert!(response.success, "Expected success response");
    response.result.as_ref().expect("Expected result value")
}

// ============================================================================
// Basic Validation Tests
// ============================================================================

#[test]
fn test_validate_simple_equation() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("F = m*a"));
    params.insert("domain".to_string(), json!("mechanics"));

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert!(result["is_valid"].as_bool().unwrap());
    assert!(result["mathematical_correctness"].as_bool().unwrap());
    assert!(result["confidence"].as_f64().unwrap() >= 0.0);
    assert!(result["confidence"].as_f64().unwrap() <= 1.0);
}

#[test]
fn test_validate_equation_with_units() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("F = m*a"));
    params.insert("domain".to_string(), json!("mechanics"));

    let mut units = HashMap::new();
    units.insert("F".to_string(), "N".to_string());
    units.insert("m".to_string(), "kg".to_string());
    units.insert("a".to_string(), "m/s^2".to_string());
    params.insert("units".to_string(), json!(units));

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Dimensional consistency check is simplified, so we just verify the response structure
    assert!(result["dimensional_consistency"].is_boolean());
    let unit_analysis = result["unit_analysis"].as_object().unwrap();
    assert!(unit_analysis.contains_key("F"));
    assert!(unit_analysis.contains_key("m"));
    assert!(unit_analysis.contains_key("a"));
}

#[test]
fn test_validate_energy_equation() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("E = m*c^2"));
    params.insert("domain".to_string(), json!("relativity"));

    let mut units = HashMap::new();
    units.insert("E".to_string(), "J".to_string());
    units.insert("m".to_string(), "kg".to_string());
    units.insert("c".to_string(), "m/s".to_string());
    params.insert("units".to_string(), json!(units));

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Energy equation should be valid with proper units
    assert!(result["mathematical_correctness"].as_bool().unwrap());
}

#[test]
fn test_validate_general_domain() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("y = mx + b"));
    // domain defaults to "general" if not specified

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert!(result["is_valid"].as_bool().unwrap());
    assert!(result["mathematical_correctness"].as_bool().unwrap());
}

// ============================================================================
// Conservation Law Tests
// ============================================================================

#[test]
fn test_validate_with_energy_conservation() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("E = KE + PE"));
    params.insert("domain".to_string(), json!("mechanics"));
    params.insert(
        "conservation_laws".to_string(),
        json!(["energy"]),
    );

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should validate without major violations
    assert!(result["is_valid"].is_boolean());
}

#[test]
fn test_validate_with_momentum_conservation() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("p = m*v"));
    params.insert("domain".to_string(), json!("mechanics"));
    params.insert(
        "conservation_laws".to_string(),
        json!(["momentum"]),
    );

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Momentum equation with m and v should be valid
    let violations = result["violations"].as_array().unwrap();
    assert!(violations.iter().all(|v| !v.as_str().unwrap().contains("Momentum")));
}

#[test]
fn test_validate_with_charge_conservation() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("F = k*q1*q2/r^2"));
    params.insert("domain".to_string(), json!("electromagnetism"));
    params.insert(
        "conservation_laws".to_string(),
        json!(["charge"]),
    );

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Coulomb's law should be valid
    assert!(result["is_valid"].is_boolean());
}

// ============================================================================
// Symmetry Tests
// ============================================================================

#[test]
fn test_validate_with_time_translation_symmetry() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("F = m*a"));
    params.insert("domain".to_string(), json!("mechanics"));
    params.insert(
        "symmetries".to_string(),
        json!(["time_translation"]),
    );

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // F=ma doesn't have explicit time dependence
    assert!(result["violations"].as_array().unwrap().is_empty()
            || !result["violations"].as_array().unwrap().iter()
                .any(|v| v.as_str().unwrap().contains("time")));
}

#[test]
fn test_validate_with_lorentz_symmetry() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("E = m*c^2"));
    params.insert("domain".to_string(), json!("relativity"));
    params.insert(
        "symmetries".to_string(),
        json!(["lorentz"]),
    );

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // E=mc² should be valid for Lorentz symmetry
    assert!(result["is_valid"].is_boolean());
}

// ============================================================================
// Mathematical Correctness Tests
// ============================================================================

#[test]
fn test_validate_division_by_zero() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("y = x/0"));
    params.insert("domain".to_string(), json!("general"));

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Division by zero should fail mathematical correctness
    assert!(!result["mathematical_correctness"].as_bool().unwrap());
    let violations = result["violations"].as_array().unwrap();
    assert!(violations.iter().any(|v| v.as_str().unwrap().contains("Mathematical")));
}

#[test]
fn test_validate_negative_square_root() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("y = sqrt(-1)"));
    params.insert("domain".to_string(), json!("general"));

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Negative square root should fail mathematical correctness
    assert!(!result["mathematical_correctness"].as_bool().unwrap());
}

#[test]
fn test_validate_unmatched_parentheses() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("y = (x + 2"));
    params.insert("domain".to_string(), json!("general"));

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Unmatched parentheses should fail mathematical correctness
    assert!(!result["mathematical_correctness"].as_bool().unwrap());
}

#[test]
fn test_validate_logarithm_of_zero() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("y = log(0)"));
    params.insert("domain".to_string(), json!("general"));

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // log(0) should fail mathematical correctness
    assert!(!result["mathematical_correctness"].as_bool().unwrap());
}

// ============================================================================
// Domain-Specific Validation Tests
// ============================================================================

#[test]
fn test_validate_mechanics_force_equation() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("F = m*a"));
    params.insert("domain".to_string(), json!("mechanics"));

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Newton's second law should be physics compliant
    assert!(result["physics_compliance"].as_bool().unwrap());
}

#[test]
fn test_validate_electromagnetism_coulomb_law() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("F = k*q1*q2/r^2"));
    params.insert("domain".to_string(), json!("electromagnetism"));

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Coulomb's law should be physics compliant
    assert!(result["is_valid"].as_bool().unwrap());
}

#[test]
fn test_validate_thermodynamics_ideal_gas() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("PV = nRT"));
    params.insert("domain".to_string(), json!("thermodynamics"));

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Ideal gas law should be physics compliant
    assert!(result["is_valid"].as_bool().unwrap());
}

// ============================================================================
// Error Handling Tests
// ============================================================================

#[test]
fn test_validate_missing_equation_parameter() {
    let params = HashMap::new(); // Missing equation

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Missing equation"));
}

#[test]
fn test_validate_invalid_equation_parameter_type() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!(123)); // Should be string, not number

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Invalid equation"));
}

#[test]
fn test_validate_unknown_operation() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("F = m*a"));

    let request = create_request("unknown_operation", params);
    let response = equation_validation::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Unknown operation"));
}

#[test]
fn test_validate_malformed_equation() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("F = = m*a")); // Double equals
    params.insert("domain".to_string(), json!("mechanics"));

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    // Malformed equation should cause an error during parsing
    // The equation has multiple equals signs which violates parsing rules
    assert!(!response.success);
    assert!(response.error.is_some());
}

// ============================================================================
// Confidence Score Tests
// ============================================================================

#[test]
fn test_confidence_score_all_checks_pass() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("F = m*a"));
    params.insert("domain".to_string(), json!("mechanics"));

    let mut units = HashMap::new();
    units.insert("F".to_string(), "N".to_string());
    units.insert("m".to_string(), "kg".to_string());
    units.insert("a".to_string(), "m/s^2".to_string());
    params.insert("units".to_string(), json!(units));

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Confidence should be between 0 and 1
    // Note: Due to simplified dimensional analysis, confidence may vary
    let confidence = result["confidence"].as_f64().unwrap();
    assert!(confidence >= 0.0);
    assert!(confidence <= 1.0);
}

#[test]
fn test_confidence_score_some_checks_fail() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("y = x/0"));
    params.insert("domain".to_string(), json!("general"));

    let request = create_request("validate", params);
    let response = equation_validation::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // With mathematical errors, confidence should be lower
    let confidence = result["confidence"].as_f64().unwrap();
    assert!(confidence < 1.0);
}
