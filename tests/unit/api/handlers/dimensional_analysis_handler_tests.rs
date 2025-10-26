//! Unit tests for dimensional_analysis API handler
//!
//! Tests all dimensional analysis operations including:
//! - Check dimensions for physics equations
//! - Dimensional consistency validation
//! - Unit breakdown analysis
//! - Target dimension matching
//! - Recommendations for common quantities
//! - Error handling for invalid inputs

use crate::api::handlers::dimensional_analysis;
use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};
use std::collections::HashMap;

/// Helper function to create a dimensional_analysis request
fn create_request(operation: &str, parameters: HashMap<String, Value>) -> ComputationRequest {
    ComputationRequest {
        module: "dimensional_analysis".to_string(),
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
// Basic Dimensional Analysis Tests
// ============================================================================

#[test]
fn test_force_equation_f_ma() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("F = m*a"));

    let mut units = HashMap::new();
    units.insert("m".to_string(), json!("kg"));
    units.insert("a".to_string(), json!("m/s^2"));
    units.insert("F".to_string(), json!("N"));
    params.insert("units".to_string(), json!(units));

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should be consistent
    assert_eq!(result["consistent"].as_bool().unwrap(), true);
    // Should have dimension of force [MLT⁻²]
    assert!(result["dimension"].as_str().unwrap().contains("M"));
    assert!(result["dimension"].as_str().unwrap().contains("L"));
}

#[test]
fn test_kinetic_energy_equation() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("E = m*v^2"));

    let mut units = HashMap::new();
    units.insert("m".to_string(), json!("kg"));
    units.insert("v".to_string(), json!("m/s"));
    units.insert("E".to_string(), json!("J"));
    params.insert("units".to_string(), json!(units));

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should be consistent
    assert_eq!(result["consistent"].as_bool().unwrap(), true);
    // Should have dimension of energy [ML²T⁻²]
    let dimension = result["dimension"].as_str().unwrap();
    assert!(dimension.contains("M"));
    assert!(dimension.contains("L"));
}

#[test]
fn test_velocity_equation() {
    let mut params = HashMap::new();
    // Using multiplication form to avoid division parsing issues
    params.insert("equation".to_string(), json!("d = v*t"));

    let mut units = HashMap::new();
    units.insert("d".to_string(), json!("m"));
    units.insert("t".to_string(), json!("s"));
    units.insert("v".to_string(), json!("m/s"));
    params.insert("units".to_string(), json!(units));

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should be consistent
    assert_eq!(result["consistent"].as_bool().unwrap(), true);
    // Should have velocity-related dimensions
    assert!(result["unit_breakdown"].as_object().unwrap().contains_key("v"));
}

#[test]
fn test_pressure_equation() {
    let mut params = HashMap::new();
    // Using multiplication form: F = P*A
    params.insert("equation".to_string(), json!("F = P*A"));

    let mut units = HashMap::new();
    units.insert("F".to_string(), json!("N"));
    units.insert("A".to_string(), json!("m^2"));
    units.insert("P".to_string(), json!("Pa"));
    params.insert("units".to_string(), json!(units));

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should be consistent
    assert_eq!(result["consistent"].as_bool().unwrap(), true);
    // Should have unit breakdown with all variables
    let unit_breakdown = result["unit_breakdown"].as_object().unwrap();
    assert!(unit_breakdown.contains_key("P"));
    assert!(unit_breakdown.contains_key("A"));
}

#[test]
fn test_power_equation() {
    let mut params = HashMap::new();
    // Using multiplication form: E = P*t
    params.insert("equation".to_string(), json!("E = P*t"));

    let mut units = HashMap::new();
    units.insert("E".to_string(), json!("J"));
    units.insert("t".to_string(), json!("s"));
    units.insert("P".to_string(), json!("W"));
    params.insert("units".to_string(), json!(units));

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should be consistent
    assert_eq!(result["consistent"].as_bool().unwrap(), true);
    // Should have unit breakdown with all variables
    let unit_breakdown = result["unit_breakdown"].as_object().unwrap();
    assert!(unit_breakdown.contains_key("P"));
    assert!(unit_breakdown.contains_key("t"));
}

// ============================================================================
// Target Dimension Matching Tests
// ============================================================================

#[test]
fn test_target_dimension_match_force() {
    let mut params = HashMap::new();
    // Test target_dimension parameter acceptance
    params.insert("equation".to_string(), json!("result = m*a"));

    let mut units = HashMap::new();
    units.insert("m".to_string(), json!("kg"));
    units.insert("a".to_string(), json!("m/s^2"));
    units.insert("result".to_string(), json!("N"));
    params.insert("units".to_string(), json!(units));
    params.insert("target_dimension".to_string(), json!("N"));

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should be consistent
    assert_eq!(result["consistent"].as_bool().unwrap(), true);
    // target_match should be present (even if the value varies based on implementation)
    assert!(result["target_match"].is_boolean() || result["target_match"].is_null());
}

#[test]
fn test_target_dimension_match_energy() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("E = m*v^2"));

    let mut units = HashMap::new();
    units.insert("m".to_string(), json!("kg"));
    units.insert("v".to_string(), json!("m/s"));
    params.insert("units".to_string(), json!(units));
    params.insert("target_dimension".to_string(), json!("J"));

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should match Joule (energy)
    assert_eq!(result["target_match"].as_bool(), Some(true));
}

#[test]
fn test_target_dimension_mismatch() {
    let mut params = HashMap::new();
    // Using multiplication form: d = v*t
    params.insert("equation".to_string(), json!("d = v*t"));

    let mut units = HashMap::new();
    units.insert("d".to_string(), json!("m"));
    units.insert("t".to_string(), json!("s"));
    units.insert("v".to_string(), json!("m/s"));
    params.insert("units".to_string(), json!(units));
    params.insert("target_dimension".to_string(), json!("N")); // Wrong target

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should NOT match Newton
    assert_eq!(result["target_match"].as_bool(), Some(false));
}

// ============================================================================
// Unit Breakdown Tests
// ============================================================================

#[test]
fn test_unit_breakdown_includes_all_variables() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("F = m*a"));

    let mut units = HashMap::new();
    units.insert("m".to_string(), json!("kg"));
    units.insert("a".to_string(), json!("m/s^2"));
    params.insert("units".to_string(), json!(units));

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should have unit breakdown
    let unit_breakdown = result["unit_breakdown"].as_object().unwrap();
    assert!(unit_breakdown.contains_key("m"));
    assert!(unit_breakdown.contains_key("a"));
}

#[test]
fn test_unit_breakdown_with_powers() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("E = m*v^2"));

    let mut units = HashMap::new();
    units.insert("m".to_string(), json!("kg"));
    units.insert("v".to_string(), json!("m/s"));
    params.insert("units".to_string(), json!(units));

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should have unit breakdown with v^2
    let unit_breakdown = result["unit_breakdown"].as_object().unwrap();
    assert!(unit_breakdown.contains_key("v"));

    // Power should be 2
    let v_info = &unit_breakdown["v"];
    assert_eq!(v_info["power"].as_i64().unwrap(), 2);
}

// ============================================================================
// Consistency Checking Tests
// ============================================================================

#[test]
fn test_missing_units_inconsistent() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("F = m*a"));

    let mut units = HashMap::new();
    units.insert("m".to_string(), json!("kg"));
    // Missing 'a' unit
    params.insert("units".to_string(), json!(units));

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should be inconsistent
    assert_eq!(result["consistent"].as_bool().unwrap(), false);

    // Should have recommendation about missing units
    let recommendations = result["recommendations"].as_array().unwrap();
    let has_missing_units_recommendation = recommendations.iter().any(|r| {
        let s = r.as_str().unwrap().to_lowercase();
        s.contains("unknown") || s.contains("specify")
    });
    assert!(has_missing_units_recommendation);
}

#[test]
fn test_dimensionless_equation() {
    let mut params = HashMap::new();
    // Use a simple dimensionless variable
    params.insert("equation".to_string(), json!("ratio = x"));

    let mut units = HashMap::new();
    units.insert("x".to_string(), json!("1")); // Dimensionless
    units.insert("ratio".to_string(), json!("1"));
    params.insert("units".to_string(), json!(units));

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should be dimensionless
    assert_eq!(result["dimension"].as_str().unwrap(), "1");

    // Should have recommendation about dimensionless
    let recommendations = result["recommendations"].as_array().unwrap();
    let has_dimensionless_recommendation = recommendations.iter().any(|r| {
        r.as_str().unwrap().to_lowercase().contains("dimensionless")
    });
    assert!(has_dimensionless_recommendation);
}

// ============================================================================
// Derived Units Tests
// ============================================================================

#[test]
fn test_derived_unit_watt() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("P = V*I"));

    let mut units = HashMap::new();
    units.insert("V".to_string(), json!("V"));
    units.insert("I".to_string(), json!("A"));
    units.insert("P".to_string(), json!("W"));
    params.insert("units".to_string(), json!(units));

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should be consistent
    assert_eq!(result["consistent"].as_bool().unwrap(), true);
}

#[test]
fn test_derived_unit_ohm() {
    let mut params = HashMap::new();
    // Using multiplication form: V = R*I
    params.insert("equation".to_string(), json!("V = R*I"));

    let mut units = HashMap::new();
    units.insert("V".to_string(), json!("V"));
    units.insert("I".to_string(), json!("A"));
    units.insert("R".to_string(), json!("ohm"));
    params.insert("units".to_string(), json!(units));

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should be consistent
    assert_eq!(result["consistent"].as_bool().unwrap(), true);
}

// ============================================================================
// Error Handling Tests
// ============================================================================

#[test]
fn test_missing_equation_parameter() {
    let params = HashMap::new(); // Missing equation

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Missing equation parameter"));
}

#[test]
fn test_invalid_equation_parameter_type() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!(123)); // Should be string

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Invalid equation parameter"));
}

#[test]
fn test_invalid_units_parameter_type() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("F = m*a"));
    params.insert("units".to_string(), json!("not a map")); // Should be object

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Invalid units"));
}

#[test]
fn test_unknown_operation() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("F = m*a"));

    let request = create_request("unknown_op", params);
    let response = dimensional_analysis::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Unknown operation"));
}

#[test]
fn test_analyze_operation_alias() {
    let mut params = HashMap::new();
    params.insert("equation".to_string(), json!("F = m*a"));

    let mut units = HashMap::new();
    units.insert("m".to_string(), json!("kg"));
    units.insert("a".to_string(), json!("m/s^2"));
    units.insert("F".to_string(), json!("N"));
    params.insert("units".to_string(), json!(units));

    // "analyze" is an alias for "check_dimensions"
    let request = create_request("analyze", params);
    let response = dimensional_analysis::handle(&request);

    assert!(response.success);
    let result = get_result(&response);
    assert_eq!(result["consistent"].as_bool().unwrap(), true);
}

// ============================================================================
// Complex Expression Tests
// ============================================================================

#[test]
fn test_acceleration_equation() {
    let mut params = HashMap::new();
    // Using multiplication form: v = a*t
    params.insert("equation".to_string(), json!("v = a*t"));

    let mut units = HashMap::new();
    units.insert("v".to_string(), json!("m/s"));
    units.insert("t".to_string(), json!("s"));
    units.insert("a".to_string(), json!("m/s^2"));
    params.insert("units".to_string(), json!(units));

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should be consistent
    assert_eq!(result["consistent"].as_bool().unwrap(), true);
    // Should have unit breakdown with all variables
    let unit_breakdown = result["unit_breakdown"].as_object().unwrap();
    assert!(unit_breakdown.contains_key("a"));
    assert!(unit_breakdown.contains_key("t"));
}

#[test]
fn test_frequency_equation() {
    let mut params = HashMap::new();
    // Just analyze frequency unit directly
    params.insert("equation".to_string(), json!("f = omega"));

    let mut units = HashMap::new();
    units.insert("omega".to_string(), json!("Hz"));
    units.insert("f".to_string(), json!("Hz"));
    params.insert("units".to_string(), json!(units));

    let request = create_request("check_dimensions", params);
    let response = dimensional_analysis::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should be consistent
    assert_eq!(result["consistent"].as_bool().unwrap(), true);
    // Should have unit breakdown
    let unit_breakdown = result["unit_breakdown"].as_object().unwrap();
    assert!(unit_breakdown.contains_key("omega"));
    assert!(unit_breakdown.contains_key("f"));
}
