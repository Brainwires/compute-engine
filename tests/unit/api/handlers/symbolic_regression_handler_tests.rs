//! Unit tests for symbolic regression API handler
//!
//! Tests all symbolic regression operations including:
//! - Equation discovery (fit/discover operations)
//! - Parameter validation (x, y data arrays)
//! - Optional parameters (max_complexity, domain, units)
//! - Response structure validation
//! - Error handling for invalid inputs
//! - Physics-informed constraints
//! - Multiple domain types

use crate::api::handlers::symbolic_regression;
use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};
use std::collections::HashMap;

/// Helper function to create a symbolic regression request
fn create_request(operation: &str, parameters: HashMap<String, Value>) -> ComputationRequest {
    ComputationRequest {
        module: "symbolic_regression".to_string(),
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
// Fit/Discover Operation Tests
// ============================================================================

#[test]
fn test_fit_linear_data_success() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0, 4.0, 5.0]));
    params.insert("y".to_string(), json!([2.0, 4.0, 6.0, 8.0, 10.0]));

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Verify result structure
    assert!(result["candidates"].is_array());
    assert!(result["constraints"].is_object());
    assert!(result["best_fitness"].is_number());
    assert!(result["convergence_data"].is_array());
}

#[test]
fn test_discover_quadratic_data_success() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0, 4.0, 5.0]));
    params.insert("y".to_string(), json!([1.0, 4.0, 9.0, 16.0, 25.0]));

    let request = create_request("discover", params);
    let response = symbolic_regression::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert!(result["candidates"].is_array());
    assert!(result["best_fitness"].is_number());
}

#[test]
fn test_fit_candidates_structure() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert("y".to_string(), json!([1.0, 2.0, 3.0]));

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let candidates = result["candidates"].as_array().unwrap();
    if !candidates.is_empty() {
        let first_candidate = &candidates[0];
        assert!(first_candidate["expression"].is_string());
        assert!(first_candidate["complexity"].is_number());
        assert!(first_candidate["fitness"].is_number());
        assert!(first_candidate["r_squared"].is_number());
        assert!(first_candidate["mse"].is_number());
        assert!(first_candidate["physical_validity"].is_number());
    }
}

#[test]
fn test_fit_constraints_structure() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert("y".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert("domain".to_string(), json!("mechanics"));

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let constraints = &result["constraints"];
    assert!(constraints["units"].is_object());
    assert!(constraints["conservation_laws"].is_array());
    assert!(constraints["symmetries"].is_array());
    assert!(constraints["dimensional_consistency"].is_boolean());
    assert!(constraints["physical_plausibility"].is_boolean());
}

#[test]
fn test_fit_convergence_data_structure() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert("y".to_string(), json!([1.0, 2.0, 3.0]));

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let convergence_data = result["convergence_data"].as_array().unwrap();
    assert!(!convergence_data.is_empty());
    // All convergence values should be numbers
    for val in convergence_data {
        assert!(val.is_number());
    }
}

// ============================================================================
// Parameter Validation Tests
// ============================================================================

#[test]
fn test_missing_x_parameter_error() {
    let mut params = HashMap::new();
    params.insert("y".to_string(), json!([1.0, 2.0, 3.0]));

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Missing x parameter"));
}

#[test]
fn test_missing_y_parameter_error() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0]));

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Missing y parameter"));
}

#[test]
fn test_invalid_x_data_type_error() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!("not_an_array"));
    params.insert("y".to_string(), json!([1.0, 2.0, 3.0]));

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Invalid x data"));
}

#[test]
fn test_invalid_y_data_type_error() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert("y".to_string(), json!("not_an_array"));

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Invalid y data"));
}

#[test]
fn test_empty_x_array() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([]));
    params.insert("y".to_string(), json!([1.0]));

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    // Should handle gracefully (may succeed with no candidates or fail)
    // The implementation should handle this case
    assert!(response.success || response.error.is_some());
}

#[test]
fn test_empty_y_array() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0]));
    params.insert("y".to_string(), json!([]));

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    // Should handle gracefully
    assert!(response.success || response.error.is_some());
}

// ============================================================================
// Optional Parameters Tests
// ============================================================================

#[test]
fn test_custom_max_complexity() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert("y".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert("max_complexity".to_string(), json!(3));

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Verify candidates respect complexity limit
    let candidates = result["candidates"].as_array().unwrap();
    for candidate in candidates {
        let complexity = candidate["complexity"].as_u64().unwrap();
        assert!(complexity <= 3);
    }
}

#[test]
fn test_default_max_complexity() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert("y".to_string(), json!([1.0, 2.0, 3.0]));
    // No max_complexity parameter - should default to 5

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let candidates = result["candidates"].as_array().unwrap();
    for candidate in candidates {
        let complexity = candidate["complexity"].as_u64().unwrap();
        assert!(complexity <= 5); // Default is 5
    }
}

#[test]
fn test_mechanics_domain() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0, 4.0]));
    params.insert("y".to_string(), json!([1.0, 4.0, 9.0, 16.0]));
    params.insert("domain".to_string(), json!("mechanics"));

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let constraints = &result["constraints"];
    let conservation_laws = constraints["conservation_laws"].as_array().unwrap();
    // Mechanics domain should include energy and momentum conservation
    assert!(conservation_laws.len() >= 2);
}

#[test]
fn test_electromagnetism_domain() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert("y".to_string(), json!([1.0, 0.25, 0.111]));
    params.insert("domain".to_string(), json!("electromagnetism"));

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let constraints = &result["constraints"];
    let conservation_laws = constraints["conservation_laws"].as_array().unwrap();
    // EM domain should include charge and energy conservation
    assert!(conservation_laws.len() >= 2);
}

#[test]
fn test_general_domain_default() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert("y".to_string(), json!([1.0, 2.0, 3.0]));
    // No domain parameter - should default to "general"

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let constraints = &result["constraints"];
    assert!(constraints["conservation_laws"].is_array());
}

#[test]
fn test_units_parameter() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert("y".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert(
        "units".to_string(),
        json!({"x": "meters", "y": "seconds"}),
    );

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let constraints = &result["constraints"];
    let units = constraints["units"].as_object().unwrap();
    assert_eq!(units.get("x").unwrap().as_str().unwrap(), "meters");
    assert_eq!(units.get("y").unwrap().as_str().unwrap(), "seconds");
}

// ============================================================================
// Unknown Operation Tests
// ============================================================================

#[test]
fn test_unknown_operation_error() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert("y".to_string(), json!([1.0, 2.0, 3.0]));

    let request = create_request("unknown_operation", params);
    let response = symbolic_regression::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response
        .error
        .unwrap()
        .contains("Unknown operation: unknown_operation"));
}

// ============================================================================
// Response Structure Tests
// ============================================================================

#[test]
fn test_response_module_field() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert("y".to_string(), json!([1.0, 2.0, 3.0]));

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    assert_eq!(response.module, "symbolic_regression");
}

#[test]
fn test_response_operation_field() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert("y".to_string(), json!([1.0, 2.0, 3.0]));

    let request = create_request("discover", params);
    let response = symbolic_regression::handle(&request);

    assert_eq!(response.operation, "discover");
}

#[test]
fn test_best_fitness_non_negative() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert("y".to_string(), json!([1.0, 2.0, 3.0]));

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let best_fitness = result["best_fitness"].as_f64().unwrap();
    assert!(best_fitness >= 0.0);
}

#[test]
fn test_physical_validity_in_range() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!([1.0, 2.0, 3.0, 4.0, 5.0]));
    params.insert("y".to_string(), json!([2.0, 4.0, 6.0, 8.0, 10.0]));
    params.insert("domain".to_string(), json!("mechanics"));

    let request = create_request("fit", params);
    let response = symbolic_regression::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let candidates = result["candidates"].as_array().unwrap();
    for candidate in candidates {
        let physical_validity = candidate["physical_validity"].as_f64().unwrap();
        assert!(physical_validity >= 0.0 && physical_validity <= 2.0);
    }
}
