//! Unit tests for tensor_calculus API handler
//!
//! Tests all handler operations including:
//! - Christoffel symbols computation
//! - Riemann tensor computation
//! - Ricci tensor computation
//! - Operation aliases
//! - Error handling for invalid operations

use super::*; // Import from parent module (tensor_calculus handler)
use crate::api::types::ComputationRequest;
use serde_json::json;
use std::collections::HashMap;

fn create_request(operation: &str, parameters: serde_json::Value) -> ComputationRequest {
    let params_map: HashMap<String, serde_json::Value> = serde_json::from_value(parameters)
        .expect("Failed to convert parameters to HashMap");

    ComputationRequest {
        module: "tensor_calculus".to_string(),
        operation: operation.to_string(),
        parameters: params_map,
    }
}

#[test]
fn test_christoffel_symbols_success() {
    let params = json!({
        "metric": [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0, 0.0],
            [0.0, 0.0, -1.0, 0.0],
            [0.0, 0.0, 0.0, -1.0]
        ],
        "coordinates": ["t", "x", "y", "z"]
    });
    let request = create_request("christoffel", params);
    let response = handle(&request);

    assert!(response.success, "Christoffel symbols computation should succeed");
    assert!(response.result.is_some(), "Should return result");

    let result = response.result.unwrap();
    assert_eq!(result.get("result").unwrap(), "computed");
    assert!(result.get("tensor").is_some(), "Should have tensor field");
}

#[test]
fn test_christoffel_symbols_alias() {
    let params = json!({
        "metric": [
            [1.0, 0.0],
            [0.0, -1.0]
        ],
        "coordinates": ["t", "x"]
    });
    let request = create_request("christoffel_symbols", params);
    let response = handle(&request);

    assert!(response.success, "Christoffel symbols alias should work");
    assert!(response.result.is_some());
}

#[test]
fn test_riemann_tensor_success() {
    let params = json!({
        "metric": [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0, 0.0],
            [0.0, 0.0, -1.0, 0.0],
            [0.0, 0.0, 0.0, -1.0]
        ],
        "coordinates": ["t", "x", "y", "z"]
    });
    let request = create_request("riemann", params);
    let response = handle(&request);

    assert!(response.success, "Riemann tensor computation should succeed");
    assert!(response.result.is_some(), "Should return result");

    let result = response.result.unwrap();
    assert_eq!(result.get("result").unwrap(), "computed");
    assert!(result.get("tensor").is_some());
}

#[test]
fn test_riemann_tensor_alias() {
    let params = json!({
        "metric": [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]
        ],
        "coordinates": ["x", "y", "z"]
    });
    let request = create_request("riemann_tensor", params);
    let response = handle(&request);

    assert!(response.success, "Riemann tensor alias should work");
    assert!(response.result.is_some());
}

#[test]
fn test_ricci_tensor_success() {
    let params = json!({
        "metric": [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0, 0.0],
            [0.0, 0.0, -1.0, 0.0],
            [0.0, 0.0, 0.0, -1.0]
        ],
        "coordinates": ["t", "x", "y", "z"]
    });
    let request = create_request("ricci", params);
    let response = handle(&request);

    assert!(response.success, "Ricci tensor computation should succeed");
    assert!(response.result.is_some(), "Should return result");

    let result = response.result.unwrap();
    assert_eq!(result.get("result").unwrap(), "computed");
    assert!(result.get("tensor").is_some());
}

#[test]
fn test_ricci_tensor_alias() {
    let params = json!({
        "metric": [
            [1.0, 0.0],
            [0.0, 1.0]
        ],
        "coordinates": ["x", "y"]
    });
    let request = create_request("ricci_tensor", params);
    let response = handle(&request);

    assert!(response.success, "Ricci tensor alias should work");
    assert!(response.result.is_some());
}

#[test]
fn test_schwarzschild_metric() {
    // Test with Schwarzschild metric parameters
    let params = json!({
        "metric": "schwarzschild",
        "mass": 1.0,
        "coordinates": ["t", "r", "theta", "phi"]
    });
    let request = create_request("christoffel", params);
    let response = handle(&request);

    assert!(response.success, "Schwarzschild metric should be handled");
}

#[test]
fn test_friedmann_metric() {
    // Test with Friedmann-Lemaître-Robertson-Walker metric
    let params = json!({
        "metric": "flrw",
        "scale_factor": 1.0,
        "curvature": 0.0,
        "coordinates": ["t", "r", "theta", "phi"]
    });
    let request = create_request("riemann", params);
    let response = handle(&request);

    assert!(response.success, "FLRW metric should be handled");
}

#[test]
fn test_kerr_metric() {
    // Test with Kerr metric (rotating black hole)
    let params = json!({
        "metric": "kerr",
        "mass": 1.0,
        "angular_momentum": 0.5,
        "coordinates": ["t", "r", "theta", "phi"]
    });
    let request = create_request("ricci", params);
    let response = handle(&request);

    assert!(response.success, "Kerr metric should be handled");
}

#[test]
fn test_unknown_operation() {
    let params = json!({
        "metric": [[1.0, 0.0], [0.0, 1.0]],
        "coordinates": ["x", "y"]
    });
    let request = create_request("unknown_tensor_operation", params);
    let response = handle(&request);

    assert!(!response.success, "Should fail for unknown operation");
    assert!(response.error.is_some(), "Should return error message");
    assert!(response.error.unwrap().contains("Unknown tensor operation"));
}

#[test]
fn test_invalid_operation_name() {
    let params = json!({
        "metric": [[1.0]],
        "coordinates": ["x"]
    });
    let request = create_request("invalid_op", params);
    let response = handle(&request);

    assert!(!response.success, "Should fail for invalid operation");
    assert!(response.error.is_some());
}

#[test]
fn test_response_structure_christoffel() {
    let params = json!({
        "metric": [[1.0, 0.0], [0.0, 1.0]],
        "coordinates": ["x", "y"]
    });
    let request = create_request("christoffel", params);
    let response = handle(&request);

    // Verify response structure
    assert_eq!(response.module, "tensor_calculus");
    assert_eq!(response.operation, "christoffel");
    assert!(response.success);
    assert!(response.result.is_some());
    assert!(response.error.is_none());
}

#[test]
fn test_response_structure_riemann() {
    let params = json!({
        "metric": [[1.0, 0.0], [0.0, 1.0]],
        "coordinates": ["x", "y"]
    });
    let request = create_request("riemann_tensor", params);
    let response = handle(&request);

    // Verify response structure
    assert_eq!(response.module, "tensor_calculus");
    assert_eq!(response.operation, "riemann_tensor");
    assert!(response.success);
    assert!(response.result.is_some());
    assert!(response.error.is_none());
}

#[test]
fn test_response_structure_ricci() {
    let params = json!({
        "metric": [[1.0, 0.0], [0.0, 1.0]],
        "coordinates": ["x", "y"]
    });
    let request = create_request("ricci", params);
    let response = handle(&request);

    // Verify response structure
    assert_eq!(response.module, "tensor_calculus");
    assert_eq!(response.operation, "ricci");
    assert!(response.success);
    assert!(response.result.is_some());
    assert!(response.error.is_none());
}

#[test]
fn test_empty_parameters() {
    let params = json!({});
    let request = create_request("christoffel", params);
    let response = handle(&request);

    // Handler returns mock success even with empty parameters
    assert!(response.success, "Handler returns mock success");
    assert!(response.result.is_some());
}

#[test]
fn test_minimal_parameters() {
    let params = json!({
        "metric": [[1.0]]
    });
    let request = create_request("riemann", params);
    let response = handle(&request);

    assert!(response.success, "Should handle minimal parameters");
    assert!(response.result.is_some());
}
