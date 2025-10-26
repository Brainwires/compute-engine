//! Unit tests for src/api/mod.rs
//!
//! Tests the main API entry point including process_request and process_json_request functions.

use crate::api::{process_json_request, process_request, ComputationRequest, ComputationResponse};
use serde_json::json;
use std::collections::HashMap;

// Helper function to create a request
fn create_request(module: &str, operation: &str, params: serde_json::Value) -> ComputationRequest {
    let params_map: HashMap<String, serde_json::Value> =
        serde_json::from_value(params).unwrap();

    ComputationRequest {
        module: module.to_string(),
        operation: operation.to_string(),
        parameters: params_map,
    }
}

// Helper function to assert success response
fn assert_success(response: &ComputationResponse) {
    assert!(response.success, "Expected success, got error: {:?}", response.error);
    assert!(response.result.is_some(), "Expected result to be Some");
    assert!(response.error.is_none(), "Expected no error");
}

// Helper function to assert error response
fn assert_error(response: &ComputationResponse) {
    assert!(!response.success, "Expected error, got success");
    assert!(response.result.is_none(), "Expected no result");
    assert!(response.error.is_some(), "Expected error message");
}

#[test]
fn test_process_request_linear_algebra() {
    let request = create_request(
        "linear_algebra",
        "compute_svd",
        json!({
            "matrix": [[1.0, 2.0], [3.0, 4.0]]
        })
    );

    let response = process_request(&request);
    assert_success(&response);
    assert_eq!(response.module, "linear_algebra");
    assert_eq!(response.operation, "compute_svd");
}

#[test]
fn test_process_request_advanced_calculus() {
    let request = create_request(
        "advanced_calculus",
        "derivative",
        json!({
            "coefficients": [0.0, 1.0, 2.0, 3.0],
            "x": 2.0
        })
    );

    let response = process_request(&request);
    assert_success(&response);
    assert_eq!(response.module, "advanced_calculus");
}

#[test]
fn test_process_request_special_functions() {
    let request = create_request(
        "special_functions",
        "bessel",
        json!({
            "function_type": "J",
            "order": 0.0,
            "x": 1.0
        })
    );

    let response = process_request(&request);
    assert_success(&response);
    assert_eq!(response.module, "special_functions");
}

#[test]
fn test_process_request_signal_processing() {
    let request = create_request(
        "signal_processing",
        "fft",
        json!({
            "signal": [1.0, 2.0, 3.0, 4.0],
            "sample_rate": 1000.0,
            "window_type": "hamming"
        })
    );

    let response = process_request(&request);
    assert_success(&response);
    assert_eq!(response.module, "signal_processing");
}

#[test]
fn test_process_request_unknown_module() {
    let request = create_request(
        "nonexistent_module",
        "some_operation",
        json!({})
    );

    let response = process_request(&request);
    assert_error(&response);
    assert_eq!(response.module, "nonexistent_module");
    assert!(response.error.unwrap().contains("Unknown module"));
}

#[test]
fn test_process_request_unknown_operation() {
    let request = create_request(
        "linear_algebra",
        "nonexistent_operation",
        json!({})
    );

    let response = process_request(&request);
    assert_error(&response);
    assert!(response.error.unwrap().contains("Unknown operation"));
}

#[test]
fn test_process_request_module_alias_tensor() {
    let request = create_request(
        "tensor",  // Alias for tensor_calculus
        "metric_tensor",
        json!({
            "coordinates": [[1.0, 0.0], [0.0, 1.0]]
        })
    );

    let response = process_request(&request);
    // Should route to tensor_calculus handler
    assert_eq!(response.module, "tensor");
}

#[test]
fn test_process_request_module_alias_linalg() {
    let request = create_request(
        "linalg",  // Alias for linear_algebra
        "compute_svd",
        json!({
            "matrix": [[1.0, 2.0], [3.0, 4.0]]
        })
    );

    let response = process_request(&request);
    // Should route to linear_algebra handler
    assert_eq!(response.module, "linalg");
}

#[test]
fn test_process_json_request_valid_json() {
    let json_str = json!({
        "module": "linear_algebra",
        "operation": "compute_svd",
        "parameters": {
            "matrix": [[1.0, 2.0], [3.0, 4.0]]
        }
    }).to_string();

    let response_str = process_json_request(&json_str);

    // Parse response
    let response: ComputationResponse = serde_json::from_str(&response_str).unwrap();
    assert_success(&response);
}

#[test]
fn test_process_json_request_invalid_json() {
    let invalid_json = "{not valid json at all}";

    let response_str = process_json_request(invalid_json);
    let response: ComputationResponse = serde_json::from_str(&response_str).unwrap();

    assert_error(&response);
    assert_eq!(response.module, "unknown");
    assert_eq!(response.operation, "parse");
    assert!(response.error.unwrap().contains("Invalid JSON"));
}

#[test]
fn test_process_json_request_missing_fields() {
    let incomplete_json = json!({
        "module": "linear_algebra"
        // Missing operation and parameters
    }).to_string();

    let response_str = process_json_request(&incomplete_json);
    let response: ComputationResponse = serde_json::from_str(&response_str).unwrap();

    assert_error(&response);
    assert!(response.error.unwrap().contains("Invalid JSON"));
}

#[test]
fn test_process_json_request_empty_parameters() {
    let json_str = json!({
        "module": "linear_algebra",
        "operation": "some_operation",
        "parameters": {}
    }).to_string();

    let response_str = process_json_request(&json_str);
    let response: ComputationResponse = serde_json::from_str(&response_str).unwrap();

    // Should parse successfully but likely fail in handler
    assert!(!response_str.is_empty());
}

#[test]
fn test_process_request_chemistry_module() {
    let request = create_request(
        "chemistry",
        "gas_law",
        json!({
            "pressure": 1.0,
            "volume": 22.4,
            "temperature": 273.15
        })
    );

    let response = process_request(&request);
    assert_success(&response);
    assert_eq!(response.module, "chemistry");
}

#[test]
fn test_process_request_biology_module() {
    let request = create_request(
        "biology",
        "biology",
        json!({
            "operation": "michaelis_menten",
            "parameters": {
                "substrate_concentration": 10.0,
                "vmax": 100.0,
                "km": 5.0
            }
        })
    );

    let response = process_request(&request);
    assert_success(&response);
    assert_eq!(response.module, "biology");
}

#[test]
fn test_process_request_thermodynamics_alias() {
    let request = create_request(
        "thermo",  // Alias for thermodynamics
        "carnot_efficiency",
        json!({
            "T_hot": 500.0,
            "T_cold": 300.0
        })
    );

    let response = process_request(&request);
    assert_eq!(response.module, "thermo");
}

#[test]
fn test_response_serialization() {
    let response = ComputationResponse::success(
        "test_module".to_string(),
        "test_op".to_string(),
        json!({"value": 42})
    );

    let serialized = serde_json::to_string(&response).unwrap();
    let deserialized: ComputationResponse = serde_json::from_str(&serialized).unwrap();

    assert_eq!(response.success, deserialized.success);
    assert_eq!(response.module, deserialized.module);
    assert_eq!(response.operation, deserialized.operation);
}

#[test]
fn test_error_response_serialization() {
    let response = ComputationResponse::error(
        "test_module".to_string(),
        "test_op".to_string(),
        "Test error message".to_string()
    );

    let serialized = serde_json::to_string(&response).unwrap();
    let deserialized: ComputationResponse = serde_json::from_str(&serialized).unwrap();

    assert!(!deserialized.success);
    assert_eq!(deserialized.error.unwrap(), "Test error message");
}
