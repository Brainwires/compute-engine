//! Unit tests for quantum_physics API handler
//!
//! Tests all quantum physics operations including:
//! - Wavefunction calculations
//! - Quantum operator operations
//! - Quantum entanglement
//! - Error handling for invalid inputs and unknown operations

use crate::api::handlers::quantum_physics;
use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};
use std::collections::HashMap;

/// Helper function to create a quantum_physics request
fn create_request(operation: &str, parameters: HashMap<String, Value>) -> ComputationRequest {
    ComputationRequest {
        module: "quantum_physics".to_string(),
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
// Wavefunction Tests
// ============================================================================

#[test]
fn test_wavefunction_basic() {
    let mut params = HashMap::new();
    params.insert("state".to_string(), json!("psi"));
    params.insert("position".to_string(), json!(0.5));

    let request = create_request("wavefunction", params);
    let response = quantum_physics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Mock handler returns computed result
    assert_eq!(result["result"], "computed");
    assert_eq!(result["value"], 0.0);
}

#[test]
fn test_wavefunction_with_time() {
    let mut params = HashMap::new();
    params.insert("state".to_string(), json!("psi"));
    params.insert("position".to_string(), json!(1.0));
    params.insert("time".to_string(), json!(0.1));

    let request = create_request("wavefunction", params);
    let response = quantum_physics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["result"], "computed");
}

#[test]
fn test_wavefunction_minimal_params() {
    let params = HashMap::new(); // Empty parameters should still work (mock handler)

    let request = create_request("wavefunction", params);
    let response = quantum_physics::handle(&request);

    // Mock implementation accepts any parameters
    assert!(response.success);
}

#[test]
fn test_wavefunction_response_structure() {
    let params = HashMap::new();
    let request = create_request("wavefunction", params);
    let response = quantum_physics::handle(&request);

    // Verify response structure
    assert!(response.success);
    assert_eq!(response.module, "quantum_physics");
    assert_eq!(response.operation, "wavefunction");
    assert!(response.result.is_some());
    assert!(response.error.is_none());
}

// ============================================================================
// Operator Tests
// ============================================================================

#[test]
fn test_operator_single() {
    let mut params = HashMap::new();
    params.insert("operator_type".to_string(), json!("hamiltonian"));
    params.insert("state".to_string(), json!("psi"));

    let request = create_request("operator", params);
    let response = quantum_physics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["result"], "computed");
    assert_eq!(result["value"], 0.0);
}

#[test]
fn test_operator_momentum() {
    let mut params = HashMap::new();
    params.insert("operator_type".to_string(), json!("momentum"));
    params.insert("dimension".to_string(), json!(3));

    let request = create_request("operator", params);
    let response = quantum_physics::handle(&request);

    assert!(response.success);
}

#[test]
fn test_operator_position() {
    let mut params = HashMap::new();
    params.insert("operator_type".to_string(), json!("position"));
    params.insert("coordinate".to_string(), json!("x"));

    let request = create_request("operator", params);
    let response = quantum_physics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["result"], "computed");
}

// ============================================================================
// Operators (Plural) Tests
// ============================================================================

#[test]
fn test_operators_multiple() {
    let mut params = HashMap::new();
    params.insert(
        "operators".to_string(),
        json!(["hamiltonian", "momentum", "position"]),
    );
    params.insert("state".to_string(), json!("psi"));

    let request = create_request("operators", params);
    let response = quantum_physics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["result"], "computed");
}

#[test]
fn test_operators_commutator() {
    let mut params = HashMap::new();
    params.insert("operators".to_string(), json!(["position", "momentum"]));
    params.insert("operation".to_string(), json!("commutator"));

    let request = create_request("operators", params);
    let response = quantum_physics::handle(&request);

    assert!(response.success);
}

#[test]
fn test_operators_empty_list() {
    let mut params = HashMap::new();
    params.insert("operators".to_string(), json!([]));

    let request = create_request("operators", params);
    let response = quantum_physics::handle(&request);

    // Mock implementation accepts any parameters
    assert!(response.success);
}

// ============================================================================
// Entanglement Tests
// ============================================================================

#[test]
fn test_entanglement_two_particles() {
    let mut params = HashMap::new();
    params.insert("particles".to_string(), json!(2));
    params.insert("state_type".to_string(), json!("bell"));

    let request = create_request("entanglement", params);
    let response = quantum_physics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["result"], "computed");
    assert_eq!(result["value"], 0.0);
}

#[test]
fn test_entanglement_bell_state() {
    let mut params = HashMap::new();
    params.insert("particles".to_string(), json!(2));
    params.insert("state_type".to_string(), json!("bell"));
    params.insert("bell_state".to_string(), json!("phi_plus"));

    let request = create_request("entanglement", params);
    let response = quantum_physics::handle(&request);

    assert!(response.success);
}

#[test]
fn test_entanglement_measure() {
    let mut params = HashMap::new();
    params.insert("particles".to_string(), json!(2));
    params.insert("measurement".to_string(), json!("entanglement_entropy"));

    let request = create_request("entanglement", params);
    let response = quantum_physics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert_eq!(result["result"], "computed");
}

#[test]
fn test_entanglement_ghz_state() {
    let mut params = HashMap::new();
    params.insert("particles".to_string(), json!(3));
    params.insert("state_type".to_string(), json!("ghz"));

    let request = create_request("entanglement", params);
    let response = quantum_physics::handle(&request);

    assert!(response.success);
}

// ============================================================================
// Error Handling Tests
// ============================================================================

#[test]
fn test_unknown_operation() {
    let params = HashMap::new();
    let request = create_request("invalid_operation", params);
    let response = quantum_physics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response
        .error
        .as_ref()
        .unwrap()
        .contains("Unknown quantum operation"));
}

#[test]
fn test_unknown_operation_with_params() {
    let mut params = HashMap::new();
    params.insert("some_param".to_string(), json!(42));

    let request = create_request("nonexistent", params);
    let response = quantum_physics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
}

#[test]
fn test_error_response_structure() {
    let params = HashMap::new();
    let request = create_request("bad_operation", params);
    let response = quantum_physics::handle(&request);

    // Verify error response structure
    assert!(!response.success);
    assert_eq!(response.module, "quantum_physics");
    assert_eq!(response.operation, "bad_operation");
    assert!(response.result.is_none());
    assert!(response.error.is_some());
}

// ============================================================================
// Integration Tests
// ============================================================================

#[test]
fn test_all_supported_operations() {
    let operations = vec!["wavefunction", "operator", "operators", "entanglement"];

    for op in operations {
        let params = HashMap::new();
        let request = create_request(op, params);
        let response = quantum_physics::handle(&request);

        assert!(
            response.success,
            "Operation '{}' should succeed",
            op
        );
        assert_eq!(response.module, "quantum_physics");
        assert_eq!(response.operation, op);
    }
}

#[test]
fn test_case_sensitivity() {
    // Test that operations are case-sensitive
    let params = HashMap::new();

    let request = create_request("Wavefunction", params.clone());
    let response = quantum_physics::handle(&request);

    // Should fail because operation names are case-sensitive
    assert!(!response.success);
}
