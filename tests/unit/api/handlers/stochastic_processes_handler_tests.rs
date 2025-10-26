//! Unit tests for stochastic_processes API handler

use super::handle; // handle is in the parent module (stochastic_processes)
use crate::api::types::ComputationRequest;
use serde_json::json;
use std::collections::HashMap;

// Helper function to create a request
fn create_request(operation: &str, parameters: serde_json::Value) -> ComputationRequest {
    let params_map = if let serde_json::Value::Object(map) = parameters {
        map.into_iter().collect()
    } else {
        HashMap::new()
    };

    ComputationRequest {
        module: "stochastic_processes".to_string(),
        operation: operation.to_string(),
        parameters: params_map,
    }
}

// ============================================================================
// Brownian Motion Tests
// ============================================================================

#[test]
fn test_brownian_motion_success() {
    let request = create_request(
        "brownian_motion",
        json!({
            "time_steps": 100,
            "dt": 0.01,
            "initial_value": 0.0,
            "drift": 0.05,
            "volatility": 0.2
        }),
    );
    let response = handle(&request);

    assert!(response.success, "Brownian motion should succeed");
    assert!(response.error.is_none(), "Should have no error");

    let result = response.result.unwrap();
    assert!(result.get("path").is_some(), "Should have path");

    let path = result["path"].as_array().unwrap();
    assert_eq!(path.len(), 101, "Path should have 101 points (100 steps + initial)");

    // Verify initial value
    assert_eq!(path[0].as_f64().unwrap(), 0.0, "Initial value should be 0.0");
}

#[test]
fn test_brownian_motion_non_zero_initial() {
    let request = create_request(
        "brownian_motion",
        json!({
            "time_steps": 50,
            "dt": 0.01,
            "initial_value": 10.0,
            "drift": 0.0,
            "volatility": 1.0
        }),
    );
    let response = handle(&request);

    assert!(response.success, "Should succeed with non-zero initial value");

    let result = response.result.unwrap();
    let path = result["path"].as_array().unwrap();
    assert_eq!(path[0].as_f64().unwrap(), 10.0, "Initial value should be 10.0");
}

#[test]
fn test_brownian_motion_missing_parameters() {
    let request = create_request("brownian_motion", json!({}));
    let response = handle(&request);

    assert!(!response.success, "Should fail with missing parameters");
    assert!(response.error.is_some(), "Should have error message");
    assert!(
        response.error.unwrap().contains("Invalid parameters"),
        "Error should mention invalid parameters"
    );
}

#[test]
fn test_brownian_motion_invalid_parameter_type() {
    let request = create_request(
        "brownian_motion",
        json!({
            "time_steps": "invalid",
            "dt": 0.01,
            "initial_value": 0.0,
            "drift": 0.0,
            "volatility": 1.0
        }),
    );
    let response = handle(&request);

    assert!(!response.success, "Should fail with invalid parameter type");
    assert!(response.error.is_some(), "Should have error message");
}

#[test]
fn test_brownian_motion_small_dt() {
    let request = create_request(
        "brownian_motion",
        json!({
            "time_steps": 10,
            "dt": 0.001,
            "initial_value": 0.0,
            "drift": 0.0,
            "volatility": 0.1
        }),
    );
    let response = handle(&request);

    assert!(response.success, "Should succeed with small dt");
}

// ============================================================================
// Markov Chain Tests
// ============================================================================

#[test]
fn test_markov_chain_success() {
    let request = create_request(
        "markov_chain",
        json!({
            "states": ["A", "B", "C"],
            "transition_matrix": [
                [0.7, 0.2, 0.1],
                [0.3, 0.4, 0.3],
                [0.1, 0.3, 0.6]
            ],
            "initial_state": "A",
            "steps": 10
        }),
    );
    let response = handle(&request);

    assert!(response.success, "Markov chain should succeed");
    assert!(response.error.is_none(), "Should have no error");

    let result = response.result.unwrap();
    assert!(result.get("chain").is_some(), "Should have chain");

    let chain = result["chain"].as_array().unwrap();
    assert_eq!(chain.len(), 11, "Chain should have 11 states (10 steps + initial)");

    // Verify initial state
    assert_eq!(chain[0].as_str().unwrap(), "A", "Initial state should be A");

    // Verify all states are valid
    for state in chain {
        let state_str = state.as_str().unwrap();
        assert!(
            state_str == "A" || state_str == "B" || state_str == "C",
            "All states should be valid"
        );
    }
}

#[test]
fn test_markov_chain_two_state() {
    let request = create_request(
        "markov_chain",
        json!({
            "states": ["On", "Off"],
            "transition_matrix": [
                [0.8, 0.2],
                [0.3, 0.7]
            ],
            "initial_state": "On",
            "steps": 20
        }),
    );
    let response = handle(&request);

    assert!(response.success, "Two-state Markov chain should succeed");
}

#[test]
fn test_markov_chain_missing_parameters() {
    let request = create_request("markov_chain", json!({}));
    let response = handle(&request);

    assert!(!response.success, "Should fail with missing parameters");
    assert!(response.error.is_some(), "Should have error message");
}

#[test]
fn test_markov_chain_invalid_transition_matrix_dimensions() {
    let request = create_request(
        "markov_chain",
        json!({
            "states": ["A", "B", "C"],
            "transition_matrix": [
                [0.5, 0.5],  // Wrong size!
                [0.3, 0.4, 0.3],
                [0.1, 0.3, 0.6]
            ],
            "initial_state": "A",
            "steps": 10
        }),
    );
    let response = handle(&request);

    assert!(!response.success, "Should fail with mismatched dimensions");
    assert!(response.error.is_some(), "Should have error message");
}

#[test]
fn test_markov_chain_invalid_initial_state() {
    let request = create_request(
        "markov_chain",
        json!({
            "states": ["A", "B", "C"],
            "transition_matrix": [
                [0.7, 0.2, 0.1],
                [0.3, 0.4, 0.3],
                [0.1, 0.3, 0.6]
            ],
            "initial_state": "D",  // Invalid state
            "steps": 10
        }),
    );
    let response = handle(&request);

    assert!(!response.success, "Should fail with invalid initial state");
}

#[test]
fn test_markov_chain_invalid_probabilities_sum() {
    let request = create_request(
        "markov_chain",
        json!({
            "states": ["A", "B"],
            "transition_matrix": [
                [0.6, 0.5],  // Sum > 1!
                [0.3, 0.7]
            ],
            "initial_state": "A",
            "steps": 5
        }),
    );
    let response = handle(&request);

    assert!(!response.success, "Should fail when probabilities don't sum to 1");
}

// ============================================================================
// Poisson Process Tests
// ============================================================================

#[test]
fn test_poisson_process_success() {
    let request = create_request(
        "poisson_process",
        json!({
            "rate": 2.0,
            "time_horizon": 10.0
        }),
    );
    let response = handle(&request);

    assert!(response.success, "Poisson process should succeed");
    assert!(response.error.is_none(), "Should have no error");

    let result = response.result.unwrap();
    assert!(result.get("events").is_some(), "Should have events");

    let events = result["events"].as_array().unwrap();

    // Verify all events are within time horizon
    for event in events {
        let time = event.as_f64().unwrap();
        assert!(time >= 0.0 && time < 10.0, "Events should be within time horizon");
    }

    // Verify events are sorted
    for i in 1..events.len() {
        let prev = events[i - 1].as_f64().unwrap();
        let curr = events[i].as_f64().unwrap();
        assert!(curr > prev, "Events should be sorted in time");
    }
}

#[test]
fn test_poisson_process_high_rate() {
    let request = create_request(
        "poisson_process",
        json!({
            "rate": 10.0,
            "time_horizon": 5.0
        }),
    );
    let response = handle(&request);

    assert!(response.success, "Should succeed with high rate");

    let result = response.result.unwrap();
    let events = result["events"].as_array().unwrap();
    // With high rate, we expect many events
    assert!(events.len() > 10, "High rate should generate many events");
}

#[test]
fn test_poisson_process_low_rate() {
    let request = create_request(
        "poisson_process",
        json!({
            "rate": 0.1,
            "time_horizon": 5.0
        }),
    );
    let response = handle(&request);

    assert!(response.success, "Should succeed with low rate");
}

#[test]
fn test_poisson_process_missing_rate() {
    let request = create_request(
        "poisson_process",
        json!({
            "time_horizon": 10.0
        }),
    );
    let response = handle(&request);

    assert!(!response.success, "Should fail without rate");
    assert!(response.error.is_some(), "Should have error message");
    assert_eq!(
        response.error.unwrap(),
        "Missing rate parameter",
        "Should report missing rate"
    );
}

#[test]
fn test_poisson_process_missing_time_horizon() {
    let request = create_request(
        "poisson_process",
        json!({
            "rate": 2.0
        }),
    );
    let response = handle(&request);

    assert!(!response.success, "Should fail without time_horizon");
    assert!(response.error.is_some(), "Should have error message");
    assert_eq!(
        response.error.unwrap(),
        "Missing time_horizon parameter",
        "Should report missing time_horizon"
    );
}

#[test]
fn test_poisson_process_invalid_rate_type() {
    let request = create_request(
        "poisson_process",
        json!({
            "rate": "invalid",
            "time_horizon": 10.0
        }),
    );
    let response = handle(&request);

    assert!(!response.success, "Should fail with invalid rate type");
    assert!(response.error.is_some(), "Should have error message");
    assert_eq!(
        response.error.unwrap(),
        "Invalid rate parameter",
        "Should report invalid rate"
    );
}

#[test]
fn test_poisson_process_invalid_time_horizon_type() {
    let request = create_request(
        "poisson_process",
        json!({
            "rate": 2.0,
            "time_horizon": "invalid"
        }),
    );
    let response = handle(&request);

    assert!(!response.success, "Should fail with invalid time_horizon type");
    assert!(response.error.is_some(), "Should have error message");
    assert_eq!(
        response.error.unwrap(),
        "Invalid time_horizon parameter",
        "Should report invalid time_horizon"
    );
}

// ============================================================================
// Unknown Operation Tests
// ============================================================================

#[test]
fn test_unknown_operation() {
    let request = create_request("unknown_operation", json!({}));
    let response = handle(&request);

    assert!(!response.success, "Should fail with unknown operation");
    assert!(response.error.is_some(), "Should have error message");
    assert!(
        response.error.unwrap().contains("Unknown operation"),
        "Error should mention unknown operation"
    );
}

#[test]
fn test_typo_in_operation_name() {
    let request = create_request("brownian_motin", json!({})); // Typo
    let response = handle(&request);

    assert!(!response.success, "Should fail with typo in operation name");
    assert!(
        response.error.unwrap().contains("Unknown operation"),
        "Should report unknown operation"
    );
}
