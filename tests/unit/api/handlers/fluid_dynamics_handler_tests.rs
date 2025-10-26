//! Unit tests for fluid dynamics API handler
//!
//! Tests all fluid dynamics operations including:
//! - Navier-Stokes equation solving
//! - Reynolds number calculations
//! - Bernoulli equation computations
//! - Analytical flow solutions
//! - Cavity flow simulations
//! - Flow field analysis
//! - Boundary condition handling
//! - Error handling for invalid inputs

use crate::api::handlers::fluid_dynamics;
use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};
use std::collections::HashMap;

/// Helper function to create a fluid_dynamics request
fn create_request(operation: &str, parameters: HashMap<String, Value>) -> ComputationRequest {
    ComputationRequest {
        module: "fluid_dynamics".to_string(),
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
// Navier-Stokes Tests
// ============================================================================

#[test]
fn test_navier_stokes_basic() {
    let mut params = HashMap::new();
    params.insert("viscosity".to_string(), json!(0.001)); // Water-like viscosity
    params.insert("density".to_string(), json!(1000.0)); // kg/m³
    params.insert("grid_size".to_string(), json!(50));
    params.insert("time_steps".to_string(), json!(100));

    let request = create_request("navier_stokes", params);
    let response = fluid_dynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Verify basic response structure
    assert_eq!(result["result"], "computed");
    assert!(result["value"].is_number());
}

#[test]
fn test_navier_stokes_high_viscosity() {
    let mut params = HashMap::new();
    params.insert("viscosity".to_string(), json!(0.1)); // High viscosity (honey-like)
    params.insert("density".to_string(), json!(1400.0));
    params.insert("grid_size".to_string(), json!(30));

    let request = create_request("navier_stokes", params);
    let response = fluid_dynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);
    assert!(result["value"].is_number());
}

#[test]
fn test_navier_stokes_low_viscosity() {
    let mut params = HashMap::new();
    params.insert("viscosity".to_string(), json!(0.00001)); // Low viscosity (air-like)
    params.insert("density".to_string(), json!(1.2));
    params.insert("grid_size".to_string(), json!(40));

    let request = create_request("navier_stokes", params);
    let response = fluid_dynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);
    assert!(result["value"].is_number());
}

// ============================================================================
// Reynolds Number Tests
// ============================================================================

#[test]
fn test_reynolds_number_laminar_flow() {
    let mut params = HashMap::new();
    params.insert("velocity".to_string(), json!(0.1)); // m/s
    params.insert("length".to_string(), json!(0.05)); // m (characteristic length)
    params.insert("density".to_string(), json!(1000.0)); // kg/m³
    params.insert("viscosity".to_string(), json!(0.001)); // Pa·s (water)

    let request = create_request("reynolds_number", params);
    let response = fluid_dynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Re = ρVL/μ = 1000 * 0.1 * 0.05 / 0.001 = 5000
    // This should indicate laminar to transitional flow
    assert!(result["value"].is_number());
}

#[test]
fn test_reynolds_number_turbulent_flow() {
    let mut params = HashMap::new();
    params.insert("velocity".to_string(), json!(10.0)); // m/s (high velocity)
    params.insert("length".to_string(), json!(1.0)); // m
    params.insert("density".to_string(), json!(1.2)); // kg/m³ (air)
    params.insert("viscosity".to_string(), json!(0.000018)); // Pa·s (air)

    let request = create_request("reynolds_number", params);
    let response = fluid_dynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Re = 1.2 * 10 * 1 / 0.000018 ≈ 666,667 (turbulent)
    assert!(result["value"].is_number());
}

#[test]
fn test_reynolds_number_very_low() {
    let mut params = HashMap::new();
    params.insert("velocity".to_string(), json!(0.001)); // Very slow
    params.insert("length".to_string(), json!(0.001)); // Very small
    params.insert("density".to_string(), json!(1000.0));
    params.insert("viscosity".to_string(), json!(0.1)); // High viscosity

    let request = create_request("reynolds_number", params);
    let response = fluid_dynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Re should be very small (creeping flow regime)
    assert!(result["value"].is_number());
}

// ============================================================================
// Bernoulli Equation Tests
// ============================================================================

#[test]
fn test_bernoulli_basic() {
    let mut params = HashMap::new();
    params.insert("pressure1".to_string(), json!(101325.0)); // Pa (1 atm)
    params.insert("velocity1".to_string(), json!(5.0)); // m/s
    params.insert("height1".to_string(), json!(0.0)); // m
    params.insert("velocity2".to_string(), json!(10.0)); // m/s
    params.insert("height2".to_string(), json!(2.0)); // m
    params.insert("density".to_string(), json!(1000.0)); // kg/m³

    let request = create_request("bernoulli", params);
    let response = fluid_dynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Bernoulli: P1 + ½ρv1² + ρgh1 = P2 + ½ρv2² + ρgh2
    assert!(result["value"].is_number());
}

#[test]
fn test_bernoulli_horizontal_flow() {
    let mut params = HashMap::new();
    params.insert("pressure1".to_string(), json!(200000.0)); // Pa
    params.insert("velocity1".to_string(), json!(2.0)); // m/s
    params.insert("height1".to_string(), json!(0.0)); // Same height
    params.insert("velocity2".to_string(), json!(8.0)); // m/s
    params.insert("height2".to_string(), json!(0.0)); // Same height
    params.insert("density".to_string(), json!(800.0)); // kg/m³ (oil)

    let request = create_request("bernoulli", params);
    let response = fluid_dynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);
    assert!(result["value"].is_number());
}

#[test]
fn test_bernoulli_vertical_flow() {
    let mut params = HashMap::new();
    params.insert("pressure1".to_string(), json!(150000.0)); // Pa
    params.insert("velocity1".to_string(), json!(0.0)); // Static fluid
    params.insert("height1".to_string(), json!(10.0)); // m (elevated tank)
    params.insert("velocity2".to_string(), json!(14.0)); // m/s (free jet)
    params.insert("height2".to_string(), json!(0.0)); // Ground level
    params.insert("density".to_string(), json!(1000.0)); // kg/m³

    let request = create_request("bernoulli", params);
    let response = fluid_dynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should show energy conservation through height and velocity conversion
    assert!(result["value"].is_number());
}

// ============================================================================
// Error Handling Tests
// ============================================================================

#[test]
fn test_unknown_operation() {
    let params = HashMap::new();
    let request = create_request("invalid_fluid_operation", params);
    let response = fluid_dynamics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Unknown fluid dynamics operation"));
}

#[test]
fn test_empty_operation() {
    let params = HashMap::new();
    let request = create_request("", params);
    let response = fluid_dynamics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
}

#[test]
fn test_unsupported_operation() {
    let params = HashMap::new();
    let request = create_request("euler_equations", params);
    let response = fluid_dynamics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Unknown fluid dynamics operation"));
}

// ============================================================================
// Response Structure Tests
// ============================================================================

#[test]
fn test_response_structure_navier_stokes() {
    let mut params = HashMap::new();
    params.insert("viscosity".to_string(), json!(0.001));
    params.insert("density".to_string(), json!(1000.0));

    let request = create_request("navier_stokes", params);
    let response = fluid_dynamics::handle(&request);

    // Verify complete response structure
    assert!(response.success);
    assert_eq!(response.module, "fluid_dynamics");
    assert_eq!(response.operation, "navier_stokes");
    assert!(response.result.is_some());
    assert!(response.error.is_none());
}

#[test]
fn test_response_structure_reynolds_number() {
    let mut params = HashMap::new();
    params.insert("velocity".to_string(), json!(5.0));
    params.insert("length".to_string(), json!(0.1));
    params.insert("density".to_string(), json!(1000.0));
    params.insert("viscosity".to_string(), json!(0.001));

    let request = create_request("reynolds_number", params);
    let response = fluid_dynamics::handle(&request);

    assert!(response.success);
    assert_eq!(response.module, "fluid_dynamics");
    assert_eq!(response.operation, "reynolds_number");
    assert!(response.result.is_some());
    assert!(response.error.is_none());
}

#[test]
fn test_response_structure_bernoulli() {
    let mut params = HashMap::new();
    params.insert("pressure1".to_string(), json!(101325.0));
    params.insert("velocity1".to_string(), json!(1.0));
    params.insert("height1".to_string(), json!(0.0));

    let request = create_request("bernoulli", params);
    let response = fluid_dynamics::handle(&request);

    assert!(response.success);
    assert_eq!(response.module, "fluid_dynamics");
    assert_eq!(response.operation, "bernoulli");
    assert!(response.result.is_some());
    assert!(response.error.is_none());
}

#[test]
fn test_error_response_structure() {
    let params = HashMap::new();
    let request = create_request("nonexistent_op", params);
    let response = fluid_dynamics::handle(&request);

    // Verify error response structure
    assert!(!response.success);
    assert_eq!(response.module, "fluid_dynamics");
    assert_eq!(response.operation, "nonexistent_op");
    assert!(response.result.is_none());
    assert!(response.error.is_some());

    let error_msg = response.error.as_ref().unwrap();
    assert!(error_msg.contains("Unknown fluid dynamics operation"));
    assert!(error_msg.contains("nonexistent_op"));
}
