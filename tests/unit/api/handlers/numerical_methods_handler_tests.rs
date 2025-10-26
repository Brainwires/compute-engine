//! Unit tests for numerical_methods API handler

use crate::api::handlers::numerical_methods::handle;
use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};
use std::collections::HashMap;

/// Helper to create a request
fn create_request(operation: &str, parameters: HashMap<String, Value>) -> ComputationRequest {
    ComputationRequest {
        module: "numerical_methods".to_string(),
        operation: operation.to_string(),
        parameters,
    }
}

/// Helper to extract result value
fn extract_result(response: ComputationResponse) -> Value {
    assert!(
        response.success,
        "Expected success but got error: {:?}",
        response.error
    );
    response.result.unwrap()
}

// ============================================================================
// ODE Solver Tests
// ============================================================================

#[test]
fn test_solve_ode_euler_method() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("euler"));
    params.insert("initial_value".to_string(), json!(1.0));
    params.insert("t_start".to_string(), json!(0.0));
    params.insert("t_end".to_string(), json!(1.0));
    params.insert("step_size".to_string(), json!(0.1));

    let request = create_request("solve_ode", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert_eq!(result["method_used"], "euler");
    assert!(result["steps_taken"].as_u64().unwrap() > 0);
    assert!(result["t_values"].is_array());
    assert!(result["y_values"].is_array());
}

#[test]
fn test_solve_ode_rk4_method() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("rk4"));
    params.insert("initial_value".to_string(), json!(1.0));
    params.insert("t_start".to_string(), json!(0.0));
    params.insert("t_end".to_string(), json!(0.5));
    params.insert("step_size".to_string(), json!(0.05));

    let request = create_request("solve_ode", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert_eq!(result["method_used"], "rk4");
    let t_values = result["t_values"].as_array().unwrap();
    let y_values = result["y_values"].as_array().unwrap();
    assert_eq!(t_values.len(), y_values.len());
}

#[test]
fn test_solve_ode_adaptive_rk45() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("adaptive_rk45"));
    params.insert("initial_value".to_string(), json!(1.0));
    params.insert("t_start".to_string(), json!(0.0));
    params.insert("t_end".to_string(), json!(1.0));
    params.insert("step_size".to_string(), json!(0.1));

    let request = create_request("solve_ode", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert_eq!(result["method_used"], "adaptive_rk45");
}

#[test]
fn test_solve_ode_invalid_method() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("invalid_method"));
    params.insert("initial_value".to_string(), json!(1.0));
    params.insert("t_start".to_string(), json!(0.0));
    params.insert("t_end".to_string(), json!(1.0));
    params.insert("step_size".to_string(), json!(0.1));

    let request = create_request("solve_ode", params);
    let response = handle(&request);

    assert!(!response.success);
    assert!(response.error.unwrap().contains("Unknown ODE method"));
}

#[test]
fn test_solve_ode_missing_parameters() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("euler"));
    // Missing other required parameters

    let request = create_request("solve_ode", params);
    let response = handle(&request);

    assert!(!response.success);
    assert!(response.error.unwrap().contains("Invalid parameters"));
}

// ============================================================================
// Root Finding Tests
// ============================================================================

#[test]
fn test_find_root_newton_method() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("newton"));
    params.insert("initial_guess".to_string(), json!(1.5));
    params.insert("tolerance".to_string(), json!(1e-6));
    params.insert("max_iterations".to_string(), json!(100));
    params.insert("function_type".to_string(), json!("polynomial"));
    params.insert("coefficients".to_string(), json!([-2.0, 0.0, 1.0])); // x^2 - 2

    let request = create_request("find_root", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert!(result["converged"].as_bool().unwrap());
    let root = result["root"].as_f64().unwrap();
    assert!((root - 1.414213).abs() < 0.001); // sqrt(2)
}

#[test]
fn test_find_root_bisection_method() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("bisection"));
    params.insert("initial_guess".to_string(), json!(1.0));
    params.insert("interval".to_string(), json!([1.0, 2.0]));
    params.insert("tolerance".to_string(), json!(1e-6));
    params.insert("max_iterations".to_string(), json!(100));
    params.insert("function_type".to_string(), json!("polynomial"));
    params.insert("coefficients".to_string(), json!([-2.0, 0.0, 1.0])); // x^2 - 2

    let request = create_request("find_root", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert!(result["converged"].as_bool().unwrap());
    let root = result["root"].as_f64().unwrap();
    assert!((root - 1.414213).abs() < 0.01);
}

#[test]
fn test_find_root_secant_method() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("secant"));
    params.insert("initial_guess".to_string(), json!(1.5));
    params.insert("tolerance".to_string(), json!(1e-6));
    params.insert("max_iterations".to_string(), json!(100));
    params.insert("function_type".to_string(), json!("polynomial"));
    params.insert("coefficients".to_string(), json!([-2.0, 0.0, 1.0]));

    let request = create_request("find_root", params);
    let response = handle(&request);
    let result = extract_result(response);

    let root = result["root"].as_f64().unwrap();
    assert!((root - 1.414213).abs() < 0.01);
}

#[test]
fn test_find_root_invalid_method() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("invalid"));
    params.insert("initial_guess".to_string(), json!(1.0));
    params.insert("tolerance".to_string(), json!(1e-6));
    params.insert("max_iterations".to_string(), json!(100));
    params.insert("function_type".to_string(), json!("polynomial"));

    let request = create_request("find_root", params);
    let response = handle(&request);

    assert!(!response.success);
    assert!(response.error.unwrap().contains("Unknown root finding method"));
}

// ============================================================================
// Integration Tests
// ============================================================================

#[test]
fn test_integrate_trapezoidal() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("trapezoidal"));
    params.insert("lower_bound".to_string(), json!(0.0));
    params.insert("upper_bound".to_string(), json!(1.0));
    params.insert("num_points".to_string(), json!(100));
    params.insert("function_type".to_string(), json!("polynomial"));
    params.insert("coefficients".to_string(), json!([0.0, 0.0, 1.0])); // x^2

    let request = create_request("integrate", params);
    let response = handle(&request);
    let result = extract_result(response);

    let integral = result["integral"].as_f64().unwrap();
    assert!((integral - 0.333333).abs() < 0.01); // Integral of x^2 from 0 to 1 is 1/3
    assert_eq!(result["method_used"], "trapezoidal");
}

#[test]
fn test_integrate_simpson() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("simpson"));
    params.insert("lower_bound".to_string(), json!(0.0));
    params.insert("upper_bound".to_string(), json!(1.0));
    params.insert("num_points".to_string(), json!(100)); // Must be even
    params.insert("function_type".to_string(), json!("polynomial"));
    params.insert("coefficients".to_string(), json!([0.0, 0.0, 1.0])); // x^2

    let request = create_request("integrate", params);
    let response = handle(&request);
    let result = extract_result(response);

    let integral = result["integral"].as_f64().unwrap();
    assert!((integral - 0.333333).abs() < 0.001);
    assert_eq!(result["method_used"], "simpson");
}

#[test]
fn test_integrate_gauss_legendre() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("gauss"));
    params.insert("lower_bound".to_string(), json!(0.0));
    params.insert("upper_bound".to_string(), json!(1.0));
    params.insert("num_points".to_string(), json!(10));
    params.insert("function_type".to_string(), json!("polynomial"));
    params.insert("coefficients".to_string(), json!([0.0, 0.0, 1.0]));

    let request = create_request("integrate", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert_eq!(result["method_used"], "gauss");
    assert!(result["integral"].is_f64());
}

#[test]
fn test_integrate_adaptive_simpson() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("adaptive_simpson"));
    params.insert("lower_bound".to_string(), json!(0.0));
    params.insert("upper_bound".to_string(), json!(1.0));
    params.insert("num_points".to_string(), json!(10));
    params.insert("function_type".to_string(), json!("polynomial"));
    params.insert("coefficients".to_string(), json!([0.0, 0.0, 1.0]));

    let request = create_request("integrate", params);
    let response = handle(&request);
    let result = extract_result(response);

    let integral = result["integral"].as_f64().unwrap();
    assert!((integral - 0.333333).abs() < 1e-6); // More accurate with adaptive
}

// ============================================================================
// Interpolation Tests
// ============================================================================

#[test]
fn test_interpolate_linear() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("linear"));
    params.insert("x_values".to_string(), json!([0.0, 1.0, 2.0]));
    params.insert("y_values".to_string(), json!([0.0, 1.0, 4.0]));
    params.insert("interpolate_at".to_string(), json!([0.5, 1.5]));

    let request = create_request("interpolate", params);
    let response = handle(&request);
    let result = extract_result(response);

    let values = result["interpolated_values"].as_array().unwrap();
    assert_eq!(values.len(), 2);
    assert!((values[0].as_f64().unwrap() - 0.5).abs() < 0.01);
    assert!((values[1].as_f64().unwrap() - 2.5).abs() < 0.01);
}

#[test]
fn test_interpolate_lagrange() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("lagrange"));
    params.insert("x_values".to_string(), json!([0.0, 1.0, 2.0]));
    params.insert("y_values".to_string(), json!([0.0, 1.0, 4.0]));
    params.insert("interpolate_at".to_string(), json!([0.5]));

    let request = create_request("interpolate", params);
    let response = handle(&request);
    let result = extract_result(response);

    let values = result["interpolated_values"].as_array().unwrap();
    assert_eq!(values.len(), 1);
}

#[test]
fn test_interpolate_mismatched_arrays() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("linear"));
    params.insert("x_values".to_string(), json!([0.0, 1.0]));
    params.insert("y_values".to_string(), json!([0.0, 1.0, 2.0])); // Mismatched length

    params.insert("interpolate_at".to_string(), json!([0.5]));

    let request = create_request("interpolate", params);
    let response = handle(&request);

    assert!(!response.success);
    assert!(response.error.unwrap().contains("same length"));
}

// ============================================================================
// Linear System Tests
// ============================================================================

#[test]
fn test_solve_linear_system_gauss() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("gauss"));
    params.insert(
        "matrix".to_string(),
        json!([[2.0, 1.0], [1.0, 3.0]]),
    );
    params.insert("rhs".to_string(), json!([5.0, 6.0]));

    let request = create_request("solve_linear_system", params);
    let response = handle(&request);
    let result = extract_result(response);

    let solution = result["solution"].as_array().unwrap();
    assert_eq!(solution.len(), 2);
    // Solution to: 2x + y = 5, x + 3y = 6
    // Solving: x = (15-6)/5 = 1.8, y = (12-5)/5 = 1.4
    // Actually: x = 1.8, y = 1.4 (verify: 2*1.8 + 1.4 = 5.0, 1.8 + 3*1.4 = 6.0)
    assert!((solution[0].as_f64().unwrap() - 1.8).abs() < 0.1);
    assert!((solution[1].as_f64().unwrap() - 1.4).abs() < 0.1);
}

#[test]
fn test_solve_linear_system_jacobi() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("jacobi"));
    params.insert(
        "matrix".to_string(),
        json!([[4.0, 1.0], [1.0, 3.0]]),
    );
    params.insert("rhs".to_string(), json!([1.0, 2.0]));
    params.insert("tolerance".to_string(), json!(1e-6));
    params.insert("max_iterations".to_string(), json!(100));

    let request = create_request("solve_linear_system", params);
    let response = handle(&request);
    let result = extract_result(response);

    let solution = result["solution"].as_array().unwrap();
    assert_eq!(solution.len(), 2);
    assert!(result["iterations"].is_number());
}

#[test]
fn test_solve_linear_system_mismatched_dimensions() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("gauss"));
    params.insert(
        "matrix".to_string(),
        json!([[2.0, 1.0], [1.0, 3.0]]),
    );
    params.insert("rhs".to_string(), json!([5.0, 6.0, 7.0])); // Wrong size

    let request = create_request("solve_linear_system", params);
    let response = handle(&request);

    assert!(!response.success);
    assert!(response.error.unwrap().contains("don't match"));
}

// ============================================================================
// Differentiation Tests
// ============================================================================

#[test]
fn test_differentiate_forward() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("forward"));
    params.insert("x_values".to_string(), json!([0.0, 1.0, 2.0, 3.0]));
    params.insert("y_values".to_string(), json!([0.0, 1.0, 4.0, 9.0])); // y = x^2
    params.insert("order".to_string(), json!(1));

    let request = create_request("differentiate", params);
    let response = handle(&request);
    let result = extract_result(response);

    let derivatives = result["derivatives"].as_array().unwrap();
    assert_eq!(derivatives.len(), 4);
    assert_eq!(result["method_used"], "forward");
}

#[test]
fn test_differentiate_backward() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("backward"));
    params.insert("x_values".to_string(), json!([0.0, 1.0, 2.0, 3.0]));
    params.insert("y_values".to_string(), json!([0.0, 1.0, 4.0, 9.0]));
    params.insert("order".to_string(), json!(1));

    let request = create_request("differentiate", params);
    let response = handle(&request);
    let result = extract_result(response);

    let derivatives = result["derivatives"].as_array().unwrap();
    assert_eq!(derivatives.len(), 4);
    assert_eq!(result["method_used"], "backward");
}

#[test]
fn test_differentiate_central() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("central"));
    params.insert("x_values".to_string(), json!([0.0, 1.0, 2.0, 3.0]));
    params.insert("y_values".to_string(), json!([0.0, 1.0, 4.0, 9.0]));
    params.insert("order".to_string(), json!(1));

    let request = create_request("differentiate", params);
    let response = handle(&request);
    let result = extract_result(response);

    let derivatives = result["derivatives"].as_array().unwrap();
    assert_eq!(derivatives.len(), 4);
    assert_eq!(result["method_used"], "central");
}

#[test]
fn test_differentiate_insufficient_points() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("forward"));
    params.insert("x_values".to_string(), json!([0.0])); // Only 1 point
    params.insert("y_values".to_string(), json!([0.0]));
    params.insert("order".to_string(), json!(1));

    let request = create_request("differentiate", params);
    let response = handle(&request);

    assert!(!response.success);
    assert!(response.error.unwrap().contains("at least 2 points"));
}

// ============================================================================
// PDE Solver Tests
// ============================================================================

#[test]
fn test_solve_pde_heat_equation() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("heat_equation"));
    params.insert("boundary_conditions".to_string(), json!([0.0, 0.0]));
    params.insert("initial_conditions".to_string(), json!([0.0, 1.0, 1.0, 0.0]));
    params.insert("spatial_steps".to_string(), json!(10));
    params.insert("time_steps".to_string(), json!(20));
    params.insert("dx".to_string(), json!(0.1));
    params.insert("dt".to_string(), json!(0.01));

    let request = create_request("solve_pde", params);
    let response = handle(&request);
    let result = extract_result(response);

    let solution = result["solution"].as_array().unwrap();
    assert_eq!(solution.len(), 20); // time_steps
    let final_state = result["final_state"].as_array().unwrap();
    assert_eq!(final_state.len(), 10); // spatial_steps
}

#[test]
fn test_solve_pde_finite_difference() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("finite_difference"));
    params.insert("boundary_conditions".to_string(), json!([0.0, 0.0]));
    params.insert("initial_conditions".to_string(), json!([0.0, 1.0, 0.0]));
    params.insert("spatial_steps".to_string(), json!(5));
    params.insert("time_steps".to_string(), json!(10));
    params.insert("dx".to_string(), json!(0.1));
    params.insert("dt".to_string(), json!(0.01));

    let request = create_request("solve_pde", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert!(result["solution"].is_array());
    assert!(result["final_state"].is_array());
}

// ============================================================================
// Error Handling Tests
// ============================================================================

#[test]
fn test_unknown_operation() {
    let params = HashMap::new();
    let request = create_request("unknown_operation", params);
    let response = handle(&request);

    assert!(!response.success);
    assert!(response
        .error
        .unwrap()
        .contains("Unknown numerical_methods operation"));
}

#[test]
fn test_malformed_json_parameters() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("euler"));
    params.insert("initial_value".to_string(), json!("not_a_number")); // Wrong type

    let request = create_request("solve_ode", params);
    let response = handle(&request);

    assert!(!response.success);
    assert!(response.error.unwrap().contains("Invalid parameters"));
}
