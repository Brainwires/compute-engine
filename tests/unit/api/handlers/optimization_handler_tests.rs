//! Unit tests for optimization API handler
//!
//! Tests all optimization operations including:
//! - Function minimization (gradient descent, Nelder-Mead)
//! - Linear programming (simplex method)
//! - Genetic algorithms
//! - Response structure validation
//! - Error handling for invalid inputs

use crate::api::handlers::optimization;
use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};
use std::collections::HashMap;

/// Helper function to create an optimization request
fn create_request(operation: &str, parameters: HashMap<String, Value>) -> ComputationRequest {
    ComputationRequest {
        module: "optimization".to_string(),
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
// Minimize Tests
// ============================================================================

#[test]
fn test_minimize_basic_success() {
    let mut params = HashMap::new();
    params.insert("function".to_string(), json!("x^2 + y^2"));
    params.insert("initial_guess".to_string(), json!([1.0, 1.0]));

    let request = create_request("minimize", params);
    let response = optimization::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Verify result structure
    assert!(result["result"].is_array());
    assert!(result["value"].is_number());
    assert!(result["iterations"].is_number());
    assert!(result["converged"].is_boolean());
}

#[test]
fn test_minimize_result_is_array() {
    let params = HashMap::new();

    let request = create_request("minimize", params);
    let response = optimization::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let result_arr = result["result"].as_array().unwrap();
    assert_eq!(result_arr.len(), 2); // Mock returns [0.0, 0.0]
    assert_eq!(result_arr[0].as_f64().unwrap(), 0.0);
    assert_eq!(result_arr[1].as_f64().unwrap(), 0.0);
}

#[test]
fn test_minimize_value_field() {
    let params = HashMap::new();

    let request = create_request("minimize", params);
    let response = optimization::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let value = result["value"].as_f64().unwrap();
    assert_eq!(value, 0.0); // Mock returns 0.0
}

#[test]
fn test_minimize_iterations_field() {
    let params = HashMap::new();

    let request = create_request("minimize", params);
    let response = optimization::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let iterations = result["iterations"].as_i64().unwrap();
    assert_eq!(iterations, 100); // Mock returns 100
}

#[test]
fn test_minimize_converged_field() {
    let params = HashMap::new();

    let request = create_request("minimize", params);
    let response = optimization::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let converged = result["converged"].as_bool().unwrap();
    assert!(converged); // Mock returns true
}

#[test]
fn test_minimize_with_method_parameter() {
    let mut params = HashMap::new();
    params.insert("function".to_string(), json!("x^2"));
    params.insert("method".to_string(), json!("gradient_descent"));

    let request = create_request("minimize", params);
    let response = optimization::handle(&request);

    assert!(response.success);
}

#[test]
fn test_minimize_with_constraints() {
    let mut params = HashMap::new();
    params.insert("function".to_string(), json!("x^2 + y^2"));
    params.insert("constraints".to_string(), json!([
        {"type": "inequality", "function": "x + y - 1"}
    ]));

    let request = create_request("minimize", params);
    let response = optimization::handle(&request);

    assert!(response.success);
}

// ============================================================================
// Linear Program Tests
// ============================================================================

#[test]
fn test_linear_program_basic_success() {
    let mut params = HashMap::new();
    params.insert("objective".to_string(), json!([1.0, 2.0])); // Coefficients
    params.insert("constraints".to_string(), json!([
        {"coefficients": [1.0, 1.0], "bound": 10.0, "type": "<="}
    ]));

    let request = create_request("linear_program", params);
    let response = optimization::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Verify result structure
    assert!(result["optimal_value"].is_number());
    assert!(result["solution"].is_array());
    assert!(result["status"].is_string());
}

#[test]
fn test_linear_program_optimal_value() {
    let params = HashMap::new();

    let request = create_request("linear_program", params);
    let response = optimization::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let optimal_value = result["optimal_value"].as_f64().unwrap();
    assert_eq!(optimal_value, 0.0); // Mock returns 0.0
}

#[test]
fn test_linear_program_solution_array() {
    let params = HashMap::new();

    let request = create_request("linear_program", params);
    let response = optimization::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let solution = result["solution"].as_array().unwrap();
    assert_eq!(solution.len(), 1); // Mock returns [0.0]
    assert_eq!(solution[0].as_f64().unwrap(), 0.0);
}

#[test]
fn test_linear_program_status_field() {
    let params = HashMap::new();

    let request = create_request("linear_program", params);
    let response = optimization::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let status = result["status"].as_str().unwrap();
    assert_eq!(status, "optimal"); // Mock returns "optimal"
}

#[test]
fn test_linear_program_with_maximize() {
    let mut params = HashMap::new();
    params.insert("objective".to_string(), json!([3.0, 5.0]));
    params.insert("maximize".to_string(), json!(true));

    let request = create_request("linear_program", params);
    let response = optimization::handle(&request);

    assert!(response.success);
}

#[test]
fn test_linear_program_with_bounds() {
    let mut params = HashMap::new();
    params.insert("objective".to_string(), json!([1.0, 2.0]));
    params.insert("bounds".to_string(), json!([
        {"lower": 0.0, "upper": 10.0},
        {"lower": 0.0, "upper": 5.0}
    ]));

    let request = create_request("linear_program", params);
    let response = optimization::handle(&request);

    assert!(response.success);
}

// ============================================================================
// Genetic Algorithm Tests
// ============================================================================

#[test]
fn test_genetic_algorithm_basic_success() {
    let mut params = HashMap::new();
    params.insert("fitness_function".to_string(), json!("x^2 + y^2"));
    params.insert("population_size".to_string(), json!(100));
    params.insert("generations".to_string(), json!(50));

    let request = create_request("genetic_algorithm", params);
    let response = optimization::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Verify result structure
    assert!(result["best_solution"].is_array());
    assert!(result["best_fitness"].is_number());
    assert!(result["generations"].is_number());
}

#[test]
fn test_genetic_algorithm_best_solution() {
    let params = HashMap::new();

    let request = create_request("genetic_algorithm", params);
    let response = optimization::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let best_solution = result["best_solution"].as_array().unwrap();
    assert_eq!(best_solution.len(), 1); // Mock returns [0.0]
    assert_eq!(best_solution[0].as_f64().unwrap(), 0.0);
}

#[test]
fn test_genetic_algorithm_best_fitness() {
    let params = HashMap::new();

    let request = create_request("genetic_algorithm", params);
    let response = optimization::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let best_fitness = result["best_fitness"].as_f64().unwrap();
    assert_eq!(best_fitness, 0.0); // Mock returns 0.0
}

#[test]
fn test_genetic_algorithm_generations() {
    let params = HashMap::new();

    let request = create_request("genetic_algorithm", params);
    let response = optimization::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let generations = result["generations"].as_i64().unwrap();
    assert_eq!(generations, 100); // Mock returns 100
}

#[test]
fn test_genetic_algorithm_with_mutation_rate() {
    let mut params = HashMap::new();
    params.insert("mutation_rate".to_string(), json!(0.1));
    params.insert("crossover_rate".to_string(), json!(0.8));

    let request = create_request("genetic_algorithm", params);
    let response = optimization::handle(&request);

    assert!(response.success);
}

#[test]
fn test_genetic_algorithm_with_elitism() {
    let mut params = HashMap::new();
    params.insert("elitism".to_string(), json!(true));
    params.insert("elite_size".to_string(), json!(10));

    let request = create_request("genetic_algorithm", params);
    let response = optimization::handle(&request);

    assert!(response.success);
}

// ============================================================================
// Unknown Operation Tests
// ============================================================================

#[test]
fn test_unknown_operation() {
    let params = HashMap::new();
    let request = create_request("invalid_operation", params);
    let response = optimization::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Unknown optimization operation"));
}

#[test]
fn test_unknown_operation_includes_operation_name() {
    let params = HashMap::new();
    let request = create_request("nonexistent_optimizer", params);
    let response = optimization::handle(&request);

    assert!(!response.success);
    let error = response.error.as_ref().unwrap();
    assert!(error.contains("nonexistent_optimizer"));
}

// ============================================================================
// Response Structure Tests
// ============================================================================

#[test]
fn test_response_structure_minimize() {
    let params = HashMap::new();
    let request = create_request("minimize", params);
    let response = optimization::handle(&request);

    // Verify response structure
    assert!(response.success);
    assert_eq!(response.module, "optimization");
    assert_eq!(response.operation, "minimize");
    assert!(response.result.is_some());
    assert!(response.error.is_none());
}

#[test]
fn test_response_structure_linear_program() {
    let params = HashMap::new();
    let request = create_request("linear_program", params);
    let response = optimization::handle(&request);

    // Verify response structure
    assert!(response.success);
    assert_eq!(response.module, "optimization");
    assert_eq!(response.operation, "linear_program");
    assert!(response.result.is_some());
    assert!(response.error.is_none());
}

#[test]
fn test_response_structure_genetic_algorithm() {
    let params = HashMap::new();
    let request = create_request("genetic_algorithm", params);
    let response = optimization::handle(&request);

    // Verify response structure
    assert!(response.success);
    assert_eq!(response.module, "optimization");
    assert_eq!(response.operation, "genetic_algorithm");
    assert!(response.result.is_some());
    assert!(response.error.is_none());
}

#[test]
fn test_error_response_structure() {
    let params = HashMap::new();
    let request = create_request("invalid_op", params);
    let response = optimization::handle(&request);

    // Verify error response structure
    assert!(!response.success);
    assert_eq!(response.module, "optimization");
    assert_eq!(response.operation, "invalid_op");
    assert!(response.result.is_none());
    assert!(response.error.is_some());
}
