//! Unit tests for function_approximator API handler

use super::*; // Import from parent module (function_approximator handler)
use crate::api::types::ComputationRequest;
use serde_json::json;
use std::collections::HashMap;

fn make_request(operation: &str, parameters: serde_json::Value) -> ComputationRequest {
    let params_map: HashMap<String, serde_json::Value> = serde_json::from_value(parameters)
        .expect("Failed to convert parameters to HashMap");

    ComputationRequest {
        module: "function_approximator".to_string(),
        operation: operation.to_string(),
        parameters: params_map,
    }
}

// ========== Basic Functionality Tests ==========

#[test]
fn test_approximate_linear_function() {
    // Test approximating y = 2x + 1
    let inputs: Vec<f64> = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let outputs: Vec<f64> = inputs.iter().map(|x| 2.0 * x + 1.0).collect();

    let request = make_request(
        "approximate",
        json!({
            "inputs": inputs,
            "outputs": outputs,
            "population_size": 100,
            "generations": 20
        }),
    );

    let response = handle(&request);
    assert!(response.success, "Response should be successful");
    assert!(response.result.is_some(), "Result should be present");

    let result = response.result.unwrap();
    assert!(result.get("best_functions").is_some(), "Best functions should be present");
    assert!(result.get("generations_run").is_some(), "Generations run should be present");
    assert!(result.get("final_population_stats").is_some(), "Population stats should be present");

    // Check that we got 20 generations as requested
    let generations = result["generations_run"].as_u64().unwrap();
    assert_eq!(generations, 20, "Should run exactly 20 generations");
}

#[test]
fn test_discover_function_alias() {
    // Test that "discover_function" works as an alias for "approximate"
    let inputs = vec![1.0, 2.0, 3.0];
    let outputs = vec![1.0, 4.0, 9.0]; // y = x^2

    let request = make_request(
        "discover_function",
        json!({
            "inputs": inputs,
            "outputs": outputs,
            "generations": 10
        }),
    );

    let response = handle(&request);
    assert!(response.success, "discover_function should work as an alias");
    assert!(response.result.is_some());
}

#[test]
fn test_approximate_quadratic_function() {
    // Test approximating y = x^2
    let inputs: Vec<f64> = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let outputs: Vec<f64> = inputs.iter().map(|x| x * x).collect();

    let request = make_request(
        "approximate",
        json!({
            "inputs": inputs,
            "outputs": outputs,
            "population_size": 200,
            "generations": 30
        }),
    );

    let response = handle(&request);
    assert!(response.success);

    let result = response.result.unwrap();
    let best_functions = result["best_functions"].as_array().unwrap();

    // Should return top 10 functions (or less if population is smaller)
    assert!(best_functions.len() > 0, "Should have at least one best function");
    assert!(best_functions.len() <= 10, "Should have at most 10 best functions");

    // Check structure of best function
    let best_fn = &best_functions[0];
    assert!(best_fn.get("expression").is_some(), "Should have expression");
    assert!(best_fn.get("fitness").is_some(), "Should have fitness");
    assert!(best_fn.get("complexity").is_some(), "Should have complexity");
    assert!(best_fn.get("predictions").is_some(), "Should have predictions");
}

// ========== Configuration Parameter Tests ==========

#[test]
fn test_custom_population_size() {
    let inputs = vec![1.0, 2.0, 3.0];
    let outputs = vec![2.0, 4.0, 6.0];

    let request = make_request(
        "approximate",
        json!({
            "inputs": inputs,
            "outputs": outputs,
            "population_size": 50,
            "generations": 5
        }),
    );

    let response = handle(&request);
    assert!(response.success, "Should handle custom population size");
}

#[test]
fn test_custom_generations() {
    let inputs = vec![1.0, 2.0, 3.0];
    let outputs = vec![1.0, 4.0, 9.0];

    let request = make_request(
        "approximate",
        json!({
            "inputs": inputs,
            "outputs": outputs,
            "generations": 50
        }),
    );

    let response = handle(&request);
    assert!(response.success);

    let result = response.result.unwrap();
    assert_eq!(result["generations_run"].as_u64().unwrap(), 50);
}

#[test]
fn test_custom_complexity_penalty() {
    let inputs = vec![0.0, 1.0, 2.0];
    let outputs = vec![0.0, 1.0, 4.0];

    let request = make_request(
        "approximate",
        json!({
            "inputs": inputs,
            "outputs": outputs,
            "complexity": 5.0,
            "generations": 10
        }),
    );

    let response = handle(&request);
    assert!(response.success, "Should handle custom complexity penalty");
}

#[test]
fn test_custom_death_probability() {
    let inputs = vec![1.0, 2.0, 3.0];
    let outputs = vec![2.0, 4.0, 6.0];

    let request = make_request(
        "approximate",
        json!({
            "inputs": inputs,
            "outputs": outputs,
            "death_probability": 0.3,
            "generations": 10
        }),
    );

    let response = handle(&request);
    assert!(response.success, "Should handle custom death probability");
}

#[test]
fn test_custom_diversity_factor() {
    let inputs = vec![1.0, 2.0, 3.0];
    let outputs = vec![1.0, 8.0, 27.0];

    let request = make_request(
        "approximate",
        json!({
            "inputs": inputs,
            "outputs": outputs,
            "diversity": 5,
            "generations": 10
        }),
    );

    let response = handle(&request);
    assert!(response.success, "Should handle custom diversity factor");
}

#[test]
fn test_all_custom_parameters() {
    let inputs = vec![1.0, 2.0, 3.0, 4.0];
    let outputs = vec![1.0, 4.0, 9.0, 16.0];

    let request = make_request(
        "approximate",
        json!({
            "inputs": inputs,
            "outputs": outputs,
            "population_size": 150,
            "generations": 25,
            "complexity": 3.0,
            "death_probability": 0.4,
            "diversity": 8
        }),
    );

    let response = handle(&request);
    assert!(response.success, "Should handle all custom parameters");

    let result = response.result.unwrap();
    assert_eq!(result["generations_run"].as_u64().unwrap(), 25);
}

// ========== Error Handling Tests ==========

#[test]
fn test_missing_inputs_parameter() {
    let request = make_request(
        "approximate",
        json!({
            "outputs": vec![1.0, 2.0, 3.0]
        }),
    );

    let response = handle(&request);
    assert!(!response.success, "Should fail without inputs");
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Missing inputs parameter"));
}

#[test]
fn test_missing_outputs_parameter() {
    let request = make_request(
        "approximate",
        json!({
            "inputs": vec![1.0, 2.0, 3.0]
        }),
    );

    let response = handle(&request);
    assert!(!response.success, "Should fail without outputs");
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Missing outputs parameter"));
}

#[test]
fn test_invalid_inputs_type() {
    let request = make_request(
        "approximate",
        json!({
            "inputs": "not_an_array",
            "outputs": vec![1.0, 2.0]
        }),
    );

    let response = handle(&request);
    assert!(!response.success, "Should fail with invalid inputs type");
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Invalid inputs"));
}

#[test]
fn test_invalid_outputs_type() {
    let request = make_request(
        "approximate",
        json!({
            "inputs": vec![1.0, 2.0],
            "outputs": "not_an_array"
        }),
    );

    let response = handle(&request);
    assert!(!response.success, "Should fail with invalid outputs type");
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Invalid outputs"));
}

#[test]
fn test_mismatched_input_output_lengths() {
    let request = make_request(
        "approximate",
        json!({
            "inputs": vec![1.0, 2.0, 3.0],
            "outputs": vec![1.0, 2.0]
        }),
    );

    let response = handle(&request);
    assert!(!response.success, "Should fail with mismatched lengths");
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("same length"));
}

#[test]
fn test_empty_input_arrays() {
    let request = make_request(
        "approximate",
        json!({
            "inputs": Vec::<f64>::new(),
            "outputs": Vec::<f64>::new()
        }),
    );

    let response = handle(&request);
    assert!(!response.success, "Should fail with empty arrays");
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("at least one"));
}

#[test]
fn test_unknown_operation() {
    let request = make_request(
        "invalid_operation",
        json!({
            "inputs": vec![1.0, 2.0],
            "outputs": vec![2.0, 4.0]
        }),
    );

    let response = handle(&request);
    assert!(!response.success, "Should fail with unknown operation");
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Unknown operation"));
}

// ========== Population Statistics Tests ==========

#[test]
fn test_population_stats_structure() {
    let inputs = vec![0.0, 1.0, 2.0, 3.0];
    let outputs = vec![0.0, 1.0, 4.0, 9.0];

    let request = make_request(
        "approximate",
        json!({
            "inputs": inputs,
            "outputs": outputs,
            "generations": 10
        }),
    );

    let response = handle(&request);
    assert!(response.success);

    let result = response.result.unwrap();
    let stats = &result["final_population_stats"];

    assert!(stats.get("min_fitness").is_some(), "Should have min_fitness");
    assert!(stats.get("max_fitness").is_some(), "Should have max_fitness");
    assert!(stats.get("avg_fitness").is_some(), "Should have avg_fitness");
    assert!(stats.get("fitness_distribution").is_some(), "Should have fitness_distribution");

    // Note: avg_fitness may be null if it's NaN/Infinity (valid JSON serialization behavior)
    // Just verify the structure exists
    let fitness_dist = stats["fitness_distribution"].as_array().unwrap();
    assert!(fitness_dist.len() > 0, "Should have fitness distribution values");
}

#[test]
fn test_predictions_match_input_length() {
    let inputs = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let outputs = vec![2.0, 4.0, 6.0, 8.0, 10.0];

    let request = make_request(
        "approximate",
        json!({
            "inputs": inputs.clone(),
            "outputs": outputs,
            "generations": 15
        }),
    );

    let response = handle(&request);
    assert!(response.success);

    let result = response.result.unwrap();
    let best_functions = result["best_functions"].as_array().unwrap();

    for best_fn in best_functions {
        let predictions = best_fn["predictions"].as_array().unwrap();
        assert_eq!(predictions.len(), inputs.len(),
                   "Predictions length should match inputs length");
    }
}

// ========== Edge Case Tests ==========

#[test]
fn test_single_point_approximation() {
    // Although not very useful, should handle single point
    let request = make_request(
        "approximate",
        json!({
            "inputs": vec![5.0],
            "outputs": vec![25.0],
            "generations": 5
        }),
    );

    let response = handle(&request);
    assert!(response.success, "Should handle single point");
}

#[test]
fn test_negative_values() {
    let inputs = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
    let outputs = vec![4.0, 1.0, 0.0, 1.0, 4.0]; // y = x^2

    let request = make_request(
        "approximate",
        json!({
            "inputs": inputs,
            "outputs": outputs,
            "generations": 20
        }),
    );

    let response = handle(&request);
    assert!(response.success, "Should handle negative values");
}

#[test]
fn test_very_small_population() {
    let inputs = vec![1.0, 2.0, 3.0];
    let outputs = vec![1.0, 4.0, 9.0];

    let request = make_request(
        "approximate",
        json!({
            "inputs": inputs,
            "outputs": outputs,
            "population_size": 10,
            "generations": 5
        }),
    );

    let response = handle(&request);
    assert!(response.success, "Should handle very small population");
}
