//! Unit tests for information_theory API handler
//!
//! This test file covers all operations in the information_theory handler:
//! - shannon_entropy (also aliased as "entropy")
//! - mutual_information
//! - channel_capacity
//! - huffman_coding
//! - kolmogorov_complexity
//! - conditional_entropy
//! - relative_entropy (also aliased as "kl_divergence")

use crate::api::handlers::information_theory::handle;
use crate::api::types::ComputationRequest;
use serde_json::{json, Value};
use std::collections::HashMap;

// Helper function to convert json! to HashMap
fn to_params(value: Value) -> HashMap<String, Value> {
    match value {
        Value::Object(map) => map.into_iter().collect(),
        _ => HashMap::new(),
    }
}

// ============================================================================
// SHANNON ENTROPY TESTS
// ============================================================================

#[test]
fn test_shannon_entropy_basic() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "shannon_entropy".to_string(),
        parameters: to_params(json!({
            "data": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Shannon entropy should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("entropy").is_some());
    assert!(result.get("entropy_type").is_some());
    assert_eq!(result.get("entropy_type").unwrap().as_str().unwrap(), "shannon");
}

#[test]
fn test_shannon_entropy_alias() {
    // Test that "entropy" works as an alias for "shannon_entropy"
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "entropy".to_string(),
        parameters: to_params(json!({
            "data": [1.0, 2.0, 3.0, 4.0, 5.0]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Entropy alias should succeed");
}

#[test]
fn test_shannon_entropy_with_base() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "shannon_entropy".to_string(),
        parameters: to_params(json!({
            "data": [1.0, 2.0, 3.0, 4.0, 5.0],
            "base": 10.0
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Shannon entropy with custom base should succeed");

    let result = response.result.unwrap();
    assert_eq!(result.get("base").unwrap().as_f64().unwrap(), 10.0);
}

#[test]
fn test_renyi_entropy() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "shannon_entropy".to_string(),
        parameters: to_params(json!({
            "data": [1.0, 2.0, 3.0, 4.0, 5.0],
            "entropy_type": "renyi",
            "alpha": 2.0
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Rényi entropy should succeed");

    let result = response.result.unwrap();
    assert_eq!(result.get("entropy_type").unwrap().as_str().unwrap(), "renyi");
}

#[test]
fn test_differential_entropy() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "shannon_entropy".to_string(),
        parameters: to_params(json!({
            "data": [1.0, 2.0, 3.0, 4.0, 5.0],
            "entropy_type": "differential"
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Differential entropy should succeed");

    let result = response.result.unwrap();
    assert_eq!(result.get("entropy_type").unwrap().as_str().unwrap(), "differential");
}

#[test]
fn test_entropy_empty_data() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "shannon_entropy".to_string(),
        parameters: to_params(json!({
            "data": []
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Empty data should fail");
}

#[test]
fn test_entropy_invalid_parameters() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "shannon_entropy".to_string(),
        parameters: to_params(json!({
            // Missing required "data" field
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Missing parameters should fail");
}

// ============================================================================
// MUTUAL INFORMATION TESTS
// ============================================================================

#[test]
fn test_mutual_information_basic() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "mutual_information".to_string(),
        parameters: to_params(json!({
            "x": [1.0, 2.0, 3.0, 4.0, 5.0],
            "y": [2.0, 4.0, 6.0, 8.0, 10.0]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Mutual information should succeed");

    let result = response.result.unwrap();
    assert!(result.get("mutual_information").is_some());
    assert!(result.get("normalized_mi").is_some());
}

#[test]
fn test_mutual_information_with_bins() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "mutual_information".to_string(),
        parameters: to_params(json!({
            "x": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
            "y": [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5],
            "bins": 5
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Mutual information with bins should succeed");
}

#[test]
fn test_mutual_information_mismatched_lengths() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "mutual_information".to_string(),
        parameters: to_params(json!({
            "x": [1.0, 2.0, 3.0],
            "y": [1.0, 2.0]
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Mismatched lengths should fail");
}

// ============================================================================
// CHANNEL CAPACITY TESTS
// ============================================================================

#[test]
fn test_channel_capacity_basic() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "channel_capacity".to_string(),
        parameters: to_params(json!({
            "transition_matrix": [
                [0.9, 0.1],
                [0.1, 0.9]
            ]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Channel capacity should succeed");

    let result = response.result.unwrap();
    assert!(result.get("capacity").is_some());
    assert!(result.get("optimal_input_distribution").is_some());
    assert!(result.get("achievable_rate").is_some());
}

#[test]
fn test_channel_capacity_with_input_distribution() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "channel_capacity".to_string(),
        parameters: to_params(json!({
            "transition_matrix": [
                [0.8, 0.2],
                [0.3, 0.7]
            ],
            "input_distribution": [0.6, 0.4]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Channel capacity with input distribution should succeed");
}

#[test]
fn test_channel_capacity_invalid_matrix() {
    // Matrix rows don't sum to 1
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "channel_capacity".to_string(),
        parameters: to_params(json!({
            "transition_matrix": [
                [0.5, 0.3],
                [0.1, 0.9]
            ]
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Invalid transition matrix should fail");
}

// ============================================================================
// HUFFMAN CODING TESTS
// ============================================================================

#[test]
fn test_huffman_coding_basic() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "huffman_coding".to_string(),
        parameters: to_params(json!({
            "symbols": ["A", "B", "C", "D"],
            "frequencies": [0.4, 0.3, 0.2, 0.1]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Huffman coding should succeed");

    let result = response.result.unwrap();
    assert!(result.get("codes").is_some());
    assert!(result.get("average_length").is_some());
    assert!(result.get("efficiency").is_some());
    assert!(result.get("compression_ratio").is_some());
}

#[test]
fn test_huffman_coding_mismatched_lengths() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "huffman_coding".to_string(),
        parameters: to_params(json!({
            "symbols": ["A", "B", "C"],
            "frequencies": [0.5, 0.5]
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Mismatched symbols and frequencies should fail");
}

#[test]
fn test_huffman_coding_empty() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "huffman_coding".to_string(),
        parameters: to_params(json!({
            "symbols": [],
            "frequencies": []
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Empty symbols should fail");
}

// ============================================================================
// KOLMOGOROV COMPLEXITY TESTS
// ============================================================================

#[test]
fn test_kolmogorov_complexity_basic() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "kolmogorov_complexity".to_string(),
        parameters: to_params(json!({
            "data": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Kolmogorov complexity should succeed");

    let result = response.result.unwrap();
    assert!(result.get("complexity_estimate").is_some());
    assert!(result.get("compressed_size").is_some());
    assert!(result.get("original_size").is_some());
    assert!(result.get("compression_ratio").is_some());
}

#[test]
fn test_kolmogorov_complexity_with_compressor() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "kolmogorov_complexity".to_string(),
        parameters: to_params(json!({
            "data": [1, 2, 3, 4, 5],
            "compressor": "lz77"
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Kolmogorov complexity with compressor should succeed");
}

#[test]
fn test_kolmogorov_complexity_empty() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "kolmogorov_complexity".to_string(),
        parameters: to_params(json!({
            "data": []
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Empty data should fail");
}

// ============================================================================
// CONDITIONAL ENTROPY TESTS
// ============================================================================

#[test]
fn test_conditional_entropy_basic() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "conditional_entropy".to_string(),
        parameters: to_params(json!({
            "x": [1.0, 2.0, 3.0, 4.0, 5.0],
            "y": [2.0, 4.0, 6.0, 8.0, 10.0]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Conditional entropy should succeed");

    let result = response.result.unwrap();
    assert!(result.get("conditional_entropy").is_some());
    assert!(result.get("normalized").is_some());
}

#[test]
fn test_conditional_entropy_with_bins() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "conditional_entropy".to_string(),
        parameters: to_params(json!({
            "x": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            "y": [1.5, 2.5, 3.5, 4.5, 5.5, 6.5],
            "bins": 4
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Conditional entropy with bins should succeed");
}

#[test]
fn test_conditional_entropy_mismatched_lengths() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "conditional_entropy".to_string(),
        parameters: to_params(json!({
            "x": [1.0, 2.0, 3.0],
            "y": [1.0, 2.0]
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Mismatched lengths should fail");
}

// ============================================================================
// RELATIVE ENTROPY (KL DIVERGENCE) TESTS
// ============================================================================

#[test]
fn test_relative_entropy_basic() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "relative_entropy".to_string(),
        parameters: to_params(json!({
            "p": [0.4, 0.3, 0.3],
            "q": [0.33, 0.33, 0.34]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Relative entropy should succeed");

    let result = response.result.unwrap();
    assert!(result.get("kl_divergence").is_some());
    assert!(result.get("js_divergence").is_some());
}

#[test]
fn test_kl_divergence_alias() {
    // Test that "kl_divergence" works as an alias for "relative_entropy"
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "kl_divergence".to_string(),
        parameters: to_params(json!({
            "p": [0.5, 0.5],
            "q": [0.6, 0.4]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "KL divergence alias should succeed");
}

#[test]
fn test_relative_entropy_mismatched_lengths() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "relative_entropy".to_string(),
        parameters: to_params(json!({
            "p": [0.5, 0.5],
            "q": [0.33, 0.33, 0.34]
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Mismatched distribution lengths should fail");
}

#[test]
fn test_relative_entropy_invalid_distribution() {
    // Distributions don't sum to 1
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "relative_entropy".to_string(),
        parameters: to_params(json!({
            "p": [0.5, 0.3],
            "q": [0.6, 0.4]
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Invalid probability distributions should fail");
}

#[test]
fn test_relative_entropy_zero_in_q() {
    // Q has zero probability where P is non-zero
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "relative_entropy".to_string(),
        parameters: to_params(json!({
            "p": [0.5, 0.5],
            "q": [1.0, 0.0]
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Zero in Q where P is non-zero should fail");
}

// ============================================================================
// UNKNOWN OPERATION TESTS
// ============================================================================

#[test]
fn test_unknown_operation() {
    let request = ComputationRequest {
        module: "information_theory".to_string(),
        operation: "unknown_operation".to_string(),
        parameters: to_params(json!({})),
    };

    let response = handle(&request);
    assert!(!response.success, "Unknown operation should fail");
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Unknown information_theory operation"));
}
