//! Unit tests for linear algebra API handler
//!
//! Tests all handler operations including:
//! - SVD computation
//! - Eigendecomposition
//! - PCA
//! - Matrix rank
//! - Pseudoinverse
//! - Matrix operations
//! - Matrix norms
//! - Matrix power
//! - Matrix exponential
//! - Cholesky decomposition

use super::*; // Import from parent module (linear_algebra handler)
use crate::api::types::ComputationRequest;
use serde_json::json;
use std::collections::HashMap;

fn create_request(operation: &str, parameters: serde_json::Value) -> ComputationRequest {
    let params_map: HashMap<String, serde_json::Value> = serde_json::from_value(parameters)
        .expect("Failed to convert parameters to HashMap");

    ComputationRequest {
        module: "linear_algebra".to_string(),
        operation: operation.to_string(),
        parameters: params_map,
    }
}

#[test]
fn test_compute_svd_success() {
    let params = json!({
        "matrix": [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
    });
    let request = create_request("compute_svd", params);
    let response = handle(&request);

    assert!(response.success, "SVD computation should succeed");
    assert!(response.result.is_some(), "Should return SVD result");

    let result = response.result.unwrap();
    assert!(result.get("u").is_some(), "Should have U matrix");
    assert!(result.get("singular_values").is_some(), "Should have singular values");
    assert!(result.get("v_transpose").is_some(), "Should have V^T matrix");
    assert!(result.get("rank").is_some(), "Should have rank");
    assert!(result.get("condition_number").is_some(), "Should have condition number");
}

#[test]
fn test_compute_svd_invalid_params() {
    let params = json!({
        "invalid_field": "test"
    });
    let request = create_request("compute_svd", params);
    let response = handle(&request);

    assert!(!response.success, "Should fail with invalid parameters");
    assert!(response.error.is_some(), "Should return error message");
    assert!(response.error.unwrap().contains("Invalid parameters"));
}

#[test]
fn test_eigendecomposition_success() {
    let params = json!({
        "matrix": [[4.0, 2.0], [2.0, 3.0]]
    });
    let request = create_request("compute_eigendecomposition", params);
    let response = handle(&request);

    assert!(response.success, "Eigendecomposition should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("eigenvalues").is_some(), "Should have eigenvalues");
    assert!(result.get("eigenvectors").is_some(), "Should have eigenvectors");
}

#[test]
fn test_eigendecomposition_alias() {
    // Test alias "eigendecomposition" (without "compute_" prefix)
    let params = json!({
        "matrix": [[2.0, 1.0], [1.0, 2.0]]
    });
    let request = create_request("eigendecomposition", params);
    let response = handle(&request);

    assert!(response.success, "Eigendecomposition alias should work");
}

#[test]
fn test_eigendecomposition_non_square() {
    let params = json!({
        "matrix": [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
    });
    let request = create_request("compute_eigendecomposition", params);
    let response = handle(&request);

    assert!(!response.success, "Should fail for non-square matrix");
    assert!(response.error.is_some());
}

#[test]
fn test_pca_success() {
    let params = json!({
        "matrix": [
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0]
        ]
    });
    let request = create_request("compute_pca", params);
    let response = handle(&request);

    assert!(response.success, "PCA should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("principal_components").is_some());
    assert!(result.get("explained_variance").is_some());
    assert!(result.get("cumulative_variance").is_some());
    assert!(result.get("transformed_data").is_some());
}

#[test]
fn test_pca_with_n_components() {
    let params = json!({
        "matrix": [
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0]
        ],
        "n_components": 2
    });
    let request = create_request("pca", params);
    let response = handle(&request);

    assert!(response.success, "PCA with n_components should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_matrix_rank_success() {
    let params = json!({
        "matrix": [[1.0, 2.0], [3.0, 4.0]]
    });
    let request = create_request("compute_matrix_rank", params);
    let response = handle(&request);

    assert!(response.success, "Matrix rank computation should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_matrix_rank_alias() {
    let params = json!({
        "matrix": [[1.0, 0.0], [0.0, 1.0]]
    });
    let request = create_request("matrix_rank", params);
    let response = handle(&request);

    assert!(response.success, "Matrix rank alias should work");
}

#[test]
fn test_pseudoinverse_success() {
    let params = json!({
        "matrix": [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
    });
    let request = create_request("compute_pseudoinverse", params);
    let response = handle(&request);

    assert!(response.success, "Pseudoinverse should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_pseudoinverse_alias() {
    let params = json!({
        "matrix": [[1.0, 2.0], [3.0, 4.0]]
    });
    let request = create_request("pseudoinverse", params);
    let response = handle(&request);

    assert!(response.success, "Pseudoinverse alias should work");
}

#[test]
fn test_matrix_operations_success() {
    let params = json!({
        "matrix": [[1.0, 2.0], [3.0, 4.0]]
    });
    let request = create_request("matrix_operations", params);
    let response = handle(&request);

    assert!(response.success, "Matrix operations should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("determinant").is_some());
    assert!(result.get("trace").is_some());
    assert!(result.get("frobenius_norm").is_some());
    assert!(result.get("max_norm").is_some());
}

#[test]
fn test_matrix_operations_alias() {
    let params = json!({
        "matrix": [[5.0, 6.0], [7.0, 8.0]]
    });
    let request = create_request("matrix_ops", params);
    let response = handle(&request);

    assert!(response.success, "Matrix operations alias should work");
}

#[test]
fn test_matrix_norm_frobenius() {
    let params = json!({
        "matrix": [[1.0, 2.0], [3.0, 4.0]],
        "norm_type": "frobenius"
    });
    let request = create_request("matrix_norm", params);
    let response = handle(&request);

    assert!(response.success, "Frobenius norm should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("norm").is_some());
    assert!(result.get("norm_type").is_some());
}

#[test]
fn test_matrix_norm_1() {
    let params = json!({
        "matrix": [[1.0, 2.0], [3.0, 4.0]],
        "norm_type": "1"
    });
    let request = create_request("norm", params);
    let response = handle(&request);

    assert!(response.success, "1-norm should succeed");
}

#[test]
fn test_matrix_norm_infinity() {
    let params = json!({
        "matrix": [[1.0, 2.0], [3.0, 4.0]],
        "norm_type": "inf"
    });
    let request = create_request("matrix_norm", params);
    let response = handle(&request);

    assert!(response.success, "Infinity norm should succeed");
}

#[test]
fn test_matrix_norm_invalid_type() {
    let params = json!({
        "matrix": [[1.0, 2.0], [3.0, 4.0]],
        "norm_type": "invalid"
    });
    let request = create_request("matrix_norm", params);
    let response = handle(&request);

    assert!(!response.success, "Should fail with invalid norm type");
    assert!(response.error.is_some());
}

#[test]
fn test_matrix_power_positive() {
    let params = json!({
        "matrix": [[2.0, 0.0], [0.0, 3.0]],
        "power": 3
    });
    let request = create_request("matrix_power", params);
    let response = handle(&request);

    assert!(response.success, "Matrix power should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_matrix_power_zero() {
    let params = json!({
        "matrix": [[1.0, 2.0], [3.0, 4.0]],
        "power": 0
    });
    let request = create_request("power", params);
    let response = handle(&request);

    assert!(response.success, "Matrix^0 should return identity");
}

#[test]
fn test_matrix_power_non_square() {
    let params = json!({
        "matrix": [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
        "power": 2
    });
    let request = create_request("matrix_power", params);
    let response = handle(&request);

    assert!(!response.success, "Should fail for non-square matrix");
}

#[test]
fn test_matrix_exp_default_terms() {
    let params = json!({
        "matrix": [[0.0, 1.0], [-1.0, 0.0]]
    });
    let request = create_request("matrix_exp", params);
    let response = handle(&request);

    assert!(response.success, "Matrix exponential should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_matrix_exp_custom_terms() {
    let params = json!({
        "matrix": [[1.0, 0.0], [0.0, 1.0]],
        "terms": 10
    });
    let request = create_request("exp", params);
    let response = handle(&request);

    assert!(response.success, "Matrix exponential with custom terms should succeed");
}

#[test]
fn test_cholesky_success() {
    let params = json!({
        "matrix": [[4.0, 2.0], [2.0, 3.0]]
    });
    let request = create_request("cholesky", params);
    let response = handle(&request);

    assert!(response.success, "Cholesky decomposition should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("lower").is_some(), "Should have lower triangular matrix");
}

#[test]
fn test_cholesky_alias() {
    let params = json!({
        "matrix": [[9.0, 3.0], [3.0, 2.0]]
    });
    let request = create_request("cholesky_decomposition", params);
    let response = handle(&request);

    assert!(response.success, "Cholesky decomposition alias should work");
}

#[test]
fn test_cholesky_not_positive_definite() {
    let params = json!({
        "matrix": [[1.0, 2.0], [2.0, 1.0]]
    });
    let request = create_request("cholesky", params);
    let response = handle(&request);

    assert!(!response.success, "Should fail for non-positive-definite matrix");
    assert!(response.error.is_some());
}

#[test]
fn test_unknown_operation() {
    let params = json!({
        "matrix": [[1.0, 2.0], [3.0, 4.0]]
    });
    let request = create_request("unknown_operation", params);
    let response = handle(&request);

    assert!(!response.success, "Should fail for unknown operation");
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Unknown operation"));
}

#[test]
fn test_empty_parameters() {
    let params = json!({});
    let request = create_request("compute_svd", params);
    let response = handle(&request);

    assert!(!response.success, "Should fail with empty parameters");
}

#[test]
fn test_response_structure() {
    let params = json!({
        "matrix": [[1.0, 2.0], [3.0, 4.0]]
    });
    let request = create_request("matrix_operations", params);
    let response = handle(&request);

    // Verify response structure
    assert_eq!(response.module, "linear_algebra");
    assert_eq!(response.operation, "matrix_operations");
    assert!(response.success);
    assert!(response.result.is_some());
    assert!(response.error.is_none());
}
