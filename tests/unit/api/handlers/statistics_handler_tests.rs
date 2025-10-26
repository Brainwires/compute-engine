//! Unit tests for statistics API handler
//!
//! This test file covers all operations in the statistics handler:
//! - descriptive/statistics - Descriptive statistics (mean, median, mode, std, variance, quartiles, etc.)
//! - monte_carlo/monte_carlo_integration - Monte Carlo integration
//! - mcmc/mcmc_sampling - MCMC sampling
//! - correlation - Correlation analysis (Pearson, Spearman)
//! - kl_divergence - KL and JS divergence
//! - regression, hypothesis_test, distribution - Mock operations (return success)

use crate::api::handlers::statistics::handle;
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
// DESCRIPTIVE STATISTICS TESTS
// ============================================================================

#[test]
fn test_descriptive_statistics_basic() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "descriptive".to_string(),
        parameters: to_params(json!({
            "data": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Descriptive statistics should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("mean").is_some());
    assert!(result.get("median").is_some());
    assert!(result.get("std").is_some());
    assert!(result.get("variance").is_some());
}

#[test]
fn test_statistics_alias() {
    // Test that "statistics" works as an alias for "descriptive"
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "statistics".to_string(),
        parameters: to_params(json!({
            "data": [1.0, 2.0, 3.0, 4.0, 5.0]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Statistics alias should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_statistics_with_quartiles() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "descriptive".to_string(),
        parameters: to_params(json!({
            "data": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
            "operations": ["mean", "median", "quartiles"]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Statistics with quartiles should succeed");

    let result = response.result.unwrap();
    assert!(result.get("quartiles").is_some());

    let quartiles = result.get("quartiles").unwrap();
    assert!(quartiles.get("q1").is_some());
    assert!(quartiles.get("q2").is_some());
    assert!(quartiles.get("q3").is_some());
    assert!(quartiles.get("iqr").is_some());
}

#[test]
fn test_statistics_skewness_kurtosis() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "descriptive".to_string(),
        parameters: to_params(json!({
            "data": [1.0, 2.0, 3.0, 4.0, 5.0],
            "operations": ["skewness", "kurtosis"]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Statistics with skewness/kurtosis should succeed");

    let result = response.result.unwrap();
    assert!(result.get("skewness").is_some());
    assert!(result.get("kurtosis").is_some());
}

#[test]
fn test_statistics_selective_operations() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "descriptive".to_string(),
        parameters: to_params(json!({
            "data": [1.0, 2.0, 3.0, 4.0, 5.0],
            "operations": ["mean", "min", "max"]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Selective operations should succeed");

    let result = response.result.unwrap();
    // These should be computed (not null)
    assert!(!result.get("mean").unwrap().is_null());
    assert!(!result.get("min").unwrap().is_null());
    assert!(!result.get("max").unwrap().is_null());
    // These shouldn't be computed when specific operations are requested (should be null)
    assert!(result.get("mode").unwrap().is_null());
    assert!(result.get("std").unwrap().is_null());
    assert!(result.get("quartiles").unwrap().is_null());
}

#[test]
fn test_statistics_empty_data() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "descriptive".to_string(),
        parameters: to_params(json!({
            "data": []
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Empty data should fail");
    assert!(response.error.is_some());
}

#[test]
fn test_statistics_invalid_parameters() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "descriptive".to_string(),
        parameters: to_params(json!({
            // Missing required "data" field
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Missing parameters should fail");
    assert!(response.error.is_some());
}

// ============================================================================
// MONTE CARLO INTEGRATION TESTS
// ============================================================================

#[test]
fn test_monte_carlo_integration_basic() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "monte_carlo".to_string(),
        parameters: to_params(json!({
            "function": "x^2",
            "bounds": [[0.0, 1.0]],
            "samples": 1000
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Monte Carlo integration should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("integral_estimate").is_some());
    assert!(result.get("standard_error").is_some());
    assert!(result.get("samples_used").is_some());
}

#[test]
fn test_monte_carlo_integration_alias() {
    // Test that "monte_carlo_integration" works as an alias
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "monte_carlo_integration".to_string(),
        parameters: to_params(json!({
            "function": "x^2",
            "bounds": [[0.0, 1.0]],
            "samples": 500
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Monte Carlo integration alias should succeed");
}

#[test]
fn test_monte_carlo_multidimensional() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "monte_carlo".to_string(),
        parameters: to_params(json!({
            "function": "x^2 + y^2",
            "bounds": [[0.0, 1.0], [0.0, 1.0]],
            "samples": 1000
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Multidimensional Monte Carlo should succeed");
}

#[test]
fn test_monte_carlo_empty_bounds() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "monte_carlo".to_string(),
        parameters: to_params(json!({
            "function": "x^2",
            "bounds": [],
            "samples": 1000
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Empty bounds should fail");
    assert!(response.error.is_some());
}

#[test]
fn test_monte_carlo_invalid_parameters() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "monte_carlo".to_string(),
        parameters: to_params(json!({
            // Missing required fields
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Missing parameters should fail");
    assert!(response.error.is_some());
}

// ============================================================================
// MCMC SAMPLING TESTS
// ============================================================================

#[test]
fn test_mcmc_sampling_basic() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "mcmc".to_string(),
        parameters: to_params(json!({
            "target_distribution": "gaussian",
            "dimensions": 2,
            "samples": 100,
            "burn_in": 50,
            "initial_state": [0.0, 0.0]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "MCMC sampling should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("samples").is_some());
    assert!(result.get("acceptance_rate").is_some());
    assert!(result.get("effective_samples").is_some());
}

#[test]
fn test_mcmc_sampling_alias() {
    // Test that "mcmc_sampling" works as an alias
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "mcmc_sampling".to_string(),
        parameters: to_params(json!({
            "target_distribution": "gaussian",
            "dimensions": 1,
            "samples": 50,
            "burn_in": 25,
            "initial_state": [0.0]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "MCMC sampling alias should succeed");
}

#[test]
fn test_mcmc_zero_dimensions() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "mcmc".to_string(),
        parameters: to_params(json!({
            "target_distribution": "gaussian",
            "dimensions": 0,
            "samples": 100,
            "burn_in": 50,
            "initial_state": []
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Zero dimensions should fail");
    assert!(response.error.is_some());
}

#[test]
fn test_mcmc_invalid_initial_state() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "mcmc".to_string(),
        parameters: to_params(json!({
            "target_distribution": "gaussian",
            "dimensions": 2,
            "samples": 100,
            "burn_in": 50,
            "initial_state": [0.0] // Wrong dimensions
        })),
    };

    let response = handle(&request);
    // Should succeed and correct the initial state internally
    assert!(response.success, "MCMC should correct initial state dimensions");
}

#[test]
fn test_mcmc_invalid_parameters() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "mcmc".to_string(),
        parameters: to_params(json!({
            // Missing required fields
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Missing parameters should fail");
    assert!(response.error.is_some());
}

// ============================================================================
// CORRELATION TESTS
// ============================================================================

#[test]
fn test_correlation_pearson() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "correlation".to_string(),
        parameters: to_params(json!({
            "x": [1.0, 2.0, 3.0, 4.0, 5.0],
            "y": [2.0, 4.0, 6.0, 8.0, 10.0],
            "method": "pearson"
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Pearson correlation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("correlation").is_some());
    assert_eq!(result.get("method").unwrap().as_str().unwrap(), "pearson");

    // Perfect positive correlation should be close to 1.0
    let corr = result.get("correlation").unwrap().as_f64().unwrap();
    assert!((corr - 1.0).abs() < 0.001, "Perfect correlation should be ~1.0");
}

#[test]
fn test_correlation_spearman() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "correlation".to_string(),
        parameters: to_params(json!({
            "x": [1.0, 2.0, 3.0, 4.0, 5.0],
            "y": [1.0, 2.0, 3.0, 4.0, 5.0],
            "method": "spearman"
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Spearman correlation should succeed");

    let result = response.result.unwrap();
    assert_eq!(result.get("method").unwrap().as_str().unwrap(), "spearman");
}

#[test]
fn test_correlation_mismatched_lengths() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "correlation".to_string(),
        parameters: to_params(json!({
            "x": [1.0, 2.0, 3.0],
            "y": [1.0, 2.0],
            "method": "pearson"
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Mismatched lengths should fail");
    assert!(response.error.is_some());
}

#[test]
fn test_correlation_empty_data() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "correlation".to_string(),
        parameters: to_params(json!({
            "x": [],
            "y": [],
            "method": "pearson"
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Empty data should fail");
    assert!(response.error.is_some());
}

#[test]
fn test_correlation_unsupported_method() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "correlation".to_string(),
        parameters: to_params(json!({
            "x": [1.0, 2.0, 3.0],
            "y": [1.0, 2.0, 3.0],
            "method": "kendall"
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Unsupported method should fail");
    assert!(response.error.is_some());
}

#[test]
fn test_correlation_invalid_parameters() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "correlation".to_string(),
        parameters: to_params(json!({
            // Missing required fields
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Missing parameters should fail");
    assert!(response.error.is_some());
}

// ============================================================================
// KL DIVERGENCE TESTS
// ============================================================================

#[test]
fn test_kl_divergence_basic() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "kl_divergence".to_string(),
        parameters: to_params(json!({
            "p": [0.5, 0.3, 0.2],
            "q": [0.4, 0.4, 0.2]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "KL divergence should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("kl_divergence").is_some());
    assert!(result.get("js_divergence").is_some());
}

#[test]
fn test_kl_divergence_identical_distributions() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "kl_divergence".to_string(),
        parameters: to_params(json!({
            "p": [0.5, 0.3, 0.2],
            "q": [0.5, 0.3, 0.2]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "KL divergence should succeed");

    let result = response.result.unwrap();
    let kl = result.get("kl_divergence").unwrap().as_f64().unwrap();
    let js = result.get("js_divergence").unwrap().as_f64().unwrap();

    // KL and JS should be close to 0 for identical distributions
    assert!(kl.abs() < 1e-6, "KL divergence of identical distributions should be ~0");
    assert!(js.abs() < 1e-6, "JS divergence of identical distributions should be ~0");
}

#[test]
fn test_kl_divergence_mismatched_lengths() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "kl_divergence".to_string(),
        parameters: to_params(json!({
            "p": [0.5, 0.3, 0.2],
            "q": [0.5, 0.5]
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Mismatched lengths should fail");
    assert!(response.error.is_some());
}

#[test]
fn test_kl_divergence_negative_probability() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "kl_divergence".to_string(),
        parameters: to_params(json!({
            "p": [0.5, -0.3, 0.8],
            "q": [0.4, 0.4, 0.2]
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Negative probabilities should fail");
    assert!(response.error.is_some());
}

#[test]
fn test_kl_divergence_invalid_parameters() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "kl_divergence".to_string(),
        parameters: to_params(json!({
            // Missing required fields
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Missing parameters should fail");
    assert!(response.error.is_some());
}

// ============================================================================
// MOCK OPERATIONS TESTS
// ============================================================================

#[test]
fn test_regression_mock() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "regression".to_string(),
        parameters: to_params(json!({})),
    };

    let response = handle(&request);
    assert!(response.success, "Regression mock should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert_eq!(result.get("result").unwrap().as_str().unwrap(), "computed");
}

#[test]
fn test_hypothesis_test_mock() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "hypothesis_test".to_string(),
        parameters: to_params(json!({})),
    };

    let response = handle(&request);
    assert!(response.success, "Hypothesis test mock should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_distribution_mock() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "distribution".to_string(),
        parameters: to_params(json!({})),
    };

    let response = handle(&request);
    assert!(response.success, "Distribution mock should succeed");
    assert!(response.result.is_some());
}

// ============================================================================
// ERROR HANDLING TESTS
// ============================================================================

#[test]
fn test_unknown_operation() {
    let request = ComputationRequest {
        module: "statistics".to_string(),
        operation: "unknown_operation".to_string(),
        parameters: to_params(json!({})),
    };

    let response = handle(&request);
    assert!(!response.success, "Unknown operation should fail");
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Unknown operation"));
}
