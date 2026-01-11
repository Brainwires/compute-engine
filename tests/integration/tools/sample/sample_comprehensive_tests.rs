//! Comprehensive SAMPLE tool test suite
//!
//! Tests for COMPUTE tool sampling operations including:
//! - Monte Carlo (Integration, MCMC, Metropolis-Hastings)
//! - Statistical Methods (BasicStats, HypothesisTest, Correlation)
//! - Signal Analysis

use computational_engine::create_default_dispatcher;
use computational_engine::engine::equations::*;
use computational_engine::engine::*;
use serde_json::json;
use std::collections::HashMap;

// ============================================================================
// MONTE CARLO TESTS
// ============================================================================

#[test]
fn test_monte_carlo_integration() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Sample(SamplingMethod::MonteCarlo(MonteCarloMethod::Integration)),
        data: json!({
            "num_samples": 10000
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform Monte Carlo integration: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        assert!(output.result.get("integral").is_some(), "Should have integral estimate");
    }
}

#[test]
fn test_mcmc_sampling() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Sample(SamplingMethod::MonteCarlo(MonteCarloMethod::MCMC)),
        data: json!({
            "num_samples": 1000
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform MCMC sampling: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        assert!(output.result.get("samples").is_some(), "Should have samples");
        assert!(output.result.get("mean").is_some(), "Should have mean");
    }
}

#[test]
fn test_metropolis_hastings() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Sample(SamplingMethod::MonteCarlo(MonteCarloMethod::MetropolisHastings)),
        data: json!({
            "num_samples": 1000
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform Metropolis-Hastings sampling: {:?}", result);
}

// ============================================================================
// STATISTICAL METHOD TESTS
// ============================================================================

#[test]
fn test_basic_stats() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Sample(SamplingMethod::Stats(StatisticalMethod::BasicStats)),
        data: json!({
            "data": data
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute basic stats: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        assert!(output.result.get("mean").is_some(), "Should have mean");
        assert!(output.result.get("variance").is_some(), "Should have variance");
        assert!(output.result.get("std_dev").is_some(), "Should have std_dev");
    }
}

#[test]
fn test_hypothesis_test() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Sample(SamplingMethod::Stats(StatisticalMethod::HypothesisTest)),
        data: json!({
            "data": data,
            "null_hypothesis_mean": 3.0
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform hypothesis test: {:?}", result);
}

#[test]
fn test_correlation() {
    let dispatcher = create_default_dispatcher();

    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y = vec![2.0, 4.0, 6.0, 8.0, 10.0];

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Sample(SamplingMethod::Stats(StatisticalMethod::Correlation)),
        data: json!({
            "x": x,
            "y": y
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute correlation: {:?}", result);
}

// ============================================================================
// SIGNAL ANALYSIS TESTS
// ============================================================================

#[test]
fn test_spectral_analysis() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..128).map(|i| (i as f64 * 0.1).sin()).collect();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Sample(SamplingMethod::SignalAnalysis(SignalMethod::SpectralAnalysis)),
        data: json!({
            "data": data
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform spectral analysis: {:?}", result);
}

#[test]
fn test_autocorrelation() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..64).map(|i| (i as f64 * 0.1).sin()).collect();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Sample(SamplingMethod::SignalAnalysis(SignalMethod::Autocorrelation)),
        data: json!({
            "data": data
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute autocorrelation: {:?}", result);
}

#[test]
fn test_power_spectrum() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..128).map(|i| (i as f64 * 0.1).sin()).collect();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Sample(SamplingMethod::SignalAnalysis(SignalMethod::PowerSpectrum)),
        data: json!({
            "data": data
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute power spectrum: {:?}", result);
}

// ============================================================================
// PATH GENERATION TESTS
// ============================================================================

#[test]
fn test_path_generation() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Sample(SamplingMethod::PathGeneration),
        data: json!({
            "num_steps": 100,
            "num_paths": 10
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should generate paths: {:?}", result);
}

// ============================================================================
// MOMENTS TESTS
// ============================================================================

#[test]
fn test_moments() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Sample(SamplingMethod::Moments),
        data: json!({
            "data": data
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute moments: {:?}", result);
}
