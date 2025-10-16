//! Comprehensive SAMPLE tool test suite
//!
//! Tests for all SAMPLE operations including:
//! - Monte Carlo (Integration, MCMC, Metropolis-Hastings, Gibbs)
//! - Statistical Methods (BasicStats, HypothesisTest, ANOVA, Regression, TimeSeries, Correlation)
//! - Signal Analysis (Spectral, Autocorrelation, CrossCorrelation, PowerSpectrum, Coherence, Cepstrum, PeakDetection)
//! - Path Generation and Moments

use computational_engine::create_default_dispatcher;
use computational_engine::engine::*;
use std::collections::HashMap;

// ============================================================================
// MONTE CARLO TESTS
// ============================================================================

#[test]
fn test_monte_carlo_integration() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("domain_min".to_string(), serde_json::json!(0.0));
    params.insert("domain_max".to_string(), serde_json::json!(1.0));

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::MonteCarlo(MonteCarloMethod::Integration),
        data: vec![],
        num_samples: Some(10000),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should perform Monte Carlo integration: {:?}",
        result
    );

    if let Ok(ToolResponse::Sample(output)) = result {
        let res = output.result.as_object().unwrap();
        assert!(
            res.contains_key("integral"),
            "Should have integral estimate"
        );
        assert!(res.contains_key("error"), "Should have error estimate");
    }
}

#[test]
fn test_mcmc_sampling() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_state".to_string(), serde_json::json!(0.0));
    params.insert("step_size".to_string(), serde_json::json!(1.0));

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::MonteCarlo(MonteCarloMethod::MCMC),
        data: vec![],
        num_samples: Some(1000),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform MCMC sampling: {:?}", result);

    if let Ok(ToolResponse::Sample(output)) = result {
        assert!(output.moments.is_some(), "Should have statistical moments");
        let moments = output.moments.unwrap();
        assert!(moments.contains_key("mean"), "Should have mean");
        assert!(moments.contains_key("variance"), "Should have variance");
    }
}

#[test]
fn test_metropolis_hastings() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_state".to_string(), serde_json::json!(0.0));
    params.insert("proposal_std".to_string(), serde_json::json!(1.0));

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::MonteCarlo(MonteCarloMethod::MetropolisHastings),
        data: vec![],
        num_samples: Some(1000),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should perform Metropolis-Hastings: {:?}",
        result
    );

    if let Ok(ToolResponse::Sample(output)) = result {
        let metadata = output.metadata.as_ref().unwrap();
        assert!(
            metadata.get("acceptance_rate").is_some(),
            "Should have acceptance rate"
        );
        assert_eq!(
            metadata.get("method").and_then(|v| v.as_str()).unwrap(),
            "metropolis_hastings"
        );
    }
}

#[test]
fn test_gibbs_sampling() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("correlation".to_string(), serde_json::json!(0.5));

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::MonteCarlo(MonteCarloMethod::Gibbs),
        data: vec![],
        num_samples: Some(1000),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should perform Gibbs sampling: {:?}",
        result
    );

    if let Ok(ToolResponse::Sample(output)) = result {
        let res = output.result.as_object().unwrap();
        assert!(res.contains_key("x"), "Should have x samples");
        assert!(res.contains_key("y"), "Should have y samples");
    }
}

// ============================================================================
// STOCHASTIC PATH GENERATION TESTS
// ============================================================================

#[test]
fn test_path_generation_brownian() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("t_max".to_string(), serde_json::json!(1.0));
    params.insert("initial_value".to_string(), serde_json::json!(0.0));
    params.insert("drift".to_string(), serde_json::json!(0.0));
    params.insert("volatility".to_string(), serde_json::json!(1.0));

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::PathGeneration,
        data: vec![],
        num_samples: Some(100),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should generate Brownian motion path: {:?}",
        result
    );

    if let Ok(ToolResponse::Sample(output)) = result {
        assert!(output.moments.is_some(), "Should have moments");
        let metadata = output.metadata.as_ref().unwrap();
        assert!(metadata.get("times").is_some(), "Should have time points");
    }
}

#[test]
fn test_moments_calculation() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::Moments,
        data,
        num_samples: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate moments: {:?}", result);

    if let Ok(ToolResponse::Sample(output)) = result {
        let moments = output.moments.as_ref().unwrap();
        assert!(moments.contains_key("mean"), "Should have mean");
        assert!(moments.contains_key("variance"), "Should have variance");
        assert!(moments.contains_key("std_dev"), "Should have std dev");
        assert!(moments.contains_key("skewness"), "Should have skewness");
        assert!(moments.contains_key("kurtosis"), "Should have kurtosis");

        // Verify mean calculation
        let mean = moments.get("mean").unwrap();
        assert!((mean - 5.5).abs() < 0.01, "Mean should be 5.5");
    }
}

// ============================================================================
// STATISTICAL METHODS TESTS
// ============================================================================

#[test]
fn test_basic_stats() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::Stats(StatisticalMethod::BasicStats),
        data,
        num_samples: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should compute basic statistics: {:?}",
        result
    );
}

#[test]
fn test_hypothesis_test() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
    let mut params = HashMap::new();
    params.insert("null_mean".to_string(), serde_json::json!(5.0));

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::Stats(StatisticalMethod::HypothesisTest),
        data,
        num_samples: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should perform hypothesis test: {:?}",
        result
    );

    if let Ok(ToolResponse::Sample(output)) = result {
        let res = output.result.as_object().unwrap();
        assert!(res.contains_key("t_statistic"), "Should have t-statistic");
        assert!(res.contains_key("p_value"), "Should have p-value");
    }
}

#[test]
fn test_anova() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::Stats(StatisticalMethod::ANOVA),
        data,
        num_samples: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should handle ANOVA request: {:?}", result);
}

#[test]
fn test_linear_regression() {
    let dispatcher = create_default_dispatcher();

    // Data for y = 2x + 1: [1, 3, 5, 7, 9]
    let data = vec![1.0, 3.0, 5.0, 7.0, 9.0];

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::Stats(StatisticalMethod::Regression),
        data,
        num_samples: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should perform linear regression: {:?}",
        result
    );

    if let Ok(ToolResponse::Sample(output)) = result {
        let res = output.result.as_object().unwrap();
        assert!(res.contains_key("slope"), "Should have slope");
        assert!(res.contains_key("intercept"), "Should have intercept");

        let slope = res.get("slope").and_then(|v| v.as_f64()).unwrap();
        let intercept = res.get("intercept").and_then(|v| v.as_f64()).unwrap();

        // Should be close to y = 2x + 1
        assert!((slope - 2.0).abs() < 0.01, "Slope should be ~2.0");
        assert!((intercept - 1.0).abs() < 0.01, "Intercept should be ~1.0");
    }
}

#[test]
fn test_time_series_analysis() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::Stats(StatisticalMethod::TimeSeries),
        data,
        num_samples: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should perform time series analysis: {:?}",
        result
    );

    if let Ok(ToolResponse::Sample(output)) = result {
        let res = output.result.as_object().unwrap();
        assert!(res.contains_key("trend"), "Should have trend");
    }
}

#[test]
fn test_correlation() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::Stats(StatisticalMethod::Correlation),
        data,
        num_samples: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute correlation: {:?}", result);

    if let Ok(ToolResponse::Sample(output)) = result {
        let res = output.result.as_object().unwrap();
        assert!(
            res.contains_key("correlation"),
            "Should have correlation coefficient"
        );
        assert!(res.contains_key("covariance"), "Should have covariance");
    }
}

// ============================================================================
// SIGNAL ANALYSIS TESTS
// ============================================================================

#[test]
fn test_spectral_analysis() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..64)
        .map(|i| (2.0 * std::f64::consts::PI * i as f64 / 64.0).sin())
        .collect();

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::SignalAnalysis(SignalMethod::SpectralAnalysis),
        data,
        num_samples: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should handle spectral analysis: {:?}",
        result
    );
}

#[test]
fn test_autocorrelation() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..100)
        .map(|i| (2.0 * std::f64::consts::PI * i as f64 / 20.0).sin())
        .collect();

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::SignalAnalysis(SignalMethod::Autocorrelation),
        data,
        num_samples: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should compute autocorrelation: {:?}",
        result
    );

    if let Ok(ToolResponse::Sample(output)) = result {
        let acf = output.result.as_array().unwrap();
        assert!(!acf.is_empty(), "Should have autocorrelation values");
        // First lag should be 1.0 (perfect correlation with itself)
        let first_lag = acf[0].as_f64().unwrap();
        assert!((first_lag - 1.0).abs() < 0.01, "ACF at lag 0 should be 1.0");
    }
}

#[test]
fn test_cross_correlation() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::SignalAnalysis(SignalMethod::CrossCorrelation),
        data,
        num_samples: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should handle cross-correlation: {:?}",
        result
    );
}

#[test]
fn test_power_spectrum() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..64)
        .map(|i| (2.0 * std::f64::consts::PI * i as f64 / 64.0).sin())
        .collect();

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::SignalAnalysis(SignalMethod::PowerSpectrum),
        data,
        num_samples: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should compute power spectrum: {:?}",
        result
    );

    if let Ok(ToolResponse::Sample(output)) = result {
        let res = output.result.as_object().unwrap();
        assert!(res.contains_key("total_power"), "Should have total power");
    }
}

#[test]
fn test_coherence() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::SignalAnalysis(SignalMethod::Coherence),
        data,
        num_samples: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should handle coherence: {:?}", result);
}

#[test]
fn test_cepstrum() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::SignalAnalysis(SignalMethod::Cepstrum),
        data,
        num_samples: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should handle cepstrum: {:?}", result);
}

#[test]
fn test_peak_detection() {
    let dispatcher = create_default_dispatcher();

    // Signal with peaks at indices 2, 5, 8
    let data = vec![1.0, 2.0, 5.0, 2.0, 3.0, 7.0, 3.0, 4.0, 6.0, 2.0];

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::SignalAnalysis(SignalMethod::PeakDetection),
        data,
        num_samples: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should detect peaks: {:?}", result);

    if let Ok(ToolResponse::Sample(output)) = result {
        let peaks = output.result.as_array().unwrap();
        assert!(peaks.len() >= 2, "Should find at least 2 peaks");
    }
}

// ============================================================================
// ERROR HANDLING TESTS
// ============================================================================

#[test]
fn test_moments_requires_data() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::Moments,
        data: vec![], // Empty data
        num_samples: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Moments should require non-empty data");
}

#[test]
fn test_regression_requires_multiple_points() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Sample(SampleInput {
        method: SamplingMethod::Stats(StatisticalMethod::Regression),
        data: vec![1.0], // Only 1 point
        num_samples: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_err(),
        "Regression should require at least 2 points"
    );
}

// ============================================================================
// INTEGRATION TESTS (Multiple operations)
// ============================================================================

#[test]
fn test_monte_carlo_convergence() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("domain_min".to_string(), serde_json::json!(0.0));
    params.insert("domain_max".to_string(), serde_json::json!(1.0));

    // Test with increasing sample sizes
    for num_samples in [100, 1000, 10000] {
        let request = ToolRequest::Sample(SampleInput {
            method: SamplingMethod::MonteCarlo(MonteCarloMethod::Integration),
            data: vec![],
            num_samples: Some(num_samples),
            parameters: params.clone(),
        });

        let result = dispatcher.dispatch(request);
        assert!(
            result.is_ok(),
            "Should succeed with {} samples",
            num_samples
        );
    }
}

#[test]
fn test_gibbs_correlation_effect() {
    let dispatcher = create_default_dispatcher();

    // Test with different correlations
    for rho in [0.0, 0.5, 0.9] {
        let mut params = HashMap::new();
        params.insert("correlation".to_string(), serde_json::json!(rho));

        let request = ToolRequest::Sample(SampleInput {
            method: SamplingMethod::MonteCarlo(MonteCarloMethod::Gibbs),
            data: vec![],
            num_samples: Some(500),
            parameters: params,
        });

        let result = dispatcher.dispatch(request);
        assert!(result.is_ok(), "Should sample with correlation {}", rho);
    }
}
