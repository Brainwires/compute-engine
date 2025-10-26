//! Unit tests for specialized::statistics

use crate::specialized::statistics::*;

#[test]
fn test_statistics_basic() {
    let request = StatisticsRequest {
        data: vec![1.0, 2.0, 3.0, 4.0, 5.0],
        operations: Some(vec!["mean".to_string(), "median".to_string()]),
    };

    let result = statistics(request).unwrap();
    assert!(result.mean.is_some());
    assert!((result.mean.unwrap() - 3.0).abs() < 1e-10);
}

#[test]
fn test_statistics_median_odd() {
    let request = StatisticsRequest {
        data: vec![1.0, 3.0, 2.0],
        operations: Some(vec!["median".to_string()]),
    };

    let result = statistics(request).unwrap();
    assert_eq!(result.median.unwrap(), 2.0);
}

#[test]
fn test_statistics_median_even() {
    let request = StatisticsRequest {
        data: vec![1.0, 2.0, 3.0, 4.0],
        operations: Some(vec!["median".to_string()]),
    };

    let result = statistics(request).unwrap();
    assert_eq!(result.median.unwrap(), 2.5);
}

#[test]
fn test_statistics_variance() {
    let request = StatisticsRequest {
        data: vec![2.0, 4.0, 6.0],
        operations: Some(vec!["variance".to_string()]),
    };

    let result = statistics(request).unwrap();
    assert!(result.variance.is_some());
}

#[test]
fn test_correlation_positive() {
    let request = CorrelationRequest {
        x: vec![1.0, 2.0, 3.0],
        y: vec![2.0, 4.0, 6.0],
        method: "pearson".to_string(),
    };

    let result = correlation(request).unwrap();
    assert!((result.correlation - 1.0).abs() < 1e-10); // Perfect positive correlation
}

#[test]
fn test_correlation_negative() {
    let request = CorrelationRequest {
        x: vec![1.0, 2.0, 3.0],
        y: vec![3.0, 2.0, 1.0],
        method: "pearson".to_string(),
    };

    let result = correlation(request).unwrap();
    assert!((result.correlation - (-1.0)).abs() < 1e-10); // Perfect negative correlation
}

#[test]
fn test_kl_divergence_identical() {
    let request = KLDivergenceRequest {
        p: vec![0.5, 0.5],
        q: vec![0.5, 0.5],
    };

    let result = kl_divergence(request).unwrap();
    assert!(result.kl_divergence.abs() < 1e-10); // Identical distributions have KL=0
}

#[test]
fn test_kl_divergence_different() {
    let request = KLDivergenceRequest {
        p: vec![0.9, 0.1],
        q: vec![0.5, 0.5],
    };

    let result = kl_divergence(request).unwrap();
    assert!(result.kl_divergence > 0.0); // Different distributions have KL>0
}

#[test]
fn test_monte_carlo_integration_basic() {
    let request = MonteCarloRequest {
        function: "x^2".to_string(),
        bounds: vec![(0.0, 1.0)],
        samples: 1000,
    };

    let result = monte_carlo_integration(request).unwrap();
    // Integral of x^2 from 0 to 1 is 1/3 ≈ 0.333
    assert!((result.integral_estimate - 0.333).abs() < 0.1);
}

#[test]
fn test_mcmc_sampling_basic() {
    let request = MCMCRequest {
        target_distribution: "gaussian".to_string(),
        dimensions: 1,
        samples: 100,
        burn_in: 10,
        initial_state: vec![0.0],
    };

    let result = mcmc_sampling(request).unwrap();
    assert_eq!(result.samples.len(), 100);
    assert!(result.acceptance_rate >= 0.0 && result.acceptance_rate <= 1.0);
}
