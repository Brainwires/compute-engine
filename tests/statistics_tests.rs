use computational_engine::statistics::*;

#[test]
fn test_statistics() {
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let result = statistics(StatisticsRequest {
        data,
        operations: None,
    }).unwrap();

    assert!(result.mean.is_some());
    assert!(result.median.is_some());
    assert!(result.variance.is_some());
    assert!(result.std.is_some());
}

#[test]
fn test_statistics_mean_median_mode() {
    let data = vec![1.0, 2.0, 2.0, 3.0, 4.0, 5.0];

    let result = statistics(StatisticsRequest {
        data,
        operations: Some(vec!["mean".to_string(), "median".to_string(), "mode".to_string()]),
    }).unwrap();

    assert_eq!(result.mean.unwrap(), 2.8333333333333335);
    assert_eq!(result.median.unwrap(), 2.5);
    assert!(result.mode.unwrap().contains(&2.0));
}

#[test]
fn test_statistics_quartiles() {
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];

    let result = statistics(StatisticsRequest {
        data,
        operations: Some(vec!["quartiles".to_string()]),
    }).unwrap();

    let quartiles = result.quartiles.unwrap();
    assert!(quartiles.q1 > 0.0);
    assert!(quartiles.q2 > 0.0);
    assert!(quartiles.q3 > 0.0);
    assert!(quartiles.iqr > 0.0);
}

#[test]
fn test_statistics_skewness_kurtosis() {
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let result = statistics(StatisticsRequest {
        data,
        operations: Some(vec!["skewness".to_string(), "kurtosis".to_string()]),
    }).unwrap();

    assert!(result.skewness.is_some());
    assert!(result.kurtosis.is_some());
}

#[test]
fn test_statistics_min_max_sum_count() {
    let data = vec![1.0, 5.0, 3.0, 9.0, 2.0];

    let result = statistics(StatisticsRequest {
        data: data.clone(),
        operations: Some(vec!["min".to_string(), "max".to_string(), "sum".to_string(), "count".to_string()]),
    }).unwrap();

    assert_eq!(result.min.unwrap(), 1.0);
    assert_eq!(result.max.unwrap(), 9.0);
    assert_eq!(result.sum.unwrap(), 20.0);
    assert_eq!(result.count.unwrap(), 5);
}

#[test]
fn test_statistics_variance_std() {
    let data = vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];

    let result = statistics(StatisticsRequest {
        data,
        operations: Some(vec!["variance".to_string(), "std".to_string()]),
    }).unwrap();

    assert!(result.variance.unwrap() > 0.0);
    assert!(result.std.unwrap() > 0.0);
    assert!((result.std.unwrap().powi(2) - result.variance.unwrap()).abs() < 1e-10);
}

#[test]
fn test_statistics_selective_operations() {
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let result = statistics(StatisticsRequest {
        data,
        operations: Some(vec!["mean".to_string(), "std".to_string()]),
    }).unwrap();

    assert!(result.mean.is_some());
    assert!(result.std.is_some());
    assert!(result.median.is_none());
    assert!(result.quartiles.is_none());
}

#[test]
fn test_monte_carlo() {
    let result = monte_carlo_integration(MonteCarloRequest {
        function: "x^2".to_string(),
        bounds: vec![(0.0, 1.0)],
        samples: 10000,
    }).unwrap();

    // Integral of x^2 from 0 to 1 is 1/3 ≈ 0.333
    assert!((result.integral_estimate - 0.333).abs() < 0.1);
}

#[test]
fn test_monte_carlo_2d() {
    let result = monte_carlo_integration(MonteCarloRequest {
        function: "x^2 + y^2".to_string(),
        bounds: vec![(0.0, 1.0), (0.0, 1.0)],
        samples: 10000,
    }).unwrap();

    assert!(result.integral_estimate > 0.0);
    assert!(result.standard_error > 0.0);
    assert_eq!(result.samples_used, 10000);
}

#[test]
fn test_monte_carlo_3d() {
    let result = monte_carlo_integration(MonteCarloRequest {
        function: "x^2 + y^2 + z^2".to_string(),
        bounds: vec![(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)],
        samples: 10000,
    }).unwrap();

    assert!(result.integral_estimate > 0.0);
    assert!(result.standard_error > 0.0);
}

#[test]
fn test_mcmc() {
    let result = mcmc_sampling(MCMCRequest {
        target_distribution: "normal".to_string(),
        dimensions: 1,
        samples: 1000,
        burn_in: 100,
        initial_state: vec![0.0],
    }).unwrap();

    assert_eq!(result.samples.len(), 1000);
    assert!(result.acceptance_rate > 0.0);
}

#[test]
fn test_mcmc_multidimensional() {
    let result = mcmc_sampling(MCMCRequest {
        target_distribution: "normal".to_string(),
        dimensions: 3,
        samples: 500,
        burn_in: 50,
        initial_state: vec![0.0, 0.0, 0.0],
    }).unwrap();

    assert_eq!(result.samples.len(), 500);
    assert!(result.acceptance_rate > 0.0 && result.acceptance_rate < 1.0);
    assert_eq!(result.effective_samples, 500);
}

#[test]
fn test_mcmc_different_initial_states() {
    let result1 = mcmc_sampling(MCMCRequest {
        target_distribution: "normal".to_string(),
        dimensions: 1,
        samples: 100,
        burn_in: 20,
        initial_state: vec![0.0],
    }).unwrap();

    let result2 = mcmc_sampling(MCMCRequest {
        target_distribution: "normal".to_string(),
        dimensions: 1,
        samples: 100,
        burn_in: 20,
        initial_state: vec![10.0],
    }).unwrap();

    assert_eq!(result1.samples.len(), 100);
    assert_eq!(result2.samples.len(), 100);
}

#[test]
fn test_correlation() {
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let result = correlation(CorrelationRequest {
        x,
        y,
        method: "pearson".to_string(),
    }).unwrap();

    assert!((result.correlation - 1.0).abs() < 1e-10);
}

#[test]
fn test_correlation_negative() {
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y = vec![5.0, 4.0, 3.0, 2.0, 1.0];

    let result = correlation(CorrelationRequest {
        x,
        y,
        method: "pearson".to_string(),
    }).unwrap();

    assert!((result.correlation + 1.0).abs() < 1e-10);
}

#[test]
fn test_correlation_uncorrelated() {
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y = vec![2.0, 3.0, 1.0, 5.0, 4.0];

    let result = correlation(CorrelationRequest {
        x,
        y,
        method: "pearson".to_string(),
    }).unwrap();

    // Correlation should be weak (not perfect)
    assert!(result.correlation.abs() < 1.0);
}

#[test]
fn test_correlation_spearman() {
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y = vec![1.0, 4.0, 9.0, 16.0, 25.0]; // y = x^2

    let result = correlation(CorrelationRequest {
        x,
        y,
        method: "spearman".to_string(),
    }).unwrap();

    assert!((result.correlation - 1.0).abs() < 1e-10);
}

#[test]
fn test_kl_divergence() {
    let p = vec![0.5, 0.5];
    let q = vec![0.4, 0.6];

    let result = kl_divergence(KLDivergenceRequest { p, q }).unwrap();

    assert!(result.kl_divergence >= 0.0);
}

#[test]
fn test_kl_divergence_identical() {
    let p = vec![0.25, 0.25, 0.25, 0.25];
    let q = vec![0.25, 0.25, 0.25, 0.25];

    let result = kl_divergence(KLDivergenceRequest { p, q }).unwrap();

    assert!(result.kl_divergence.abs() < 1e-10);
    assert!(result.js_divergence.abs() < 1e-10);
}

#[test]
fn test_kl_divergence_js_symmetric() {
    let p = vec![0.7, 0.3];
    let q = vec![0.3, 0.7];

    let result1 = kl_divergence(KLDivergenceRequest { p: p.clone(), q: q.clone() }).unwrap();
    let result2 = kl_divergence(KLDivergenceRequest { p: q, q: p }).unwrap();

    // Both should compute valid divergences
    assert!(result1.kl_divergence >= 0.0);
    assert!(result2.kl_divergence >= 0.0);
    // JS divergence is symmetric (with small tolerance for floating point)
    assert!((result1.js_divergence - result2.js_divergence).abs() < 1e-6);
}

#[test]
fn test_mutual_information() {
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let result = mutual_information(MutualInformationRequest {
        x,
        y,
        bins: 5,
    }).unwrap();

    assert!(result.mutual_information >= 0.0);
}

#[test]
fn test_mutual_information_independent() {
    // Create independent variables using different patterns
    let x = vec![1.0, 1.0, 2.0, 2.0, 3.0, 3.0];
    let y = vec![1.0, 2.0, 1.0, 2.0, 1.0, 2.0];

    let result = mutual_information(MutualInformationRequest {
        x,
        y,
        bins: 3,
    }).unwrap();

    // Independent variables should have low MI
    assert!(result.mutual_information < 1.0);
}

#[test]
fn test_mutual_information_bins() {
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let y = vec![2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0];

    let result_small = mutual_information(MutualInformationRequest {
        x: x.clone(),
        y: y.clone(),
        bins: 2,
    }).unwrap();

    let result_large = mutual_information(MutualInformationRequest {
        x,
        y,
        bins: 4,
    }).unwrap();

    assert!(result_small.mutual_information > 0.0);
    assert!(result_large.mutual_information > 0.0);
}
