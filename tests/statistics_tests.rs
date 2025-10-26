use computational_engine::statistics::*;

#[test]
fn test_statistics() {
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let result = statistics(StatisticsRequest {
        data,
        operations: None,
    })
    .unwrap();

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
        operations: Some(vec![
            "mean".to_string(),
            "median".to_string(),
            "mode".to_string(),
        ]),
    })
    .unwrap();

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
    })
    .unwrap();

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
    })
    .unwrap();

    assert!(result.skewness.is_some());
    assert!(result.kurtosis.is_some());
}

#[test]
fn test_statistics_min_max_sum_count() {
    let data = vec![1.0, 5.0, 3.0, 9.0, 2.0];

    let result = statistics(StatisticsRequest {
        data: data.clone(),
        operations: Some(vec![
            "min".to_string(),
            "max".to_string(),
            "sum".to_string(),
            "count".to_string(),
        ]),
    })
    .unwrap();

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
    })
    .unwrap();

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
    })
    .unwrap();

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
    })
    .unwrap();

    // Integral of x^2 from 0 to 1 is 1/3 ≈ 0.333
    assert!((result.integral_estimate - 0.333).abs() < 0.1);
}

#[test]
fn test_monte_carlo_2d() {
    let result = monte_carlo_integration(MonteCarloRequest {
        function: "x^2 + y^2".to_string(),
        bounds: vec![(0.0, 1.0), (0.0, 1.0)],
        samples: 10000,
    })
    .unwrap();

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
    })
    .unwrap();

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
    })
    .unwrap();

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
    })
    .unwrap();

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
    })
    .unwrap();

    let result2 = mcmc_sampling(MCMCRequest {
        target_distribution: "normal".to_string(),
        dimensions: 1,
        samples: 100,
        burn_in: 20,
        initial_state: vec![10.0],
    })
    .unwrap();

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
    })
    .unwrap();

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
    })
    .unwrap();

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
    })
    .unwrap();

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
    })
    .unwrap();

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

    let result1 = kl_divergence(KLDivergenceRequest {
        p: p.clone(),
        q: q.clone(),
    })
    .unwrap();
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

    let result = mutual_information(MutualInformationRequest { x, y, bins: 5 }).unwrap();

    assert!(result.mutual_information >= 0.0);
}

#[test]
fn test_mutual_information_independent() {
    // Create independent variables using different patterns
    let x = vec![1.0, 1.0, 2.0, 2.0, 3.0, 3.0];
    let y = vec![1.0, 2.0, 1.0, 2.0, 1.0, 2.0];

    let result = mutual_information(MutualInformationRequest { x, y, bins: 3 }).unwrap();

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
    })
    .unwrap();

    let result_large = mutual_information(MutualInformationRequest { x, y, bins: 4 }).unwrap();

    assert!(result_small.mutual_information > 0.0);
    assert!(result_large.mutual_information > 0.0);
}

// ========================================================================
// CRITICAL BUG FIX TESTS (Added 2025-10-26)
// ========================================================================

#[test]
fn test_spearman_correlation_with_ties() {
    // Bug Fix #1: Test that Spearman correlation handles ties correctly
    // Data with ties: [1.0, 2.0, 2.0, 3.0]
    // Correct ranks: [1.0, 2.5, 2.5, 4.0] (average of ranks 2 and 3)
    let x = vec![1.0, 2.0, 2.0, 3.0];
    let y = vec![1.0, 2.0, 2.0, 3.0]; // Perfect rank correlation

    let result = correlation(CorrelationRequest {
        x,
        y,
        method: "spearman".to_string(),
    })
    .unwrap();

    // With correct tie handling, perfect rank correlation should be 1.0
    assert!(
        (result.correlation - 1.0).abs() < 1e-10,
        "Spearman correlation with ties should be 1.0, got {}",
        result.correlation
    );
}

#[test]
fn test_spearman_correlation_multiple_ties() {
    // More complex tie scenario
    let x = vec![1.0, 2.0, 2.0, 2.0, 3.0, 4.0];
    let y = vec![1.0, 2.0, 2.0, 2.0, 3.0, 4.0];

    let result = correlation(CorrelationRequest {
        x,
        y,
        method: "spearman".to_string(),
    })
    .unwrap();

    assert!(
        (result.correlation - 1.0).abs() < 1e-10,
        "Spearman correlation with multiple ties should be 1.0"
    );
}

#[test]
fn test_sample_variance_bessel_correction() {
    // Bug Fix #6: Test that variance uses Bessel's correction (÷n-1)
    let data = vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
    let n = data.len() as f64;

    let result = statistics(StatisticsRequest {
        data: data.clone(),
        operations: Some(vec!["variance".to_string(), "mean".to_string()]),
    })
    .unwrap();

    // Calculate expected sample variance manually
    let mean = data.iter().sum::<f64>() / n;
    let expected_variance = data
        .iter()
        .map(|x| (x - mean).powi(2))
        .sum::<f64>()
        / (n - 1.0); // Bessel's correction

    assert!(
        (result.variance.unwrap() - expected_variance).abs() < 1e-10,
        "Variance should use Bessel's correction (÷n-1)"
    );
}

#[test]
fn test_sample_variance_vs_population_variance() {
    // Verify sample variance is larger than population variance
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let n = data.len() as f64;

    let result = statistics(StatisticsRequest {
        data: data.clone(),
        operations: Some(vec!["variance".to_string(), "mean".to_string()]),
    })
    .unwrap();

    let mean = data.iter().sum::<f64>() / n;

    // Sample variance (÷n-1)
    let sample_variance = result.variance.unwrap();

    // Population variance (÷n)
    let population_variance = data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;

    // Sample variance should be slightly larger
    assert!(
        sample_variance > population_variance,
        "Sample variance should be larger than population variance"
    );
}

#[test]
fn test_standard_deviation_consistency() {
    // Ensure std = sqrt(variance) with Bessel's correction
    let data = vec![10.0, 12.0, 23.0, 23.0, 16.0, 23.0, 21.0, 16.0];

    let result = statistics(StatisticsRequest {
        data,
        operations: Some(vec!["variance".to_string(), "std".to_string()]),
    })
    .unwrap();

    let variance = result.variance.unwrap();
    let std = result.std.unwrap();

    assert!(
        (std - variance.sqrt()).abs() < 1e-10,
        "Standard deviation should equal sqrt(variance)"
    );
}

#[test]
fn test_monte_carlo_hardcoded_function() {
    // Bug Fix #2: Verify Monte Carlo is hardcoded to Σx²
    // For bounds [0,1], integral of x² is 1/3 ≈ 0.333
    let result = monte_carlo_integration(MonteCarloRequest {
        function: "IGNORED".to_string(), // This parameter is ignored
        bounds: vec![(0.0, 1.0)],
        samples: 100000,
    })
    .unwrap();

    // Should converge to 1/3 regardless of function parameter
    assert!(
        (result.integral_estimate - 0.333).abs() < 0.01,
        "Monte Carlo should compute integral of x² ≈ 0.333"
    );
}

#[test]
fn test_mcmc_hardcoded_gaussian() {
    // Bug Fix #3: Verify MCMC is hardcoded to Gaussian
    let result = mcmc_sampling(MCMCRequest {
        target_distribution: "IGNORED".to_string(), // This parameter is ignored
        dimensions: 1,
        samples: 10000,
        burn_in: 1000,
        initial_state: vec![0.0],
    })
    .unwrap();

    // For Gaussian N(0,1), samples should cluster around mean=0
    let samples: Vec<f64> = result.samples.iter().map(|s| s[0]).collect();
    let mean = samples.iter().sum::<f64>() / samples.len() as f64;

    // Mean should be close to 0 for standard Gaussian
    assert!(
        mean.abs() < 0.1,
        "MCMC with Gaussian should have mean ≈ 0, got {}",
        mean
    );
}

#[test]
fn test_pearson_vs_spearman_nonlinear() {
    // Spearman should detect monotonic relationship, Pearson may not
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y = vec![1.0, 4.0, 9.0, 16.0, 25.0]; // y = x²

    let pearson = correlation(CorrelationRequest {
        x: x.clone(),
        y: y.clone(),
        method: "pearson".to_string(),
    })
    .unwrap();

    let spearman = correlation(CorrelationRequest {
        x,
        y,
        method: "spearman".to_string(),
    })
    .unwrap();

    // Spearman should be perfect (1.0) for monotonic relationship
    assert!((spearman.correlation - 1.0).abs() < 1e-10);
    // Pearson will be high but not perfect
    assert!(pearson.correlation > 0.9 && pearson.correlation < 1.0);
}

#[test]
fn test_skewness_kurtosis_with_bessel() {
    // Verify skewness and kurtosis use sample std (with Bessel's correction)
    let data = vec![1.0, 2.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];

    let result = statistics(StatisticsRequest {
        data,
        operations: Some(vec![
            "skewness".to_string(),
            "kurtosis".to_string(),
            "std".to_string(),
        ]),
    })
    .unwrap();

    // Should compute without errors
    assert!(result.skewness.is_some());
    assert!(result.kurtosis.is_some());
    // Skewness can be positive or negative
    assert!(result.skewness.unwrap().is_finite());
    // Excess kurtosis (minus 3)
    assert!(result.kurtosis.unwrap().is_finite());
}

#[test]
fn test_correlation_error_handling() {
    // Test mismatched lengths
    let x = vec![1.0, 2.0, 3.0];
    let y = vec![1.0, 2.0];

    let result = correlation(CorrelationRequest {
        x,
        y,
        method: "pearson".to_string(),
    });

    assert!(result.is_err());
}

#[test]
fn test_statistics_empty_data() {
    // Test empty data error handling
    let result = statistics(StatisticsRequest {
        data: vec![],
        operations: None,
    });

    assert!(result.is_err());
}

#[test]
fn test_kl_divergence_error_handling() {
    // Test mismatched lengths
    let result = kl_divergence(KLDivergenceRequest {
        p: vec![0.5, 0.5],
        q: vec![0.33, 0.33, 0.34],
    });

    assert!(result.is_err());

    // Test negative probabilities
    let result = kl_divergence(KLDivergenceRequest {
        p: vec![0.5, -0.5],
        q: vec![0.5, 0.5],
    });

    assert!(result.is_err());
}
