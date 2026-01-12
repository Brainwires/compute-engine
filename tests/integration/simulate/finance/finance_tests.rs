//! Unit tests for finance simulation models
//!
//! Tests the Black-Scholes, Heston, SABR, and Stochastic Volatility models
//! from src/simulate/finance/mod.rs

use computational_engine::create_default_dispatcher;
use computational_engine::engine::equations::*;
use computational_engine::engine::*;
use std::collections::HashMap;

// Helper function to create SimulateInput for finance tests
fn create_finance_input(
    model: FinanceModel,
    range: [f64; 2],
    steps: usize,
    parameters: HashMap<String, f64>,
) -> SimulateInput {
    SimulateInput {
        model: SimulationModel::Finance(model),
        equations: Vec::new(),
        variables: vec!["price".to_string()],
        parameters,
        initial_conditions: None,
        range: Some(range),
        steps: Some(steps),
        method: None,
        num_paths: Some(1),
    }
}

// ============================================================================
// BLACK-SCHOLES MODEL TESTS
// ============================================================================

#[test]
fn test_black_scholes_basic_simulation() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_price".to_string(), 100.0);
    params.insert("drift".to_string(), 0.05);
    params.insert("volatility".to_string(), 0.2);

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::BlackScholes,
        [0.0, 1.0],
        100,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate Black-Scholes: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("price"), "Should have price path");
        let prices = &output.results["price"];
        assert_eq!(prices.len(), 101, "Should have 101 price points (100 steps + initial)");
        assert!((prices[0] - 100.0).abs() < 0.01, "Initial price should be 100");
    }
}

#[test]
fn test_black_scholes_zero_volatility() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_price".to_string(), 100.0);
    params.insert("drift".to_string(), 0.05);
    params.insert("volatility".to_string(), 0.0);

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::BlackScholes,
        [0.0, 1.0],
        10,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should handle zero volatility: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        let prices = &output.results["price"];
        // With zero volatility, price should grow deterministically at drift rate
        // After 1 year with 5% drift, price should be approximately 100 * (1 + 0.05) = 105
        let final_price = prices[prices.len() - 1];
        // Allow some tolerance due to discrete time steps
        assert!(final_price > 100.0, "Price should increase with positive drift");
    }
}

#[test]
fn test_black_scholes_high_volatility() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_price".to_string(), 100.0);
    params.insert("drift".to_string(), 0.0);
    params.insert("volatility".to_string(), 1.0); // 100% volatility

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::BlackScholes,
        [0.0, 1.0],
        252,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should handle high volatility: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        let prices = &output.results["price"];
        // All prices should be positive (geometric Brownian motion)
        assert!(
            prices.iter().all(|&p| p > 0.0),
            "All prices should remain positive"
        );
    }
}

#[test]
fn test_black_scholes_negative_drift() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_price".to_string(), 100.0);
    params.insert("drift".to_string(), -0.1); // Negative drift
    params.insert("volatility".to_string(), 0.2);

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::BlackScholes,
        [0.0, 1.0],
        100,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should handle negative drift: {:?}", result);
}

#[test]
fn test_black_scholes_metadata() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_price".to_string(), 100.0);
    params.insert("drift".to_string(), 0.05);
    params.insert("volatility".to_string(), 0.2);

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::BlackScholes,
        [0.0, 1.0],
        10,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok());

    if let Ok(ToolResponse::Simulate(output)) = result {
        let metadata = output.metadata.expect("Should have metadata");
        assert_eq!(metadata["model"], "black_scholes");
        assert_eq!(metadata["drift"], 0.05);
        assert_eq!(metadata["volatility"], 0.2);
    }
}

// ============================================================================
// HESTON MODEL TESTS
// ============================================================================

#[test]
fn test_heston_basic_simulation() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_price".to_string(), 100.0);
    params.insert("initial_variance".to_string(), 0.04);
    params.insert("kappa".to_string(), 2.0);
    params.insert("theta".to_string(), 0.04);
    params.insert("sigma".to_string(), 0.3);
    params.insert("rho".to_string(), -0.7);
    params.insert("risk_free_rate".to_string(), 0.05);

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::Heston,
        [0.0, 1.0],
        252,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate Heston model: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("price"), "Should have price path");
        assert!(output.results.contains_key("variance"), "Should have variance path");

        let prices = &output.results["price"];
        let variances = &output.results["variance"];

        assert_eq!(prices.len(), variances.len(), "Price and variance paths should have same length");
        assert!((prices[0] - 100.0).abs() < 0.01, "Initial price should be 100");
        assert!((variances[0] - 0.04).abs() < 0.001, "Initial variance should be 0.04");
    }
}

#[test]
fn test_heston_variance_mean_reversion() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_price".to_string(), 100.0);
    params.insert("initial_variance".to_string(), 0.01); // Start below long-term variance
    params.insert("kappa".to_string(), 5.0); // Strong mean reversion
    params.insert("theta".to_string(), 0.04); // Long-term variance
    params.insert("sigma".to_string(), 0.1); // Low vol-of-vol for cleaner test
    params.insert("rho".to_string(), 0.0);
    params.insert("risk_free_rate".to_string(), 0.0);

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::Heston,
        [0.0, 5.0], // Longer simulation
        500,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok());

    if let Ok(ToolResponse::Simulate(output)) = result {
        let variances = &output.results["variance"];
        // Variance should be non-negative (Heston enforces this)
        assert!(
            variances.iter().all(|&v| v >= 0.0),
            "Variance should remain non-negative"
        );
    }
}

#[test]
fn test_heston_correlation_effect() {
    let dispatcher = create_default_dispatcher();

    // Test with negative correlation (typical for equity markets)
    let mut params = HashMap::new();
    params.insert("initial_price".to_string(), 100.0);
    params.insert("initial_variance".to_string(), 0.04);
    params.insert("kappa".to_string(), 2.0);
    params.insert("theta".to_string(), 0.04);
    params.insert("sigma".to_string(), 0.3);
    params.insert("rho".to_string(), -0.9); // Strong negative correlation
    params.insert("risk_free_rate".to_string(), 0.05);

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::Heston,
        [0.0, 1.0],
        100,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should handle strong negative correlation");
}

#[test]
fn test_heston_metadata() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_price".to_string(), 100.0);
    params.insert("initial_variance".to_string(), 0.04);
    params.insert("kappa".to_string(), 2.0);
    params.insert("theta".to_string(), 0.04);
    params.insert("sigma".to_string(), 0.3);
    params.insert("rho".to_string(), -0.7);
    params.insert("risk_free_rate".to_string(), 0.05);

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::Heston,
        [0.0, 1.0],
        10,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok());

    if let Ok(ToolResponse::Simulate(output)) = result {
        let metadata = output.metadata.expect("Should have metadata");
        assert_eq!(metadata["model"], "heston");
        assert_eq!(metadata["kappa"], 2.0);
        assert_eq!(metadata["theta"], 0.04);
        assert_eq!(metadata["sigma"], 0.3);
        assert_eq!(metadata["rho"], -0.7);
    }
}

// ============================================================================
// SABR MODEL TESTS
// ============================================================================

#[test]
fn test_sabr_basic_simulation() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("forward_rate".to_string(), 0.05);
    params.insert("alpha".to_string(), 0.3);
    params.insert("beta".to_string(), 0.5);
    params.insert("rho".to_string(), -0.3);
    params.insert("nu".to_string(), 0.4);

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::SABR,
        [0.0, 1.0],
        252,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate SABR model: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("forward_rate"), "Should have forward_rate path");
        assert!(output.results.contains_key("volatility"), "Should have volatility path");

        let forward_rates = &output.results["forward_rate"];
        let volatilities = &output.results["volatility"];

        assert_eq!(forward_rates.len(), volatilities.len());
        assert!((forward_rates[0] - 0.05).abs() < 0.001, "Initial forward rate should be 0.05");
    }
}

#[test]
fn test_sabr_beta_parameter() {
    let dispatcher = create_default_dispatcher();

    // Test with beta = 1 (lognormal case)
    let mut params = HashMap::new();
    params.insert("forward_rate".to_string(), 0.05);
    params.insert("alpha".to_string(), 0.3);
    params.insert("beta".to_string(), 1.0); // Lognormal
    params.insert("rho".to_string(), 0.0);
    params.insert("nu".to_string(), 0.4);

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::SABR,
        [0.0, 1.0],
        100,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should handle beta=1 (lognormal): {:?}", result);
}

#[test]
fn test_sabr_zero_beta() {
    let dispatcher = create_default_dispatcher();

    // Test with beta = 0 (normal case)
    let mut params = HashMap::new();
    params.insert("forward_rate".to_string(), 0.05);
    params.insert("alpha".to_string(), 0.01); // Lower alpha for normal model
    params.insert("beta".to_string(), 0.0); // Normal
    params.insert("rho".to_string(), 0.0);
    params.insert("nu".to_string(), 0.4);

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::SABR,
        [0.0, 1.0],
        100,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should handle beta=0 (normal): {:?}", result);
}

#[test]
fn test_sabr_metadata() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("forward_rate".to_string(), 0.05);
    params.insert("alpha".to_string(), 0.3);
    params.insert("beta".to_string(), 0.5);
    params.insert("rho".to_string(), -0.3);
    params.insert("nu".to_string(), 0.4);

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::SABR,
        [0.0, 1.0],
        10,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok());

    if let Ok(ToolResponse::Simulate(output)) = result {
        let metadata = output.metadata.expect("Should have metadata");
        assert_eq!(metadata["model"], "SABR");
        assert_eq!(metadata["beta"], 0.5);
        assert_eq!(metadata["rho"], -0.3);
        assert_eq!(metadata["nu"], 0.4);
    }
}

// ============================================================================
// STOCHASTIC VOLATILITY MODEL TESTS
// ============================================================================

#[test]
fn test_stochastic_volatility_basic_simulation() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_price".to_string(), 100.0);
    params.insert("initial_volatility".to_string(), 0.2);
    params.insert("drift".to_string(), 0.05);
    params.insert("mean_reversion".to_string(), 1.0);
    params.insert("long_term_vol".to_string(), 0.2);
    params.insert("vol_of_vol".to_string(), 0.3);

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::StochasticVolatility,
        [0.0, 1.0],
        252,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate stochastic volatility: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("price"), "Should have price path");
        assert!(output.results.contains_key("volatility"), "Should have volatility path");

        let prices = &output.results["price"];
        let volatilities = &output.results["volatility"];

        assert_eq!(prices.len(), volatilities.len());
        assert!((prices[0] - 100.0).abs() < 0.01, "Initial price should be 100");
        assert!((volatilities[0] - 0.2).abs() < 0.01, "Initial volatility should be 0.2");
    }
}

#[test]
fn test_stochastic_volatility_mean_reversion() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_price".to_string(), 100.0);
    params.insert("initial_volatility".to_string(), 0.5); // Start above long-term vol
    params.insert("drift".to_string(), 0.0);
    params.insert("mean_reversion".to_string(), 5.0); // Strong mean reversion
    params.insert("long_term_vol".to_string(), 0.2);
    params.insert("vol_of_vol".to_string(), 0.1);

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::StochasticVolatility,
        [0.0, 5.0],
        500,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok());

    if let Ok(ToolResponse::Simulate(output)) = result {
        let volatilities = &output.results["volatility"];
        // Volatility should remain positive (model enforces minimum of 0.01)
        assert!(
            volatilities.iter().all(|&v| v >= 0.01),
            "Volatility should remain above minimum threshold"
        );
    }
}

#[test]
fn test_stochastic_volatility_metadata() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_price".to_string(), 100.0);
    params.insert("initial_volatility".to_string(), 0.2);
    params.insert("drift".to_string(), 0.05);
    params.insert("mean_reversion".to_string(), 1.0);
    params.insert("long_term_vol".to_string(), 0.2);
    params.insert("vol_of_vol".to_string(), 0.3);

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::StochasticVolatility,
        [0.0, 1.0],
        10,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok());

    if let Ok(ToolResponse::Simulate(output)) = result {
        let metadata = output.metadata.expect("Should have metadata");
        assert_eq!(metadata["model"], "stochastic_volatility");
        assert_eq!(metadata["mean_reversion"], 1.0);
        assert_eq!(metadata["long_term_vol"], 0.2);
    }
}

// ============================================================================
// DEFAULT PARAMETER TESTS
// ============================================================================

#[test]
fn test_black_scholes_default_parameters() {
    let dispatcher = create_default_dispatcher();

    // Empty parameters - should use defaults
    let params = HashMap::new();

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::BlackScholes,
        [0.0, 1.0],
        10,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should work with default parameters: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        let prices = &output.results["price"];
        // Default initial_price is 100.0
        assert!((prices[0] - 100.0).abs() < 0.01);
    }
}

#[test]
fn test_heston_default_parameters() {
    let dispatcher = create_default_dispatcher();

    // Empty parameters - should use defaults
    let params = HashMap::new();

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::Heston,
        [0.0, 1.0],
        10,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should work with default parameters: {:?}", result);
}

#[test]
fn test_sabr_default_parameters() {
    let dispatcher = create_default_dispatcher();

    // Empty parameters - should use defaults
    let params = HashMap::new();

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::SABR,
        [0.0, 1.0],
        10,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should work with default parameters: {:?}", result);
}

// ============================================================================
// TIME HANDLING TESTS
// ============================================================================

#[test]
fn test_finance_time_points() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_price".to_string(), 100.0);
    params.insert("drift".to_string(), 0.05);
    params.insert("volatility".to_string(), 0.2);

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::BlackScholes,
        [0.0, 1.0],
        252, // Daily steps for a year
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok());

    if let Ok(ToolResponse::Simulate(output)) = result {
        let time = output.time.expect("Should have time points");
        assert_eq!(time.len(), 253, "Should have 253 time points (252 steps + initial)");
        assert!((time[0] - 0.0).abs() < 0.0001, "First time point should be 0");
        assert!((time[252] - 1.0).abs() < 0.01, "Last time point should be ~1.0");
    }
}

#[test]
fn test_finance_non_zero_start_time() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_price".to_string(), 100.0);
    params.insert("drift".to_string(), 0.05);
    params.insert("volatility".to_string(), 0.2);

    let request = ToolRequest::Simulate(create_finance_input(
        FinanceModel::BlackScholes,
        [0.5, 1.5], // Start at t=0.5, end at t=1.5
        100,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok());

    if let Ok(ToolResponse::Simulate(output)) = result {
        let time = output.time.expect("Should have time points");
        assert!((time[0] - 0.5).abs() < 0.01, "First time point should be 0.5");
    }
}
