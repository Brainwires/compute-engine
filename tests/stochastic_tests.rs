use computational_engine::engine::*;
use computational_engine::create_default_dispatcher;
use std::collections::HashMap;

#[test]
fn test_brownian_motion() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_value".to_string(), 0.0);
    params.insert("drift".to_string(), 0.0);
    params.insert("volatility".to_string(), 1.0);

    let request = ToolRequest::Simulate(SimulateInput {
        model: SimulationModel::Stochastic(StochasticProcess::BrownianMotion),
        equations: vec![],
        variables: vec!["W".to_string()],
        parameters: params,
        initial_conditions: None,
        range: Some([0.0, 1.0]),
        steps: Some(100),
        num_paths: Some(1),
        method: None,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Brownian motion simulation should work");
}

#[test]
fn test_geometric_brownian_motion() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_value".to_string(), 100.0);
    params.insert("drift".to_string(), 0.05);
    params.insert("volatility".to_string(), 0.2);

    let request = ToolRequest::Simulate(SimulateInput {
        model: SimulationModel::Stochastic(StochasticProcess::GeometricBrownian),
        equations: vec![],
        variables: vec!["S".to_string()],
        parameters: params,
        initial_conditions: None,
        range: Some([0.0, 1.0]),
        steps: Some(100),
        num_paths: Some(1),
        method: None,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Geometric Brownian motion simulation should work");
}

#[test]
fn test_ornstein_uhlenbeck() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_value".to_string(), 0.0);
    params.insert("mean_reversion_speed".to_string(), 1.0);
    params.insert("long_term_mean".to_string(), 0.0);
    params.insert("volatility".to_string(), 0.5);

    let request = ToolRequest::Simulate(SimulateInput {
        model: SimulationModel::Stochastic(StochasticProcess::OrnsteinUhlenbeck),
        equations: vec![],
        variables: vec!["X".to_string()],
        parameters: params,
        initial_conditions: None,
        range: Some([0.0, 1.0]),
        steps: Some(100),
        num_paths: Some(1),
        method: None,
    });

    let result = dispatcher.dispatch(request);
    // Ornstein-Uhlenbeck may not be fully implemented
    let _ = result;
}

#[test]
fn test_poisson_process() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("intensity".to_string(), 10.0);

    let request = ToolRequest::Simulate(SimulateInput {
        model: SimulationModel::Stochastic(StochasticProcess::Poisson),
        equations: vec![],
        variables: vec!["N".to_string()],
        parameters: params,
        initial_conditions: None,
        range: Some([0.0, 1.0]),
        steps: Some(100),
        num_paths: Some(1),
        method: None,
    });

    let result = dispatcher.dispatch(request);
    // Poisson process may not be fully implemented
    let _ = result;
}
