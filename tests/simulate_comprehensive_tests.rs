//! Comprehensive simulation test suite
//!
//! Tests for all SIMULATE tool operations including:
//! - Stochastic processes (Brownian, Geometric Brownian, Ornstein-Uhlenbeck, Poisson, Lévy, Jump Diffusion, Fractional Brownian, Variance Gamma)
//! - Finance models (Heston, SABR, Stochastic Volatility, Black-Scholes)
//! - Fluid dynamics (Lattice Boltzmann)
//! - Time evolution (Euler, RK4, Adaptive Step, Implicit Euler)

use computational_engine::engine::*;
use computational_engine::engine::equations::*;
use computational_engine::create_default_dispatcher;
use std::collections::HashMap;

// Helper function to create SimulateInput with required fields
fn create_simulate_input(
    model: SimulationModel,
    variables: Vec<String>,
    initial_conditions: Option<HashMap<String, f64>>,
    range: Option<[f64; 2]>,
    steps: Option<usize>,
    num_paths: Option<usize>,
    parameters: HashMap<String, f64>,
) -> SimulateInput {
    SimulateInput {
        model,
        equations: Vec::new(), // Not used by current implementation
        variables,
        parameters,
        initial_conditions,
        range,
        steps,
        method: None,
        num_paths,
    }
}

// ============================================================================
// STOCHASTIC PROCESS TESTS
// ============================================================================

#[test]
fn test_brownian_motion() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_value".to_string(), 0.0);
    params.insert("drift".to_string(), 0.0);
    params.insert("volatility".to_string(), 1.0);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::Stochastic(StochasticProcess::BrownianMotion),
        vec!["W".to_string()],
        None,
        Some([0.0, 1.0]),
        Some(100),
        Some(1),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate Brownian motion: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("W"), "Should have W path");
        assert!(output.time.is_some(), "Should have time points");
    }
}

#[test]
fn test_geometric_brownian() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_value".to_string(), 100.0);
    params.insert("drift".to_string(), 0.05);
    params.insert("volatility".to_string(), 0.2);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::Stochastic(StochasticProcess::GeometricBrownian),
        vec!["S".to_string()],
        None,
        Some([0.0, 1.0]),
        Some(252),
        Some(1),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate geometric Brownian: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("S"), "Should have S path");
    }
}

#[test]
fn test_ornstein_uhlenbeck() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("theta".to_string(), 0.5);
    params.insert("mu".to_string(), 0.0);
    params.insert("sigma".to_string(), 0.3);
    params.insert("initial_value".to_string(), 0.0);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::Stochastic(StochasticProcess::OrnsteinUhlenbeck),
        vec!["X".to_string()],
        None,
        Some([0.0, 10.0]),
        Some(1000),
        Some(1),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate Ornstein-Uhlenbeck: {:?}", result);
}

#[test]
fn test_poisson_process() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("lambda".to_string(), 2.0);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::Stochastic(StochasticProcess::Poisson),
        vec!["N".to_string()],
        None,
        Some([0.0, 10.0]),
        Some(1000),
        Some(1),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate Poisson process: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("N"), "Should have N path");
    }
}

#[test]
fn test_levy_process() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("alpha".to_string(), 1.5);
    params.insert("beta".to_string(), 0.0);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::Stochastic(StochasticProcess::Levy),
        vec!["L".to_string()],
        None,
        Some([0.0, 1.0]),
        Some(500),
        Some(1),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate Lévy process: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("L"), "Should have L path");
        let metadata = output.metadata.unwrap();
        assert_eq!(metadata["process"], "levy");
    }
}

#[test]
fn test_jump_diffusion() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_value".to_string(), 100.0);
    params.insert("drift".to_string(), 0.05);
    params.insert("volatility".to_string(), 0.2);
    params.insert("jump_intensity".to_string(), 0.1);
    params.insert("jump_mean".to_string(), 0.0);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::Stochastic(StochasticProcess::JumpDiffusion),
        vec!["S".to_string()],
        None,
        Some([0.0, 1.0]),
        Some(252),
        Some(1),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate jump diffusion: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        let metadata = output.metadata.unwrap();
        assert_eq!(metadata["process"], "jump_diffusion");
    }
}

#[test]
fn test_fractional_brownian() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("hurst".to_string(), 0.7);
    params.insert("initial_value".to_string(), 0.0);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::Stochastic(StochasticProcess::FractionalBrownian),
        vec!["B_H".to_string()],
        None,
        Some([0.0, 1.0]),
        Some(100),
        Some(1),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate fractional Brownian: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("B_H"), "Should have B_H path");
        let metadata = output.metadata.unwrap();
        assert_eq!(metadata["process"], "fractional_brownian");
        assert_eq!(metadata["hurst"], 0.7);
    }
}

#[test]
fn test_variance_gamma() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("nu".to_string(), 0.2);
    params.insert("theta".to_string(), 0.0);
    params.insert("sigma".to_string(), 0.3);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::Stochastic(StochasticProcess::VarianceGamma),
        vec!["X_VG".to_string()],
        None,
        Some([0.0, 1.0]),
        Some(500),
        Some(1),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate variance gamma: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("X_VG"), "Should have X_VG path");
        let metadata = output.metadata.unwrap();
        assert_eq!(metadata["process"], "variance_gamma");
    }
}

#[test]
fn test_mean_reverting() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("kappa".to_string(), 0.5);
    params.insert("theta".to_string(), 1.0);
    params.insert("sigma".to_string(), 0.2);
    params.insert("initial_value".to_string(), 0.5);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::Stochastic(StochasticProcess::MeanReverting),
        vec!["X".to_string()],
        None,
        Some([0.0, 10.0]),
        Some(1000),
        Some(1),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate mean-reverting process: {:?}", result);
}

// ============================================================================
// FINANCE MODEL TESTS
// ============================================================================

#[test]
fn test_heston_model() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_price".to_string(), 100.0);
    params.insert("initial_variance".to_string(), 0.04);
    params.insert("kappa".to_string(), 2.0);
    params.insert("theta".to_string(), 0.04);
    params.insert("sigma".to_string(), 0.3);
    params.insert("rho".to_string(), -0.7);
    params.insert("risk_free_rate".to_string(), 0.05);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::Finance(FinanceModel::Heston),
        vec!["price".to_string(), "variance".to_string()],
        None,
        Some([0.0, 1.0]),
        Some(252),
        Some(1),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate Heston model: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("price"), "Should have price path");
        assert!(output.results.contains_key("variance"), "Should have variance path");
        let metadata = output.metadata.unwrap();
        assert_eq!(metadata["model"], "heston");
    }
}

#[test]
fn test_sabr_model() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("forward_rate".to_string(), 0.05);
    params.insert("alpha".to_string(), 0.3);
    params.insert("beta".to_string(), 0.5);
    params.insert("rho".to_string(), -0.3);
    params.insert("nu".to_string(), 0.4);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::Finance(FinanceModel::SABR),
        vec!["forward_rate".to_string(), "volatility".to_string()],
        None,
        Some([0.0, 1.0]),
        Some(252),
        Some(1),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate SABR model: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("forward_rate"), "Should have forward_rate path");
        assert!(output.results.contains_key("volatility"), "Should have volatility path");
        let metadata = output.metadata.unwrap();
        assert_eq!(metadata["model"], "SABR");
    }
}

#[test]
fn test_stochastic_volatility() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_price".to_string(), 100.0);
    params.insert("initial_volatility".to_string(), 0.2);
    params.insert("drift".to_string(), 0.05);
    params.insert("mean_reversion".to_string(), 1.0);
    params.insert("long_term_vol".to_string(), 0.2);
    params.insert("vol_of_vol".to_string(), 0.3);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::Finance(FinanceModel::StochasticVolatility),
        vec!["price".to_string(), "volatility".to_string()],
        None,
        Some([0.0, 1.0]),
        Some(252),
        Some(1),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate stochastic volatility: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("price"), "Should have price path");
        assert!(output.results.contains_key("volatility"), "Should have volatility path");
    }
}

#[test]
fn test_black_scholes() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_price".to_string(), 100.0);
    params.insert("drift".to_string(), 0.05);
    params.insert("volatility".to_string(), 0.2);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::Finance(FinanceModel::BlackScholes),
        vec!["price".to_string()],
        None,
        Some([0.0, 1.0]),
        Some(252),
        Some(1),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate Black-Scholes: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("price"), "Should have price path");
        let metadata = output.metadata.unwrap();
        assert_eq!(metadata["model"], "black_scholes");
    }
}

// ============================================================================
// FLUID DYNAMICS TESTS
// ============================================================================

#[test]
fn test_lattice_boltzmann() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("nx".to_string(), 20.0);
    params.insert("ny".to_string(), 20.0);
    params.insert("tau".to_string(), 0.6);
    params.insert("lid_velocity".to_string(), 0.1);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::FluidDynamics(FluidSim::LatticeBotzmann),
        vec!["velocity_x".to_string(), "velocity_y".to_string()],
        None,
        None,
        Some(100),
        None,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate Lattice Boltzmann: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("velocity_x"), "Should have velocity_x");
        assert!(output.results.contains_key("velocity_y"), "Should have velocity_y");
        assert!(output.results.contains_key("density"), "Should have density");

        let metadata = output.metadata.unwrap();
        assert_eq!(metadata["method"], "lattice_boltzmann_d2q9");

        // Verify we have the right number of grid points
        let vel_x = &output.results["velocity_x"];
        assert_eq!(vel_x.len(), 20 * 20, "Should have 400 grid points (20x20)");
    }
}

// ============================================================================
// TIME EVOLUTION TESTS
// ============================================================================

#[test]
fn test_euler_method() {
    let dispatcher = create_default_dispatcher();

    let mut initial = HashMap::new();
    initial.insert("y".to_string(), 1.0);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::TimeEvolution(TimeEvolutionMethod::Euler),
        vec!["y".to_string()],
        Some(initial),
        Some([0.0, 1.0]),
        Some(100),
        None,
        HashMap::new(),
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate with Euler method: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("y"), "Should have y values");
        assert!(output.time.is_some(), "Should have time points");
    }
}

#[test]
fn test_runge_kutta_4() {
    let dispatcher = create_default_dispatcher();

    let mut initial = HashMap::new();
    initial.insert("y".to_string(), 1.0);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::TimeEvolution(TimeEvolutionMethod::RungeKutta4),
        vec!["y".to_string()],
        Some(initial),
        Some([0.0, 1.0]),
        Some(100),
        None,
        HashMap::new(),
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate with RK4: {:?}", result);
}

#[test]
fn test_adaptive_step() {
    let dispatcher = create_default_dispatcher();

    let mut initial = HashMap::new();
    initial.insert("y".to_string(), 1.0);

    let mut params = HashMap::new();
    params.insert("tolerance".to_string(), 1e-6);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::TimeEvolution(TimeEvolutionMethod::AdaptiveStep),
        vec!["y".to_string()],
        Some(initial),
        Some([0.0, 1.0]),
        Some(100),
        None,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate with adaptive step: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        let metadata = output.metadata.unwrap();
        assert_eq!(metadata["method"], "adaptive_step");
    }
}

#[test]
fn test_implicit_euler() {
    let dispatcher = create_default_dispatcher();

    let mut initial = HashMap::new();
    initial.insert("y".to_string(), 1.0);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::TimeEvolution(TimeEvolutionMethod::ImplicitEuler),
        vec!["y".to_string()],
        Some(initial),
        Some([0.0, 1.0]),
        Some(100),
        None,
        HashMap::new(),
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate with implicit Euler: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        let metadata = output.metadata.unwrap();
        assert_eq!(metadata["method"], "implicit_euler");
    }
}

// ============================================================================
// ERROR HANDLING TESTS
// ============================================================================

#[test]
fn test_ou_requires_parameters() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::Stochastic(StochasticProcess::OrnsteinUhlenbeck),
        vec!["X".to_string()],
        None,
        Some([0.0, 1.0]),
        Some(100),
        Some(1),
        HashMap::new(),
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Ornstein-Uhlenbeck should require theta parameter");
}

#[test]
fn test_time_evolution_requires_initial() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::TimeEvolution(TimeEvolutionMethod::Euler),
        vec!["y".to_string()],
        None,
        Some([0.0, 1.0]),
        Some(100),
        None,
        HashMap::new(),
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Time evolution should require initial conditions");
}

#[test]
fn test_stochastic_requires_range() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("lambda".to_string(), 1.0);

    let request = ToolRequest::Simulate(create_simulate_input(
        SimulationModel::Stochastic(StochasticProcess::Poisson),
        vec!["N".to_string()],
        None,
        None,
        Some(100),
        Some(1),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Stochastic simulation should require range");
}
