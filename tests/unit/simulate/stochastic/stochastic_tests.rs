//! Unit tests for specialized::stochastic_processes

use crate::simulate::stochastic::*;

#[test]
fn test_brownian_motion_basic() {
    let params = BrownianMotionParams {
        time_steps: 100,
        dt: 0.01,
        initial_value: 0.0,
        drift: 0.0,
        volatility: 1.0,
    };

    let result = generate_brownian_motion(params).unwrap();
    // Returns time_steps + 1 elements (initial value + 100 steps)
    assert_eq!(result.len(), 101);
    assert_eq!(result[0], 0.0); // Initial value
}

#[test]
fn test_markov_chain_basic() {
    let params = MarkovChainParams {
        states: vec!["A".to_string(), "B".to_string()],
        transition_matrix: vec![
            vec![0.7, 0.3],
            vec![0.4, 0.6],
        ],
        initial_state: "A".to_string(),
        steps: 10,
    };

    let result = simulate_markov_chain(params).unwrap();
    // Returns steps + 1 elements (initial state + 10 steps)
    assert_eq!(result.len(), 11);
    assert_eq!(result[0], "A");
}

#[test]
fn test_poisson_process_basic() {
    let result = simulate_poisson_process(5.0, 10.0).unwrap();

    // Should have some events
    assert!(!result.is_empty());

    // All event times should be within time horizon
    assert!(result.iter().all(|&t| t >= 0.0 && t <= 10.0));
}
