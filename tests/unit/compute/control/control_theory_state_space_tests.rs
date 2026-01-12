//! Unit tests for state_space module

use super::{StateSpace, siso_state_space};

#[test]
fn test_state_space_creation() {
    // Simple integrator: ẋ = u, y = x
    let sys = StateSpace::new(
        vec![vec![0.0]],  // A
        vec![vec![1.0]],  // B
        vec![vec![1.0]],  // C
        vec![vec![0.0]],  // D
    );
    assert!(sys.is_ok());
    let sys = sys.unwrap();
    assert_eq!(sys.num_states(), 1);
    assert_eq!(sys.num_inputs(), 1);
    assert_eq!(sys.num_outputs(), 1);
}

#[test]
fn test_state_space_output() {
    let mut sys = StateSpace::new(
        vec![vec![0.0]],
        vec![vec![1.0]],
        vec![vec![1.0]],
        vec![vec![0.0]],
    ).unwrap();

    // Set state and compute output
    sys.set_state(vec![2.0]).unwrap();
    let output = sys.output(&[0.0]).unwrap();
    assert_eq!(output.len(), 1);
    assert!((output[0] - 2.0).abs() < 1e-10);
}

#[test]
fn test_state_space_step() {
    let mut sys = StateSpace::new(
        vec![vec![-1.0]],  // A
        vec![vec![1.0]],   // B
        vec![vec![1.0]],   // C
        vec![vec![0.0]],   // D
    ).unwrap();

    let output = sys.step(&[1.0], 0.1);
    assert!(output.is_ok());
}

#[test]
fn test_state_space_simulate() {
    let mut sys = StateSpace::new(
        vec![vec![0.0]],
        vec![vec![1.0]],
        vec![vec![1.0]],
        vec![vec![0.0]],
    ).unwrap();

    let inputs = vec![vec![1.0], vec![1.0], vec![1.0]];
    let outputs = sys.simulate(&inputs, 0.1);
    assert!(outputs.is_ok());
    assert_eq!(outputs.unwrap().len(), 3);
}

#[test]
fn test_siso_state_space() {
    let sys = siso_state_space(
        vec![vec![-1.0]],
        vec![1.0],
        vec![1.0],
        0.0,
    );
    assert!(sys.is_ok());
    let sys = sys.unwrap();
    assert_eq!(sys.num_states(), 1);
    assert_eq!(sys.num_inputs(), 1);
    assert_eq!(sys.num_outputs(), 1);
}
