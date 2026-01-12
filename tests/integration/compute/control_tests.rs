//! Comprehensive tests for control systems compute operations
//!
//! Tests control theory operations including:
//! - Transfer functions
//! - Pole-zero analysis
//! - Bode plots
//! - Stability analysis
//! - State-space methods
//! - Time-domain analysis

use computational_engine::engine::*;
use serde_json::{json, Value};
use std::collections::HashMap;

// ============================================================================
// TRANSFER FUNCTION TESTS
// ============================================================================

#[test]
fn test_transfer_function_evaluate() {
    let mut parameters = HashMap::new();
    parameters.insert("numerator".to_string(), json!([1.0]));
    parameters.insert("denominator".to_string(), json!([1.0, 2.0, 1.0]));
    parameters.insert("operation".to_string(), json!("evaluate"));
    parameters.insert("frequency".to_string(), json!(1.0));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::TransferFunction)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Transfer function evaluation should succeed");
}

#[test]
fn test_transfer_function_series() {
    let mut parameters = HashMap::new();
    parameters.insert("numerator".to_string(), json!([1.0]));
    parameters.insert("denominator".to_string(), json!([1.0, 1.0]));
    parameters.insert("operation".to_string(), json!("series"));
    parameters.insert(
        "second_tf".to_string(),
        json!([[2.0], [1.0, 2.0]]),
    );

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::TransferFunction)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Series transfer function should succeed");
}

#[test]
fn test_transfer_function_parallel() {
    let mut parameters = HashMap::new();
    parameters.insert("numerator".to_string(), json!([1.0]));
    parameters.insert("denominator".to_string(), json!([1.0, 1.0]));
    parameters.insert("operation".to_string(), json!("parallel"));
    parameters.insert(
        "second_tf".to_string(),
        json!([[1.0], [1.0, 2.0]]),
    );

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::TransferFunction)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Parallel transfer function should succeed");
}

// ============================================================================
// POLE-ZERO ANALYSIS TESTS
// ============================================================================

#[test]
fn test_pole_zero_analysis_first_order() {
    let mut parameters = HashMap::new();
    parameters.insert("numerator".to_string(), json!([1.0]));
    parameters.insert("denominator".to_string(), json!([1.0, 1.0]));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::PoleZeroAnalysis)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Pole-zero analysis should succeed");
}

#[test]
fn test_pole_zero_analysis_second_order() {
    let mut parameters = HashMap::new();
    parameters.insert("numerator".to_string(), json!([1.0, 0.5]));
    parameters.insert("denominator".to_string(), json!([1.0, 2.0, 5.0]));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::PoleZeroAnalysis)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Second-order pole-zero analysis should succeed");
}

// ============================================================================
// BODE PLOT TESTS
// ============================================================================

#[test]
fn test_bode_plot_first_order() {
    let mut parameters = HashMap::new();
    parameters.insert("numerator".to_string(), json!([1.0]));
    parameters.insert("denominator".to_string(), json!([1.0, 1.0]));
    parameters.insert("start_frequency".to_string(), json!(0.1));
    parameters.insert("end_frequency".to_string(), json!(100.0));
    parameters.insert("points".to_string(), json!(50));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::BodePlot)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Bode plot should succeed");
}

#[test]
fn test_bode_plot_second_order() {
    let mut parameters = HashMap::new();
    parameters.insert("numerator".to_string(), json!([1.0]));
    parameters.insert("denominator".to_string(), json!([1.0, 0.1, 1.0]));
    parameters.insert("start_frequency".to_string(), json!(0.01));
    parameters.insert("end_frequency".to_string(), json!(10.0));
    parameters.insert("points".to_string(), json!(100));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::BodePlot)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Second-order Bode plot should succeed");
}

// ============================================================================
// NYQUIST PLOT TESTS
// ============================================================================

#[test]
fn test_nyquist_plot() {
    let mut parameters = HashMap::new();
    parameters.insert("numerator".to_string(), json!([1.0]));
    parameters.insert("denominator".to_string(), json!([1.0, 2.0, 1.0]));
    parameters.insert("start_frequency".to_string(), json!(0.1));
    parameters.insert("end_frequency".to_string(), json!(10.0));
    parameters.insert("points".to_string(), json!(100));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::NyquistPlot)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Nyquist plot should succeed");
}

// ============================================================================
// ROOT LOCUS TESTS
// ============================================================================

#[test]
fn test_root_locus() {
    let mut parameters = HashMap::new();
    parameters.insert("numerator".to_string(), json!([1.0]));
    parameters.insert("denominator".to_string(), json!([1.0, 3.0, 2.0]));
    parameters.insert("gain_start".to_string(), json!(0.0));
    parameters.insert("gain_end".to_string(), json!(10.0));
    parameters.insert("points".to_string(), json!(50));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::RootLocus)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Root locus should succeed");
}

// ============================================================================
// STATE-SPACE TESTS
// ============================================================================

#[test]
fn test_state_space_first_order() {
    let mut parameters = HashMap::new();
    parameters.insert("a_matrix".to_string(), json!([[-1.0]]));
    parameters.insert("b_matrix".to_string(), json!([[1.0]]));
    parameters.insert("c_matrix".to_string(), json!([[1.0]]));
    parameters.insert("d_matrix".to_string(), json!([[0.0]]));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::StateSpace)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "State-space representation should succeed");
}

#[test]
fn test_state_space_second_order() {
    let mut parameters = HashMap::new();
    parameters.insert("a_matrix".to_string(), json!([[0.0, 1.0], [-2.0, -3.0]]));
    parameters.insert("b_matrix".to_string(), json!([[0.0], [1.0]]));
    parameters.insert("c_matrix".to_string(), json!([[1.0, 0.0]]));
    parameters.insert("d_matrix".to_string(), json!([[0.0]]));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::StateSpace)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Second-order state-space should succeed");
}

// ============================================================================
// CONTROLLABILITY TESTS
// ============================================================================

#[test]
fn test_controllability_controllable() {
    let mut parameters = HashMap::new();
    parameters.insert("a_matrix".to_string(), json!([[0.0, 1.0], [-2.0, -3.0]]));
    parameters.insert("b_matrix".to_string(), json!([[0.0], [1.0]]));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::Controllability)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Controllability test should succeed");
}

#[test]
fn test_controllability_uncontrollable() {
    let mut parameters = HashMap::new();
    parameters.insert("a_matrix".to_string(), json!([[1.0, 0.0], [0.0, 2.0]]));
    parameters.insert("b_matrix".to_string(), json!([[1.0], [0.0]]));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::Controllability)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Uncontrollable system test should succeed");
}

// ============================================================================
// OBSERVABILITY TESTS
// ============================================================================

#[test]
fn test_observability_observable() {
    let mut parameters = HashMap::new();
    parameters.insert("a_matrix".to_string(), json!([[0.0, 1.0], [-2.0, -3.0]]));
    parameters.insert("c_matrix".to_string(), json!([[1.0, 0.0]]));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::Observability)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Observability test should succeed");
}

#[test]
fn test_observability_unobservable() {
    let mut parameters = HashMap::new();
    parameters.insert("a_matrix".to_string(), json!([[1.0, 0.0], [0.0, 2.0]]));
    parameters.insert("c_matrix".to_string(), json!([[1.0, 0.0]]));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::Observability)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Unobservable system test should succeed");
}

// ============================================================================
// ROUTH-HURWITZ TESTS
// ============================================================================

#[test]
fn test_routh_hurwitz_stable() {
    let mut parameters = HashMap::new();
    // s^3 + 3s^2 + 3s + 1 (stable)
    parameters.insert("characteristic_polynomial".to_string(), json!([1.0, 3.0, 3.0, 1.0]));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::RouthHurwitz)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Routh-Hurwitz stable system should succeed");
}

#[test]
fn test_routh_hurwitz_unstable() {
    let mut parameters = HashMap::new();
    // s^2 - 1 (unstable)
    parameters.insert("characteristic_polynomial".to_string(), json!([1.0, 0.0, -1.0]));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::RouthHurwitz)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Routh-Hurwitz unstable system should succeed");
}

// ============================================================================
// GAIN MARGIN TESTS
// ============================================================================

#[test]
fn test_gain_margin() {
    let mut parameters = HashMap::new();
    parameters.insert("numerator".to_string(), json!([1.0]));
    parameters.insert("denominator".to_string(), json!([1.0, 2.0, 2.0, 1.0]));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::GainMargin)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Gain margin calculation should succeed");
}

// ============================================================================
// PHASE MARGIN TESTS
// ============================================================================

#[test]
fn test_phase_margin() {
    let mut parameters = HashMap::new();
    parameters.insert("numerator".to_string(), json!([1.0]));
    parameters.insert("denominator".to_string(), json!([1.0, 2.0, 2.0, 1.0]));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::PhaseMargin)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Phase margin calculation should succeed");
}

// ============================================================================
// STEP RESPONSE TESTS
// ============================================================================

#[test]
fn test_step_response_first_order() {
    let mut parameters = HashMap::new();
    parameters.insert("numerator".to_string(), json!([1.0]));
    parameters.insert("denominator".to_string(), json!([1.0, 1.0]));
    parameters.insert("time_span".to_string(), json!(10.0));
    parameters.insert("points".to_string(), json!(100));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::StepResponse)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "First-order step response should succeed");
}

#[test]
fn test_step_response_second_order_underdamped() {
    let mut parameters = HashMap::new();
    parameters.insert("numerator".to_string(), json!([1.0]));
    parameters.insert("denominator".to_string(), json!([1.0, 0.2, 1.0]));
    parameters.insert("time_span".to_string(), json!(20.0));
    parameters.insert("points".to_string(), json!(200));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::StepResponse)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Underdamped step response should succeed"
    );
}

#[test]
fn test_step_response_second_order_critically_damped() {
    let mut parameters = HashMap::new();
    parameters.insert("numerator".to_string(), json!([1.0]));
    parameters.insert("denominator".to_string(), json!([1.0, 2.0, 1.0]));
    parameters.insert("time_span".to_string(), json!(10.0));
    parameters.insert("points".to_string(), json!(100));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::StepResponse)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Critically damped step response should succeed"
    );
}

#[test]
fn test_step_response_second_order_overdamped() {
    let mut parameters = HashMap::new();
    parameters.insert("numerator".to_string(), json!([1.0]));
    parameters.insert("denominator".to_string(), json!([1.0, 5.0, 4.0]));
    parameters.insert("time_span".to_string(), json!(5.0));
    parameters.insert("points".to_string(), json!(100));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::ControlSystems(ControlSystemsOp::StepResponse)),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Overdamped step response should succeed");
}
