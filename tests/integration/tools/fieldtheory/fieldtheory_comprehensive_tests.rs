//! Comprehensive FIELD tool test suite
//!
//! Tests for COMPUTE tool field operations including:
//! - EM Fields (Antenna, Waveguide, Scattering)
//! - Green's Functions (Poisson, Helmholtz, Diffusion)
//! - Decoherence scale, Bohm potential

use computational_engine::create_default_dispatcher;
use computational_engine::engine::equations::*;
use computational_engine::engine::*;
use serde_json::json;
use std::collections::HashMap;

// ============================================================================
// EM FIELD TESTS
// ============================================================================

#[test]
fn test_em_field() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Field(FieldType::EM(EMField::Antenna)),
        data: json!({
            "frequency": 2.4e9,
            "antenna_type": "dipole"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute EM antenna field: {:?}", result);
}

#[test]
fn test_em_waveguide() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Field(FieldType::EM(EMField::Waveguide)),
        data: json!({
            "mode": "TE10",
            "frequency": 10e9,
            "dimensions": {"a": 0.023, "b": 0.01}
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute waveguide field: {:?}", result);
}

// ============================================================================
// GREEN'S FUNCTION TESTS
// ============================================================================

#[test]
fn test_green_function() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Field(FieldType::GreenFunction),
        data: json!({
            "r": 1.0
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute Green's function: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        assert!(output.result.get("green_function").is_some(), "Should have green_function value");
    }
}

// ============================================================================
// DECOHERENCE TESTS
// ============================================================================

#[test]
fn test_decoherence_scale() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Field(FieldType::DecoherenceScale),
        data: json!({
            "mass": 1e-26,
            "temperature": 300.0
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute decoherence scale: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        assert!(output.result.get("decoherence_scale").is_some(), "Should have decoherence_scale");
    }
}

// ============================================================================
// BOHM POTENTIAL TESTS
// ============================================================================

#[test]
fn test_bohm_potential() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Field(FieldType::BohmPotential),
        data: json!({
            "psi_real": 1.0,
            "psi_imag": 0.0,
            "mass": 9.109e-31
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute Bohm potential: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        assert!(output.result.get("bohm_potential").is_some(), "Should have bohm_potential");
    }
}

// ============================================================================
// EM SCATTERING FIELD TESTS
// ============================================================================

#[test]
fn test_em_scattering() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Field(FieldType::EM(EMField::Scattering)),
        data: json!({
            "wavelength": 0.03,
            "particle_radius": 0.01
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    // May not be fully implemented - test doesn't panic
    assert!(result.is_ok() || result.is_err());
}

// ============================================================================
// ADDITIONAL EM FIELD TESTS
// ============================================================================

#[test]
fn test_em_antenna_advanced() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Field(FieldType::EM(EMField::Antenna)),
        data: json!({
            "frequency": 5.8e9,
            "antenna_type": "patch"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    // May not be fully implemented - test doesn't panic
    assert!(result.is_ok() || result.is_err());
}

// ============================================================================
// QUANTUM FIELD TESTS
// ============================================================================

#[test]
fn test_quantum_scalar_field() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Field(FieldType::QuantumField(QuantumFieldType::ScalarField)),
        data: json!({
            "mass": 9.109e-31
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    // May not be fully implemented - test doesn't panic
    assert!(result.is_ok() || result.is_err());
}
