//! Comprehensive integration test suite
//!
//! Tests for COMPUTE tool integration operations including:
//! - Symbolic integration
//! - Numeric integration (Simpson, Trapezoidal, Monte Carlo)
//! - Line/Surface/Volume integrals

use computational_engine::create_default_dispatcher;
use computational_engine::engine::equations::*;
use computational_engine::engine::*;
use serde_json::json;
use std::collections::HashMap;

// ============================================================================
// SYMBOLIC INTEGRATION TESTS
// ============================================================================

#[test]
fn test_symbolic_polynomial() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::Symbolic),
        data: json!({
            "expression": "x^2",
            "variable": "x"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should integrate polynomial: {:?}", result);
}

#[test]
fn test_symbolic_trigonometric() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::Symbolic),
        data: json!({
            "expression": "sin(x)",
            "variable": "x"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should integrate trigonometric: {:?}", result);
}

// ============================================================================
// NUMERIC INTEGRATION TESTS
// ============================================================================

#[test]
fn test_numeric_simpson() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::Numeric(NumericIntegration::Simpson)),
        data: json!({
            "expression": "x^2",
            "lower": 0.0,
            "upper": 1.0
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute Simpson integration: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        // Integral of x^2 from 0 to 1 = 1/3
        if let Some(integral) = output.result.get("integral") {
            let value = integral.as_f64().unwrap_or(0.0);
            assert!((value - 1.0/3.0).abs() < 0.01, "Integral of x^2 from 0 to 1 should be ~0.333, got {}", value);
        }
    }
}

#[test]
fn test_numeric_trapezoidal() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::Numeric(NumericIntegration::Trapezoidal)),
        data: json!({
            "expression": "x^2",
            "lower": 0.0,
            "upper": 1.0
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute trapezoidal integration: {:?}", result);
}

#[test]
fn test_monte_carlo_integration() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::MonteCarlo),
        data: json!({
            "expression": "x^2",
            "lower": 0.0,
            "upper": 1.0,
            "samples": 10000
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute Monte Carlo integration: {:?}", result);
}

// ============================================================================
// GEOMETRIC INTEGRAL TESTS
// ============================================================================

#[test]
fn test_line_integral() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::Geometric(GeometricIntegral::Line)),
        data: json!({
            "vector_field": ["x", "y"],
            "curve": "circle",
            "parameters": {"radius": 1.0}
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute line integral: {:?}", result);
}

#[test]
fn test_surface_integral() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::Geometric(GeometricIntegral::Surface)),
        data: json!({
            "function": "x^2 + y^2",
            "surface": "sphere",
            "parameters": {"radius": 1.0}
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute surface integral: {:?}", result);
}

#[test]
fn test_volume_integral() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::Geometric(GeometricIntegral::Volume)),
        data: json!({
            "function": "1",
            "region": "sphere",
            "parameters": {"radius": 1.0}
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute volume integral: {:?}", result);
}

// ============================================================================
// INTEGRAL THEOREM TESTS
// ============================================================================

#[test]
fn test_greens_theorem() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::Theorem(IntegralTheorem::Greens)),
        data: json!({
            "P": "y",
            "Q": "x",
            "region": "unit_square"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply Green's theorem: {:?}", result);
}

#[test]
fn test_stokes_theorem() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::Theorem(IntegralTheorem::Stokes)),
        data: json!({
            "vector_field": ["y", "-x", "0"],
            "surface": "hemisphere"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply Stokes' theorem: {:?}", result);
}

#[test]
fn test_divergence_theorem() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::Theorem(IntegralTheorem::Divergence)),
        data: json!({
            "vector_field": ["x", "y", "z"],
            "volume": "unit_sphere"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply divergence theorem: {:?}", result);
}

// ============================================================================
// COMPLEX ANALYSIS TESTS
// ============================================================================

#[test]
fn test_contour_integral() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::ComplexAnalysis(ComplexIntegral::Contour)),
        data: json!({
            "function": "1/z",
            "contour": "unit_circle"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute contour integral: {:?}", result);
}

#[test]
fn test_residue_integration() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::ComplexAnalysis(ComplexIntegral::Residue)),
        data: json!({
            "function": "1/(z^2 + 1)",
            "poles": [{"re": 0, "im": 1}, {"re": 0, "im": -1}]
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute residue integration: {:?}", result);
}
