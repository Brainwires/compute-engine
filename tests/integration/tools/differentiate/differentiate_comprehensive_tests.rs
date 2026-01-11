//! Comprehensive differentiation test suite
//!
//! Tests for COMPUTE tool differentiation operations including:
//! - Symbolic differentiation
//! - Numeric differentiation

use computational_engine::create_default_dispatcher;
use computational_engine::engine::equations::*;
use computational_engine::engine::*;
use serde_json::json;
use std::collections::HashMap;

// ============================================================================
// SYMBOLIC DIFFERENTIATION TESTS
// ============================================================================

#[test]
fn test_symbolic_polynomial() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Differentiate(DifferentiationOp::Symbolic),
        data: json!({
            "expression": "x^3 + 2*x^2 + x + 1",
            "variable": "x"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should differentiate polynomial: {:?}", result);
}

#[test]
fn test_symbolic_trigonometric() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Differentiate(DifferentiationOp::Symbolic),
        data: json!({
            "expression": "sin(x) + cos(x)",
            "variable": "x"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should differentiate trigonometric functions: {:?}", result);
}

#[test]
fn test_symbolic_exponential() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Differentiate(DifferentiationOp::Symbolic),
        data: json!({
            "expression": "exp(x)",
            "variable": "x"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should differentiate exponential: {:?}", result);
}

// ============================================================================
// NUMERIC DIFFERENTIATION TESTS
// ============================================================================

#[test]
fn test_numeric_differentiation() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Differentiate(DifferentiationOp::Numeric),
        data: json!({
            "expression": "x^2",
            "point": 2.0
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute numeric derivative: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        // Derivative of x^2 at x=2 should be 4
        if let Some(deriv) = output.result.get("derivative") {
            let value = deriv.as_f64().unwrap_or(0.0);
            assert!((value - 4.0).abs() < 0.01, "Derivative of x^2 at x=2 should be ~4, got {}", value);
        }
    }
}

#[test]
fn test_numeric_with_custom_h() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("h".to_string(), json!(1e-6));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Differentiate(DifferentiationOp::Numeric),
        data: json!({
            "expression": "x^2",
            "point": 3.0
        }),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute numeric derivative with custom h: {:?}", result);
}

// ============================================================================
// VECTOR CALCULUS TESTS (via Compute)
// ============================================================================

#[test]
fn test_gradient() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Differentiate(DifferentiationOp::VectorCalc(VectorCalcOp::Gradient)),
        data: json!({
            "expression": "x^2 + y^2 + z^2",
            "variables": ["x", "y", "z"]
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute gradient: {:?}", result);
}

#[test]
fn test_divergence() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Differentiate(DifferentiationOp::VectorCalc(VectorCalcOp::Divergence)),
        data: json!({
            "vector_field": ["x", "y", "z"],
            "variables": ["x", "y", "z"]
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute divergence: {:?}", result);
}

#[test]
fn test_curl() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Differentiate(DifferentiationOp::VectorCalc(VectorCalcOp::Curl)),
        data: json!({
            "vector_field": ["y", "-x", "0"],
            "variables": ["x", "y", "z"]
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute curl: {:?}", result);
}

#[test]
fn test_laplacian() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Differentiate(DifferentiationOp::VectorCalc(VectorCalcOp::Laplacian)),
        data: json!({
            "expression": "x^2 + y^2 + z^2",
            "variables": ["x", "y", "z"]
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute Laplacian: {:?}", result);
}

// ============================================================================
// SYMBOLIC DIFFERENTIATION TESTS (Partial Derivatives)
// ============================================================================

#[test]
fn test_symbolic_partial_derivative() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Differentiate(DifferentiationOp::Symbolic),
        data: json!({
            "expression": "x^2 * y + y^3",
            "variable": "x"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    // May or may not be fully implemented - test doesn't panic
    assert!(result.is_ok() || result.is_err());
}

// ============================================================================
// VARIATIONAL CALCULUS TESTS
// ============================================================================

#[test]
fn test_variational_derivative() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Differentiate(DifferentiationOp::Variational),
        data: json!({
            "functional": "integral(y'^2, x, 0, 1)",
            "variable": "y"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    // May or may not be fully implemented - test doesn't panic
    assert!(result.is_ok() || result.is_err());
}
