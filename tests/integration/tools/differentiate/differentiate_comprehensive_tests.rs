//! Comprehensive differentiation test suite
//!
//! Tests for all DIFFERENTIATE tool operations including:
//! - Symbolic differentiation
//! - Numeric differentiation
//! - Vector calculus (gradient, divergence, curl, laplacian, directional)
//! - Tensor calculus (covariant derivative, Lie derivative, exterior derivative)
//! - Variational calculus
//! - Differential forms

use computational_engine::create_default_dispatcher;
use computational_engine::engine::equations::*;
use computational_engine::engine::*;
use std::collections::HashMap;

// ============================================================================
// SYMBOLIC DIFFERENTIATION TESTS
// ============================================================================

#[test]
fn test_symbolic_polynomial() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::Symbolic,
        expression: "x^3 + 2*x^2 + x + 1".to_string(),
        variables: vec!["x".to_string()],
        order: None,
        evaluate_at: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should differentiate polynomial");
}

#[test]
fn test_symbolic_trigonometric() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::Symbolic,
        expression: "sin(x) + cos(x)".to_string(),
        variables: vec!["x".to_string()],
        order: None,
        evaluate_at: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should differentiate trigonometric functions"
    );
}

#[test]
fn test_symbolic_exponential() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::Symbolic,
        expression: "exp(x)".to_string(),
        variables: vec!["x".to_string()],
        order: None,
        evaluate_at: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should differentiate exponential");
}

// ============================================================================
// VECTOR CALCULUS TESTS
// ============================================================================

#[test]
fn test_gradient() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::VectorCalc(VectorCalcOp::Gradient),
        expression: "x^2 + y^2 + z^2".to_string(),
        variables: vec!["x".to_string(), "y".to_string(), "z".to_string()],
        order: None,
        evaluate_at: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute gradient: {:?}", result);

    if let Ok(ToolResponse::Differentiate(output)) = result {
        assert!(
            output.derivatives.contains_key("gradient"),
            "Should have gradient key"
        );
    }
}

#[test]
fn test_divergence() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("components".to_string(), serde_json::json!(["x", "y", "z"]));

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::VectorCalc(VectorCalcOp::Divergence),
        expression: "vector_field".to_string(),
        variables: vec!["x".to_string(), "y".to_string(), "z".to_string()],
        order: None,
        evaluate_at: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute divergence: {:?}", result);
}

#[test]
fn test_curl() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("components".to_string(), serde_json::json!(["y", "x", "0"]));

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::VectorCalc(VectorCalcOp::Curl),
        expression: "vector_field".to_string(),
        variables: vec!["x".to_string(), "y".to_string(), "z".to_string()],
        order: None,
        evaluate_at: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute curl: {:?}", result);
}

#[test]
fn test_laplacian() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::VectorCalc(VectorCalcOp::Laplacian),
        expression: "x^2 + y^2".to_string(),
        variables: vec!["x".to_string(), "y".to_string()],
        order: None,
        evaluate_at: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute Laplacian: {:?}", result);
}

#[test]
fn test_directional_derivative() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("direction".to_string(), serde_json::json!([1.0, 0.0, 0.0]));

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::VectorCalc(VectorCalcOp::Directional),
        expression: "x^2 + y^2 + z^2".to_string(),
        variables: vec!["x".to_string(), "y".to_string(), "z".to_string()],
        order: None,
        evaluate_at: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should compute directional derivative: {:?}",
        result
    );
}

// ============================================================================
// TENSOR CALCULUS TESTS
// ============================================================================

#[test]
fn test_covariant_derivative() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert(
        "metric".to_string(),
        serde_json::json!([[1.0, 0.0], [0.0, 1.0]]),
    );
    params.insert(
        "tensor_components".to_string(),
        serde_json::json!(["x", "y"]),
    );

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::TensorCalc(TensorDiffOp::Covariant),
        expression: "vector_field".to_string(),
        variables: vec!["x".to_string(), "y".to_string()],
        order: None,
        evaluate_at: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should compute covariant derivative: {:?}",
        result
    );
}

#[test]
fn test_lie_derivative_scalar() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("vector_field".to_string(), serde_json::json!(["y", "x"]));

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::TensorCalc(TensorDiffOp::Lie),
        expression: "x^2 + y^2".to_string(),
        variables: vec!["x".to_string(), "y".to_string()],
        order: None,
        evaluate_at: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should compute Lie derivative of scalar: {:?}",
        result
    );
}

#[test]
fn test_lie_derivative_vector() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("vector_field".to_string(), serde_json::json!(["y", "x"]));
    params.insert(
        "target_components".to_string(),
        serde_json::json!(["x", "y"]),
    );

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::TensorCalc(TensorDiffOp::Lie),
        expression: "vector_field".to_string(),
        variables: vec!["x".to_string(), "y".to_string()],
        order: None,
        evaluate_at: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should compute Lie derivative of vector: {:?}",
        result
    );
}

#[test]
fn test_exterior_derivative_0form() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("form_degree".to_string(), serde_json::json!(0));

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::TensorCalc(TensorDiffOp::ExteriorDerivative),
        expression: "x^2 + y^2".to_string(),
        variables: vec!["x".to_string(), "y".to_string()],
        order: None,
        evaluate_at: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should compute exterior derivative of 0-form: {:?}",
        result
    );
}

#[test]
fn test_exterior_derivative_1form() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("form_degree".to_string(), serde_json::json!(1));
    params.insert("form_components".to_string(), serde_json::json!(["x", "y"]));

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::TensorCalc(TensorDiffOp::ExteriorDerivative),
        expression: "omega".to_string(),
        variables: vec!["x".to_string(), "y".to_string()],
        order: None,
        evaluate_at: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should compute exterior derivative of 1-form: {:?}",
        result
    );
}

// ============================================================================
// VARIATIONAL CALCULUS TESTS
// ============================================================================

#[test]
fn test_variational_euler_lagrange() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("independent_var".to_string(), serde_json::json!("x"));
    params.insert("dependent_var".to_string(), serde_json::json!("y"));
    params.insert("derivative_var".to_string(), serde_json::json!("y_prime"));

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::Variational,
        expression: "y_prime^2 - y^2".to_string(),
        variables: vec!["x".to_string(), "y".to_string(), "y_prime".to_string()],
        order: None,
        evaluate_at: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should compute Euler-Lagrange equation: {:?}",
        result
    );

    if let Ok(ToolResponse::Differentiate(output)) = result {
        assert!(
            output.derivatives.contains_key("euler_lagrange"),
            "Should have Euler-Lagrange equation"
        );
        assert!(
            output.derivatives.contains_key("∂L/∂y"),
            "Should have ∂L/∂y"
        );
        assert!(
            output.derivatives.contains_key("∂L/∂y'"),
            "Should have ∂L/∂y'"
        );
    }
}

// ============================================================================
// DIFFERENTIAL FORMS TESTS
// ============================================================================

#[test]
fn test_differential_forms_scalar() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::DifferentialForms,
        expression: "x^2 + y^2".to_string(),
        variables: vec!["x".to_string(), "y".to_string()],
        order: None,
        evaluate_at: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should compute differential of scalar function: {:?}",
        result
    );
}

// ============================================================================
// NUMERIC DIFFERENTIATION TESTS
// ============================================================================

#[test]
fn test_numeric_differentiation_central() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert(
        "x_values".to_string(),
        serde_json::json!([0.0, 0.1, 0.2, 0.3, 0.4]),
    );
    params.insert(
        "y_values".to_string(),
        serde_json::json!([0.0, 0.01, 0.04, 0.09, 0.16]),
    );
    params.insert("method".to_string(), serde_json::json!("central"));

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::Numeric,
        expression: "f(x)".to_string(),
        variables: vec!["x".to_string()],
        order: Some(vec![1]),
        evaluate_at: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should perform numeric differentiation: {:?}",
        result
    );
}

#[test]
fn test_numeric_differentiation_forward() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert(
        "x_values".to_string(),
        serde_json::json!([0.0, 1.0, 2.0, 3.0]),
    );
    params.insert(
        "y_values".to_string(),
        serde_json::json!([0.0, 1.0, 4.0, 9.0]),
    );
    params.insert("method".to_string(), serde_json::json!("forward"));

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::Numeric,
        expression: "f(x)".to_string(),
        variables: vec!["x".to_string()],
        order: Some(vec![1]),
        evaluate_at: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should perform forward difference: {:?}",
        result
    );
}

// ============================================================================
// ERROR HANDLING TESTS
// ============================================================================

#[test]
fn test_curl_requires_3d() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("components".to_string(), serde_json::json!(["x", "y"]));

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::VectorCalc(VectorCalcOp::Curl),
        expression: "vector_field".to_string(),
        variables: vec!["x".to_string(), "y".to_string()],
        order: None,
        evaluate_at: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Curl should require 3D vector field");
}

#[test]
fn test_directional_requires_direction() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::VectorCalc(VectorCalcOp::Directional),
        expression: "x^2 + y^2".to_string(),
        variables: vec!["x".to_string(), "y".to_string()],
        order: None,
        evaluate_at: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_err(),
        "Directional derivative should require direction parameter"
    );
}

#[test]
fn test_covariant_requires_metric() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::TensorCalc(TensorDiffOp::Covariant),
        expression: "vector_field".to_string(),
        variables: vec!["x".to_string(), "y".to_string()],
        order: None,
        evaluate_at: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_err(),
        "Covariant derivative should require metric tensor"
    );
}

#[test]
fn test_lie_requires_vector_field() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::TensorCalc(TensorDiffOp::Lie),
        expression: "x^2 + y^2".to_string(),
        variables: vec!["x".to_string(), "y".to_string()],
        order: None,
        evaluate_at: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_err(),
        "Lie derivative should require vector field parameter"
    );
}
