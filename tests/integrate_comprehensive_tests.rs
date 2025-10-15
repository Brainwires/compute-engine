//! Comprehensive integration test suite
//!
//! Tests for all INTEGRATE tool operations including:
//! - Symbolic integration
//! - Numeric integration
//! - Geometric integrals (line, surface, volume, contour)
//! - Integral theorems (Green's, Stokes, Divergence, Cauchy)
//! - Complex analysis integrals (residue, Cauchy formula)

use computational_engine::engine::*;
use computational_engine::engine::equations::*;
use computational_engine::create_default_dispatcher;
use std::collections::HashMap;

// ============================================================================
// SYMBOLIC INTEGRATION TESTS
// ============================================================================

#[test]
fn test_symbolic_polynomial() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Symbolic,
        expression: "x^2".to_string(),
        variables: vec!["x".to_string()],
        limits: None,
        path: None,
        method: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should integrate polynomial: {:?}", result);

    if let Ok(ToolResponse::Integrate(output)) = result {
        assert!(output.symbolic.is_some(), "Should have symbolic result");
    }
}

#[test]
fn test_symbolic_trigonometric() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Symbolic,
        expression: "sin(x)".to_string(),
        variables: vec!["x".to_string()],
        limits: None,
        path: None,
        method: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should integrate trigonometric function: {:?}", result);
}

// ============================================================================
// NUMERIC INTEGRATION TESTS
// ============================================================================

#[test]
fn test_numeric_trapezoidal() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("coefficients".to_string(), serde_json::json!([1.0, 0.0, 1.0]));
    params.insert("function_type".to_string(), serde_json::json!("polynomial"));

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Numeric(NumericIntegration::Trapezoidal),
        expression: "x^2".to_string(),
        variables: vec!["x".to_string()],
        limits: Some(vec![[0.0, 1.0]]),
        path: None,
        method: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform trapezoidal integration: {:?}", result);
}

#[test]
fn test_numeric_simpson() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("coefficients".to_string(), serde_json::json!([1.0, 0.0, 1.0]));
    params.insert("function_type".to_string(), serde_json::json!("polynomial"));

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Numeric(NumericIntegration::Simpson),
        expression: "x^2".to_string(),
        variables: vec!["x".to_string()],
        limits: Some(vec![[0.0, 1.0]]),
        path: None,
        method: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform Simpson integration: {:?}", result);
}

// ============================================================================
// GEOMETRIC INTEGRAL TESTS
// ============================================================================

#[test]
fn test_line_integral_simple() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Geometric(GeometricIntegral::Line),
        expression: "t".to_string(),
        variables: vec!["t".to_string()],
        limits: Some(vec![[0.0, 1.0]]),
        path: None,
        method: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute line integral: {:?}", result);
}

#[test]
fn test_line_integral_parametric() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("path".to_string(), serde_json::json!(["t", "t^2"]));

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Geometric(GeometricIntegral::Line),
        expression: "x + y".to_string(),
        variables: vec!["x".to_string(), "y".to_string()],
        limits: Some(vec![[0.0, 1.0]]),
        path: None,
        method: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute parametric line integral: {:?}", result);
}

#[test]
fn test_surface_integral() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Geometric(GeometricIntegral::Surface),
        expression: "1".to_string(),
        variables: vec!["u".to_string(), "v".to_string()],
        limits: Some(vec![[0.0, 1.0], [0.0, 1.0]]),
        path: None,
        method: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute surface integral: {:?}", result);
}

#[test]
fn test_volume_integral() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Geometric(GeometricIntegral::Volume),
        expression: "1".to_string(),
        variables: vec!["x".to_string(), "y".to_string(), "z".to_string()],
        limits: Some(vec![[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]),
        path: None,
        method: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute volume integral: {:?}", result);

    if let Ok(ToolResponse::Integrate(output)) = result {
        if let Some(val) = output.result.as_f64() {
            assert!((val - 1.0).abs() < 0.1, "Volume should be approximately 1.0");
        }
    }
}

#[test]
fn test_contour_integral_circle() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("contour".to_string(), serde_json::json!("circle"));
    params.insert("radius".to_string(), serde_json::json!(1.0));

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Geometric(GeometricIntegral::Contour),
        expression: "1".to_string(),
        variables: vec!["z".to_string()],
        limits: None,
        path: None,
        method: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute contour integral: {:?}", result);
}

// ============================================================================
// INTEGRAL THEOREM TESTS
// ============================================================================

#[test]
fn test_greens_theorem() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("components".to_string(), serde_json::json!(["y", "x"]));

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Theorem(IntegralTheorem::Greens),
        expression: "vector_field".to_string(),
        variables: vec!["x".to_string(), "y".to_string()],
        limits: Some(vec![[0.0, 1.0], [0.0, 1.0]]),
        path: None,
        method: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply Green's theorem: {:?}", result);

    if let Ok(ToolResponse::Integrate(output)) = result {
        assert!(output.symbolic.is_some(), "Should have symbolic representation");
    }
}

#[test]
fn test_stokes_theorem() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("components".to_string(), serde_json::json!(["y", "x", "0"]));

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Theorem(IntegralTheorem::Stokes),
        expression: "vector_field".to_string(),
        variables: vec!["x".to_string(), "y".to_string(), "z".to_string()],
        limits: Some(vec![[0.0, 1.0], [0.0, 1.0]]),
        path: None,
        method: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply Stokes' theorem: {:?}", result);
}

#[test]
fn test_divergence_theorem() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("components".to_string(), serde_json::json!(["x", "y", "z"]));

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Theorem(IntegralTheorem::Divergence),
        expression: "vector_field".to_string(),
        variables: vec!["x".to_string(), "y".to_string(), "z".to_string()],
        limits: Some(vec![[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]),
        path: None,
        method: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply Divergence theorem: {:?}", result);

    if let Ok(ToolResponse::Integrate(output)) = result {
        assert!(output.symbolic.is_some(), "Should have symbolic representation");
    }
}

#[test]
fn test_cauchy_theorem_analytic() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("analytic".to_string(), serde_json::json!(true));

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Theorem(IntegralTheorem::CauchyIntegral),
        expression: "z^2".to_string(),
        variables: vec!["z".to_string()],
        limits: None,
        path: None,
        method: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply Cauchy's theorem: {:?}", result);

    if let Ok(ToolResponse::Integrate(output)) = result {
        // Analytic function should have integral = 0
        if let Some(val) = output.result.as_f64() {
            assert_eq!(val, 0.0, "Cauchy theorem: analytic function integral should be 0");
        }
    }
}

// ============================================================================
// COMPLEX ANALYSIS TESTS
// ============================================================================

#[test]
fn test_residue_theorem() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("residues".to_string(), serde_json::json!([1.0, 2.0]));

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::ComplexAnalysis(ComplexIntegral::Residue),
        expression: "f(z)".to_string(),
        variables: vec!["z".to_string()],
        limits: None,
        path: None,
        method: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute residue integral: {:?}", result);

    if let Ok(ToolResponse::Integrate(output)) = result {
        assert!(output.symbolic.is_some(), "Should have symbolic representation");
    }
}

#[test]
fn test_cauchy_formula() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("point_value".to_string(), serde_json::json!(3.14));

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::ComplexAnalysis(ComplexIntegral::Cauchy),
        expression: "f(z)".to_string(),
        variables: vec!["z".to_string()],
        limits: None,
        path: None,
        method: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply Cauchy formula: {:?}", result);
}

// ============================================================================
// MONTE CARLO INTEGRATION TESTS
// ============================================================================

#[test]
fn test_monte_carlo_integration() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::MonteCarlo,
        expression: "1".to_string(),
        variables: vec!["x".to_string()],
        limits: Some(vec![[0.0, 1.0]]),
        path: None,
        method: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform Monte Carlo integration: {:?}", result);

    if let Ok(ToolResponse::Integrate(output)) = result {
        assert!(output.error_estimate.is_some(), "Should have error estimate");
    }
}

// ============================================================================
// ERROR HANDLING TESTS
// ============================================================================

#[test]
fn test_greens_requires_components() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Theorem(IntegralTheorem::Greens),
        expression: "vector_field".to_string(),
        variables: vec!["x".to_string(), "y".to_string()],
        limits: Some(vec![[0.0, 1.0], [0.0, 1.0]]),
        path: None,
        method: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Green's theorem should require components parameter");
}

#[test]
fn test_stokes_requires_3d() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("components".to_string(), serde_json::json!(["x", "y"]));

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Theorem(IntegralTheorem::Stokes),
        expression: "vector_field".to_string(),
        variables: vec!["x".to_string(), "y".to_string()],
        limits: Some(vec![[0.0, 1.0], [0.0, 1.0]]),
        path: None,
        method: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Stokes' theorem should require 3D vector field");
}

#[test]
fn test_volume_requires_3_limits() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Geometric(GeometricIntegral::Volume),
        expression: "1".to_string(),
        variables: vec!["x".to_string(), "y".to_string()],
        limits: Some(vec![[0.0, 1.0], [0.0, 1.0]]),
        path: None,
        method: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Volume integral should require three limits");
}
