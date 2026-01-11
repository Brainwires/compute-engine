//! Comprehensive OPTIMIZE tool test suite
//!
//! Tests for SOLVE tool optimization operations including:
//! - Curve Fitting (Polynomial, Exponential, Logarithmic)
//! - Interpolation
//! - Minimization

use computational_engine::create_default_dispatcher;
use computational_engine::engine::equations::*;
use computational_engine::engine::*;
use serde_json::json;
use std::collections::HashMap;

// ============================================================================
// CURVE FITTING TESTS (via Solve)
// ============================================================================

#[test]
fn test_fit_polynomial() {
    let dispatcher = create_default_dispatcher();

    // Linear data: y = 2x + 1
    let mut params = HashMap::new();
    params.insert("x_data".to_string(), json!([1.0, 2.0, 3.0, 4.0, 5.0]));
    params.insert("y_data".to_string(), json!([3.0, 5.0, 7.0, 9.0, 11.0]));
    params.insert("degree".to_string(), json!(1));

    let request = ToolRequest::Solve(SolveInput {
        equation_type: EquationType::Optimize(OptimizationMethod::Fit(FitMethod::Polynomial)),
        equations: vec![],
        variables: None,
        initial_guess: None,
        boundary_conditions: None,
        domain: None,
        method: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should fit polynomial: {:?}", result);
}

#[test]
fn test_fit_exponential() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("x_data".to_string(), json!([0.0, 1.0, 2.0]));
    params.insert("y_data".to_string(), json!([1.0, 2.7, 7.4]));

    let request = ToolRequest::Solve(SolveInput {
        equation_type: EquationType::Optimize(OptimizationMethod::Fit(FitMethod::Exponential)),
        equations: vec![],
        variables: None,
        initial_guess: None,
        boundary_conditions: None,
        domain: None,
        method: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should fit exponential: {:?}", result);
}

// ============================================================================
// INTERPOLATION TESTS
// ============================================================================

#[test]
fn test_interpolation_linear() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("x_data".to_string(), json!([0.0, 1.0, 2.0]));
    params.insert("y_data".to_string(), json!([0.0, 1.0, 4.0]));

    let request = ToolRequest::Solve(SolveInput {
        equation_type: EquationType::Optimize(OptimizationMethod::Interpolation(InterpolationMethod::Linear)),
        equations: vec![],
        variables: None,
        initial_guess: None,
        boundary_conditions: None,
        domain: None,
        method: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform linear interpolation: {:?}", result);
}

#[test]
fn test_interpolation_spline() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("x_data".to_string(), json!([0.0, 1.0, 2.0, 3.0]));
    params.insert("y_data".to_string(), json!([0.0, 1.0, 0.0, 1.0]));

    let request = ToolRequest::Solve(SolveInput {
        equation_type: EquationType::Optimize(OptimizationMethod::Interpolation(InterpolationMethod::Spline)),
        equations: vec![],
        variables: None,
        initial_guess: None,
        boundary_conditions: None,
        domain: None,
        method: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform spline interpolation: {:?}", result);
}

// ============================================================================
// MINIMIZATION TESTS
// ============================================================================

#[test]
fn test_minimize_conjugate_gradient() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Solve(SolveInput {
        equation_type: EquationType::Optimize(OptimizationMethod::Minimize(MinimizationMethod::ConjugateGradient)),
        equations: vec!["x^2 + y^2".to_string()],
        variables: Some(vec!["x".to_string(), "y".to_string()]),
        initial_guess: None,
        boundary_conditions: None,
        domain: None,
        method: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should minimize with conjugate gradient: {:?}", result);
}

#[test]
fn test_minimize_bfgs() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Solve(SolveInput {
        equation_type: EquationType::Optimize(OptimizationMethod::Minimize(MinimizationMethod::BFGS)),
        equations: vec!["(x-1)^2 + (y-2)^2".to_string()],
        variables: Some(vec!["x".to_string(), "y".to_string()]),
        initial_guess: None,
        boundary_conditions: None,
        domain: None,
        method: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should minimize with BFGS: {:?}", result);
}

// ============================================================================
// SYMBOLIC REGRESSION TESTS
// ============================================================================

#[test]
fn test_symbolic_regression() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("x_data".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert("y_data".to_string(), json!([2.0, 4.0, 6.0]));

    let request = ToolRequest::Solve(SolveInput {
        equation_type: EquationType::Optimize(OptimizationMethod::SymbolicRegression),
        equations: vec![],
        variables: None,
        initial_guess: None,
        boundary_conditions: None,
        domain: None,
        method: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform symbolic regression: {:?}", result);
}
