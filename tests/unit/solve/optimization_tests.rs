//! Unit tests for solve::optimization

use crate::solve::optimization::*;
use std::collections::HashMap;

// Simple test functions
fn quadratic(x: &[f64]) -> f64 {
    // f(x, y) = x^2 + y^2 (minimum at origin)
    x.iter().map(|&v| v * v).sum()
}

fn quadratic_gradient(x: &[f64]) -> Vec<f64> {
    // gradient = 2*x
    x.iter().map(|&v| 2.0 * v).collect()
}

#[test]
fn test_gradient_descent_basic() {
    let options = OptimizationOptions {
        max_iterations: 100,
        tolerance: 1e-6,
        step_size: 0.1,
        maximize: false,
    };

    let result = gradient_descent(quadratic, quadratic_gradient, vec![5.0, 5.0], options).unwrap();

    assert!(result.success);
    assert!(result.convergence);
    assert!(result.objective_value < 0.1); // Should be close to 0
}

#[test]
fn test_gradient_descent_minimization() {
    let options = OptimizationOptions {
        max_iterations: 200,
        tolerance: 1e-8,
        step_size: 0.05,
        maximize: false,
    };

    let result = gradient_descent(quadratic, quadratic_gradient, vec![10.0], options).unwrap();

    assert!(result.success);
    assert!(result.objective_value < 0.01);
}

#[test]
fn test_nelder_mead_basic() {
    let options = OptimizationOptions {
        max_iterations: 100,
        tolerance: 1e-6,
        step_size: 1.0,
        maximize: false,
    };

    let result = nelder_mead(quadratic, vec![5.0, 5.0], options).unwrap();

    assert!(result.success);
    assert!(result.objective_value < 1.0); // Should converge close to 0
}

#[test]
fn test_curve_fitting_linear() {
    let request = CurveFitRequest {
        x_data: vec![1.0, 2.0, 3.0, 4.0, 5.0],
        y_data: vec![2.0, 4.0, 6.0, 8.0, 10.0],
        model: "linear".to_string(),
    };

    let result = curve_fitting(request).unwrap();

    assert_eq!(result.coefficients.len(), 2); // slope and intercept
    assert!(result.r_squared > 0.95); // Good fit
}

#[test]
fn test_curve_fitting_quadratic() {
    let request = CurveFitRequest {
        x_data: vec![0.0, 1.0, 2.0, 3.0],
        y_data: vec![0.0, 1.0, 4.0, 9.0], // y = x^2
        model: "quadratic".to_string(),
    };

    let result = curve_fitting(request).unwrap();

    assert_eq!(result.coefficients.len(), 3); // a, b, c for ax^2 + bx + c
    assert!(result.r_squared > 0.99);
}

#[test]
fn test_curve_fitting_exponential() {
    let request = CurveFitRequest {
        x_data: vec![0.0, 1.0, 2.0, 3.0],
        y_data: vec![1.0, 2.718, 7.389, 20.086], // y ≈ e^x
        model: "exponential".to_string(),
    };

    let result = curve_fitting(request).unwrap();

    assert_eq!(result.coefficients.len(), 2); // a, b for y = a*e^(bx)
    assert!(result.r_squared > 0.95);
}

#[test]
fn test_curve_fitting_power() {
    let request = CurveFitRequest {
        x_data: vec![1.0, 2.0, 3.0, 4.0],
        y_data: vec![1.0, 4.0, 9.0, 16.0], // y = x^2
        model: "power".to_string(),
    };

    let result = curve_fitting(request).unwrap();

    assert_eq!(result.coefficients.len(), 2); // a, b for y = a*x^b
    assert!(result.r_squared > 0.99);
}

#[test]
fn test_sensitivity_analysis_basic() {
    // Simple objective function: f(x, y) = x^2 + y^2
    fn objective(params: &HashMap<String, f64>) -> f64 {
        let x = params.get("x").unwrap_or(&0.0);
        let y = params.get("y").unwrap_or(&0.0);
        x * x + y * y
    }

    let mut parameters = HashMap::new();
    parameters.insert("x".to_string(), 1.0);
    parameters.insert("y".to_string(), 2.0);

    let request = SensitivityRequest {
        parameters,
        perturbation: 0.01,
    };

    let result = sensitivity_analysis(objective, request).unwrap();

    assert!(result.sensitivities.contains_key("x"));
    assert!(result.sensitivities.contains_key("y"));
    assert!(result.total_sensitivity >= 0.0);
}

#[test]
fn test_optimization_options_default() {
    let options = OptimizationOptions::default();

    assert_eq!(options.max_iterations, 1000);
    assert_eq!(options.tolerance, 1e-6);
    assert_eq!(options.step_size, 0.01);
    assert!(!options.maximize);
}
