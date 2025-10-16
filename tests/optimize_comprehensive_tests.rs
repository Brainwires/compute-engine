//! Comprehensive OPTIMIZE tool test suite
//!
//! Tests for all OPTIMIZE operations including:
//! - Curve Fitting (Polynomial, Exponential, Logarithmic, PowerLaw, Rational, Trigonometric, Custom)
//! - Interpolation (Linear, Polynomial, Spline, Cubic)
//! - Minimization (ConjugateGradient, BFGS, LevenbergMarquardt)
//! - Symbolic Regression (Evolutionary algorithm)
//! - Auto Model Selection (AIC, BIC, AICc, R²)

use computational_engine::create_default_dispatcher;
use computational_engine::engine::*;
use std::collections::HashMap;

// ============================================================================
// CURVE FITTING TESTS
// ============================================================================

#[test]
fn test_fit_polynomial() {
    let dispatcher = create_default_dispatcher();

    // Linear data: y = 2x + 1
    let x_data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y_data = vec![3.0, 5.0, 7.0, 9.0, 11.0];

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Fit(FitMethod::Polynomial),
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should fit polynomial: {:?}", result);

    if let Ok(ToolResponse::Optimize(output)) = result {
        assert!(output.function.is_some(), "Should have fitted function");
        assert!(output.r_squared.is_some(), "Should have R²");
    }
}

#[test]
fn test_fit_exponential() {
    let dispatcher = create_default_dispatcher();

    let x_data = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let y_data = vec![1.0, 2.7, 7.4, 20.1, 54.6]; // Roughly exp(x)

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Fit(FitMethod::Exponential),
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should fit exponential: {:?}", result);
}

#[test]
fn test_fit_logarithmic() {
    let dispatcher = create_default_dispatcher();

    let x_data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y_data = vec![0.0, 0.7, 1.1, 1.4, 1.6]; // Roughly ln(x)

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Fit(FitMethod::Logarithmic),
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should fit logarithmic: {:?}", result);
}

#[test]
fn test_fit_power_law() {
    let dispatcher = create_default_dispatcher();

    let x_data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y_data = vec![1.0, 4.0, 9.0, 16.0, 25.0]; // y = x^2

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Fit(FitMethod::PowerLaw),
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should fit power law: {:?}", result);
}

#[test]
fn test_fit_trigonometric() {
    let dispatcher = create_default_dispatcher();

    // Sine wave data
    let x_data: Vec<f64> = (0..10).map(|i| i as f64 * 0.5).collect();
    let y_data: Vec<f64> = x_data.iter().map(|x| x.sin()).collect();

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Fit(FitMethod::Trigonometric),
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should fit trigonometric: {:?}", result);
}

#[test]
fn test_fit_rational() {
    let dispatcher = create_default_dispatcher();

    let x_data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y_data = vec![0.5, 0.67, 0.75, 0.8, 0.83]; // Roughly x/(1+x)

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Fit(FitMethod::Rational),
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should fit rational function: {:?}", result);
}

#[test]
fn test_fit_custom() {
    let dispatcher = create_default_dispatcher();

    let x_data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y_data = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let mut params = HashMap::new();
    params.insert("function_form".to_string(), serde_json::json!("a*x + b"));

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Fit(FitMethod::Custom),
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should handle custom fit: {:?}", result);
}

// ============================================================================
// INTERPOLATION TESTS
// ============================================================================

#[test]
fn test_interpolation_linear() {
    let dispatcher = create_default_dispatcher();

    let x_data = vec![0.0, 1.0, 2.0, 3.0];
    let y_data = vec![0.0, 1.0, 4.0, 9.0];

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Interpolation(InterpolationMethod::Linear),
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should perform linear interpolation: {:?}",
        result
    );

    if let Ok(ToolResponse::Optimize(output)) = result {
        assert_eq!(
            output.r_squared,
            Some(1.0),
            "Interpolation should have R²=1"
        );
    }
}

#[test]
fn test_interpolation_polynomial() {
    let dispatcher = create_default_dispatcher();

    let x_data = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let y_data = vec![0.0, 1.0, 4.0, 9.0, 16.0];

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Interpolation(InterpolationMethod::Polynomial),
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should perform polynomial interpolation: {:?}",
        result
    );
}

#[test]
fn test_interpolation_spline() {
    let dispatcher = create_default_dispatcher();

    let x_data = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let y_data = vec![0.0, 1.0, 0.0, -1.0, 0.0];

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Interpolation(InterpolationMethod::Spline),
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should perform spline interpolation: {:?}",
        result
    );

    if let Ok(ToolResponse::Optimize(output)) = result {
        let metadata = output.metadata.as_ref().unwrap();
        assert_eq!(
            metadata.get("method").and_then(|v| v.as_str()).unwrap(),
            "cubic_spline"
        );
    }
}

#[test]
fn test_interpolation_cubic() {
    let dispatcher = create_default_dispatcher();

    let x_data = vec![0.0, 1.0, 2.0, 3.0];
    let y_data = vec![0.0, 1.0, 8.0, 27.0]; // y = x^3

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Interpolation(InterpolationMethod::Cubic),
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should perform cubic interpolation: {:?}",
        result
    );
}

// ============================================================================
// MINIMIZATION TESTS
// ============================================================================

#[test]
fn test_minimize_conjugate_gradient() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Minimize(MinimizationMethod::ConjugateGradient),
        data: None,
        objective: Some("x^2 + y^2".to_string()),
        initial_guess: Some(vec![1.0, 1.0]),
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should handle conjugate gradient: {:?}",
        result
    );

    if let Ok(ToolResponse::Optimize(output)) = result {
        assert!(output.convergence.is_some(), "Should have convergence info");
    }
}

#[test]
fn test_minimize_bfgs() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Minimize(MinimizationMethod::BFGS),
        data: None,
        objective: Some("x^2".to_string()),
        initial_guess: Some(vec![5.0]),
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should handle BFGS: {:?}", result);
}

#[test]
fn test_minimize_levenberg_marquardt() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Minimize(MinimizationMethod::LevenbergMarquardt),
        data: None,
        objective: Some("nonlinear least squares".to_string()),
        initial_guess: Some(vec![1.0]),
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should handle Levenberg-Marquardt: {:?}",
        result
    );

    if let Ok(ToolResponse::Optimize(output)) = result {
        let metadata = output.metadata.as_ref().unwrap();
        assert_eq!(
            metadata.get("method").and_then(|v| v.as_str()).unwrap(),
            "levenberg_marquardt"
        );
    }
}

#[test]
fn test_minimize_gradient_descent_not_supported() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Minimize(MinimizationMethod::GradientDescent),
        data: None,
        objective: Some("x^2".to_string()),
        initial_guess: Some(vec![1.0]),
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_err(),
        "Gradient descent should return error (not supported via JSON API)"
    );
}

#[test]
fn test_minimize_nelder_mead_not_supported() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Minimize(MinimizationMethod::NelderMead),
        data: None,
        objective: Some("x^2".to_string()),
        initial_guess: Some(vec![1.0]),
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_err(),
        "Nelder-Mead should return error (not supported via JSON API)"
    );
}

// ============================================================================
// SYMBOLIC REGRESSION TESTS
// ============================================================================

#[test]
fn test_symbolic_regression_simple() {
    let dispatcher = create_default_dispatcher();

    // Simple linear relationship: y = 2x
    let x_data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y_data = vec![2.0, 4.0, 6.0, 8.0, 10.0];

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::SymbolicRegression,
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should perform symbolic regression: {:?}",
        result
    );

    if let Ok(ToolResponse::Optimize(output)) = result {
        assert!(output.function.is_some(), "Should discover an equation");
        assert!(output.convergence.is_some(), "Should have convergence data");
    }
}

#[test]
fn test_symbolic_regression_with_domain() {
    let dispatcher = create_default_dispatcher();

    let x_data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y_data = vec![1.0, 4.0, 9.0, 16.0, 25.0]; // y = x^2

    let mut params = HashMap::new();
    params.insert("domain".to_string(), serde_json::json!("physics"));
    params.insert("max_complexity".to_string(), serde_json::json!(5));

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::SymbolicRegression,
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should perform symbolic regression with physics domain: {:?}",
        result
    );
}

// ============================================================================
// AUTO MODEL SELECTION TESTS
// ============================================================================

#[test]
fn test_auto_selection_aicc() {
    let dispatcher = create_default_dispatcher();

    // Data that could fit multiple models
    let x_data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y_data = vec![2.0, 4.0, 6.0, 8.0, 10.0]; // Linear

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Auto {
            criteria: SelectionCriteria::Aicc,
            candidates: vec![],
        },
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should auto-select model using AICc: {:?}",
        result
    );

    if let Ok(ToolResponse::Optimize(output)) = result {
        let metadata = output.metadata.as_ref().unwrap();
        assert_eq!(
            metadata
                .get("auto_selection")
                .and_then(|v| v.as_bool())
                .unwrap(),
            true
        );
        assert_eq!(
            metadata.get("criteria").and_then(|v| v.as_str()).unwrap(),
            "AICc"
        );
    }
}

#[test]
fn test_auto_selection_aic() {
    let dispatcher = create_default_dispatcher();

    let x_data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y_data = vec![1.0, 4.0, 9.0, 16.0, 25.0]; // Quadratic

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Auto {
            criteria: SelectionCriteria::Aic,
            candidates: vec![],
        },
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should auto-select model using AIC: {:?}",
        result
    );
}

#[test]
fn test_auto_selection_bic() {
    let dispatcher = create_default_dispatcher();

    let x_data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y_data = vec![1.0, 2.0, 4.0, 8.0, 16.0]; // Exponential-ish

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Auto {
            criteria: SelectionCriteria::Bic,
            candidates: vec![],
        },
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should auto-select model using BIC: {:?}",
        result
    );
}

#[test]
fn test_auto_selection_r_squared() {
    let dispatcher = create_default_dispatcher();

    let x_data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y_data = vec![3.0, 5.0, 7.0, 9.0, 11.0]; // y = 2x + 1

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Auto {
            criteria: SelectionCriteria::RSquared,
            candidates: vec![],
        },
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should auto-select model using R²: {:?}",
        result
    );
}

#[test]
fn test_auto_selection_with_candidates() {
    let dispatcher = create_default_dispatcher();

    let x_data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y_data = vec![2.0, 4.0, 6.0, 8.0, 10.0];

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Auto {
            criteria: SelectionCriteria::Aicc,
            candidates: vec!["linear".to_string(), "exponential".to_string()],
        },
        data: Some((x_data, y_data)),
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Should auto-select from specified candidates: {:?}",
        result
    );
}

// ============================================================================
// DIMENSIONAL ANALYSIS TESTS
// ============================================================================

#[test]
fn test_dimensional_analysis_not_implemented() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::DimensionalAnalysis(DimAnalysisMethod::BuckinghamPi),
        data: None,
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_err(),
        "Dimensional analysis should return error (not yet implemented)"
    );
}

// ============================================================================
// ERROR HANDLING TESTS
// ============================================================================

#[test]
fn test_fit_requires_data() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Fit(FitMethod::Polynomial),
        data: None, // Missing data
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Fitting should require data");
}

#[test]
fn test_interpolation_requires_minimum_points() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Interpolation(InterpolationMethod::Linear),
        data: Some((vec![1.0], vec![1.0])), // Only 1 point
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_err(),
        "Linear interpolation should require at least 2 points"
    );
}

#[test]
fn test_spline_requires_minimum_points() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Optimize(OptimizeInput {
        method: OptimizationMethod::Interpolation(InterpolationMethod::Spline),
        data: Some((vec![1.0, 2.0], vec![1.0, 2.0])), // Only 2 points
        objective: None,
        initial_guess: None,
        constraints: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_err(),
        "Spline interpolation should require at least 3 points"
    );
}
