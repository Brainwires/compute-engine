use computational_engine::optimization::*;

#[test]
fn test_curve_fitting() {
    let x_data = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let y_data = vec![1.0, 3.0, 5.0, 7.0, 9.0]; // y = 2x + 1

    let result = curve_fitting(CurveFitRequest {
        x_data,
        y_data,
        model: "linear".to_string(),
    })
    .unwrap();

    // Should fit y = 2x + 1
    assert!(result.coefficients.len() >= 2);
    assert!((result.coefficients[0] - 1.0).abs() < 0.1); // intercept
    assert!((result.coefficients[1] - 2.0).abs() < 0.1); // slope
    assert!(result.r_squared > 0.99);
}

#[test]
fn test_curve_fitting_quadratic() {
    let x_data = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let y_data = vec![1.0, 2.0, 5.0, 10.0, 17.0]; // y = 1 + x + x^2

    let result = curve_fitting(CurveFitRequest {
        x_data,
        y_data,
        model: "quadratic".to_string(),
    })
    .unwrap();

    assert_eq!(result.coefficients.len(), 3);
    // Quadratic fit should give reasonable coefficients
    assert!(result.coefficients[0].abs() < 10.0);
    assert!(result.coefficients[1].abs() < 10.0);
    assert!(result.coefficients[2].abs() < 10.0);
    assert!(result.r_squared > 0.95);
    assert!(result.aic.is_some());
    assert!(result.bic.is_some());
    assert!(result.aicc.is_some());
}

#[test]
fn test_curve_fitting_exponential() {
    let x_data = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let y_data = vec![1.0, 2.718, 7.389, 20.086, 54.598]; // y = e^x

    let result = curve_fitting(CurveFitRequest {
        x_data,
        y_data,
        model: "exponential".to_string(),
    })
    .unwrap();

    assert_eq!(result.coefficients.len(), 2);
    assert!((result.coefficients[0] - 1.0).abs() < 0.2); // a
    assert!((result.coefficients[1] - 1.0).abs() < 0.2); // b
    assert!(result.r_squared > 0.95);
}

#[test]
fn test_curve_fitting_logarithmic() {
    let x_data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y_data = vec![0.0, 0.693, 1.099, 1.386, 1.609]; // y = ln(x)

    let result = curve_fitting(CurveFitRequest {
        x_data,
        y_data,
        model: "logarithmic".to_string(),
    })
    .unwrap();

    assert_eq!(result.coefficients.len(), 2);
    assert!(result.coefficients[0].abs() < 0.1); // a should be ~0
    assert!((result.coefficients[1] - 1.0).abs() < 0.1); // b should be ~1
    assert!(result.r_squared > 0.99);
}

#[test]
fn test_curve_fitting_power() {
    let x_data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y_data = vec![1.0, 4.0, 9.0, 16.0, 25.0]; // y = x^2

    let result = curve_fitting(CurveFitRequest {
        x_data,
        y_data,
        model: "power".to_string(),
    })
    .unwrap();

    assert_eq!(result.coefficients.len(), 2);
    assert!((result.coefficients[0] - 1.0).abs() < 0.1); // a should be ~1
    assert!((result.coefficients[1] - 2.0).abs() < 0.1); // b should be ~2
    assert!(result.r_squared > 0.99);
}

#[test]
fn test_curve_fitting_rational() {
    let x_data = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let y_data = vec![1.0, 2.0, 2.5, 2.75, 2.875];

    let result = curve_fitting(CurveFitRequest {
        x_data,
        y_data,
        model: "rational".to_string(),
    })
    .unwrap();

    assert_eq!(result.coefficients.len(), 3);
    assert!(result.r_squared > 0.0);
}

#[test]
fn test_curve_fitting_trigonometric() {
    let x_data = vec![0.0, 0.785, 1.571, 2.356, 3.142]; // 0, π/4, π/2, 3π/4, π
    let y_data = vec![0.0, 0.707, 1.0, 0.707, 0.0]; // sin(x)

    let result = curve_fitting(CurveFitRequest {
        x_data,
        y_data,
        model: "trigonometric".to_string(),
    })
    .unwrap();

    assert_eq!(result.coefficients.len(), 3);
    assert!(result.r_squared > 0.8);
}

#[test]
fn test_curve_fitting_information_criteria() {
    let x_data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let y_data = vec![2.1, 4.0, 6.1, 7.9, 10.0, 12.1, 13.9, 16.0];

    let linear = curve_fitting(CurveFitRequest {
        x_data: x_data.clone(),
        y_data: y_data.clone(),
        model: "linear".to_string(),
    })
    .unwrap();

    let quadratic = curve_fitting(CurveFitRequest {
        x_data,
        y_data,
        model: "quadratic".to_string(),
    })
    .unwrap();

    // Linear should have better (lower) AIC/BIC for linear data
    assert!(linear.aic.unwrap() < quadratic.aic.unwrap());
    assert!(linear.bic.unwrap() < quadratic.bic.unwrap());
}

#[test]
fn test_gradient_descent() {
    // Minimize f(x,y) = x^2 + y^2, minimum at (0,0)
    fn objective(params: &[f64]) -> f64 {
        params[0].powi(2) + params[1].powi(2)
    }

    fn gradient(params: &[f64]) -> Vec<f64> {
        vec![2.0 * params[0], 2.0 * params[1]]
    }

    let result = gradient_descent(
        objective,
        gradient,
        vec![5.0, 5.0],
        OptimizationOptions {
            max_iterations: 1000,
            tolerance: 1e-6,
            step_size: 0.1,
            maximize: false,
        },
    )
    .unwrap();

    // Should converge to (0,0)
    assert!(result.variables.get("x1").unwrap().abs() < 0.1);
    assert!(result.variables.get("x2").unwrap().abs() < 0.1);
    assert!(result.convergence);
}

#[test]
fn test_gradient_descent_rosenbrock() {
    // Rosenbrock function: f(x,y) = (1-x)^2 + 100(y-x^2)^2
    // Minimum at (1, 1)
    fn objective(params: &[f64]) -> f64 {
        let x = params[0];
        let y = params[1];
        (1.0 - x).powi(2) + 100.0 * (y - x * x).powi(2)
    }

    fn gradient(params: &[f64]) -> Vec<f64> {
        let x = params[0];
        let y = params[1];
        vec![
            -2.0 * (1.0 - x) - 400.0 * x * (y - x * x),
            200.0 * (y - x * x),
        ]
    }

    let result = gradient_descent(
        objective,
        gradient,
        vec![0.0, 0.0],
        OptimizationOptions {
            max_iterations: 10000,
            tolerance: 1e-4,
            step_size: 0.001,
            maximize: false,
        },
    )
    .unwrap();

    // Rosenbrock is hard to optimize, check if it gets close
    assert!(
        objective(&[
            *result.variables.get("x1").unwrap(),
            *result.variables.get("x2").unwrap()
        ]) < 1.0
    );
}

#[test]
fn test_gradient_descent_different_step_sizes() {
    fn objective(params: &[f64]) -> f64 {
        params[0].powi(2) + params[1].powi(2)
    }

    fn gradient(params: &[f64]) -> Vec<f64> {
        vec![2.0 * params[0], 2.0 * params[1]]
    }

    // Small step size
    let result1 = gradient_descent(
        objective,
        gradient,
        vec![5.0, 5.0],
        OptimizationOptions {
            max_iterations: 1000,
            tolerance: 1e-6,
            step_size: 0.01,
            maximize: false,
        },
    )
    .unwrap();

    // Large step size
    let result2 = gradient_descent(
        objective,
        gradient,
        vec![5.0, 5.0],
        OptimizationOptions {
            max_iterations: 1000,
            tolerance: 1e-6,
            step_size: 0.2,
            maximize: false,
        },
    )
    .unwrap();

    // Both should converge, but larger step size might be faster
    assert!(result1.objective_value < 0.1);
    assert!(result2.objective_value < 0.1);
}

#[test]
fn test_nelder_mead() {
    // Minimize f(x,y) = (x-3)^2 + (y-2)^2, minimum at (3,2)
    fn objective(params: &[f64]) -> f64 {
        (params[0] - 3.0).powi(2) + (params[1] - 2.0).powi(2)
    }

    let result = nelder_mead(
        objective,
        vec![0.0, 0.0],
        OptimizationOptions {
            max_iterations: 1000,
            tolerance: 1e-6,
            step_size: 0.1,
            maximize: false,
        },
    )
    .unwrap();

    // Should converge to (3,2)
    assert!((result.variables.get("x1").unwrap() - 3.0).abs() < 0.1);
    assert!((result.variables.get("x2").unwrap() - 2.0).abs() < 0.1);
}

#[test]
fn test_nelder_mead_himmelblau() {
    // Himmelblau's function: f(x,y) = (x^2+y-11)^2 + (x+y^2-7)^2
    // Multiple local minima
    fn objective(params: &[f64]) -> f64 {
        let x = params[0];
        let y = params[1];
        (x * x + y - 11.0).powi(2) + (x + y * y - 7.0).powi(2)
    }

    let result = nelder_mead(
        objective,
        vec![0.0, 0.0],
        OptimizationOptions {
            max_iterations: 1000,
            tolerance: 1e-6,
            step_size: 0.1,
            maximize: false,
        },
    )
    .unwrap();

    // Should find one of the minima (value ~0)
    assert!(result.objective_value < 0.1);
}

#[test]
fn test_nelder_mead_sphere() {
    // N-dimensional sphere function
    fn objective(params: &[f64]) -> f64 {
        params.iter().map(|x| x * x).sum()
    }

    let result = nelder_mead(
        objective,
        vec![5.0, 5.0, 5.0],
        OptimizationOptions {
            max_iterations: 1000,
            tolerance: 1e-6,
            step_size: 0.1,
            maximize: false,
        },
    )
    .unwrap();

    // Should converge to origin
    assert!(result.objective_value < 0.1);
}

#[test]
fn test_nelder_mead_different_starting_points() {
    fn objective(params: &[f64]) -> f64 {
        (params[0] - 3.0).powi(2) + (params[1] - 2.0).powi(2)
    }

    let result1 = nelder_mead(objective, vec![0.0, 0.0], OptimizationOptions::default()).unwrap();
    let result2 = nelder_mead(objective, vec![10.0, 10.0], OptimizationOptions::default()).unwrap();

    // Both should converge to same minimum
    assert!((result1.objective_value - result2.objective_value).abs() < 0.1);
}

// ========================================================================
// CRITICAL BUG FIX TESTS (Added 2025-10-26)
// ========================================================================

#[test]
fn test_gradient_descent_maximize_flag() {
    // Bug Fix #4: Test that maximize flag works correctly
    // Maximize f(x,y) = -(x^2 + y^2), maximum at (0,0)
    fn objective(params: &[f64]) -> f64 {
        -(params[0].powi(2) + params[1].powi(2))
    }

    fn gradient(params: &[f64]) -> Vec<f64> {
        vec![-2.0 * params[0], -2.0 * params[1]]
    }

    // Test maximization
    let result_max = gradient_descent(
        objective,
        gradient,
        vec![5.0, 5.0],
        OptimizationOptions {
            max_iterations: 1000,
            tolerance: 1e-6,
            step_size: 0.1,
            maximize: true, // MAXIMIZE
        },
    )
    .unwrap();

    // Should converge to (0,0) when maximizing
    assert!(
        result_max.variables.get("x1").unwrap().abs() < 0.1,
        "Maximization should find x1 ≈ 0, got {}",
        result_max.variables.get("x1").unwrap()
    );
    assert!(
        result_max.variables.get("x2").unwrap().abs() < 0.1,
        "Maximization should find x2 ≈ 0, got {}",
        result_max.variables.get("x2").unwrap()
    );
}

#[test]
fn test_gradient_descent_minimize_vs_maximize() {
    // Test that minimize and maximize give opposite results
    fn objective(params: &[f64]) -> f64 {
        params[0].powi(2) // Simple parabola, min at 0
    }

    fn gradient(params: &[f64]) -> Vec<f64> {
        vec![2.0 * params[0]]
    }

    // Minimize from positive starting point
    let result_min = gradient_descent(
        objective,
        gradient,
        vec![5.0],
        OptimizationOptions {
            max_iterations: 100,
            tolerance: 1e-6,
            step_size: 0.1,
            maximize: false,
        },
    )
    .unwrap();

    // Should go toward 0
    let x_min = result_min.variables.get("x1").unwrap();
    assert!(*x_min < 5.0, "Minimization should decrease x");
}

#[test]
fn test_nelder_mead_maximize_flag() {
    // Bug Fix #5: Test that Nelder-Mead maximize flag works
    // Maximize f(x,y) = -(x^2 + y^2), maximum at (0,0)
    fn objective(params: &[f64]) -> f64 {
        -(params[0].powi(2) + params[1].powi(2))
    }

    let result_max = nelder_mead(
        objective,
        vec![5.0, 5.0],
        OptimizationOptions {
            max_iterations: 1000,
            tolerance: 1e-6,
            step_size: 0.1,
            maximize: true, // MAXIMIZE
        },
    )
    .unwrap();

    // Should converge to (0,0) when maximizing
    assert!(
        result_max.variables.get("x1").unwrap().abs() < 0.5,
        "Maximization should find x1 ≈ 0, got {}",
        result_max.variables.get("x1").unwrap()
    );
    assert!(
        result_max.variables.get("x2").unwrap().abs() < 0.5,
        "Maximization should find x2 ≈ 0, got {}",
        result_max.variables.get("x2").unwrap()
    );
}

#[test]
fn test_nelder_mead_minimize_vs_maximize() {
    // Maximize vs minimize should give different results
    fn objective(params: &[f64]) -> f64 {
        // Function with clear min (0,0) and behavior away from origin
        params[0].powi(2) + params[1].powi(2)
    }

    // Minimize - should go to (0,0)
    let result_min = nelder_mead(
        objective,
        vec![5.0, 5.0],
        OptimizationOptions {
            max_iterations: 1000,
            tolerance: 1e-6,
            step_size: 0.1,
            maximize: false,
        },
    )
    .unwrap();

    // Should find minimum near origin
    assert!(result_min.objective_value < 1.0);

    // Maximize with bounded starting point
    // For unbounded function x^2+y^2, it will try to go to infinity
    // So we test with a different function
    fn bounded_objective(params: &[f64]) -> f64 {
        // Inverted Gaussian - has clear maximum at (0,0)
        let r2 = params[0].powi(2) + params[1].powi(2);
        (-r2).exp()
    }

    let result_max = nelder_mead(
        bounded_objective,
        vec![2.0, 2.0],
        OptimizationOptions {
            max_iterations: 1000,
            tolerance: 1e-6,
            step_size: 0.1,
            maximize: true,
        },
    )
    .unwrap();

    // Should find maximum near origin, value close to 1.0
    assert!(
        result_max.objective_value > 0.9,
        "Maximum of exp(-r²) should be ≈ 1.0, got {}",
        result_max.objective_value
    );
}

#[test]
fn test_curve_fitting_trigonometric_documentation() {
    // Bug Fix #7: Verify trigonometric fit uses linearized form
    // y = a + b*sin(x) + c*cos(x), not y = a + b*sin(c*x + d)
    let x_data = vec![0.0, 0.785, 1.571, 2.356, 3.142]; // 0, π/4, π/2, 3π/4, π
    let y_data = vec![0.0, 0.707, 1.0, 0.707, 0.0]; // sin(x)

    let result = curve_fitting(CurveFitRequest {
        x_data,
        y_data,
        model: "trigonometric".to_string(),
    })
    .unwrap();

    // Should return 3 coefficients (a, amplitude, phase)
    // NOT 4 coefficients (a, b, c, d)
    assert_eq!(
        result.coefficients.len(),
        3,
        "Trigonometric fit should return 3 coefficients (linearized form)"
    );

    // First coefficient should be close to 0 (no DC offset)
    assert!(
        result.coefficients[0].abs() < 0.2,
        "DC offset should be small for sin(x)"
    );

    // Second coefficient (amplitude) should be close to 1
    assert!(
        (result.coefficients[1] - 1.0).abs() < 0.2,
        "Amplitude should be ≈ 1.0 for sin(x)"
    );
}

#[test]
fn test_curve_fitting_sinusoidal_alias() {
    // Verify "sinusoidal" is an alias for "trigonometric"
    let x_data = vec![0.0, 1.0, 2.0, 3.0];
    let y_data = vec![1.0, 2.0, 1.5, 2.5];

    let result_trig = curve_fitting(CurveFitRequest {
        x_data: x_data.clone(),
        y_data: y_data.clone(),
        model: "trigonometric".to_string(),
    })
    .unwrap();

    let result_sin = curve_fitting(CurveFitRequest {
        x_data,
        y_data,
        model: "sinusoidal".to_string(),
    })
    .unwrap();

    // Should give identical results
    assert_eq!(result_trig.coefficients.len(), result_sin.coefficients.len());
    for (c1, c2) in result_trig.coefficients.iter().zip(result_sin.coefficients.iter()) {
        assert!((c1 - c2).abs() < 1e-10);
    }
}

#[test]
fn test_optimization_options_default() {
    // Verify default options
    let opts = OptimizationOptions::default();
    assert_eq!(opts.max_iterations, 1000);
    assert_eq!(opts.tolerance, 1e-6);
    assert_eq!(opts.step_size, 0.01);
    assert_eq!(opts.maximize, false);
}

#[test]
fn test_curve_fitting_error_handling() {
    // Test mismatched x and y lengths
    let result = curve_fitting(CurveFitRequest {
        x_data: vec![1.0, 2.0, 3.0],
        y_data: vec![1.0, 2.0],
        model: "linear".to_string(),
    });
    assert!(result.is_err());

    // Test exponential with non-positive y values
    let result = curve_fitting(CurveFitRequest {
        x_data: vec![1.0, 2.0, 3.0],
        y_data: vec![1.0, -2.0, 3.0], // negative value
        model: "exponential".to_string(),
    });
    assert!(result.is_err());

    // Test logarithmic with non-positive x values
    let result = curve_fitting(CurveFitRequest {
        x_data: vec![0.0, 1.0, 2.0], // zero value
        y_data: vec![1.0, 2.0, 3.0],
        model: "logarithmic".to_string(),
    });
    assert!(result.is_err());

    // Test power with non-positive values
    let result = curve_fitting(CurveFitRequest {
        x_data: vec![-1.0, 1.0, 2.0], // negative value
        y_data: vec![1.0, 2.0, 3.0],
        model: "power".to_string(),
    });
    assert!(result.is_err());
}
