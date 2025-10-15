use computational_engine::optimization::*;

#[test]
fn test_curve_fitting() {
    let x_data = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let y_data = vec![1.0, 3.0, 5.0, 7.0, 9.0]; // y = 2x + 1

    let result = curve_fitting(CurveFitRequest {
        x_data,
        y_data,
        model: "linear".to_string(),
    }).unwrap();

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
    }).unwrap();

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
    }).unwrap();

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
    }).unwrap();

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
    }).unwrap();

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
    }).unwrap();

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
    }).unwrap();

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
    }).unwrap();

    let quadratic = curve_fitting(CurveFitRequest {
        x_data,
        y_data,
        model: "quadratic".to_string(),
    }).unwrap();

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
    ).unwrap();

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
    ).unwrap();

    // Rosenbrock is hard to optimize, check if it gets close
    assert!(objective(&[*result.variables.get("x1").unwrap(), *result.variables.get("x2").unwrap()]) < 1.0);
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
    ).unwrap();

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
    ).unwrap();

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
    ).unwrap();

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
    ).unwrap();

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
    ).unwrap();

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
