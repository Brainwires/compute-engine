// Unit tests for tools::numerical_methods::mod
use computational_engine::tools::numerical_methods::mod::*;

use super::*;

    #[test]
    fn test_bisection() {
        let result = find_root(RootFindingRequest {
            method: "bisection".to_string(),
            initial_guess: 0.0,
            interval: Some((1.0, 2.0)),
            tolerance: 1e-6,
            max_iterations: 100,
            function_type: "polynomial".to_string(),
            coefficients: Some(vec![-2.0, 0.0, 1.0]), // x^2 - 2 = 0
        })
        .unwrap();

        assert!((result.root - 1.41421356).abs() < 1e-5);
    }

    #[test]
    fn test_trapezoidal_integration() {
        let result = integrate(IntegrationRequest {
            method: "trapezoidal".to_string(),
            lower_bound: 0.0,
            upper_bound: 1.0,
            num_points: 100,
            function_type: "polynomial".to_string(),
            coefficients: Some(vec![0.0, 0.0, 1.0]), // x^2
        })
        .unwrap();

        // Integral of x^2 from 0 to 1 is 1/3
        assert!((result.integral - 0.333333).abs() < 0.01);
    }

    #[test]
    fn test_newton_method() {
        let result = find_root(RootFindingRequest {
            method: "newton".to_string(),
            initial_guess: 1.5,
            interval: None,
            tolerance: 1e-6,
            max_iterations: 100,
            function_type: "polynomial".to_string(),
            coefficients: Some(vec![-2.0, 0.0, 1.0]), // x^2 - 2 = 0
        })
        .unwrap();

        assert!((result.root - 1.41421356).abs() < 1e-5);
        assert!(result.converged);
    }

    #[test]
    fn test_secant_method() {
        let result = find_root(RootFindingRequest {
            method: "secant".to_string(),
            initial_guess: 1.5,
            interval: None,
            tolerance: 1e-6,
            max_iterations: 100,
            function_type: "polynomial".to_string(),
            coefficients: Some(vec![-2.0, 0.0, 1.0]),
        })
        .unwrap();

        assert!((result.root - 1.41421356).abs() < 1e-3);
    }

    #[test]
    fn test_simpson_integration() {
        let result = integrate(IntegrationRequest {
            method: "simpson".to_string(),
            lower_bound: 0.0,
            upper_bound: 1.0,
            num_points: 100,
            function_type: "polynomial".to_string(),
            coefficients: Some(vec![0.0, 0.0, 1.0]), // x^2
        })
        .unwrap();

        assert!((result.integral - 0.333333).abs() < 0.001);
    }

    #[test]
    fn test_gauss_integration() {
        let result = integrate(IntegrationRequest {
            method: "gauss".to_string(),
            lower_bound: 0.0,
            upper_bound: 1.0,
            num_points: 2,
            function_type: "polynomial".to_string(),
            coefficients: Some(vec![0.0, 0.0, 1.0]),
        })
        .unwrap();

        assert!((result.integral - 0.333333).abs() < 0.1);
    }

    #[test]
    fn test_simpson_odd_intervals() {
        let result = integrate(IntegrationRequest {
            method: "simpson".to_string(),
            lower_bound: 0.0,
            upper_bound: 1.0,
            num_points: 99,
            function_type: "polynomial".to_string(),
            coefficients: Some(vec![0.0, 0.0, 1.0]),
        });

        assert!(result.is_err());
    }

    #[test]
    fn test_linear_interpolation() {
        let result = interpolate(InterpolationRequest {
            method: "linear".to_string(),
            x_values: vec![0.0, 1.0, 2.0],
            y_values: vec![0.0, 1.0, 4.0],
            interpolate_at: vec![0.5, 1.5],
        })
        .unwrap();

        assert!((result.interpolated_values[0] - 0.5).abs() < 1e-6);
        assert!((result.interpolated_values[1] - 2.5).abs() < 1e-6);
    }

    #[test]
    fn test_lagrange_interpolation() {
        let result = interpolate(InterpolationRequest {
            method: "lagrange".to_string(),
            x_values: vec![0.0, 1.0, 2.0],
            y_values: vec![0.0, 1.0, 4.0],
            interpolate_at: vec![0.5],
        })
        .unwrap();

        // Should match polynomial x^2
        assert!((result.interpolated_values[0] - 0.25).abs() < 1e-6);
    }

    #[test]
    fn test_interpolation_mismatched_lengths() {
        let result = interpolate(InterpolationRequest {
            method: "linear".to_string(),
            x_values: vec![0.0, 1.0],
            y_values: vec![0.0],
            interpolate_at: vec![0.5],
        });

        assert!(result.is_err());
    }

    #[test]
    fn test_gauss_elimination() {
        let result = solve_linear_system(LinearSystemRequest {
            method: "gauss".to_string(),
            matrix: vec![vec![2.0, 1.0], vec![1.0, 3.0]],
            rhs: vec![5.0, 6.0],
            tolerance: None,
            max_iterations: None,
        })
        .unwrap();

        assert!((result.solution[0] - 1.8).abs() < 1e-6);
        assert!((result.solution[1] - 1.4).abs() < 1e-6);
    }

    #[test]
    fn test_jacobi_iteration() {
        let result = solve_linear_system(LinearSystemRequest {
            method: "jacobi".to_string(),
            matrix: vec![vec![4.0, 1.0], vec![1.0, 3.0]],
            rhs: vec![1.0, 2.0],
            tolerance: Some(1e-6),
            max_iterations: Some(100),
        })
        .unwrap();

        assert!(result.solution.len() == 2);
        assert!(result.iterations.is_some());
    }

    #[test]
    fn test_linear_system_dimension_mismatch() {
        let result = solve_linear_system(LinearSystemRequest {
            method: "gauss".to_string(),
            matrix: vec![vec![1.0, 2.0]],
            rhs: vec![1.0, 2.0],
            tolerance: None,
            max_iterations: None,
        });

        assert!(result.is_err());
    }

    #[test]
    fn test_forward_differentiation() {
        let result = differentiate(DifferentiationRequest {
            method: "forward".to_string(),
            x_values: vec![0.0, 1.0, 2.0, 3.0],
            y_values: vec![0.0, 1.0, 4.0, 9.0], // x^2
            order: 1,
        })
        .unwrap();

        assert!(result.derivatives.len() == 4);
        assert!((result.derivatives[0] - 1.0).abs() < 0.5);
    }

    #[test]
    fn test_backward_differentiation() {
        let result = differentiate(DifferentiationRequest {
            method: "backward".to_string(),
            x_values: vec![0.0, 1.0, 2.0, 3.0],
            y_values: vec![0.0, 1.0, 4.0, 9.0],
            order: 1,
        })
        .unwrap();

        assert!(result.derivatives.len() == 4);
    }

    #[test]
    fn test_central_differentiation() {
        let result = differentiate(DifferentiationRequest {
            method: "central".to_string(),
            x_values: vec![0.0, 1.0, 2.0, 3.0],
            y_values: vec![0.0, 1.0, 4.0, 9.0],
            order: 1,
        })
        .unwrap();

        assert!(result.derivatives.len() == 4);
        assert!((result.derivatives[1] - 2.0).abs() < 0.5);
    }

    #[test]
    fn test_differentiation_too_few_points() {
        let result = differentiate(DifferentiationRequest {
            method: "forward".to_string(),
            x_values: vec![0.0],
            y_values: vec![0.0],
            order: 1,
        });

        assert!(result.is_err());
    }

    #[test]
    fn test_euler_ode() {
        let result = solve_ode(ODESolverRequest {
            method: "euler".to_string(),
            initial_value: 1.0,
            t_start: 0.0,
            t_end: 1.0,
            step_size: 0.1,
            derivative_expression: None,
        })
        .unwrap();

        assert!(result.t_values.len() > 0);
        assert_eq!(result.t_values.len(), result.y_values.len());
        assert!(result.steps_taken > 0);
    }

    #[test]
    fn test_rk4_ode() {
        let result = solve_ode(ODESolverRequest {
            method: "rk4".to_string(),
            initial_value: 1.0,
            t_start: 0.0,
            t_end: 1.0,
            step_size: 0.1,
            derivative_expression: None,
        })
        .unwrap();

        assert!(result.method_used == "rk4");
        assert!(result.y_values.len() > 0);
    }

    #[test]
    fn test_adaptive_rk45_ode() {
        let result = solve_ode(ODESolverRequest {
            method: "rk45".to_string(),
            initial_value: 1.0,
            t_start: 0.0,
            t_end: 1.0,
            step_size: 0.1,
            derivative_expression: None,
        })
        .unwrap();

        assert!(result.y_values.len() > 0);
    }

    #[test]
    fn test_pde_heat_equation() {
        let result = solve_pde(PDESolverRequest {
            method: "heat_equation".to_string(),
            boundary_conditions: vec![0.0, 0.0],
            initial_conditions: vec![1.0, 1.0, 1.0, 1.0, 1.0],
            spatial_steps: 5,
            time_steps: 10,
            dx: 0.1,
            dt: 0.01,
        })
        .unwrap();

        assert_eq!(result.solution.len(), 10);
        assert_eq!(result.final_state.len(), 5);
    }

    #[test]
    fn test_bisection_same_sign() {
        let result = find_root(RootFindingRequest {
            method: "bisection".to_string(),
            initial_guess: 0.0,
            interval: Some((2.0, 3.0)),
            tolerance: 1e-6,
            max_iterations: 100,
            function_type: "polynomial".to_string(),
            coefficients: Some(vec![-2.0, 0.0, 1.0]),
        });

        assert!(result.is_err());
    }

    #[test]
    fn test_newton_max_iterations() {
        let result = find_root(RootFindingRequest {
            method: "newton".to_string(),
            initial_guess: 10.0,
            interval: None,
            tolerance: 1e-10,
            max_iterations: 5,
            function_type: "polynomial".to_string(),
            coefficients: Some(vec![-2.0, 0.0, 1.0]),
        })
        .unwrap();

        assert!(!result.converged || result.iterations <= 5);
    }

    #[test]
    fn test_integration_constant_function() {
        let result = integrate(IntegrationRequest {
            method: "trapezoidal".to_string(),
            lower_bound: 0.0,
            upper_bound: 1.0,
            num_points: 10,
            function_type: "polynomial".to_string(),
            coefficients: Some(vec![1.0]), // f(x) = 1
        })
        .unwrap();

        assert!((result.integral - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_ode_zero_steps() {
        let result = solve_ode(ODESolverRequest {
            method: "euler".to_string(),
            initial_value: 1.0,
            t_start: 0.0,
            t_end: 0.0,
            step_size: 0.1,
            derivative_expression: None,
        })
        .unwrap();

        assert_eq!(result.t_values.len(), 1);
    }

    #[test]
    fn test_invalid_ode_method() {
        let result = solve_ode(ODESolverRequest {
            method: "invalid_method".to_string(),
            initial_value: 1.0,
            t_start: 0.0,
            t_end: 1.0,
            step_size: 0.1,
            derivative_expression: None,
        });

        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_root_finding_method() {
        let result = find_root(RootFindingRequest {
            method: "invalid".to_string(),
            initial_guess: 1.0,
            interval: None,
            tolerance: 1e-6,
            max_iterations: 100,
            function_type: "polynomial".to_string(),
            coefficients: None,
        });

        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_integration_method() {
        let result = integrate(IntegrationRequest {
            method: "invalid".to_string(),
            lower_bound: 0.0,
            upper_bound: 1.0,
            num_points: 10,
            function_type: "polynomial".to_string(),
            coefficients: None,
        });

        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_interpolation_method() {
        let result = interpolate(InterpolationRequest {
            method: "invalid".to_string(),
            x_values: vec![0.0, 1.0],
            y_values: vec![0.0, 1.0],
            interpolate_at: vec![0.5],
        });

        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_linear_system_method() {
        let result = solve_linear_system(LinearSystemRequest {
            method: "invalid".to_string(),
            matrix: vec![vec![1.0]],
            rhs: vec![1.0],
            tolerance: None,
            max_iterations: None,
        });

        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_differentiation_method() {
        let result = differentiate(DifferentiationRequest {
            method: "invalid".to_string(),
            x_values: vec![0.0, 1.0],
            y_values: vec![0.0, 1.0],
            order: 1,
        });

        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_pde_method() {
        let result = solve_pde(PDESolverRequest {
            method: "invalid".to_string(),
            boundary_conditions: vec![0.0, 0.0],
            initial_conditions: vec![1.0],
            spatial_steps: 5,
            time_steps: 10,
            dx: 0.1,
            dt: 0.01,
        });

        assert!(result.is_err());
    }

    #[test]
    fn test_adaptive_integration() {
        let result = integrate(IntegrationRequest {
            method: "adaptive".to_string(),
            lower_bound: 0.0,
            upper_bound: 1.0,
            num_points: 100, // Not used in adaptive, but required by struct
            function_type: "polynomial".to_string(),
            coefficients: Some(vec![0.0, 0.0, 1.0]), // x^2
        })
        .unwrap();

        // Integral of x^2 from 0 to 1 is 1/3 ≈ 0.333333
        assert!((result.integral - 0.333333).abs() < 1e-6);
        assert_eq!(result.method_used, "adaptive");
    }

    #[test]
    fn test_adaptive_integration_cubic() {
        // Test with a cubic function: x^3
        let result = integrate(IntegrationRequest {
            method: "adaptive".to_string(),
            lower_bound: 0.0,
            upper_bound: 2.0,
            num_points: 100,
            function_type: "polynomial".to_string(),
            coefficients: Some(vec![0.0, 0.0, 0.0, 1.0]), // x^3
        })
        .unwrap();

        // Integral of x^3 from 0 to 2 is 2^4/4 = 4
        assert!((result.integral - 4.0).abs() < 1e-6);
        assert_eq!(result.method_used, "adaptive");
    }
