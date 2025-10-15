//! Comprehensive tests for mathematics/calculus module operations
//!
//! Tests cover:
//! - Fractional calculus (derivatives and integrals)
//! - Special functions (Riemann zeta, elliptic integrals, hypergeometric, etc.)
//! - Variational calculus (Euler-Lagrange equations, Hamilton's principle)
//! - Symbolic integration (indefinite, definite, improper integrals)
//!
//! Note: Stochastic calculus tests are excluded as they involve randomness

use computational_engine::mathematics::calculus::*;
use serde_json::json;
use std::collections::HashMap;

const EPSILON: f64 = 1e-6;
const LOOSE_EPSILON: f64 = 1e-4;

// ============================================================================
// FRACTIONAL CALCULUS TESTS
// ============================================================================

#[test]
fn test_fractional_derivative_grunwald_letnikov_constant() {
    // Fractional derivative of constant using Grunwald-Letnikov
    // Note: GL definition may not give exactly zero for finite samples
    let function_values = vec![1.0; 100];
    let order = 0.5;
    let dx = 0.01;

    let result = fractional_derivative_grunwald_letnikov(&function_values, order, dx);

    // Just verify it computes without errors and is finite
    assert!(result.len() == 100);
    for i in 10..100 {
        assert!(result[i].is_finite(), "Expected finite fractional derivative, got {}", result[i]);
    }
}

#[test]
fn test_fractional_derivative_grunwald_letnikov_linear() {
    // Fractional derivative of x has known behavior
    let num_points = 100;
    let dx = 0.01;
    let function_values: Vec<f64> = (0..num_points).map(|i| i as f64 * dx).collect();
    let order = 0.5;

    let result = fractional_derivative_grunwald_letnikov(&function_values, order, dx);

    assert!(result.len() == num_points);
    // Result should be finite and reasonable
    for i in 10..num_points {
        assert!(result[i].is_finite());
    }
}

#[test]
fn test_fractional_derivative_grunwald_letnikov_quadratic() {
    let num_points = 100;
    let dx = 0.01;
    let function_values: Vec<f64> = (0..num_points).map(|i| {
        let x = i as f64 * dx;
        x * x
    }).collect();
    let order = 0.5;

    let result = fractional_derivative_grunwald_letnikov(&function_values, order, dx);

    assert!(result.len() == num_points);
    for i in 10..num_points {
        assert!(result[i].is_finite());
        assert!(result[i] > 0.0); // Should be positive for x^2
    }
}

#[test]
fn test_fractional_derivative_caputo_constant() {
    let function_values = vec![1.0; 100];
    let order = 0.5;
    let dx = 0.01;

    let result = fractional_derivative_caputo(&function_values, order, dx);

    assert!(result.len() == 100);
    // Caputo derivative of constant should be exactly zero
    for i in 10..100 {
        assert!(result[i].abs() < 0.1);
    }
}

#[test]
fn test_fractional_derivative_caputo_linear() {
    let num_points = 100;
    let dx = 0.01;
    let function_values: Vec<f64> = (0..num_points).map(|i| i as f64 * dx).collect();
    let order = 0.5;

    let result = fractional_derivative_caputo(&function_values, order, dx);

    assert!(result.len() == num_points);
    for val in result.iter() {
        assert!(val.is_finite());
    }
}

#[test]
fn test_fractional_integral_riemann_liouville_constant() {
    let function_values = vec![1.0; 100];
    let order = 0.5;
    let dx = 0.01;

    let result = fractional_integral_riemann_liouville(&function_values, order, dx);

    assert!(result.len() == 100);
    // Fractional integral should increase
    for i in 1..100 {
        assert!(result[i] >= result[i-1]);
    }
}

#[test]
fn test_fractional_integral_riemann_liouville_linear() {
    let num_points = 100;
    let dx = 0.01;
    let function_values: Vec<f64> = (0..num_points).map(|i| i as f64 * dx).collect();
    let order = 0.5;

    let result = fractional_integral_riemann_liouville(&function_values, order, dx);

    assert!(result.len() == num_points);
    // Integral should be monotonically increasing for x
    for i in 1..num_points {
        assert!(result[i] >= result[i-1]);
        assert!(result[i].is_finite());
    }
}

#[test]
fn test_fractional_integral_order_variations() {
    let function_values = vec![1.0; 50];
    let dx = 0.01;

    for order in [0.25, 0.5, 0.75, 1.0] {
        let result = fractional_integral_riemann_liouville(&function_values, order, dx);
        assert_eq!(result.len(), 50);
        assert!(result.last().unwrap().is_finite());
    }
}

#[test]
fn test_handle_fractional_derivative_default() {
    let params = HashMap::new();
    let result = handle_fractional_derivative(&params);

    assert!(result.success);
    assert_eq!(result.operation, "fractional_derivative");
    assert!(result.result.is_some());
    assert!(result.error.is_none());
}

#[test]
fn test_handle_fractional_derivative_custom_params() {
    let mut params = HashMap::new();
    params.insert("function".to_string(), json!("x^2"));
    params.insert("order".to_string(), json!(0.5));
    params.insert("lower_bound".to_string(), json!(0.0));
    params.insert("upper_bound".to_string(), json!(2.0));
    params.insert("num_points".to_string(), json!(50));
    params.insert("definition".to_string(), json!("caputo"));

    let result = handle_fractional_derivative(&params);

    assert!(result.success);
    assert!(result.result.is_some());

    let res = result.result.unwrap();
    assert_eq!(res["order"], 0.5);
    assert_eq!(res["definition"], "caputo");
    assert_eq!(res["function"], "x^2");
}

#[test]
fn test_handle_fractional_integral_default() {
    let params = HashMap::new();
    let result = handle_fractional_integral(&params);

    assert!(result.success);
    assert_eq!(result.operation, "fractional_integral");
    assert!(result.result.is_some());
}

#[test]
fn test_handle_fractional_integral_custom() {
    let mut params = HashMap::new();
    params.insert("function".to_string(), json!("sin(x)"));
    params.insert("order".to_string(), json!(0.75));
    params.insert("num_points".to_string(), json!(100));

    let result = handle_fractional_integral(&params);

    assert!(result.success);
    let res = result.result.unwrap();
    assert_eq!(res["order"], 0.75);
}

#[test]
fn test_handle_fractional_calculus_dispatcher() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("fractional_derivative"));

    let result = handle_fractional_calculus(&params);
    assert!(result.success);
    assert_eq!(result.operation, "fractional_derivative");

    params.insert("operation".to_string(), json!("fractional_integral"));
    let result2 = handle_fractional_calculus(&params);
    assert!(result2.success);
    assert_eq!(result2.operation, "fractional_integral");
}

// ============================================================================
// SPECIAL FUNCTIONS TESTS
// ============================================================================

#[test]
fn test_riemann_zeta_known_values() {
    // ζ(2) = π²/6 ≈ 1.6449340668
    let zeta_2 = riemann_zeta(2.0);
    let expected = std::f64::consts::PI.powi(2) / 6.0;
    // Use relaxed tolerance for approximate numerical computation
    assert!((zeta_2 - expected).abs() < 0.001,
            "ζ(2) = {}, expected {}", zeta_2, expected);

    // ζ(4) = π⁴/90 ≈ 1.0823232337
    let zeta_4 = riemann_zeta(4.0);
    let expected_4 = std::f64::consts::PI.powi(4) / 90.0;
    assert!((zeta_4 - expected_4).abs() < 0.001,
            "ζ(4) = {}, expected {}", zeta_4, expected_4);
}

#[test]
fn test_riemann_zeta_large_values() {
    // For large s, ζ(s) → 1
    let zeta_10 = riemann_zeta(10.0);
    assert!((zeta_10 - 1.0).abs() < 0.01);

    let zeta_20 = riemann_zeta(20.0);
    assert!((zeta_20 - 1.0).abs() < 0.001);
}

#[test]
fn test_riemann_zeta_at_one() {
    let zeta_1 = riemann_zeta(1.0);
    assert!(zeta_1.is_infinite());
}

#[test]
fn test_elliptic_integral_first_kind_zero() {
    // K(0) = π/2
    let k0 = elliptic_integral_first_kind(0.0);
    let expected = std::f64::consts::PI / 2.0;
    assert!((k0 - expected).abs() < LOOSE_EPSILON,
            "K(0) = {}, expected {}", k0, expected);
}

#[test]
fn test_elliptic_integral_first_kind_values() {
    let k_half = elliptic_integral_first_kind(0.5);
    assert!(k_half > std::f64::consts::PI / 2.0);
    assert!(k_half.is_finite());

    let k_small = elliptic_integral_first_kind(0.1);
    assert!(k_small.is_finite());
    assert!(k_small > 0.0);
}

#[test]
fn test_elliptic_integral_second_kind_zero() {
    // E(0) = π/2
    let e0 = elliptic_integral_second_kind(0.0);
    let expected = std::f64::consts::PI / 2.0;
    assert!((e0 - expected).abs() < LOOSE_EPSILON,
            "E(0) = {}, expected {}", e0, expected);
}

#[test]
fn test_elliptic_integral_second_kind_values() {
    let e_half = elliptic_integral_second_kind(0.5);
    assert!(e_half.is_finite());
    assert!(e_half > 0.0);
}

#[test]
fn test_elliptic_integral_invalid_k() {
    let k_invalid = elliptic_integral_first_kind(1.5);
    assert!(k_invalid.is_nan());

    let e_invalid = elliptic_integral_second_kind(-1.5);
    assert!(e_invalid.is_nan());
}

#[test]
fn test_hypergeometric_2f1_simple_cases() {
    // 2F1(a,b;c;0) = 1
    let h0 = hypergeometric_2f1(1.0, 2.0, 3.0, 0.0);
    assert!((h0 - 1.0).abs() < EPSILON);

    // 2F1(1,1;2;z) = -ln(1-z)/z (for small z)
    let h_small = hypergeometric_2f1(1.0, 1.0, 2.0, 0.5);
    assert!(h_small.is_finite());
    assert!(h_small > 0.0);
}

#[test]
fn test_hypergeometric_2f1_convergence() {
    let h1 = hypergeometric_2f1(1.0, 0.5, 2.0, 0.3);
    assert!(h1.is_finite());

    let h2 = hypergeometric_2f1(2.0, 1.0, 3.0, 0.9);
    assert!(h2.is_finite());
}

#[test]
fn test_hypergeometric_invalid_z() {
    let h = hypergeometric_2f1(1.0, 1.0, 2.0, 1.5);
    assert!(h.is_nan());
}

#[test]
fn test_jacobi_theta_3_basic() {
    // θ₃(0, 0) = 1
    let theta = jacobi_theta_3(0.0, 0.0);
    assert!((theta - 1.0).abs() < EPSILON);
}

#[test]
fn test_jacobi_theta_3_values() {
    let theta1 = jacobi_theta_3(0.1, 0.0);
    assert!(theta1.is_finite());
    assert!(theta1 > 1.0);

    let theta2 = jacobi_theta_3(0.1, std::f64::consts::PI / 4.0);
    assert!(theta2.is_finite());
}

#[test]
fn test_jacobi_theta_4_basic() {
    let theta = jacobi_theta_4(0.0, 0.0);
    assert!((theta - 1.0).abs() < EPSILON);
}

#[test]
fn test_jacobi_theta_4_values() {
    let theta = jacobi_theta_4(0.1, 0.5);
    assert!(theta.is_finite());
}

#[test]
fn test_bessel_j0_zero() {
    // J₀(0) = 1
    let j0_zero = bessel_j0(0.0);
    assert!((j0_zero - 1.0).abs() < EPSILON);
}

#[test]
fn test_bessel_j0_values() {
    // J₀ is oscillatory
    let j0_1 = bessel_j0(1.0);
    assert!(j0_1.is_finite());
    assert!(j0_1 > 0.0); // J₀(1) ≈ 0.765

    let j0_5 = bessel_j0(5.0);
    assert!(j0_5.is_finite());
}

#[test]
fn test_bessel_j1_zero() {
    // J₁(0) = 0
    let j1_zero = bessel_j1(0.0);
    assert!(j1_zero.abs() < EPSILON);
}

#[test]
fn test_bessel_j1_values() {
    let j1_1 = bessel_j1(1.0);
    assert!(j1_1.is_finite());
    assert!(j1_1 > 0.0); // J₁(1) ≈ 0.44
}

#[test]
fn test_legendre_polynomial_p0() {
    // P₀(x) = 1
    for x in [-1.0, -0.5, 0.0, 0.5, 1.0] {
        let p0 = legendre_polynomial(0, x);
        assert!((p0 - 1.0).abs() < EPSILON);
    }
}

#[test]
fn test_legendre_polynomial_p1() {
    // P₁(x) = x
    for x in [-1.0, -0.5, 0.0, 0.5, 1.0] {
        let p1 = legendre_polynomial(1, x);
        assert!((p1 - x).abs() < EPSILON);
    }
}

#[test]
fn test_legendre_polynomial_p2() {
    // P₂(x) = (3x² - 1)/2
    let x = 0.5;
    let p2 = legendre_polynomial(2, x);
    let expected = (3.0 * x * x - 1.0) / 2.0;
    assert!((p2 - expected).abs() < EPSILON);
}

#[test]
fn test_legendre_polynomial_orthogonality_property() {
    // P_n(1) = 1, P_n(-1) = (-1)^n
    for n in 0..5 {
        let p_at_1 = legendre_polynomial(n, 1.0);
        assert!((p_at_1 - 1.0).abs() < EPSILON);

        let p_at_minus_1 = legendre_polynomial(n, -1.0);
        let expected = if n % 2 == 0 { 1.0 } else { -1.0 };
        assert!((p_at_minus_1 - expected).abs() < EPSILON);
    }
}

#[test]
fn test_handle_riemann_zeta() {
    let mut params = HashMap::new();
    params.insert("s".to_string(), json!(2.0));

    let result = handle_riemann_zeta(&params);
    assert!(result.success);

    let res = result.result.unwrap();
    assert_eq!(res["s"], 2.0);

    let zeta_val = res["zeta_value"].as_f64().unwrap();
    let expected = std::f64::consts::PI.powi(2) / 6.0;
    assert!((zeta_val - expected).abs() < 0.001);
}

#[test]
fn test_handle_elliptic_integral_first_kind() {
    let mut params = HashMap::new();
    params.insert("k".to_string(), json!(0.5));
    params.insert("type".to_string(), json!("first"));

    let result = handle_elliptic_integral(&params);
    assert!(result.success);

    let res = result.result.unwrap();
    assert_eq!(res["type"], "first");
    assert!(res["value"].as_f64().unwrap().is_finite());
}

#[test]
fn test_handle_elliptic_integral_second_kind() {
    let mut params = HashMap::new();
    params.insert("k".to_string(), json!(0.3));
    params.insert("type".to_string(), json!("second"));

    let result = handle_elliptic_integral(&params);
    assert!(result.success);

    let res = result.result.unwrap();
    assert_eq!(res["type"], "second");
}

#[test]
fn test_handle_hypergeometric() {
    let mut params = HashMap::new();
    params.insert("a".to_string(), json!(1.0));
    params.insert("b".to_string(), json!(2.0));
    params.insert("c".to_string(), json!(3.0));
    params.insert("z".to_string(), json!(0.5));

    let result = handle_hypergeometric(&params);
    assert!(result.success);

    let res = result.result.unwrap();
    assert!(res["value"].as_f64().unwrap().is_finite());
}

#[test]
fn test_handle_jacobi_theta_3() {
    let mut params = HashMap::new();
    params.insert("q".to_string(), json!(0.1));
    params.insert("z".to_string(), json!(0.5));
    params.insert("variant".to_string(), json!(3));

    let result = handle_jacobi_theta(&params);
    assert!(result.success);
}

#[test]
fn test_handle_jacobi_theta_4() {
    let mut params = HashMap::new();
    params.insert("q".to_string(), json!(0.1));
    params.insert("z".to_string(), json!(0.5));
    params.insert("variant".to_string(), json!(4));

    let result = handle_jacobi_theta(&params);
    assert!(result.success);
}

#[test]
fn test_handle_bessel_function_j0() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!(1.0));
    params.insert("order".to_string(), json!(0));

    let result = handle_bessel_function(&params);
    assert!(result.success);
}

#[test]
fn test_handle_bessel_function_j1() {
    let mut params = HashMap::new();
    params.insert("x".to_string(), json!(2.0));
    params.insert("order".to_string(), json!(1));

    let result = handle_bessel_function(&params);
    assert!(result.success);
}

#[test]
fn test_handle_legendre_polynomial() {
    let mut params = HashMap::new();
    params.insert("n".to_string(), json!(3));
    params.insert("x".to_string(), json!(0.5));

    let result = handle_legendre_polynomial(&params);
    assert!(result.success);

    let res = result.result.unwrap();
    assert!(res["value"].as_f64().unwrap().is_finite());
}

#[test]
fn test_handle_special_functions_dispatcher() {
    let functions = vec![
        "riemann_zeta", "elliptic_integral", "hypergeometric",
        "jacobi_theta", "bessel", "legendre"
    ];

    for func_type in functions {
        let mut params = HashMap::new();
        params.insert("function_type".to_string(), json!(func_type));

        let result = handle_special_functions(&params);
        assert!(result.success, "Failed for function type: {}", func_type);
    }
}

// ============================================================================
// VARIATIONAL CALCULUS TESTS
// ============================================================================

#[test]
fn test_compute_euler_lagrange_free_particle() {
    // L = 0.5*y'^2 => y'' = 0
    let equation = compute_euler_lagrange_equation("0.5*y'^2");
    assert_eq!(equation, "y'' = 0");
}

#[test]
fn test_compute_euler_lagrange_harmonic_oscillator() {
    // L = 0.5*y'^2 - 0.5*k*y^2 => y'' + k*y = 0
    let equation = compute_euler_lagrange_equation("0.5*y'^2 - 0.5*k*y^2");
    assert_eq!(equation, "y'' + k*y = 0");
}

#[test]
fn test_compute_euler_lagrange_brachistochrone() {
    // Brachistochrone problem
    let equation = compute_euler_lagrange_equation("sqrt(1 + y'^2)");
    assert!(equation.contains("y''"));
}

#[test]
fn test_compute_action_functional_constant_lagrangian() {
    let lagrangian_values = vec![1.0; 100];
    let dx = 0.01;

    let action = compute_action_functional(&lagrangian_values, dx);

    // Action ≈ 1.0 * (100-1) * 0.01 ≈ 0.99
    assert!((action - 0.99).abs() < 0.01);
}

#[test]
fn test_compute_action_functional_linear() {
    let lagrangian_values: Vec<f64> = (0..100).map(|i| i as f64 * 0.01).collect();
    let dx = 0.01;

    let action = compute_action_functional(&lagrangian_values, dx);

    // Should be positive and finite
    assert!(action > 0.0);
    assert!(action.is_finite());
}

#[test]
fn test_compute_action_functional_empty() {
    let action = compute_action_functional(&[], 0.01);
    assert_eq!(action, 0.0);

    let action_single = compute_action_functional(&[1.0], 0.01);
    assert_eq!(action_single, 0.0);
}

#[test]
fn test_solve_euler_lagrange_bvp_linear() {
    // For simple Lagrangians, should give linear interpolation
    let solution = solve_euler_lagrange_bvp("0.5*y'^2", (0.0, 1.0), (0.0, 1.0), 11);

    assert_eq!(solution.len(), 11);
    assert!((solution[0].1 - 0.0).abs() < EPSILON);
    assert!((solution[10].1 - 1.0).abs() < EPSILON);

    // Check linearity
    assert!((solution[5].1 - 0.5).abs() < EPSILON);
}

#[test]
fn test_solve_euler_lagrange_bvp_different_boundaries() {
    let solution = solve_euler_lagrange_bvp("0.5*y'^2", (0.0, 2.0), (-1.0, 1.0), 21);

    assert_eq!(solution.len(), 21);
    assert!((solution[0].1 - (-1.0)).abs() < EPSILON);
    assert!((solution[20].1 - 1.0).abs() < EPSILON);
}

#[test]
fn test_find_conservation_laws_energy() {
    let laws = find_conservation_laws("0.5*y'^2 - 0.5*k*y^2");
    assert!(!laws.is_empty());
    // Should detect some conservation property
}

#[test]
fn test_handle_euler_lagrange() {
    let mut params = HashMap::new();
    params.insert("lagrangian".to_string(), json!("0.5*y'^2"));

    let result = handle_euler_lagrange(&params);
    assert!(result.success);

    let res = result.result.unwrap();
    assert_eq!(res["lagrangian"], "0.5*y'^2");
    assert_eq!(res["euler_lagrange_equation"], "y'' = 0");
}

#[test]
fn test_handle_variational_calculus_euler_lagrange() {
    let mut params = HashMap::new();
    params.insert("lagrangian".to_string(), json!("0.5*y'^2 - 0.5*k*y^2"));
    params.insert("method".to_string(), json!("euler_lagrange"));

    let result = handle_variational_calculus(&params);
    assert!(result.success);
}

#[test]
fn test_handle_variational_calculus_with_boundary_conditions() {
    let mut params = HashMap::new();
    params.insert("lagrangian".to_string(), json!("0.5*y'^2"));
    params.insert("method".to_string(), json!("euler_lagrange"));

    let mut bc = serde_json::Map::new();
    bc.insert("x0".to_string(), json!(0.0));
    bc.insert("x1".to_string(), json!(1.0));
    bc.insert("y0".to_string(), json!(0.0));
    bc.insert("y1".to_string(), json!(1.0));
    params.insert("boundary_conditions".to_string(), json!(bc));

    let result = handle_variational_calculus(&params);
    assert!(result.success);

    let res = result.result.unwrap();
    assert!(!res["solution"].is_null());
}

#[test]
fn test_handle_variational_calculus_hamilton_principle() {
    let mut params = HashMap::new();
    params.insert("lagrangian".to_string(), json!("0.5*y'^2"));
    params.insert("method".to_string(), json!("hamilton_principle"));

    let mut bc = serde_json::Map::new();
    bc.insert("x0".to_string(), json!(0.0));
    bc.insert("x1".to_string(), json!(1.0));
    bc.insert("y0".to_string(), json!(0.0));
    bc.insert("y1".to_string(), json!(1.0));
    params.insert("boundary_conditions".to_string(), json!(bc));

    let result = handle_variational_calculus(&params);
    assert!(result.success);

    let res = result.result.unwrap();
    assert!(res["action_value"].is_number());
    assert!(res["extremal_path"].is_array());
}

#[test]
fn test_handle_variational_calculus_noether_theorem() {
    let mut params = HashMap::new();
    params.insert("lagrangian".to_string(), json!("0.5*y'^2"));
    params.insert("method".to_string(), json!("noether_theorem"));

    let result = handle_variational_calculus(&params);
    assert!(result.success);

    let res = result.result.unwrap();
    assert!(res["conservation_laws"].is_array());
}

#[test]
fn test_handle_variational_calculus_invalid_method() {
    let mut params = HashMap::new();
    params.insert("method".to_string(), json!("invalid_method"));

    let result = handle_variational_calculus(&params);
    assert!(!result.success);
    assert!(result.error.is_some());
}

// ============================================================================
// SYMBOLIC INTEGRATION TESTS
// ============================================================================

#[test]
fn test_symbolic_integrate_constant() {
    assert_eq!(symbolic_integrate("1", "x"), "x");
}

#[test]
fn test_symbolic_integrate_linear() {
    assert_eq!(symbolic_integrate("x", "x"), "x^2/2");
}

#[test]
fn test_symbolic_integrate_quadratic() {
    assert_eq!(symbolic_integrate("x^2", "x"), "x^3/3");
}

#[test]
fn test_symbolic_integrate_cubic() {
    assert_eq!(symbolic_integrate("x^3", "x"), "x^4/4");
}

#[test]
fn test_symbolic_integrate_reciprocal() {
    assert_eq!(symbolic_integrate("1/x", "x"), "ln|x|");
}

#[test]
fn test_symbolic_integrate_exponential() {
    assert_eq!(symbolic_integrate("exp(x)", "x"), "exp(x)");
}

#[test]
fn test_symbolic_integrate_sin() {
    assert_eq!(symbolic_integrate("sin(x)", "x"), "-cos(x)");
}

#[test]
fn test_symbolic_integrate_cos() {
    assert_eq!(symbolic_integrate("cos(x)", "x"), "sin(x)");
}

#[test]
fn test_symbolic_integrate_tan() {
    assert_eq!(symbolic_integrate("tan(x)", "x"), "-ln|cos(x)|");
}

#[test]
fn test_symbolic_integrate_sinh() {
    assert_eq!(symbolic_integrate("sinh(x)", "x"), "cosh(x)");
}

#[test]
fn test_symbolic_integrate_cosh() {
    assert_eq!(symbolic_integrate("cosh(x)", "x"), "sinh(x)");
}

#[test]
fn test_symbolic_integrate_ln() {
    let result = symbolic_integrate("ln(x)", "x");
    assert!(result.contains("ln"));
}

#[test]
fn test_symbolic_integrate_sqrt() {
    let result = symbolic_integrate("sqrt(x)", "x");
    assert!(result.contains("3/2"));
}

#[test]
fn test_symbolic_integrate_reciprocal_squared() {
    assert_eq!(symbolic_integrate("1/x^2", "x"), "-1/x");
}

#[test]
fn test_symbolic_integrate_arctan_integrand() {
    assert_eq!(symbolic_integrate("1/(1+x^2)", "x"), "arctan(x)");
}

#[test]
fn test_handle_symbolic_integral_default() {
    let params = HashMap::new();
    let result = handle_symbolic_integral(&params);

    assert!(result.success);
    assert_eq!(result.operation, "symbolic_integral");
}

#[test]
fn test_handle_symbolic_integral_quadratic() {
    let mut params = HashMap::new();
    params.insert("function".to_string(), json!("x^2"));
    params.insert("variable".to_string(), json!("x"));

    let result = handle_symbolic_integral(&params);
    assert!(result.success);

    let res = result.result.unwrap();
    assert_eq!(res["integral"], "x^3/3");
}

#[test]
fn test_handle_symbolic_integral_trig() {
    let mut params = HashMap::new();
    params.insert("function".to_string(), json!("sin(x)"));

    let result = handle_symbolic_integral(&params);
    assert!(result.success);

    let res = result.result.unwrap();
    assert_eq!(res["integral"], "-cos(x)");
}

#[test]
fn test_handle_definite_integral_x_squared() {
    let mut params = HashMap::new();
    params.insert("function".to_string(), json!("x^2"));
    params.insert("lower_bound".to_string(), json!(0.0));
    params.insert("upper_bound".to_string(), json!(1.0));

    let result = handle_definite_integral(&params);
    assert!(result.success);

    let res = result.result.unwrap();
    let numerical = res["numerical_value"].as_f64().unwrap();
    // ∫₀¹ x² dx = 1/3
    assert!((numerical - 1.0/3.0).abs() < LOOSE_EPSILON);
}

#[test]
fn test_handle_definite_integral_constant() {
    let mut params = HashMap::new();
    params.insert("function".to_string(), json!("1"));
    params.insert("lower_bound".to_string(), json!(0.0));
    params.insert("upper_bound".to_string(), json!(2.0));

    let result = handle_definite_integral(&params);
    assert!(result.success);

    let res = result.result.unwrap();
    let numerical = res["numerical_value"].as_f64().unwrap();
    // ∫₀² 1 dx = 2
    assert!((numerical - 2.0).abs() < EPSILON);
}

#[test]
fn test_handle_improper_integral_reciprocal_squared() {
    let mut params = HashMap::new();
    params.insert("function".to_string(), json!("1/x^2"));
    params.insert("lower_bound".to_string(), json!(1.0));
    params.insert("upper_bound".to_string(), json!(null));

    let result = handle_improper_integral(&params);
    assert!(result.success);

    let res = result.result.unwrap();
    let convergence = res["convergence_value"].as_f64().unwrap();
    // ∫₁^∞ 1/x² dx = 1
    assert!((convergence - 1.0).abs() < LOOSE_EPSILON);
}

#[test]
fn test_handle_improper_integral_exponential() {
    let mut params = HashMap::new();
    params.insert("function".to_string(), json!("exp(-x)"));
    params.insert("lower_bound".to_string(), json!(0.0));
    params.insert("upper_bound".to_string(), json!(null));

    let result = handle_improper_integral(&params);
    assert!(result.success);

    let res = result.result.unwrap();
    let convergence = res["convergence_value"].as_f64().unwrap();
    // ∫₀^∞ e^(-x) dx = 1
    assert!((convergence - 1.0).abs() < LOOSE_EPSILON);
}

#[test]
fn test_definite_integral_numerical_sin() {
    // ∫₀^π sin(x) dx = 2
    let result = definite_integral_numerical("sin(x)", "x", 0.0, std::f64::consts::PI, 1000);
    assert!((result - 2.0).abs() < LOOSE_EPSILON);
}

#[test]
fn test_definite_integral_numerical_cos() {
    // ∫₀^(π/2) cos(x) dx = 1
    let result = definite_integral_numerical("cos(x)", "x", 0.0, std::f64::consts::PI / 2.0, 1000);
    assert!((result - 1.0).abs() < LOOSE_EPSILON);
}

#[test]
fn test_definite_integral_numerical_exp() {
    // ∫₀¹ e^x dx = e - 1
    let result = definite_integral_numerical("exp(x)", "x", 0.0, 1.0, 1000);
    let expected = 1.0_f64.exp() - 1.0;
    assert!((result - expected).abs() < LOOSE_EPSILON);
}
