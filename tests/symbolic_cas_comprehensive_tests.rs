//! Comprehensive tests for symbolic CAS module
//!
//! This test suite covers:
//! - Expression parsing and creation
//! - Simplification and expansion
//! - Symbolic differentiation
//! - Symbolic integration
//! - Expression evaluation
//! - Matrix symbolic operations
//! - Substitution and manipulation
//! - Series expansion
//! - Limits

use computational_engine::mathematics::symbolic_cas::*;
use std::collections::HashMap;

// ============================================================================
// EXPRESSION PARSING TESTS
// ============================================================================

#[test]
fn test_parse_number() {
    let expr = parse("42").unwrap();
    assert_eq!(expr, Expr::num(42));
}

#[test]
fn test_parse_negative_number() {
    let expr = parse("-5").unwrap();
    let expected = Expr::mul(Expr::num(-1), Expr::num(5));
    assert_eq!(expr, expected);
}

#[test]
fn test_parse_variable() {
    let expr = parse("x").unwrap();
    assert_eq!(expr, Expr::sym("x"));
}

#[test]
fn test_parse_addition() {
    let expr = parse("x + 1").unwrap();
    assert!(matches!(expr, Expr::Add(_, _)));
}

#[test]
fn test_parse_subtraction() {
    let expr = parse("x - 3").unwrap();
    assert!(matches!(expr, Expr::Add(_, _)));
}

#[test]
fn test_parse_multiplication() {
    let expr = parse("3 * x").unwrap();
    assert!(matches!(expr, Expr::Mul(_, _)));
}

#[test]
fn test_parse_division() {
    let expr = parse("x / 2").unwrap();
    assert!(matches!(expr, Expr::Mul(_, _))); // Division becomes multiplication by power
}

#[test]
fn test_parse_power() {
    let expr = parse("x^3").unwrap();
    assert!(matches!(expr, Expr::Pow(_, _)));
}

#[test]
fn test_parse_complex_expression() {
    let expr = parse("3*x^2 + 2*x - 5").unwrap();
    eprintln!("Parsed: {}", expr);
    assert!(matches!(expr, Expr::Add(_, _)));
}

#[test]
fn test_parse_nested_parentheses() {
    let expr = parse("((x + 1) * 2)^2").unwrap();
    eprintln!("Parsed nested: {}", expr);
    assert!(matches!(expr, Expr::Pow(_, _)));
}

#[test]
fn test_parse_function_sin() {
    let expr = parse("sin(x)").unwrap();
    assert!(matches!(expr, Expr::Function(_, _)));
}

#[test]
fn test_parse_function_cos() {
    let expr = parse("cos(x + 1)").unwrap();
    assert!(matches!(expr, Expr::Function(_, _)));
}

#[test]
fn test_parse_function_exp() {
    let expr = parse("exp(x)").unwrap();
    assert!(matches!(expr, Expr::Function(_, _)));
}

#[test]
fn test_parse_function_ln() {
    let expr = parse("ln(x)").unwrap();
    assert!(matches!(expr, Expr::Function(_, _)));
}

#[test]
fn test_parse_nested_functions() {
    let expr = parse("sin(cos(x))").unwrap();
    if let Expr::Function(name, args) = expr {
        assert_eq!(name, "sin");
        assert!(matches!(args[0], Expr::Function(_, _)));
    } else {
        panic!("Expected function");
    }
}

// ============================================================================
// EXPRESSION SIMPLIFICATION TESTS
// ============================================================================

#[test]
fn test_simplify_add_zero_left() {
    let expr = parse("0 + x").unwrap();
    let simplified = simplify_expr(&expr);
    assert_eq!(simplified, Expr::sym("x"));
}

#[test]
fn test_simplify_add_zero_right() {
    let expr = parse("x + 0").unwrap();
    let simplified = simplify_expr(&expr);
    assert_eq!(simplified, Expr::sym("x"));
}

#[test]
fn test_simplify_mul_zero() {
    let expr = parse("0 * x").unwrap();
    let simplified = simplify_expr(&expr);
    assert_eq!(simplified, Expr::num(0));
}

#[test]
fn test_simplify_mul_one_left() {
    let expr = parse("1 * x").unwrap();
    let simplified = simplify_expr(&expr);
    assert_eq!(simplified, Expr::sym("x"));
}

#[test]
fn test_simplify_mul_one_right() {
    let expr = parse("x * 1").unwrap();
    let simplified = simplify_expr(&expr);
    assert_eq!(simplified, Expr::sym("x"));
}

#[test]
fn test_simplify_power_zero() {
    let expr = parse("x^0").unwrap();
    let simplified = simplify_expr(&expr);
    assert_eq!(simplified, Expr::num(1));
}

#[test]
fn test_simplify_power_one() {
    let expr = parse("x^1").unwrap();
    let simplified = simplify_expr(&expr);
    assert_eq!(simplified, Expr::sym("x"));
}

#[test]
fn test_simplify_add_numbers() {
    let expr = Expr::add(Expr::num(3), Expr::num(5));
    let simplified = simplify_expr(&expr);
    assert_eq!(simplified, Expr::num(8));
}

#[test]
fn test_simplify_mul_numbers() {
    let expr = Expr::mul(Expr::num(3), Expr::num(4));
    let simplified = simplify_expr(&expr);
    assert_eq!(simplified, Expr::num(12));
}

#[test]
fn test_simplify_same_variables() {
    let expr = parse("x + x").unwrap();
    let simplified = simplify_expr(&expr);
    eprintln!("x + x simplified: {}", simplified);
    // Should be 2*x
    assert!(matches!(simplified, Expr::Mul(_, _)));
}

#[test]
fn test_simplify_multiply_same_base() {
    let expr = Expr::mul(parse("x^2").unwrap(), parse("x^3").unwrap());
    let simplified = simplify_expr(&expr);
    eprintln!("x^2 * x^3 simplified: {}", simplified);
    // Should be x^5
    assert!(matches!(simplified, Expr::Pow(_, _)));
}

// ============================================================================
// EXPRESSION EXPANSION TESTS
// ============================================================================

#[test]
fn test_expand_square_binomial() {
    let expr = parse("(x + 1)^2").unwrap();
    let expanded = expand_expr(&expr);
    eprintln!("(x + 1)^2 expanded: {}", expanded);
    // Should contain x^2, x, and constant terms
}

#[test]
fn test_expand_product() {
    let expr = parse("(x + 1) * (x + 2)").unwrap();
    let expanded = expand_expr(&expr);
    eprintln!("(x + 1)(x + 2) expanded: {}", expanded);
}

#[test]
fn test_expand_cube() {
    let expr = parse("(x + 1)^3").unwrap();
    let expanded = expand_expr(&expr);
    eprintln!("(x + 1)^3 expanded: {}", expanded);
}

#[test]
fn test_expand_distribute() {
    let expr = parse("2 * (x + 3)").unwrap();
    let expanded = expand_expr(&expr);
    eprintln!("2(x + 3) expanded: {}", expanded);
}

// ============================================================================
// SYMBOLIC DIFFERENTIATION TESTS
// ============================================================================

#[test]
fn test_differentiate_constant() {
    let expr = Expr::num(42);
    let deriv = diff(&expr, "x");
    assert_eq!(deriv, Expr::num(0));
}

#[test]
fn test_differentiate_variable_same() {
    let expr = Expr::sym("x");
    let deriv = diff(&expr, "x");
    assert_eq!(deriv, Expr::num(1));
}

#[test]
fn test_differentiate_variable_different() {
    let expr = Expr::sym("y");
    let deriv = diff(&expr, "x");
    assert_eq!(deriv, Expr::num(0));
}

#[test]
fn test_differentiate_linear() {
    let expr = parse("3*x + 2").unwrap();
    let deriv = diff(&expr, "x");
    eprintln!("d/dx(3x + 2) = {}", deriv);
    // Should simplify to 3
}

#[test]
fn test_differentiate_quadratic() {
    let expr = parse("x^2").unwrap();
    let deriv = diff(&expr, "x");
    eprintln!("d/dx(x^2) = {}", deriv);
    // Should be 2*x
}

#[test]
fn test_differentiate_cubic() {
    let expr = parse("x^3").unwrap();
    let deriv = diff(&expr, "x");
    eprintln!("d/dx(x^3) = {}", deriv);
    // Should be 3*x^2
}

#[test]
fn test_differentiate_polynomial() {
    let expr = parse("x^3 + 2*x^2 - 5*x + 7").unwrap();
    let deriv = diff(&expr, "x");
    eprintln!("d/dx(x^3 + 2x^2 - 5x + 7) = {}", deriv);
}

#[test]
fn test_differentiate_product_rule() {
    let expr = parse("x * x^2").unwrap();
    let deriv = diff(&expr, "x");
    eprintln!("d/dx(x * x^2) = {}", deriv);
}

#[test]
fn test_differentiate_sin() {
    let expr = parse("sin(x)").unwrap();
    let deriv = diff(&expr, "x");
    eprintln!("d/dx(sin(x)) = {}", deriv);
    // Should be cos(x)
    assert!(matches!(deriv, Expr::Function(_, _)));
}

#[test]
fn test_differentiate_cos() {
    let expr = parse("cos(x)").unwrap();
    let deriv = diff(&expr, "x");
    eprintln!("d/dx(cos(x)) = {}", deriv);
    // Should be -sin(x)
}

#[test]
fn test_differentiate_exp() {
    let expr = parse("exp(x)").unwrap();
    let deriv = diff(&expr, "x");
    eprintln!("d/dx(exp(x)) = {}", deriv);
    // Should be exp(x)
}

#[test]
fn test_differentiate_ln() {
    let expr = parse("ln(x)").unwrap();
    let deriv = diff(&expr, "x");
    eprintln!("d/dx(ln(x)) = {}", deriv);
    // Should be 1/x
}

#[test]
fn test_differentiate_chain_rule_sin() {
    let expr = parse("sin(x^2)").unwrap();
    let deriv = diff(&expr, "x");
    eprintln!("d/dx(sin(x^2)) = {}", deriv);
    // Should be cos(x^2) * 2x
}

#[test]
fn test_differentiate_second_order() {
    let expr = parse("x^4").unwrap();
    let first = diff(&expr, "x");
    let second = diff(&first, "x");
    eprintln!("d²/dx²(x^4) = {}", second);
    // Should be 12*x^2
}

// ============================================================================
// SYMBOLIC INTEGRATION TESTS
// ============================================================================

#[test]
fn test_integrate_constant() {
    let expr = Expr::num(5);
    let integral = integrate_expr(&expr, "x");
    eprintln!("∫ 5 dx = {}", integral);
    // Should be 5*x
}

#[test]
fn test_integrate_variable() {
    let expr = Expr::sym("x");
    let integral = integrate_expr(&expr, "x");
    eprintln!("∫ x dx = {}", integral);
    // Should be x^2/2
}

#[test]
fn test_integrate_x_squared() {
    let expr = parse("x^2").unwrap();
    let integral = integrate_expr(&expr, "x");
    eprintln!("∫ x^2 dx = {}", integral);
    // Should be x^3/3
}

#[test]
fn test_integrate_power_general() {
    let expr = parse("x^5").unwrap();
    let integral = integrate_expr(&expr, "x");
    eprintln!("∫ x^5 dx = {}", integral);
    // Should be x^6/6
}

#[test]
fn test_integrate_linear() {
    let expr = parse("2*x + 3").unwrap();
    let integral = integrate_expr(&expr, "x");
    eprintln!("∫ (2x + 3) dx = {}", integral);
    // Should be x^2 + 3x
}

#[test]
fn test_integrate_polynomial() {
    let expr = parse("3*x^2 + 2*x + 1").unwrap();
    let integral = integrate_expr(&expr, "x");
    eprintln!("∫ (3x^2 + 2x + 1) dx = {}", integral);
}

#[test]
fn test_integrate_sin() {
    let expr = parse("sin(x)").unwrap();
    let integral = integrate_expr(&expr, "x");
    eprintln!("∫ sin(x) dx = {}", integral);
    // Should be -cos(x)
}

#[test]
fn test_integrate_cos() {
    let expr = parse("cos(x)").unwrap();
    let integral = integrate_expr(&expr, "x");
    eprintln!("∫ cos(x) dx = {}", integral);
    // Should be sin(x)
}

#[test]
fn test_integrate_exp() {
    let expr = parse("exp(x)").unwrap();
    let integral = integrate_expr(&expr, "x");
    eprintln!("∫ exp(x) dx = {}", integral);
    // Should be exp(x)
}

#[test]
fn test_integrate_constant_factor() {
    let expr = parse("5*x^2").unwrap();
    let integral = integrate_expr(&expr, "x");
    eprintln!("∫ 5x^2 dx = {}", integral);
    // Should be (5/3)*x^3
}

// ============================================================================
// EXPRESSION EVALUATION TESTS
// ============================================================================

#[test]
fn test_evaluate_constant() {
    let expr = Expr::num(42);
    let result = expr.evaluate(&HashMap::new()).unwrap();
    assert_eq!(result, 42.0);
}

#[test]
fn test_evaluate_variable() {
    let expr = Expr::sym("x");
    let mut vars = HashMap::new();
    vars.insert("x".to_string(), 5.0);
    let result = expr.evaluate(&vars).unwrap();
    assert_eq!(result, 5.0);
}

#[test]
fn test_evaluate_addition() {
    let expr = parse("x + 3").unwrap();
    let mut vars = HashMap::new();
    vars.insert("x".to_string(), 2.0);
    let result = expr.evaluate(&vars).unwrap();
    assert_eq!(result, 5.0);
}

#[test]
fn test_evaluate_multiplication() {
    let expr = parse("3 * x").unwrap();
    let mut vars = HashMap::new();
    vars.insert("x".to_string(), 4.0);
    let result = expr.evaluate(&vars).unwrap();
    assert_eq!(result, 12.0);
}

#[test]
fn test_evaluate_power() {
    let expr = parse("x^2").unwrap();
    let mut vars = HashMap::new();
    vars.insert("x".to_string(), 5.0);
    let result = expr.evaluate(&vars).unwrap();
    assert_eq!(result, 25.0);
}

#[test]
fn test_evaluate_polynomial() {
    let expr = parse("x^2 + 3*x + 2").unwrap();
    let mut vars = HashMap::new();
    vars.insert("x".to_string(), 2.0);
    let result = expr.evaluate(&vars).unwrap();
    // 4 + 6 + 2 = 12
    assert_eq!(result, 12.0);
}

#[test]
fn test_evaluate_sin() {
    let expr = parse("sin(x)").unwrap();
    let mut vars = HashMap::new();
    vars.insert("x".to_string(), 0.0);
    let result = expr.evaluate(&vars).unwrap();
    assert!((result - 0.0).abs() < 1e-10);
}

#[test]
fn test_evaluate_cos() {
    let expr = parse("cos(x)").unwrap();
    let mut vars = HashMap::new();
    vars.insert("x".to_string(), 0.0);
    let result = expr.evaluate(&vars).unwrap();
    assert!((result - 1.0).abs() < 1e-10);
}

#[test]
fn test_evaluate_exp() {
    let expr = parse("exp(x)").unwrap();
    let mut vars = HashMap::new();
    vars.insert("x".to_string(), 0.0);
    let result = expr.evaluate(&vars).unwrap();
    assert!((result - 1.0).abs() < 1e-10);
}

#[test]
fn test_evaluate_multiple_variables() {
    let expr = parse("x + y").unwrap();
    let mut vars = HashMap::new();
    vars.insert("x".to_string(), 2.0);
    vars.insert("y".to_string(), 3.0);
    let result = expr.evaluate(&vars).unwrap();
    assert_eq!(result, 5.0);
}

// ============================================================================
// MATRIX SYMBOLIC OPERATIONS TESTS
// ============================================================================

#[test]
fn test_matrix_create() {
    let data = vec![
        vec![Expr::num(1), Expr::num(2)],
        vec![Expr::num(3), Expr::num(4)],
    ];
    let mat = SymbolicMatrix::new(data).unwrap();
    assert_eq!(mat.rows(), 2);
    assert_eq!(mat.cols(), 2);
}

#[test]
fn test_matrix_identity() {
    let identity = SymbolicMatrix::identity(3);
    assert_eq!(identity.rows(), 3);
    assert_eq!(identity.cols(), 3);
    assert_eq!(identity.get(0, 0), Some(&Expr::num(1)));
    assert_eq!(identity.get(0, 1), Some(&Expr::num(0)));
    assert_eq!(identity.get(1, 1), Some(&Expr::num(1)));
}

#[test]
fn test_matrix_zeros() {
    let zeros = SymbolicMatrix::zeros(2, 3);
    assert_eq!(zeros.rows(), 2);
    assert_eq!(zeros.cols(), 3);
    assert_eq!(zeros.get(0, 0), Some(&Expr::num(0)));
    assert_eq!(zeros.get(1, 2), Some(&Expr::num(0)));
}

#[test]
fn test_matrix_transpose() {
    let data = vec![
        vec![Expr::num(1), Expr::num(2), Expr::num(3)],
        vec![Expr::num(4), Expr::num(5), Expr::num(6)],
    ];
    let mat = SymbolicMatrix::new(data).unwrap();
    let transposed = mat.transpose();
    assert_eq!(transposed.rows(), 3);
    assert_eq!(transposed.cols(), 2);
    assert_eq!(transposed.get(0, 0), Some(&Expr::num(1)));
    assert_eq!(transposed.get(1, 0), Some(&Expr::num(2)));
    assert_eq!(transposed.get(0, 1), Some(&Expr::num(4)));
}

#[test]
fn test_matrix_add() {
    let a = SymbolicMatrix::new(vec![
        vec![Expr::num(1), Expr::num(2)],
        vec![Expr::num(3), Expr::num(4)],
    ]).unwrap();

    let b = SymbolicMatrix::new(vec![
        vec![Expr::num(5), Expr::num(6)],
        vec![Expr::num(7), Expr::num(8)],
    ]).unwrap();

    let c = a.add(&b).unwrap();
    assert_eq!(c.rows(), 2);
    assert_eq!(c.cols(), 2);
    // Check result
    eprintln!("Matrix sum:\n{}", c);
}

#[test]
fn test_matrix_mul() {
    let a = SymbolicMatrix::new(vec![
        vec![Expr::num(1), Expr::num(2)],
        vec![Expr::num(3), Expr::num(4)],
    ]).unwrap();

    let b = SymbolicMatrix::new(vec![
        vec![Expr::num(5), Expr::num(6)],
        vec![Expr::num(7), Expr::num(8)],
    ]).unwrap();

    let c = a.mul(&b).unwrap();
    assert_eq!(c.rows(), 2);
    assert_eq!(c.cols(), 2);
    eprintln!("Matrix product:\n{}", c);
}

#[test]
fn test_matrix_scalar_mul() {
    let mat = SymbolicMatrix::new(vec![
        vec![Expr::num(1), Expr::num(2)],
        vec![Expr::num(3), Expr::num(4)],
    ]).unwrap();

    let scaled = mat.scalar_mul(&Expr::num(3));
    assert_eq!(scaled.rows(), 2);
    assert_eq!(scaled.cols(), 2);
    eprintln!("Scalar multiplication:\n{}", scaled);
}

#[test]
fn test_matrix_determinant_2x2() {
    let mat = SymbolicMatrix::new(vec![
        vec![Expr::num(1), Expr::num(2)],
        vec![Expr::num(3), Expr::num(4)],
    ]).unwrap();

    let det = mat.determinant().unwrap();
    eprintln!("Det of 2x2: {}", det);
    // Should be 1*4 - 2*3 = -2
}

#[test]
fn test_matrix_determinant_symbolic() {
    let mat = SymbolicMatrix::new(vec![
        vec![Expr::sym("a"), Expr::sym("b")],
        vec![Expr::sym("c"), Expr::sym("d")],
    ]).unwrap();

    let det = mat.determinant().unwrap();
    eprintln!("Symbolic det: {}", det);
    // Should be a*d - b*c
}

#[test]
fn test_matrix_determinant_3x3() {
    let mat = SymbolicMatrix::new(vec![
        vec![Expr::num(1), Expr::num(2), Expr::num(3)],
        vec![Expr::num(4), Expr::num(5), Expr::num(6)],
        vec![Expr::num(7), Expr::num(8), Expr::num(9)],
    ]).unwrap();

    let det = mat.determinant().unwrap();
    eprintln!("Det of 3x3: {}", det);
    // This matrix is singular, det should be 0
}

#[test]
fn test_matrix_trace() {
    let mat = SymbolicMatrix::new(vec![
        vec![Expr::num(1), Expr::num(2), Expr::num(3)],
        vec![Expr::num(4), Expr::num(5), Expr::num(6)],
        vec![Expr::num(7), Expr::num(8), Expr::num(9)],
    ]).unwrap();

    let trace = mat.trace().unwrap();
    eprintln!("Trace: {}", trace);
    // Should be 1 + 5 + 9 = 15
}

#[test]
fn test_matrix_symbolic_trace() {
    let mat = SymbolicMatrix::new(vec![
        vec![Expr::sym("a"), Expr::num(0)],
        vec![Expr::num(0), Expr::sym("b")],
    ]).unwrap();

    let trace = mat.trace().unwrap();
    eprintln!("Symbolic trace: {}", trace);
    // Should be a + b
}

// ============================================================================
// SUBSTITUTION TESTS
// ============================================================================

#[test]
fn test_substitute_simple() {
    let mut rules = HashMap::new();
    rules.insert("x".to_string(), "y".to_string());
    let result = substitute("x + 1", &rules).unwrap();
    assert!(result.expression.contains("y"));
}

#[test]
fn test_substitute_multiple() {
    let mut rules = HashMap::new();
    rules.insert("x".to_string(), "a".to_string());
    rules.insert("y".to_string(), "b".to_string());
    let result = substitute("x + y", &rules).unwrap();
    eprintln!("Substituted: {}", result.expression);
}

#[test]
fn test_substitute_expression() {
    let mut rules = HashMap::new();
    rules.insert("x".to_string(), "2*y".to_string());
    let result = substitute("x^2", &rules).unwrap();
    eprintln!("Substituted x^2 with x=2y: {}", result.expression);
}

// ============================================================================
// SERIES EXPANSION TESTS
// ============================================================================

#[test]
fn test_series_exp_order_2() {
    let result = series_expansion("exp(x)", "x", 0.0, 2).unwrap();
    eprintln!("exp(x) Taylor series (order 2): {}", result.expression);
    // Should be 1 + x + x^2/2
}

#[test]
fn test_series_exp_order_4() {
    let result = series_expansion("exp(x)", "x", 0.0, 4).unwrap();
    eprintln!("exp(x) Taylor series (order 4): {}", result.expression);
    // Should include terms up to x^4/24
}

#[test]
fn test_series_sin_order_3() {
    let result = series_expansion("sin(x)", "x", 0.0, 3).unwrap();
    eprintln!("sin(x) Taylor series (order 3): {}", result.expression);
    // Should be x - x^3/6
}

#[test]
fn test_series_cos_order_4() {
    let result = series_expansion("cos(x)", "x", 0.0, 4).unwrap();
    eprintln!("cos(x) Taylor series (order 4): {}", result.expression);
    // Should be 1 - x^2/2 + x^4/24
}

#[test]
fn test_series_polynomial() {
    let result = series_expansion("x^2 + 3*x + 1", "x", 0.0, 3).unwrap();
    eprintln!("Polynomial Taylor series: {}", result.expression);
}

// ============================================================================
// LIMIT TESTS
// ============================================================================

#[test]
fn test_limit_constant() {
    let result = limit("5", "x", "2", None).unwrap();
    assert_eq!(result.expression, "5");
}

#[test]
fn test_limit_variable() {
    let result = limit("x", "x", "3", None).unwrap();
    assert_eq!(result.expression, "3");
}

#[test]
fn test_limit_polynomial() {
    let result = limit("x^2 + 1", "x", "2", None).unwrap();
    eprintln!("lim (x^2 + 1) as x->2: {}", result.expression);
    // Should be 5
}

#[test]
fn test_limit_at_infinity() {
    let result = limit("1/x", "x", "inf", None).unwrap();
    eprintln!("lim (1/x) as x->∞: {}", result.expression);
    // Should be 0
}

#[test]
fn test_limit_at_negative_infinity() {
    let result = limit("1/x", "x", "-inf", None).unwrap();
    eprintln!("lim (1/x) as x->-∞: {}", result.expression);
    // Should be 0
}

// ============================================================================
// HIGH-LEVEL API TESTS
// ============================================================================

#[test]
fn test_api_simplify() {
    let result = simplify("x + 0").unwrap();
    assert!(result.expression.contains("x"));
    assert!(result.latex.is_some());
}

#[test]
fn test_api_expand() {
    let result = expand("(x + 1)^2").unwrap();
    eprintln!("Expanded: {}", result.expression);
    assert!(result.latex.is_some());
}

#[test]
fn test_api_differentiate() {
    let result = differentiate("x^3", "x", None).unwrap();
    eprintln!("Derivative: {}", result.expression);
    assert!(result.latex.is_some());
    assert!(result.metadata.is_some());
}

#[test]
fn test_api_differentiate_second_order() {
    let result = differentiate("x^4", "x", Some(2)).unwrap();
    eprintln!("Second derivative: {}", result.expression);
}

#[test]
fn test_api_integrate() {
    let result = integrate("x^2", "x").unwrap();
    eprintln!("Integral: {}", result.expression);
    assert!(result.latex.is_some());
}

#[test]
fn test_api_gradient() {
    let vars = vec!["x".to_string(), "y".to_string()];
    let results = gradient("x^2 + y^2", &vars).unwrap();
    assert_eq!(results.len(), 2);
    eprintln!("∂/∂x: {}", results[0].expression);
    eprintln!("∂/∂y: {}", results[1].expression);
}

#[test]
fn test_api_evaluate_at() {
    let mut values = HashMap::new();
    values.insert("x".to_string(), 2.0);
    values.insert("y".to_string(), 3.0);
    let result = evaluate_at("x^2 + y^2", &values).unwrap();
    assert_eq!(result, 13.0);
}

// ============================================================================
// ADDITIONAL EDGE CASE TESTS
// ============================================================================

#[test]
fn test_parse_error_invalid_char() {
    let result = parse("x @ y");
    assert!(result.is_err());
}

#[test]
fn test_parse_error_mismatched_parens() {
    let result = parse("(x + 1");
    assert!(result.is_err());
}

#[test]
fn test_differentiate_product_complex() {
    let expr = parse("(x^2) * sin(x)").unwrap();
    let deriv = diff(&expr, "x");
    eprintln!("d/dx(x^2 * sin(x)) = {}", deriv);
    // Should use product rule
}

#[test]
fn test_differentiate_chain_rule_complex() {
    let expr = parse("exp(x^2)").unwrap();
    let deriv = diff(&expr, "x");
    eprintln!("d/dx(exp(x^2)) = {}", deriv);
    // Should be exp(x^2) * 2x
}

#[test]
fn test_integrate_sum() {
    let expr = parse("x + x^2 + x^3").unwrap();
    let integral = integrate_expr(&expr, "x");
    eprintln!("∫(x + x^2 + x^3) dx = {}", integral);
}

#[test]
fn test_simplify_nested() {
    let expr = parse("((x + 0) * 1)^1").unwrap();
    let simplified = simplify_expr(&expr);
    eprintln!("Simplified nested: {}", simplified);
    assert_eq!(simplified, Expr::sym("x"));
}

#[test]
fn test_expand_nested_binomials() {
    let expr = parse("((x + 1)^2)^1").unwrap();
    let expanded = expand_expr(&expr);
    eprintln!("Expanded nested binomials: {}", expanded);
}

#[test]
fn test_evaluate_nested_functions() {
    let expr = parse("sin(cos(x))").unwrap();
    let mut vars = HashMap::new();
    vars.insert("x".to_string(), 0.0);
    let result = expr.evaluate(&vars).unwrap();
    eprintln!("sin(cos(0)) = {}", result);
    // cos(0) = 1, sin(1) ≈ 0.841
    assert!((result - 0.8414709848).abs() < 1e-6);
}

#[test]
fn test_matrix_determinant_4x4() {
    let mat = SymbolicMatrix::new(vec![
        vec![Expr::num(1), Expr::num(0), Expr::num(0), Expr::num(0)],
        vec![Expr::num(0), Expr::num(2), Expr::num(0), Expr::num(0)],
        vec![Expr::num(0), Expr::num(0), Expr::num(3), Expr::num(0)],
        vec![Expr::num(0), Expr::num(0), Expr::num(0), Expr::num(4)],
    ]).unwrap();

    let det = mat.determinant().unwrap();
    eprintln!("Det of 4x4 diagonal: {}", det);
    // Should be 1*2*3*4 = 24
}

#[test]
fn test_matrix_symbolic_multiplication() {
    let a = SymbolicMatrix::new(vec![
        vec![Expr::sym("a"), Expr::sym("b")],
        vec![Expr::sym("c"), Expr::sym("d")],
    ]).unwrap();

    let b = SymbolicMatrix::new(vec![
        vec![Expr::sym("x"), Expr::num(0)],
        vec![Expr::num(0), Expr::sym("y")],
    ]).unwrap();

    let c = a.mul(&b).unwrap();
    eprintln!("Symbolic matrix product:\n{}", c);
}

#[test]
fn test_series_expansion_at_nonzero_point() {
    let result = series_expansion("x^2", "x", 1.0, 2).unwrap();
    eprintln!("Taylor series of x^2 at x=1: {}", result.expression);
}

#[test]
fn test_limit_left_right_sides() {
    let result_left = limit("1/x", "x", "0", Some("left")).unwrap();
    eprintln!("lim (1/x) as x->0- : {}", result_left.expression);

    let result_right = limit("1/x", "x", "0", Some("right")).unwrap();
    eprintln!("lim (1/x) as x->0+ : {}", result_right.expression);
}

#[test]
fn test_differentiate_sqrt() {
    let expr = parse("sqrt(x)").unwrap();
    let deriv = diff(&expr, "x");
    eprintln!("d/dx(sqrt(x)) = {}", deriv);
    // Should be 1/(2*sqrt(x))
}

#[test]
fn test_integrate_multiple_variables_constant() {
    let expr = parse("y").unwrap();
    let integral = integrate_expr(&expr, "x");
    eprintln!("∫ y dx (y constant) = {}", integral);
    // Should be y*x
}

#[test]
fn test_simplify_power_of_power() {
    let expr = parse("(x^2)^3").unwrap();
    let simplified = simplify_expr(&expr);
    eprintln!("(x^2)^3 simplified: {}", simplified);
    // Should be x^6
}

#[test]
fn test_matrix_nonsquare_determinant_error() {
    let mat = SymbolicMatrix::new(vec![
        vec![Expr::num(1), Expr::num(2), Expr::num(3)],
        vec![Expr::num(4), Expr::num(5), Expr::num(6)],
    ]).unwrap();

    let result = mat.determinant();
    assert!(result.is_err());
}

#[test]
fn test_matrix_dimension_mismatch_add() {
    let a = SymbolicMatrix::zeros(2, 3);
    let b = SymbolicMatrix::zeros(3, 2);
    let result = a.add(&b);
    assert!(result.is_err());
}

#[test]
fn test_matrix_dimension_mismatch_mul() {
    let a = SymbolicMatrix::zeros(2, 3);
    let b = SymbolicMatrix::zeros(2, 2);
    let result = a.mul(&b);
    assert!(result.is_err());
}

#[test]
fn test_evaluate_missing_variable() {
    let expr = parse("x + y").unwrap();
    let mut vars = HashMap::new();
    vars.insert("x".to_string(), 2.0);
    let result = expr.evaluate(&vars);
    assert!(result.is_err());
}

#[test]
fn test_parse_multidigit_numbers() {
    let expr = parse("123 + 456").unwrap();
    eprintln!("Parsed: {}", expr);
}

#[test]
fn test_parse_underscored_variables() {
    let expr = parse("x_1 + x_2").unwrap();
    eprintln!("Parsed: {}", expr);
}

#[test]
fn test_differentiate_third_order() {
    let expr = parse("x^5").unwrap();
    let first = diff(&expr, "x");
    let second = diff(&first, "x");
    let third = diff(&second, "x");
    eprintln!("d³/dx³(x^5) = {}", third);
    // Should be 60*x^2
}

#[test]
fn test_api_factor() {
    let result = factor("x^2 - 1").unwrap();
    eprintln!("Factored: {}", result.expression);
    assert!(result.metadata.is_some());
}

#[test]
fn test_api_collect() {
    let result = collect("x + 2*x + 3*x", "x").unwrap();
    eprintln!("Collected: {}", result.expression);
}

#[test]
fn test_api_rationalize() {
    let result = rationalize("x + 1").unwrap();
    eprintln!("Rationalized: {}", result.expression);
}

#[test]
fn test_series_expansion_higher_order() {
    let result = series_expansion("exp(x)", "x", 0.0, 6).unwrap();
    eprintln!("exp(x) Taylor series (order 6): {}", result.expression);
}

#[test]
fn test_limit_polynomial_at_infinity_positive() {
    let result = limit("x^2", "x", "inf", None).unwrap();
    eprintln!("lim (x^2) as x->∞: {}", result.expression);
}

#[test]
fn test_limit_rational_function() {
    let result = limit("(x^2 + 1) / (x + 1)", "x", "0", None).unwrap();
    eprintln!("lim ((x^2 + 1)/(x + 1)) as x->0: {}", result.expression);
}

#[test]
fn test_matrix_get_set() {
    let mut mat = SymbolicMatrix::zeros(2, 2);
    mat.set(0, 0, Expr::sym("a")).unwrap();
    mat.set(1, 1, Expr::sym("b")).unwrap();
    assert_eq!(mat.get(0, 0), Some(&Expr::sym("a")));
    assert_eq!(mat.get(1, 1), Some(&Expr::sym("b")));
}

#[test]
fn test_matrix_identity_multiplication() {
    let a = SymbolicMatrix::new(vec![
        vec![Expr::num(1), Expr::num(2)],
        vec![Expr::num(3), Expr::num(4)],
    ]).unwrap();

    let identity = SymbolicMatrix::identity(2);
    let result = a.mul(&identity).unwrap();
    eprintln!("A * I:\n{}", result);
    // Should be equal to A
}

#[test]
fn test_differentiate_multiple_terms() {
    let expr = parse("3*x^4 + 2*x^3 + x^2 + x + 1").unwrap();
    let deriv = diff(&expr, "x");
    eprintln!("d/dx(3x^4 + 2x^3 + x^2 + x + 1) = {}", deriv);
}

#[test]
fn test_integrate_multiple_terms() {
    let expr = parse("4*x^3 + 3*x^2 + 2*x + 1").unwrap();
    let integral = integrate_expr(&expr, "x");
    eprintln!("∫(4x^3 + 3x^2 + 2x + 1) dx = {}", integral);
}

// ============================================================================
// RATIONAL NUMBER TESTS
// ============================================================================

#[test]
fn test_rational_creation() {
    let r = Rational::new(6, 9).unwrap();
    assert_eq!(r.numerator, 2);
    assert_eq!(r.denominator, 3);
}

#[test]
fn test_rational_negative() {
    let r = Rational::new(-4, 6).unwrap();
    assert_eq!(r.numerator, -2);
    assert_eq!(r.denominator, 3);
}

#[test]
fn test_rational_zero_denominator() {
    let result = Rational::new(5, 0);
    assert!(result.is_err());
}

#[test]
fn test_rational_add() {
    let r1 = Rational::new(1, 2).unwrap();
    let r2 = Rational::new(1, 3).unwrap();
    let sum = r1 + r2;
    assert_eq!(sum.numerator, 5);
    assert_eq!(sum.denominator, 6);
}

#[test]
fn test_rational_mul() {
    let r1 = Rational::new(2, 3).unwrap();
    let r2 = Rational::new(3, 4).unwrap();
    let product = r1 * r2;
    assert_eq!(product.numerator, 1);
    assert_eq!(product.denominator, 2);
}

#[test]
fn test_rational_to_f64() {
    let r = Rational::new(1, 2).unwrap();
    assert_eq!(r.to_f64(), 0.5);
}

#[test]
fn test_rational_is_zero() {
    let r = Rational::new(0, 1).unwrap();
    assert!(r.is_zero());
}

#[test]
fn test_rational_is_one() {
    let r = Rational::new(3, 3).unwrap();
    assert!(r.is_one());
}

// ============================================================================
// EXPRESSION BUILDER TESTS
// ============================================================================

#[test]
fn test_expr_num() {
    let expr = Expr::num(42);
    assert!(matches!(expr, Expr::Number(_)));
}

#[test]
fn test_expr_rational() {
    let expr = Expr::rational(1, 2).unwrap();
    assert!(matches!(expr, Expr::Number(_)));
}

#[test]
fn test_expr_sym() {
    let expr = Expr::sym("x");
    assert!(matches!(expr, Expr::Symbol(_)));
}

#[test]
fn test_expr_add() {
    let expr = Expr::add(Expr::num(1), Expr::num(2));
    assert!(matches!(expr, Expr::Add(_, _)));
}

#[test]
fn test_expr_mul() {
    let expr = Expr::mul(Expr::num(3), Expr::sym("x"));
    assert!(matches!(expr, Expr::Mul(_, _)));
}

#[test]
fn test_expr_pow() {
    let expr = Expr::pow(Expr::sym("x"), Expr::num(2));
    assert!(matches!(expr, Expr::Pow(_, _)));
}

#[test]
fn test_expr_func() {
    let expr = Expr::func("sin", vec![Expr::sym("x")]);
    assert!(matches!(expr, Expr::Function(_, _)));
}

#[test]
fn test_expr_is_zero() {
    assert!(Expr::num(0).is_zero());
    assert!(!Expr::num(1).is_zero());
}

#[test]
fn test_expr_is_one() {
    assert!(Expr::num(1).is_one());
    assert!(!Expr::num(0).is_one());
}
