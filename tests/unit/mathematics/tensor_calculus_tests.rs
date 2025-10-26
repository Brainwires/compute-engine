//! Unit tests for mathematics::tensor_calculus

use crate::mathematics::tensor_calculus::symbolic::SymbolicExpr;
use crate::mathematics::tensor_calculus::tensor::*;

// ==================== Symbolic Expression Tests ====================

#[test]
fn test_symbolic_parse_constant() {
    let expr = SymbolicExpr::parse("42").unwrap();
    assert!(matches!(expr, SymbolicExpr::Constant(42.0)));
}

#[test]
fn test_symbolic_parse_zero() {
    let expr = SymbolicExpr::parse("0").unwrap();
    assert!(matches!(expr, SymbolicExpr::Zero));
}

#[test]
fn test_symbolic_parse_one() {
    let expr = SymbolicExpr::parse("1").unwrap();
    assert!(matches!(expr, SymbolicExpr::One));
}

#[test]
fn test_symbolic_parse_variable() {
    let expr = SymbolicExpr::parse("x").unwrap();
    assert!(matches!(expr, SymbolicExpr::Variable(_)));
}

#[test]
fn test_symbolic_is_zero() {
    assert!(SymbolicExpr::Zero.is_zero());
    assert!(!SymbolicExpr::One.is_zero());
    assert!(!SymbolicExpr::Constant(5.0).is_zero());
}

#[test]
fn test_symbolic_simplify_zero_addition() {
    let expr = SymbolicExpr::Add(
        Box::new(SymbolicExpr::Zero),
        Box::new(SymbolicExpr::Constant(5.0)),
    );
    let simplified = expr.simplify();
    assert!(matches!(simplified, SymbolicExpr::Constant(5.0)));
}

#[test]
fn test_symbolic_simplify_zero_multiplication() {
    let expr = SymbolicExpr::Multiply(
        Box::new(SymbolicExpr::Zero),
        Box::new(SymbolicExpr::Variable("x".to_string())),
    );
    let simplified = expr.simplify();
    assert!(simplified.is_zero());
}

#[test]
fn test_symbolic_simplify_one_multiplication() {
    let expr = SymbolicExpr::Multiply(
        Box::new(SymbolicExpr::One),
        Box::new(SymbolicExpr::Constant(5.0)),
    );
    let simplified = expr.simplify();
    assert!(matches!(simplified, SymbolicExpr::Constant(5.0)));
}

#[test]
fn test_symbolic_derivative_constant() {
    let expr = SymbolicExpr::Constant(42.0);
    let deriv = expr.derivative("x");
    assert!(deriv.is_zero());
}

#[test]
fn test_symbolic_derivative_variable() {
    let expr = SymbolicExpr::Variable("x".to_string());
    let deriv = expr.derivative("x");
    // d/dx(x) = 1
    assert!(matches!(deriv, SymbolicExpr::One));
}

#[test]
fn test_symbolic_derivative_different_variable() {
    let expr = SymbolicExpr::Variable("y".to_string());
    let deriv = expr.derivative("x");
    // d/dx(y) = 0
    assert!(deriv.is_zero());
}

#[test]
fn test_symbolic_display() {
    let expr = SymbolicExpr::Constant(5.0);
    assert_eq!(expr.to_string(), "5");

    let var = SymbolicExpr::Variable("x".to_string());
    assert_eq!(var.to_string(), "x");
}

// ==================== Metric Tensor Tests ====================

#[test]
fn test_parse_flat_metric_2d() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string()],
    ];
    let coords = vec!["x".to_string(), "y".to_string()];

    let result = parse_metric_tensor(metric_strings, &coords);
    assert!(result.is_ok());

    let metric = result.unwrap();
    assert_eq!(metric.len(), 2);
    assert_eq!(metric[0].len(), 2);
}

#[test]
fn test_parse_minkowski_metric() {
    let metric_strings = vec![
        vec!["-1".to_string(), "0".to_string(), "0".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string(), "0".to_string(), "0".to_string()],
        vec!["0".to_string(), "0".to_string(), "1".to_string(), "0".to_string()],
        vec!["0".to_string(), "0".to_string(), "0".to_string(), "1".to_string()],
    ];
    let coords = vec!["t".to_string(), "x".to_string(), "y".to_string(), "z".to_string()];

    let result = parse_metric_tensor(metric_strings, &coords);
    assert!(result.is_ok());

    let metric = result.unwrap();
    assert_eq!(metric.len(), 4);
}

#[test]
fn test_parse_metric_non_square() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string(), "0".to_string()], // Different length
    ];
    let coords = vec!["x".to_string(), "y".to_string()];

    let result = parse_metric_tensor(metric_strings, &coords);
    assert!(result.is_err());
}

// ==================== Christoffel Symbol Tests ====================

#[test]
fn test_christoffel_flat_metric() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string()],
    ];
    let coords = vec!["x".to_string(), "y".to_string()];

    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();
    let result = calculate_christoffel_symbols(&metric, &coords).unwrap();

    // Flat metric should have all zero Christoffel symbols
    assert_eq!(result.symbols.len(), 0); // No non-zero components
    assert_eq!(result.dimension, 2);
}

#[test]
fn test_christoffel_minkowski() {
    let metric_strings = vec![
        vec!["-1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string()],
    ];
    let coords = vec!["t".to_string(), "x".to_string()];

    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();
    let result = calculate_christoffel_symbols(&metric, &coords).unwrap();

    // Minkowski metric is constant, all Christoffel symbols are zero
    assert_eq!(result.symbols.len(), 0);
}

// ==================== Riemann Tensor Tests ====================

#[test]
fn test_riemann_flat_metric() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string()],
    ];
    let coords = vec!["x".to_string(), "y".to_string()];

    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();
    let result = calculate_riemann_tensor(&metric, &coords).unwrap();

    // Flat metric has zero Riemann curvature
    assert_eq!(result.components.len(), 0); // All components are zero
    assert_eq!(result.dimension, 2);
}

// ==================== Ricci Tensor Tests ====================

#[test]
fn test_ricci_flat_metric() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string()],
    ];
    let coords = vec!["x".to_string(), "y".to_string()];

    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();
    let result = calculate_ricci_tensor(&metric, &coords).unwrap();

    // Flat metric has zero Ricci tensor
    assert_eq!(result.components.len(), 0);
    assert_eq!(result.dimension, 2);
}

// ==================== Ricci Scalar Tests ====================

#[test]
fn test_ricci_scalar_flat_metric() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string()],
    ];
    let coords = vec!["x".to_string(), "y".to_string()];

    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();
    let result = calculate_ricci_scalar(&metric, &coords).unwrap();

    // Flat metric has zero Ricci scalar
    assert_eq!(result.expression, "0");
}

// ==================== Einstein Tensor Tests ====================

#[test]
fn test_einstein_tensor_flat_metric() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string()],
    ];
    let coords = vec!["x".to_string(), "y".to_string()];

    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();
    let result = calculate_einstein_tensor(&metric, &coords).unwrap();

    // Flat metric has zero Einstein tensor
    assert_eq!(result.components.len(), 0);
    assert_eq!(result.dimension, 2);
}

// ==================== Helper Function Tests ====================

#[test]
fn test_tensor_component_creation() {
    let component = TensorComponent {
        indices: vec![0, 1, 2],
        expression: "x + y".to_string(),
    };

    assert_eq!(component.indices.len(), 3);
    assert_eq!(component.expression, "x + y");
}

#[test]
fn test_christoffel_result_creation() {
    let result = ChristoffelResult {
        symbols: vec![],
        dimension: 4,
    };

    assert_eq!(result.dimension, 4);
    assert_eq!(result.symbols.len(), 0);
}

#[test]
fn test_riemann_result_creation() {
    let result = RiemannResult {
        components: vec![],
        dimension: 3,
    };

    assert_eq!(result.dimension, 3);
    assert_eq!(result.components.len(), 0);
}
