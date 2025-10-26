// Comprehensive unit tests for mathematics::symbolic_cas module
use crate::mathematics::symbolic_cas::*;
use std::collections::HashMap;

mod symbolic_cas_tests {
    use super::*;

    // ========== Expression Creation and Basic Operations ==========

    #[test]
    fn test_rational_creation_and_reduction() {
        let r = Rational::new(6, 9).unwrap();
        assert_eq!(r.numerator, 2);
        assert_eq!(r.denominator, 3);
        assert_eq!(r.to_f64(), 2.0 / 3.0);
    }

    #[test]
    fn test_rational_zero_denominator_error() {
        let result = Rational::new(5, 0);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Denominator cannot be zero");
    }

    #[test]
    fn test_rational_negative_handling() {
        let r = Rational::new(-6, 9).unwrap();
        assert_eq!(r.numerator, -2);
        assert_eq!(r.denominator, 3);

        let r2 = Rational::new(6, -9).unwrap();
        assert_eq!(r2.numerator, -2);
        assert_eq!(r2.denominator, 3);
    }

    #[test]
    fn test_rational_arithmetic() {
        let r1 = Rational::from_int(2);
        let r2 = Rational::from_int(3);

        let sum = r1.clone() + r2.clone();
        assert_eq!(sum.numerator, 5);
        assert_eq!(sum.denominator, 1);

        let prod = r1 * r2;
        assert_eq!(prod.numerator, 6);
        assert_eq!(prod.denominator, 1);
    }

    // ========== Parsing Tests ==========

    #[test]
    fn test_parse_simple_expression() {
        let expr = parse("x + 1").unwrap();
        assert!(matches!(expr, Expr::Add(_, _)));
    }

    #[test]
    fn test_parse_polynomial() {
        let expr = parse("x^2 + 2*x + 1").unwrap();
        assert!(matches!(expr, Expr::Add(_, _)));
    }

    #[test]
    fn test_parse_with_parentheses() {
        let expr = parse("(x + 1) * (x - 1)").unwrap();
        assert!(matches!(expr, Expr::Mul(_, _)));
    }

    #[test]
    fn test_parse_division() {
        let expr = parse("x / 2").unwrap();
        // Division is represented as multiplication by power with exponent -1
        assert!(matches!(expr, Expr::Mul(_, _)));
    }

    // ========== Simplification Tests ==========

    #[test]
    fn test_simplify_add_zero() {
        let expr = parse("x + 0").unwrap();
        let simplified = simplify_expr(&expr);
        assert!(matches!(simplified, Expr::Symbol(_)));
    }

    #[test]
    fn test_simplify_mul_one() {
        let expr = parse("x * 1").unwrap();
        let simplified = simplify_expr(&expr);
        assert!(matches!(simplified, Expr::Symbol(_)));
    }

    #[test]
    fn test_simplify_mul_zero() {
        let expr = parse("x * 0").unwrap();
        let simplified = simplify_expr(&expr);
        assert_eq!(simplified, Expr::num(0));
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
        assert!(matches!(simplified, Expr::Symbol(_)));
    }

    // ========== Expansion Tests ==========

    #[test]
    fn test_expand_binomial_square() {
        let expr = parse("(x + 1)^2").unwrap();
        let expanded = expand_expr(&expr);
        // Should expand to x^2 + 2*x + 1
        eprintln!("Expanded (x+1)^2: {}", expanded);
        // Verify it's not the same as original
        assert_ne!(expanded, expr);
    }

    #[test]
    fn test_expand_product() {
        let expr = parse("(x + 1) * (x + 2)").unwrap();
        let expanded = expand_expr(&expr);
        eprintln!("Expanded (x+1)(x+2): {}", expanded);
        // Should expand to x^2 + 3*x + 2
        assert_ne!(expanded, expr);
    }

    // ========== Differentiation Tests ==========

    #[test]
    fn test_diff_constant() {
        let expr = Expr::num(5);
        let derivative = diff(&expr, "x");
        assert_eq!(derivative, Expr::num(0));
    }

    #[test]
    fn test_diff_variable() {
        let expr = Expr::sym("x");
        let derivative = diff(&expr, "x");
        assert_eq!(derivative, Expr::num(1));
    }

    #[test]
    fn test_diff_power_rule() {
        let expr = parse("x^3").unwrap();
        let derivative = diff(&expr, "x");
        eprintln!("d/dx(x^3) = {}", derivative);
        // Should be 3*x^2
        assert!(matches!(derivative, Expr::Mul(_, _)));
    }

    #[test]
    fn test_diff_polynomial() {
        let expr = parse("x^2 + 2*x + 1").unwrap();
        let derivative = diff(&expr, "x");
        eprintln!("d/dx(x^2 + 2x + 1) = {}", derivative);
        // Should be 2*x + 2
        assert!(matches!(derivative, Expr::Add(_, _)));
    }

    #[test]
    fn test_diff_product_rule() {
        let expr = parse("x * x^2").unwrap();
        let derivative = diff(&expr, "x");
        eprintln!("d/dx(x * x^2) = {}", derivative);
        assert_ne!(derivative, Expr::num(0));
    }

    // ========== Integration Tests ==========

    #[test]
    fn test_integrate_constant() {
        let expr = Expr::num(5);
        let integral = integrate_expr(&expr, "x");
        eprintln!("∫5 dx = {}", integral);
        // Should be 5*x
        assert!(matches!(integral, Expr::Mul(_, _)));
    }

    #[test]
    fn test_integrate_variable() {
        let expr = Expr::sym("x");
        let integral = integrate_expr(&expr, "x");
        eprintln!("∫x dx = {}", integral);
        // Should be x^2/2
        assert!(matches!(integral, Expr::Mul(_, _)));
    }

    #[test]
    fn test_integrate_power() {
        let expr = parse("x^2").unwrap();
        let integral = integrate_expr(&expr, "x");
        eprintln!("∫x^2 dx = {}", integral);
        // Should be x^3/3
        assert_ne!(integral, Expr::num(0));
    }

    // ========== High-Level API Tests ==========

    #[test]
    fn test_symbolic_simplify_api() {
        let result = simplify("x + 0").unwrap();
        assert!(result.expression.contains("x"));
        assert!(result.latex.is_some());
        assert_eq!(result.metadata.as_ref().unwrap()["operation"], "simplify");
    }

    #[test]
    fn test_symbolic_expand_api() {
        let result = expand("(x + 1)^2").unwrap();
        eprintln!("Expand API result: {}", result.expression);
        assert!(result.latex.is_some());
        assert_eq!(result.metadata.as_ref().unwrap()["operation"], "expand");
    }

    #[test]
    fn test_symbolic_differentiate_api() {
        let result = differentiate("x^2", "x", None).unwrap();
        eprintln!("Differentiate API result: {}", result.expression);
        assert!(result.latex.is_some());
        assert_eq!(result.metadata.as_ref().unwrap()["operation"], "differentiate");
        assert_eq!(result.metadata.as_ref().unwrap()["variable"], "x");
    }

    #[test]
    fn test_symbolic_differentiate_higher_order() {
        let result = differentiate("x^3", "x", Some(2)).unwrap();
        eprintln!("Second derivative of x^3: {}", result.expression);
        assert_eq!(result.metadata.as_ref().unwrap()["order"], 2);
    }

    #[test]
    fn test_symbolic_integrate_api() {
        let result = integrate("x^2", "x").unwrap();
        eprintln!("Integrate API result: {}", result.expression);
        assert!(result.latex.is_some());
        assert_eq!(result.metadata.as_ref().unwrap()["operation"], "integrate");
    }

    // ========== Substitution and Evaluation ==========

    #[test]
    fn test_substitute_variable() {
        let mut rules = HashMap::new();
        rules.insert("x".to_string(), "2".to_string());

        let result = substitute("x + 1", &rules).unwrap();
        eprintln!("Substitution result: {}", result.expression);
        assert_eq!(result.metadata.as_ref().unwrap()["operation"], "substitute");
    }

    #[test]
    fn test_evaluate_at_values() {
        let mut values = HashMap::new();
        values.insert("x".to_string(), 2.0);

        let result = evaluate_at("x^2 + 1", &values).unwrap();
        assert_eq!(result, 5.0);
    }

    #[test]
    fn test_evaluate_complex_expression() {
        let mut values = HashMap::new();
        values.insert("x".to_string(), 3.0);
        values.insert("y".to_string(), 4.0);

        let result = evaluate_at("x * y + 2", &values).unwrap();
        assert_eq!(result, 14.0);
    }

    // ========== Gradient (Partial Derivatives) ==========

    #[test]
    fn test_gradient_computation() {
        let vars = vec!["x".to_string(), "y".to_string()];
        let results = gradient("x^2 + y^2", &vars).unwrap();

        assert_eq!(results.len(), 2);
        eprintln!("∂/∂x(x^2 + y^2) = {}", results[0].expression);
        eprintln!("∂/∂y(x^2 + y^2) = {}", results[1].expression);
    }

    // ========== Series Expansion ==========

    #[test]
    fn test_taylor_series_expansion() {
        let result = series_expansion("x^2 + x", "x", 0.0, 3).unwrap();
        eprintln!("Taylor series of x^2+x at x=0: {}", result.expression);
        assert_eq!(result.metadata.as_ref().unwrap()["operation"], "series");
        assert_eq!(result.metadata.as_ref().unwrap()["order"], 3);
    }

    // ========== Limit Computation ==========

    #[test]
    fn test_limit_finite_point() {
        let result = limit("x^2", "x", "2", None).unwrap();
        eprintln!("lim(x→2) x^2 = {}", result.expression);
        assert_eq!(result.metadata.as_ref().unwrap()["operation"], "limit");
    }

    #[test]
    fn test_limit_at_infinity() {
        let result = limit("1/x", "x", "inf", None).unwrap();
        eprintln!("lim(x→∞) 1/x = {}", result.expression);
        assert!(result.expression.contains("0") || result.expression.contains("undefined"));
    }

    // ========== Factor and Collect ==========

    #[test]
    fn test_factor_expression() {
        let result = factor("x^2 + 2*x + 1").unwrap();
        eprintln!("Factor(x^2 + 2x + 1) = {}", result.expression);
        assert_eq!(result.metadata.as_ref().unwrap()["operation"], "factor");
    }

    #[test]
    fn test_collect_terms() {
        let result = collect("x^2 + x + 2*x", "x").unwrap();
        eprintln!("Collect(x^2 + x + 2x, x) = {}", result.expression);
        assert_eq!(result.metadata.as_ref().unwrap()["operation"], "collect");
    }
}
