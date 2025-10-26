// Unit tests for mathematics::symbolic_cas::differentiate
use computational_engine::mathematics::symbolic_cas::differentiate::*;

use super::*;
    use crate::mathematics::symbolic_cas::parser::parse;

    #[test]
    fn test_diff_constant() {
        let expr = Expr::num(5);
        let diff = differentiate(&expr, "x");
        assert_eq!(diff, Expr::num(0));
    }

    #[test]
    fn test_diff_variable() {
        let expr = Expr::sym("x");
        let diff = differentiate(&expr, "x");
        assert_eq!(diff, Expr::num(1));
    }

    #[test]
    fn test_diff_power() {
        let expr = parse("x^2").unwrap();
        let diff = differentiate(&expr, "x");
        println!("d/dx(x^2) = {}", diff);
    }

    #[test]
    fn test_diff_polynomial() {
        let expr = parse("x^3 + 2*x^2 + x").unwrap();
        let diff = differentiate(&expr, "x");
        println!("d/dx(x^3 + 2*x^2 + x) = {}", diff);
    }
