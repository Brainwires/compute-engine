// Unit tests for mathematics::symbolic_cas::parser
use computational_engine::analyze::symbolic::parser::*;

use super::*;

    #[test]
    fn test_parse_simple() {
        let expr = parse("x + 1").unwrap();
        assert!(matches!(expr, Expr::Add(_, _)));
    }

    #[test]
    fn test_parse_multiplication() {
        let expr = parse("2 * x").unwrap();
        assert!(matches!(expr, Expr::Mul(_, _)));
    }

    #[test]
    fn test_parse_power() {
        let expr = parse("x^2").unwrap();
        assert!(matches!(expr, Expr::Pow(_, _)));
    }

    #[test]
    fn test_parse_function() {
        let expr = parse("sin(x)").unwrap();
        assert!(matches!(expr, Expr::Function(_, _)));
    }

    #[test]
    fn test_parse_complex() {
        let expr = parse("(x + 1)^2").unwrap();
        println!("{}", expr);
    }
