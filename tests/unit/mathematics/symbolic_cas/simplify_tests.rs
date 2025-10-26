// Unit tests for mathematics::symbolic_cas::simplify
use computational_engine::mathematics::symbolic_cas::simplify::*;

use super::*;
    use crate::mathematics::symbolic_cas::parser::parse;

    #[test]
    fn test_simplify_add_zero() {
        let expr = parse("x + 0").unwrap();
        let simplified = simplify(&expr);
        assert!(matches!(simplified, Expr::Symbol(_)));
    }

    #[test]
    fn test_simplify_mul_one() {
        let expr = parse("x * 1").unwrap();
        let simplified = simplify(&expr);
        assert!(matches!(simplified, Expr::Symbol(_)));
    }

    #[test]
    fn test_expand_square() {
        let expr = parse("(x + 1)^2").unwrap();
        let expanded = expand(&expr);
        println!("Expanded: {}", expanded);
    }
