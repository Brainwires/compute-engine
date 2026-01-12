// Unit tests for mathematics::symbolic_cas::expr
use computational_engine::analyze::symbolic::expr::*;

use super::*;

    #[test]
    fn test_rational() {
        let r = Rational::new(6, 9).unwrap();
        assert_eq!(r.numerator, 2);
        assert_eq!(r.denominator, 3);
    }

    #[test]
    fn test_rational_zero_denominator() {
        let result = Rational::new(5, 0);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Denominator cannot be zero");
    }

    #[test]
    fn test_expr_display() {
        let expr = Expr::add(Expr::sym("x"), Expr::num(1));
        assert_eq!(format!("{}", expr), "(x + 1)");
    }

    #[test]
    fn test_evaluate() {
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), 2.0);

        let expr = Expr::add(Expr::sym("x"), Expr::num(3));
        assert_eq!(expr.evaluate(&vars).unwrap(), 5.0);
    }
