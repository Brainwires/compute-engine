// Unit tests for mathematics::tensor_calculus::symbolic
use computational_engine::mathematics::tensor_calculus::symbolic::*;

use super::*;

    #[test]
    fn test_parse_simple() {
        assert_eq!(
            SymbolicExpr::parse("x").unwrap(),
            SymbolicExpr::Variable("x".to_string())
        );
        assert_eq!(
            SymbolicExpr::parse("42").unwrap(),
            SymbolicExpr::Constant(42.0)
        );
        assert_eq!(SymbolicExpr::parse("0").unwrap(), SymbolicExpr::Zero);
        assert_eq!(SymbolicExpr::parse("1").unwrap(), SymbolicExpr::One);
    }

    #[test]
    fn test_parse_power() {
        match SymbolicExpr::parse("r^2").unwrap() {
            SymbolicExpr::Power(base, exp) => {
                assert_eq!(*base, SymbolicExpr::Variable("r".to_string()));
                assert_eq!(*exp, SymbolicExpr::Constant(2.0));
            }
            _ => panic!("Expected power expression"),
        }
    }

    #[test]
    fn test_simplify() {
        let expr = SymbolicExpr::Add(
            Box::new(SymbolicExpr::Constant(2.0)),
            Box::new(SymbolicExpr::Constant(3.0)),
        );
        assert_eq!(expr.simplify(), SymbolicExpr::Constant(5.0));
    }

    #[test]
    fn test_derivative() {
        let expr = SymbolicExpr::Variable("x".to_string());
        assert_eq!(expr.derivative("x"), SymbolicExpr::One);
        assert_eq!(expr.derivative("y"), SymbolicExpr::Zero);
    }
