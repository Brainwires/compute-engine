// Unit tests for mathematics::symbolic_cas::mod
use computational_engine::mathematics::symbolic_cas::mod::*;

use super::*;

    #[test]
    fn test_expand() {
        let result = expand("(x + 1)^2").unwrap();
        assert!(result.expression.len() > 0);
    }

    #[test]
    fn test_differentiate() {
        let result = differentiate("x^2 + 3*x", "x", None).unwrap();
        assert!(result.expression.len() > 0);
    }

    #[test]
    fn test_evaluate() {
        let mut values = HashMap::new();
        values.insert("x".to_string(), 2.0);
        let result = evaluate_at("x^2 + 3*x + 2", &values).unwrap();
        assert!((result - 12.0).abs() < 1e-5);
    }

    #[test]
    fn test_substitute() {
        let mut rules = HashMap::new();
        rules.insert("x".to_string(), "y".to_string());
        let result = substitute("x + 1", &rules).unwrap();
        assert!(result.expression.contains("y"));
    }
