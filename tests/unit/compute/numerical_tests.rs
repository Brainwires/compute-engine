// Unit tests for mathematics::numerical
use computational_engine::compute::numerical::*;

use super::*;

    #[test]
    fn test_evaluate_numeric() {
        let expr = Expr::add(Expr::mul(Expr::num(2), Expr::sym("x")), Expr::num(3));

        let mut values = HashMap::new();
        values.insert("x".to_string(), 5.0);

        let result = evaluate_numeric(&expr, &values).unwrap();
        assert!((result - 13.0).abs() < 1e-10);
    }

    #[test]
    fn test_matrix_exponential() {
        // Test exp(0) = I
        let zero_matrix = SymbolicMatrix::zeros(2, 2);
        let result = matrix_exponential(&zero_matrix, &HashMap::new(), Some(10)).unwrap();

        // Should be identity
        assert!((result[0][0] - 1.0).abs() < 1e-10);
        assert!((result[1][1] - 1.0).abs() < 1e-10);
        assert!(result[0][1].abs() < 1e-10);
    }

    #[test]
    fn test_eigenvalues_numeric() {
        // Test 2x2 matrix with known eigenvalues
        let matrix = SymbolicMatrix::new(vec![
            vec![Expr::num(2), Expr::num(0)],
            vec![Expr::num(0), Expr::num(3)],
        ])
        .unwrap();

        let eigenvalues = eigenvalues_numeric(&matrix, &HashMap::new(), Some(50)).unwrap();

        println!("Eigenvalues: {:?}", eigenvalues);

        // Should be approximately [2, 3]
        assert!(eigenvalues.iter().any(|&ev| (ev - 2.0).abs() < 1e-5));
        assert!(eigenvalues.iter().any(|&ev| (ev - 3.0).abs() < 1e-5));
    }

    #[test]
    fn test_matrix_power() {
        let matrix = SymbolicMatrix::new(vec![
            vec![Expr::num(2), Expr::num(0)],
            vec![Expr::num(0), Expr::num(3)],
        ])
        .unwrap();

        let result = matrix_power_numeric(&matrix, 2, &HashMap::new()).unwrap();

        // [[2, 0], [0, 3]]^2 = [[4, 0], [0, 9]]
        assert!((result[0][0] - 4.0).abs() < 1e-10);
        assert!((result[1][1] - 9.0).abs() < 1e-10);
    }

    #[test]
    fn test_evaluate_functions() {
        let mut values = HashMap::new();
        values.insert("x".to_string(), std::f64::consts::PI / 2.0);

        let sin_expr = Expr::func("sin", vec![Expr::sym("x")]);
        let result = evaluate_numeric(&sin_expr, &values).unwrap();

        assert!((result - 1.0).abs() < 1e-10);
    }
