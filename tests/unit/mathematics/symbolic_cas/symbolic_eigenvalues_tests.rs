// Unit tests for mathematics::symbolic_cas::symbolic_eigenvalues
use computational_engine::analyze::symbolic::symbolic_eigenvalues::*;

use super::*;

    #[test]
    fn test_characteristic_polynomial() {
        let mat = SymbolicMatrix::new(vec![
            vec![Expr::num(1), Expr::num(2)],
            vec![Expr::num(3), Expr::num(4)],
        ])
        .unwrap();

        let char_poly = characteristic_polynomial(&mat).unwrap();
        println!("Characteristic polynomial: {}", char_poly);
    }

    #[test]
    fn test_eigenvalues_2x2_numeric() {
        let mat = SymbolicMatrix::new(vec![
            vec![Expr::num(3), Expr::num(1)],
            vec![Expr::num(1), Expr::num(3)],
        ])
        .unwrap();

        let eigenvalues = eigenvalues_2x2(&mat).unwrap();
        assert_eq!(eigenvalues.len(), 2);
        println!("Eigenvalue 1: {}", eigenvalues[0]);
        println!("Eigenvalue 2: {}", eigenvalues[1]);
    }

    #[test]
    fn test_eigenvalues_2x2_symbolic() {
        let mat = SymbolicMatrix::new(vec![
            vec![Expr::sym("a"), Expr::sym("b")],
            vec![Expr::sym("b"), Expr::sym("a")],
        ])
        .unwrap();

        let eigenvalues = eigenvalues_2x2(&mat).unwrap();
        println!("Symbolic eigenvalue 1: {}", eigenvalues[0]);
        println!("Symbolic eigenvalue 2: {}", eigenvalues[1]);
    }

    #[test]
    fn test_matrix_inverse() {
        let mat = SymbolicMatrix::new(vec![
            vec![Expr::num(1), Expr::num(2)],
            vec![Expr::num(3), Expr::num(4)],
        ])
        .unwrap();

        let inv = matrix_inverse(&mat).unwrap();
        println!("Inverse matrix:\n{}", inv);
    }

    #[test]
    fn test_symbolic_inverse() {
        let mat = SymbolicMatrix::new(vec![
            vec![Expr::sym("a"), Expr::sym("b")],
            vec![Expr::sym("c"), Expr::sym("d")],
        ])
        .unwrap();

        let inv = matrix_inverse(&mat).unwrap();
        println!("Symbolic inverse:\n{}", inv);
    }
