// Unit tests for mathematics::symbolic_cas::symbolic_matrix
use computational_engine::mathematics::symbolic_cas::symbolic_matrix::*;

use super::*;
    use crate::mathematics::symbolic_cas::Expr;

    #[test]
    fn test_identity_matrix() {
        let identity = SymbolicMatrix::identity(3);
        assert_eq!(identity.rows(), 3);
        assert_eq!(identity.cols(), 3);
        assert_eq!(identity.get(0, 0), Some(&Expr::num(1)));
        assert_eq!(identity.get(0, 1), Some(&Expr::num(0)));
    }

    #[test]
    fn test_transpose() {
        let data = vec![
            vec![Expr::num(1), Expr::num(2)],
            vec![Expr::num(3), Expr::num(4)],
        ];
        let mat = SymbolicMatrix::new(data).unwrap();
        let transposed = mat.transpose();

        assert_eq!(transposed.rows(), 2);
        assert_eq!(transposed.cols(), 2);
        assert_eq!(transposed.get(0, 1), Some(&Expr::num(3)));
    }

    #[test]
    fn test_matrix_mul() {
        let a = SymbolicMatrix::new(vec![
            vec![Expr::num(1), Expr::num(2)],
            vec![Expr::num(3), Expr::num(4)],
        ])
        .unwrap();

        let b = SymbolicMatrix::new(vec![
            vec![Expr::num(5), Expr::num(6)],
            vec![Expr::num(7), Expr::num(8)],
        ])
        .unwrap();

        let c = a.mul(&b).unwrap();
        assert_eq!(c.rows(), 2);
        assert_eq!(c.cols(), 2);
    }

    #[test]
    fn test_determinant_2x2() {
        let mat = SymbolicMatrix::new(vec![
            vec![Expr::num(1), Expr::num(2)],
            vec![Expr::num(3), Expr::num(4)],
        ])
        .unwrap();

        let det = mat.determinant().unwrap();
        // det = 1*4 - 2*3 = -2
        println!("Determinant: {}", det);
    }

    #[test]
    fn test_symbolic_determinant() {
        let mat = SymbolicMatrix::new(vec![
            vec![Expr::sym("a"), Expr::sym("b")],
            vec![Expr::sym("c"), Expr::sym("d")],
        ])
        .unwrap();

        let det = mat.determinant().unwrap();
        // det = a*d - b*c
        println!("Symbolic determinant: {}", det);
    }
