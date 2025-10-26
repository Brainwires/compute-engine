// Unit tests for mathematics::symbolic_cas::symbolic_tensor
use computational_engine::mathematics::symbolic_cas::symbolic_tensor::*;

use super::*;

    #[test]
    fn test_scalar_tensor() {
        let scalar = SymbolicTensor::scalar(Expr::num(42));
        assert_eq!(scalar.rank(), 0);
        assert_eq!(scalar.get(&[]).unwrap(), &Expr::num(42));
    }

    #[test]
    fn test_vector_tensor() {
        let vector = SymbolicTensor::vector(
            vec![Expr::sym("x"), Expr::sym("y"), Expr::sym("z")],
            IndexType::Contravariant,
        )
        .unwrap();

        assert_eq!(vector.rank(), 1);
        assert_eq!(vector.get(&[0]).unwrap(), &Expr::sym("x"));
        assert_eq!(vector.get(&[2]).unwrap(), &Expr::sym("z"));
    }

    #[test]
    fn test_matrix_to_tensor() {
        let mat = SymbolicMatrix::new(vec![
            vec![Expr::num(1), Expr::num(2)],
            vec![Expr::num(3), Expr::num(4)],
        ])
        .unwrap();

        let tensor =
            SymbolicTensor::from_matrix(&mat, [IndexType::Contravariant, IndexType::Covariant]);
        assert_eq!(tensor.rank(), 2);
        assert_eq!(tensor.get(&[0, 1]).unwrap(), &Expr::num(2));
    }

    #[test]
    fn test_tensor_contraction() {
        let mat = SymbolicMatrix::new(vec![
            vec![Expr::sym("a"), Expr::sym("b")],
            vec![Expr::sym("c"), Expr::sym("d")],
        ])
        .unwrap();

        let tensor =
            SymbolicTensor::from_matrix(&mat, [IndexType::Contravariant, IndexType::Covariant]);
        let trace = tensor.contract(0, 1).unwrap();

        assert_eq!(trace.rank(), 0);
        println!("Trace: {}", trace);
    }

    #[test]
    fn test_outer_product() {
        let v1 = SymbolicTensor::vector(
            vec![Expr::sym("a"), Expr::sym("b")],
            IndexType::Contravariant,
        )
        .unwrap();

        let v2 = SymbolicTensor::vector(vec![Expr::sym("c"), Expr::sym("d")], IndexType::Covariant)
            .unwrap();

        let tensor = v1.outer_product(&v2).unwrap();
        assert_eq!(tensor.rank(), 2);
        println!("Outer product:\n{}", tensor);
    }

    #[test]
    fn test_raise_lower_index() {
        // Metric tensor (Minkowski for 2D)
        let metric = SymbolicMatrix::new(vec![
            vec![Expr::num(1), Expr::num(0)],
            vec![Expr::num(0), Expr::num(-1)],
        ])
        .unwrap();

        let covariant_vector = SymbolicTensor::vector(
            vec![Expr::sym("v_0"), Expr::sym("v_1")],
            IndexType::Covariant,
        )
        .unwrap();

        let contravariant = covariant_vector.raise_index(&metric, 0).unwrap();
        println!("Raised vector:\n{}", contravariant);
    }
