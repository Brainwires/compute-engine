// Unit tests for mathematics::symbolic_cas::christoffel
use computational_engine::mathematics::symbolic_cas::christoffel::*;

use super::*;
    use crate::mathematics::symbolic_cas::{euclidean_metric, minkowski_2d};

    #[test]
    fn test_flat_space_christoffel() {
        // In flat Euclidean space, all Christoffel symbols should be zero
        let metric = euclidean_metric(2);
        let coords = vec!["x", "y"];

        let christoffel = christoffel_symbols(&metric, &coords).unwrap();

        // Check that all components are zero
        let data = christoffel.data();
        for component in data {
            assert_eq!(component, &Expr::num(0));
        }
    }

    #[test]
    fn test_minkowski_christoffel() {
        // Minkowski spacetime in standard coordinates also has zero Christoffel symbols
        let metric = minkowski_2d(true).unwrap();
        let coords = vec!["t", "x"];

        let christoffel = christoffel_symbols(&metric, &coords).unwrap();

        println!("Minkowski Christoffel symbols:");
        println!("{}", christoffel);

        // All should be zero for Minkowski in standard coordinates
        let data = christoffel.data();
        for component in data {
            assert_eq!(component, &Expr::num(0));
        }
    }

    #[test]
    fn test_christoffel_first_kind() {
        let metric = euclidean_metric(2);
        let coords = vec!["x", "y"];

        let christoffel = christoffel_first_kind(&metric, &coords).unwrap();

        // All should be zero for flat space
        let data = christoffel.data();
        for component in data {
            assert_eq!(component, &Expr::num(0));
        }
    }

    #[test]
    fn test_schwarzschild_christoffel() {
        use crate::mathematics::symbolic_cas::schwarzschild_metric;

        let metric = schwarzschild_metric().unwrap();
        let coords = vec!["t", "r", "θ", "φ"];

        let christoffel = christoffel_symbols(&metric, &coords).unwrap();

        println!("Schwarzschild Christoffel symbols (rank-3 tensor):");
        println!("{}", christoffel);

        // Should have non-zero components for curved spacetime
        assert_eq!(christoffel.rank(), 3);
    }

    #[test]
    fn test_flat_riemann() {
        // Flat space should have zero Riemann tensor
        let metric = euclidean_metric(2);
        let coords = vec!["x", "y"];

        let riemann = riemann_tensor(&metric, &coords).unwrap();

        println!("Euclidean Riemann tensor:");
        println!("{}", riemann);

        // All components should be zero
        assert_eq!(riemann.rank(), 4);
        for component in riemann.data() {
            assert_eq!(component, &Expr::num(0));
        }
    }

    #[test]
    fn test_flat_ricci() {
        // Flat space should have zero Ricci tensor
        let metric = euclidean_metric(2);
        let coords = vec!["x", "y"];

        let ricci = ricci_tensor(&metric, &coords).unwrap();

        println!("Euclidean Ricci tensor:");
        println!("{}", ricci);

        // All components should be zero
        for i in 0..2 {
            for j in 0..2 {
                assert_eq!(ricci.get(i, j).unwrap(), &Expr::num(0));
            }
        }
    }

    #[test]
    fn test_ricci_scalar_flat() {
        // Flat space should have zero Ricci scalar
        let metric = euclidean_metric(2);
        let coords = vec!["x", "y"];

        let scalar = ricci_scalar(&metric, &coords).unwrap();

        println!("Euclidean Ricci scalar: {}", scalar);
        assert_eq!(scalar, Expr::num(0));
    }

    #[test]
    fn test_einstein_tensor_flat() {
        // Flat space should have zero Einstein tensor
        let metric = euclidean_metric(2);
        let coords = vec!["x", "y"];

        let einstein = einstein_tensor(&metric, &coords).unwrap();

        println!("Euclidean Einstein tensor:");
        println!("{}", einstein);

        // All components should be zero
        for i in 0..2 {
            for j in 0..2 {
                assert_eq!(einstein.get(i, j).unwrap(), &Expr::num(0));
            }
        }
    }
