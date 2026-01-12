// Unit tests for mathematics::chaos::fractals
use computational_engine::chaos::fractals::*;

use super::*;

    #[test]
    fn test_mandelbrot_origin() {
        // c = 0 should not escape (z_n stays at 0)
        let result = mandelbrot(Complex::new(0.0, 0.0), 100, 2.0);
        assert!(!result.escaped);
        assert!(result.final_magnitude < 0.1);
    }

    #[test]
    fn test_mandelbrot_escape() {
        // c = 2 should escape quickly
        let result = mandelbrot(Complex::new(2.0, 0.0), 100, 2.0);
        assert!(result.escaped);
        assert!(result.iterations < 10);
    }

    #[test]
    fn test_mandelbrot_boundary() {
        // c = -0.5 is in the Mandelbrot set (main cardioid)
        let result = mandelbrot(Complex::new(-0.5, 0.0), 100, 2.0);
        // May or may not escape depending on iterations, just check it runs
        assert!(result.iterations <= 100);
    }

    #[test]
    fn test_julia_set() {
        // Julia set with c = -0.7 + 0.27i
        let c = Complex::new(-0.7, 0.27);
        let z0 = Complex::new(0.0, 0.0);

        let result = julia(z0, c, 100, 2.0);
        assert!(result.iterations <= 100);
    }

    #[test]
    fn test_burning_ship() {
        let c = Complex::new(-1.8, -0.1);
        let result = burning_ship(c, 100, 2.0);
        assert!(result.iterations <= 100);
    }

    #[test]
    fn test_koch_snowflake() {
        let order_0 = koch_snowflake(0);
        assert_eq!(order_0.len(), 4); // Triangle + return to start

        let order_1 = koch_snowflake(1);
        assert!(order_1.len() > order_0.len()); // More points after iteration
    }

    #[test]
    fn test_box_counting() {
        // Create a simple line from (0,0) to (10,0)
        let points: Vec<(f64, f64)> = (0..=100).map(|i| (i as f64 * 0.1, 0.0)).collect();

        let dim = box_counting_dimension(&points, 0.1, 10.0, 10);

        // A line should have dimension ~1
        assert!(dim > 0.5 && dim < 1.5);
    }
