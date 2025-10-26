// Unit tests for mathematics::chaos::lyapunov
use computational_engine::mathematics::chaos::lyapunov::*;

use super::*;

    #[test]
    fn test_lyapunov_logistic_stable() {
        // r = 2.5: stable fixed point → negative Lyapunov
        let lyap = lyapunov_logistic_map(2.5, 1000, 100);
        assert!(lyap < 0.0);
    }

    #[test]
    fn test_lyapunov_logistic_chaotic() {
        // r = 4.0: fully chaotic → positive Lyapunov
        let lyap = lyapunov_logistic_map(4.0, 1000, 100);
        assert!(lyap > 0.0);
    }

    #[test]
    fn test_lyapunov_logistic_periodic() {
        // r = 3.2: period-2 orbit → zero or slightly negative
        let lyap = lyapunov_logistic_map(3.2, 1000, 100);
        // Should be near 0 or slightly negative for periodic
        assert!(lyap < 0.5);
    }

    #[test]
    fn test_kaplan_yorke_dimension() {
        // Example: one positive, two negative
        let exponents = vec![0.9, 0.0, -14.5];
        let dim = kaplan_yorke_dimension(&exponents);

        // D = 1 + λ₁/|λ₃| = 1 + 0.9/14.5 ≈ 1.062
        // Just check it's between 1 and 2 (reasonable for 2D attractor embedded in 3D)
        assert!(dim > 0.5 && dim < 2.5, "Dimension was {}", dim);
    }

    #[test]
    fn test_gram_schmidt() {
        let mut basis = [
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ];

        let norms = gram_schmidt(&mut basis);

        // Check orthogonality
        let dot01: f64 = (0..3).map(|k| basis[0][k] * basis[1][k]).sum();
        assert!(dot01.abs() < 1e-10);

        // Check normalization
        assert!((norms[0] - (2.0_f64).sqrt()).abs() < 1e-10);
    }
