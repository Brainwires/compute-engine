// Unit tests for specialized::control_theory::transfer_function
use computational_engine::specialized::control_theory::transfer_function::*;

use super::*;

    #[test]
    fn test_transfer_function_creation() {
        // G(s) = 1 / (s + 1)
        let tf = TransferFunction::new(vec![1.0], vec![1.0, 1.0]);
        assert!(tf.is_ok());
    }

    #[test]
    fn test_dc_gain() {
        // G(s) = 2 / (s + 1) → DC gain = 2
        let tf = TransferFunction::new(vec![2.0], vec![1.0, 1.0]).unwrap();
        assert!((tf.dc_gain() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_frequency_response() {
        // G(s) = 1 / (s + 1) at s = j*0
        let tf = TransferFunction::new(vec![1.0], vec![1.0, 1.0]).unwrap();
        let h = tf.frequency_response(0.0);
        assert!((h.re - 1.0).abs() < 1e-10);
        assert!(h.im.abs() < 1e-10);
    }

    #[test]
    fn test_poles_linear() {
        // G(s) = 1 / (s + 2) → pole at s = -2
        let tf = TransferFunction::new(vec![1.0], vec![1.0, 2.0]).unwrap();
        let poles = tf.poles();
        assert_eq!(poles.len(), 1);
        assert!((poles[0].re + 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_stability() {
        // Stable: pole at -1
        let stable = TransferFunction::new(vec![1.0], vec![1.0, 1.0]).unwrap();
        assert!(stable.is_stable());

        // Unstable: pole at +1
        let unstable = TransferFunction::new(vec![1.0], vec![1.0, -1.0]).unwrap();
        assert!(!unstable.is_stable());
    }

    #[test]
    fn test_series_connection() {
        // G1(s) = 1, G2(s) = 2 → G1*G2 = 2
        let g1 = TransferFunction::new(vec![1.0], vec![1.0]).unwrap();
        let g2 = TransferFunction::new(vec![2.0], vec![1.0]).unwrap();
        let series = g1.series(&g2);
        assert_eq!(series.dc_gain(), 2.0);
    }

    #[test]
    fn test_convolve() {
        // (x + 1) * (x + 2) = x² + 3x + 2
        let result = convolve(&[1.0, 1.0], &[1.0, 2.0]);
        assert_eq!(result, vec![1.0, 3.0, 2.0]);
    }
