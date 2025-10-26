// Unit tests for electrical::filter_design
use computational_engine::electrical::filter_design::*;

use super::*;

    #[test]
    fn test_rc_lowpass_design() {
        let (r, c) = design_rc_lowpass(1000.0, Some(10_000.0), None);
        assert_eq!(r, 10_000.0);
        // f_c = 1/(2πRC) => C = 1/(2π*1000*10000) ≈ 15.9 nF
        assert!((c - 15.9e-9).abs() < 1e-9);
    }

    #[test]
    fn test_butterworth() {
        let (r1, r2, c1, c2) = design_butterworth_lowpass(1000.0, 10_000.0);
        assert_eq!(r1, r2); // Butterworth has equal resistors
        assert_eq!(c1, c2); // And equal capacitors
    }
