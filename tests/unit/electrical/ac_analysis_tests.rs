// Unit tests for electrical::ac_analysis
use computational_engine::electrical::ac_analysis::*;

use super::*;

    #[test]
    fn test_rms_conversions() {
        assert!((peak_to_rms(10.0) - 7.071).abs() < 0.001);
        assert!((rms_to_peak(7.071) - 10.0).abs() < 0.01);
    }

    #[test]
    fn test_reactance() {
        let xc = capacitive_reactance(60.0, 10e-6);
        assert!((xc - 265.26).abs() < 0.1);

        let xl = inductive_reactance(60.0, 0.1);
        assert!((xl - 37.699).abs() < 0.01);
    }

    #[test]
    fn test_impedance() {
        let z_r = resistor_impedance(100.0);
        assert_eq!(z_r.re, 100.0);
        assert_eq!(z_r.im, 0.0);

        let z_c = capacitor_impedance(60.0, 10e-6);
        assert!(z_c.re.abs() < 0.01);
        assert!(z_c.im < 0.0); // Capacitive reactance is negative
    }

    #[test]
    fn test_resonance() {
        let f0 = resonant_frequency(0.1, 10e-6);
        assert!((f0 - 159.15).abs() < 0.1);

        let q = quality_factor(10.0, 0.1, 10e-6);
        assert!(q > 0.0);
    }

    #[test]
    fn test_time_constants() {
        let tau_rc = rc_time_constant(1000.0, 10e-6);
        assert_eq!(tau_rc, 0.01); // 10ms

        let tau_rl = rl_time_constant(0.1, 10.0);
        assert_eq!(tau_rl, 0.01); // 10ms
    }

    #[test]
    fn test_transfer_functions() {
        let cutoff = 1000.0;
        let h_low = lowpass_transfer_function(cutoff, cutoff);
        assert!((h_low - 0.707).abs() < 0.01); // -3dB point

        let h_high = highpass_transfer_function(cutoff, cutoff);
        assert!((h_high - 0.707).abs() < 0.01);
    }

    #[test]
    fn test_db_conversions() {
        assert!((magnitude_to_db(1.0) - 0.0).abs() < 0.001);
        assert!((magnitude_to_db(0.707) - (-3.0)).abs() < 0.1);
        assert!((db_to_magnitude(0.0) - 1.0).abs() < 0.001);
    }
