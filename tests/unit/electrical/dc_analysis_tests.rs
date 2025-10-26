// Unit tests for electrical::dc_analysis
use computational_engine::electrical::dc_analysis::*;

use super::*;

    #[test]
    fn test_ohms_law() {
        assert_eq!(ohms_law_voltage(2.0, 10.0), 20.0);
        assert_eq!(ohms_law_current(20.0, 10.0), 2.0);
        assert_eq!(ohms_law_resistance(20.0, 2.0), 10.0);
    }

    #[test]
    fn test_series_parallel() {
        assert_eq!(series_resistance(&[10.0, 20.0, 30.0]), 60.0);
        assert_eq!(parallel_resistance(&[10.0, 10.0]), 5.0);
        assert!((parallel_resistance(&[10.0, 20.0, 30.0]) - 5.454545).abs() < 0.001);
    }

    #[test]
    fn test_dividers() {
        assert_eq!(voltage_divider(12.0, 10.0, 5.0), 4.0);
        assert_eq!(current_divider(6.0, 10.0, 5.0), 2.0);
    }

    #[test]
    fn test_power() {
        assert_eq!(power_dissipation(12.0, 2.0), 24.0);
        assert_eq!(power_from_voltage_resistance(12.0, 6.0), 24.0);
        assert_eq!(power_from_current_resistance(2.0, 6.0), 24.0);
    }

    #[test]
    fn test_energy_storage() {
        assert!((capacitor_energy(1e-6, 10.0) - 5e-5).abs() < 1e-10); // 1μF at 10V
        assert!((inductor_energy(1e-3, 2.0) - 2e-3).abs() < 1e-10); // 1mH with 2A
    }
