// Unit tests for electrical::power_analysis
use computational_engine::electrical::power_analysis::*;

use super::*;

    #[test]
    fn test_power_calculations() {
        let p = real_power(120.0, 10.0, 0.8);
        assert_eq!(p, 960.0); // 960W

        let s = apparent_power(120.0, 10.0);
        assert_eq!(s, 1200.0); // 1200VA

        let pf = power_factor_from_powers(960.0, 1200.0);
        assert_eq!(pf, 0.8);
    }

    #[test]
    fn test_three_phase() {
        let p = three_phase_power(480.0, 100.0, 1.0);
        assert!((p - 83138.0).abs() < 1.0);
    }

    #[test]
    fn test_voltage_conversions() {
        let v_phase = line_to_phase_voltage(480.0);
        assert!((v_phase - 277.0).abs() < 1.0);

        let v_line = phase_to_line_voltage(277.0);
        assert!((v_line - 480.0).abs() < 1.0);
    }
