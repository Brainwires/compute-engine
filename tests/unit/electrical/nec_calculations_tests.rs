// Unit tests for electrical::nec_calculations
use computational_engine::electrical::nec_calculations::*;

use super::*;

    #[test]
    fn test_temperature_correction() {
        let ampacity = temperature_corrected_ampacity(100.0, 40.0, 75);
        assert_eq!(ampacity, 88.0); // 100A * 0.88
    }

    #[test]
    fn test_conduit_fill() {
        let ampacity = conduit_fill_ampacity(100.0, 5);
        assert_eq!(ampacity, 80.0); // 100A * 0.8 (4-6 conductors)
    }

    #[test]
    fn test_voltage_drop() {
        // 20A load, 100ft run, 12 AWG copper (1.93 Ω/1000ft), single-phase
        let vd = voltage_drop(20.0, 100.0, 1.93, false);
        assert!((vd - 7.72).abs() < 0.01); // 2 * 20 * (1.93 * 100/1000) = 7.72V

        let percent = percent_voltage_drop(vd, 120.0);
        assert!((percent - 6.43).abs() < 0.01); // ~6.43%
    }

    #[test]
    fn test_branch_circuit() {
        let min_ampacity = branch_circuit_ampacity(20.0, 5.0);
        assert_eq!(min_ampacity, 30.0); // (20 * 1.25) + 5
    }

    #[test]
    fn test_service_size() {
        let load_va = 20_000.0;
        let min_amps = minimum_service_size(load_va, 240.0);
        assert!((min_amps - 83.33).abs() < 0.1);

        let standard = standard_service_size(min_amps);
        assert_eq!(standard, 100); // Rounds up to 100A
    }

    #[test]
    fn test_motor_calculations() {
        let fla = motor_full_load_current(5.0, 230, false);
        assert_eq!(fla, 25.0); // 5 HP * 5.0

        let ocpd = motor_ocpd_size(25.0, "ac");
        assert_eq!(ocpd, 62.5); // 25A * 2.5
    }
