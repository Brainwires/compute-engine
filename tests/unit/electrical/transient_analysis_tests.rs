// Unit tests for electrical::transient_analysis
use computational_engine::electrical::transient_analysis::*;

use super::*;

    #[test]
    fn test_rc_transients() {
        let tau = 0.01; // 10ms
        let v_s = 12.0;

        // At t = τ, voltage should be ≈ 63.2% of final
        let v_at_tau = rc_charging_voltage(v_s, tau, tau);
        assert!((v_at_tau / v_s - 0.632).abs() < 0.01);

        // At t = 5τ, voltage should be ≈ 99.3% of final
        let v_at_5tau = rc_charging_voltage(v_s, 5.0 * tau, tau);
        assert!((v_at_5tau / v_s - 0.993).abs() < 0.01);

        // Discharging from 12V
        let v_discharge = rc_discharging_voltage(12.0, tau, tau);
        assert!((v_discharge - 4.416).abs() < 0.01); // 12 * e^-1
    }

    #[test]
    fn test_rl_transients() {
        let tau = 0.01; // 10ms
        let v_s = 12.0;
        let r = 10.0;

        let i_at_tau = rl_rising_current(v_s, r, tau, tau);
        let i_final = v_s / r;
        assert!((i_at_tau / i_final - 0.632).abs() < 0.01);
    }

    #[test]
    fn test_settling_time() {
        let tau = 0.01;
        let t_95 = rc_settling_time(tau, 95.0);
        assert!((t_95 / tau - 3.0).abs() < 0.1); // ~3τ for 95%

        let t_99 = rc_settling_time(tau, 99.0);
        assert!((t_99 / tau - 4.6).abs() < 0.1); // ~4.6τ for 99%
    }

    #[test]
    fn test_overshoot() {
        let os = percent_overshoot(0.5); // ζ = 0.5
        assert!((os - 16.3).abs() < 0.5); // ~16% overshoot

        let os2 = percent_overshoot(0.707); // ζ = 1/√2
        assert!((os2 - 4.3).abs() < 0.5); // ~4% overshoot
    }
