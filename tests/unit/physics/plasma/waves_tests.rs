// Unit tests for physics::plasma::waves
use computational_engine::compute::physics::plasma::waves::*;

use super::*;
    use crate::compute::physics::plasma::{plasma_frequency, electron_cyclotron_frequency, C};

    #[test]
    fn test_langmuir_wave() {
        let params = PlasmaParams::tokamak();
        let k = 1e3; // m⁻¹

        let wave = langmuir_wave(k, &params);

        // Frequency should be close to plasma frequency
        let omega_pe = plasma_frequency(params.n_e);
        assert!(wave.omega >= omega_pe);
        assert!(wave.v_phase > 0.0);
    }

    #[test]
    fn test_alfven_wave() {
        let params = PlasmaParams::tokamak();
        let k = 1e3;

        let wave = alfven_wave(k, &params);

        // Alfvén wave is non-dispersive
        assert!((wave.v_phase - wave.v_group).abs() < 1e-6);
        assert!(wave.omega > 0.0);
    }

    #[test]
    fn test_ion_acoustic_wave() {
        let params = PlasmaParams::tokamak();
        let k = 1e3;

        let wave = ion_acoustic_wave(k, &params);

        // Ion acoustic wave is also non-dispersive
        assert!((wave.v_phase - wave.v_group).abs() / wave.v_phase < 0.01);
        assert!(wave.omega > 0.0);
    }

    #[test]
    fn test_hybrid_frequencies() {
        let params = PlasmaParams::tokamak();

        let omega_uh = upper_hybrid_frequency(&params);
        let omega_lh = lower_hybrid_frequency(&params);

        let omega_pe = plasma_frequency(params.n_e);
        let omega_ce = electron_cyclotron_frequency(params.b_field);

        // Upper hybrid should be higher than both
        assert!(omega_uh > omega_pe);
        assert!(omega_uh > omega_ce);

        // Lower hybrid should be lower
        assert!(omega_lh < omega_uh);
    }

    #[test]
    fn test_two_stream_instability() {
        let n_b = 1e18;
        let n_0 = 1e20;

        let gamma = two_stream_growth_rate(n_b, n_0);

        // Growth rate should be positive
        assert!(gamma > 0.0);
        assert!(gamma.is_finite());
    }

    #[test]
    fn test_weibel_instability() {
        let params = PlasmaParams::tokamak();
        let t_perp = 20000.0; // Hot perpendicular
        let t_parallel = 10000.0;

        let gamma = weibel_growth_rate(&params, t_perp, t_parallel);

        // Should be unstable (positive growth rate)
        assert!(gamma > 0.0);
    }

    #[test]
    fn test_weibel_stable() {
        let params = PlasmaParams::tokamak();
        let t_perp = 10000.0;
        let t_parallel = 20000.0; // Hot parallel

        let gamma = weibel_growth_rate(&params, t_perp, t_parallel);

        // Should be stable (zero growth rate)
        assert_eq!(gamma, 0.0);
    }

    #[test]
    fn test_em_wave_dispersion() {
        let params = PlasmaParams::tokamak();
        let k = 1e6;

        let wave = em_wave_in_plasma(k, &params);

        // At high k, should approach light speed
        assert!(wave.v_phase.is_finite());
        assert!(wave.v_group.is_finite());
        assert!(wave.v_group < C);
    }

    #[test]
    fn test_cutoff_frequency() {
        let params = PlasmaParams::tokamak();
        let f_cutoff = cutoff_frequency(&params);

        // Should equal plasma frequency
        let omega_pe = plasma_frequency(params.n_e);
        assert_eq!(f_cutoff, omega_pe);
    }
