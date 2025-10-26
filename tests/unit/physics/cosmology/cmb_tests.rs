// Unit tests for physics::cosmology::cmb
use computational_engine::physics::cosmology::cmb::*;

use super::*;

    #[test]
    fn test_planck_spectrum() {
        let t_cmb = 2.7255;

        // Peak should be at microwave frequencies (~160 GHz)
        let freq_peak = cmb_peak_frequency(t_cmb);
        assert!(freq_peak > 1e11 && freq_peak < 2e11);

        // Spectrum at peak should be non-zero
        let intensity = planck_spectrum(freq_peak, t_cmb);
        assert!(intensity > 0.0);
        assert!(intensity.is_finite());
    }

    #[test]
    fn test_cmb_peak_frequency() {
        let params = CosmologyParams::planck_2018();
        let freq = cmb_peak_frequency(params.t_cmb);

        // Should be ~160 GHz
        assert!(freq > 1.4e11 && freq < 1.8e11);
    }

    #[test]
    fn test_photon_number_density() {
        let params = CosmologyParams::planck_2018();
        let n_gamma = photon_number_density(params.t_cmb);

        // CMB photon density should be ~400 photons/cm³ = 4e8 photons/m³
        assert!(n_gamma > 1e8 && n_gamma < 1e9);
    }

    #[test]
    fn test_cmb_energy_density() {
        let params = CosmologyParams::planck_2018();
        let rho_gamma = cmb_energy_density(params.t_cmb);

        // Should be finite and positive
        assert!(rho_gamma > 0.0);
        assert!(rho_gamma.is_finite());
    }

    #[test]
    fn test_recombination() {
        let params = CosmologyParams::planck_2018();
        let recomb = recombination(&params);

        // Recombination at z ~ 1100
        assert!((recomb.redshift - 1100.0).abs() < 1.0);

        // Temperature ~ 3000 K
        assert!(recomb.temperature > 2900.0 && recomb.temperature < 3100.0);

        // Scale factor ~ 1/1101
        assert!(recomb.scale_factor > 0.0009 && recomb.scale_factor < 0.001);
    }

    #[test]
    fn test_matter_radiation_equality() {
        let params = CosmologyParams::planck_2018();
        let z_eq = matter_radiation_equality(&params);

        // Matter-radiation equality at z ~ 3400
        assert!(z_eq > 3000.0 && z_eq < 4000.0);
    }

    #[test]
    fn test_sound_horizon() {
        let params = CosmologyParams::planck_2018();
        let r_s = sound_horizon_recombination(&params);

        // Sound horizon should be positive and finite
        assert!(r_s > 0.0);
        assert!(r_s.is_finite());
    }

    #[test]
    fn test_first_acoustic_peak() {
        let params = CosmologyParams::planck_2018();
        let theta = first_acoustic_peak_angle(&params);

        // First peak should be a positive, finite angle
        assert!(theta > 0.0);
        assert!(theta.is_finite());
    }

    #[test]
    fn test_sachs_wolfe() {
        let potential = 1e-5; // Typical gravitational potential
        let dt_t = sachs_wolfe_amplitude(potential);

        // Temperature fluctuation should be tiny
        assert!(dt_t.abs() < 1e-10);
    }
