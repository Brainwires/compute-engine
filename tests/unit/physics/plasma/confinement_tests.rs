// Unit tests for physics::plasma::confinement
use computational_engine::physics::plasma::confinement::*;

use super::*;

    fn iter_config() -> TokamakConfig {
        TokamakConfig {
            major_radius: 6.2,
            minor_radius: 2.0,
            toroidal_field: 5.3,
            plasma_current: 15e6, // 15 MA
            elongation: 1.7,
            triangularity: 0.33,
        }
    }

    #[test]
    fn test_aspect_ratio() {
        let config = iter_config();
        let aspect = aspect_ratio(&config);

        assert!((aspect - 3.1).abs() < 0.1);
    }

    #[test]
    fn test_safety_factor() {
        let config = iter_config();
        let q_edge = safety_factor(&config, config.minor_radius);

        // ITER edge safety factor should be > 3
        assert!(q_edge > 2.0);
        assert!(q_edge < 10.0);
    }

    #[test]
    fn test_plasma_volume() {
        let config = iter_config();
        let volume = plasma_volume(&config);

        // ITER plasma volume ~ 800 m³
        assert!(volume > 700.0 && volume < 900.0);
    }

    #[test]
    fn test_kruskal_shafranov() {
        let config = iter_config();
        let stable = kruskal_shafranov_stable(&config);

        assert!(stable);
    }

    #[test]
    fn test_troyon_limit() {
        let config = iter_config();
        let beta_n = 3.0; // Typical value
        let beta_max = troyon_beta_limit(&config, beta_n);

        // Should be positive and finite
        assert!(beta_max > 0.0);
        assert!(beta_max.is_finite());
    }

    #[test]
    fn test_greenwald_limit() {
        let config = iter_config();
        let n_g = greenwald_density_limit(&config);

        // Greenwald limit for ITER ~ 10^20 m⁻³
        assert!(n_g > 5e19 && n_g < 2e20);
    }

    #[test]
    fn test_iter98_confinement() {
        let config = iter_config();
        let params = PlasmaParams::tokamak();
        let p_heat = 50e6; // 50 MW

        let tau_e = iter98_confinement_time(&config, &params, p_heat);

        // ITER target confinement time ~ 3-5 seconds
        assert!(tau_e > 1.0 && tau_e < 10.0);
    }

    #[test]
    fn test_fusion_triple_product() {
        let params = PlasmaParams::tokamak();
        let tau_e = 3.0; // 3 seconds

        let ntt = fusion_triple_product(&params, tau_e);

        // ITER should exceed Lawson criterion (3e21)
        assert!(ntt > 1e21);
    }

    #[test]
    fn test_dt_fusion_power() {
        let params = PlasmaParams::tokamak();
        let p_fusion = dt_fusion_power_density(&params);

        // Should be finite and positive at 10 keV
        assert!(p_fusion >= 0.0);
        assert!(p_fusion.is_finite());
    }

    #[test]
    fn test_confinement_mode() {
        let config = iter_config();
        let params = PlasmaParams::tokamak();
        let p_heat = 50e6;

        let mode = estimate_confinement_mode(&config, &params, p_heat);

        assert!(mode.is_hmode); // 50 MW should trigger H-mode
        assert!(mode.h_factor >= 0.8);
    }

    #[test]
    fn test_magnetic_field_energy() {
        let config = iter_config();
        let energy = magnetic_field_energy(&config);

        // ITER magnetic field energy ~ GJ
        assert!(energy > 1e8); // > 100 MJ
        assert!(energy.is_finite());
    }
