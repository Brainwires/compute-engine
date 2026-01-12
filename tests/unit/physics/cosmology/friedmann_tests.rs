// Unit tests for physics::cosmology::friedmann
use computational_engine::compute::physics::cosmology::friedmann::*;

use super::*;

    #[test]
    fn test_hubble_parameter_today() {
        let params = CosmologyParams::planck_2018();
        let h_today = hubble_parameter(Redshift(0.0), &params);
        let h0 = params.hubble_constant_si();

        // Should equal H0 today
        assert!((h_today - h0).abs() < 1e-10);
    }

    #[test]
    fn test_hubble_parameter_evolution() {
        let params = CosmologyParams::planck_2018();

        let h_today = hubble_parameter(Redshift(0.0), &params);
        let h_early = hubble_parameter(Redshift(1000.0), &params);

        // H(z) should be larger in the early universe
        assert!(h_early > h_today);
    }

    #[test]
    fn test_deceleration_parameter() {
        let params = CosmologyParams::planck_2018();

        let q_today = deceleration_parameter(Redshift(0.0), &params);

        // With dark energy, universe is accelerating (q < 0)
        assert!(q_today < 0.0);
    }

    #[test]
    fn test_lookback_time() {
        let params = CosmologyParams::planck_2018();

        let t_lookback = lookback_time(Redshift(1.0), &params);

        // Lookback time to z=1 should be several billion years
        assert!(t_lookback > 1e17); // >3 billion years in seconds
        assert!(t_lookback.is_finite());
    }

    #[test]
    fn test_comoving_distance() {
        let params = CosmologyParams::planck_2018();

        let d_c = comoving_distance(Redshift(1.0), &params);

        // Comoving distance to z=1 should be a few Gpc
        assert!(d_c > 0.0);
        assert!(d_c.is_finite());
    }

    #[test]
    fn test_luminosity_distance() {
        let params = CosmologyParams::planck_2018();

        let z = Redshift(1.0);
        let d_c = comoving_distance(z, &params);
        let d_l = luminosity_distance(z, &params);

        // D_L = (1+z) D_c
        assert!((d_l - 2.0 * d_c).abs() < 1e-6);
    }

    #[test]
    fn test_angular_diameter_distance() {
        let params = CosmologyParams::planck_2018();

        let z = Redshift(1.0);
        let d_c = comoving_distance(z, &params);
        let d_a = angular_diameter_distance(z, &params);

        // D_A = D_c / (1+z)
        assert!((d_a * 2.0 - d_c).abs() < 1e-6);
    }

    #[test]
    fn test_distance_relationship() {
        let params = CosmologyParams::planck_2018();
        let z = Redshift(2.0);

        let d_l = luminosity_distance(z, &params);
        let d_a = angular_diameter_distance(z, &params);

        // D_L = (1+z)² D_A
        let expected_ratio = (1.0 + z.0).powi(2);
        let actual_ratio = d_l / d_a;

        assert!((actual_ratio - expected_ratio).abs() < 1e-6);
    }

    #[test]
    fn test_evolve_universe() {
        let params = CosmologyParams::planck_2018();
        let evolution = evolve_universe(10.0, 100, &params);

        assert_eq!(evolution.redshifts.len(), 101);
        assert_eq!(evolution.scale_factors.len(), 101);

        // Scale factor should decrease going back in time
        assert!(evolution.scale_factors[0] < *evolution.scale_factors.last().unwrap());

        // Hubble parameter should increase going back in time
        assert!(evolution.hubble_params[0] > *evolution.hubble_params.last().unwrap());
    }
