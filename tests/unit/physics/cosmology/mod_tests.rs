// Unit tests for physics::cosmology::mod
use computational_engine::physics::cosmology::mod::*;

use super::*;

    #[test]
    fn test_planck_parameters() {
        let params = CosmologyParams::planck_2018();

        assert!((params.omega_m - 0.315).abs() < 0.001);
        assert!((params.omega_lambda - 0.685).abs() < 0.001);

        // Total should be close to 1 (flat universe)
        let omega_total = params.omega_total();
        assert!((omega_total - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_critical_density() {
        let params = CosmologyParams::planck_2018();
        let rho_c = params.critical_density();

        // Critical density should be ~10^-26 kg/m³
        assert!(rho_c > 1e-27 && rho_c < 1e-25);
    }

    #[test]
    fn test_age_of_universe() {
        let params = CosmologyParams::planck_2018();
        let age = params.age_of_universe();

        // Age should be positive and finite (billions of years)
        assert!(age > 1e17); // More than 3 billion years
        assert!(age.is_finite());
    }

    #[test]
    fn test_redshift_scale_factor() {
        let z0 = Redshift(0.0);
        assert_eq!(z0.to_scale_factor(), 1.0);

        let z1 = Redshift(1.0);
        assert_eq!(z1.to_scale_factor(), 0.5);

        let z_recomb = Redshift(1100.0); // Recombination
        let a = z_recomb.to_scale_factor();
        assert!((a - 1.0/1101.0).abs() < 1e-6);
    }

    #[test]
    fn test_cmb_temperature_evolution() {
        let params = CosmologyParams::planck_2018();

        let z_today = Redshift(0.0);
        let t_today = z_today.temperature(&params);
        assert!((t_today - 2.7255).abs() < 0.001);

        let z_recomb = Redshift(1100.0);
        let t_recomb = z_recomb.temperature(&params);
        // At recombination, T ~ 3000 K
        assert!(t_recomb > 2900.0 && t_recomb < 3100.0);
    }

    #[test]
    fn test_omega_k() {
        let params = CosmologyParams::planck_2018();
        let omega_k = params.omega_k();

        // Should be very close to 0 (flat universe)
        assert!(omega_k.abs() < 0.01);
    }
