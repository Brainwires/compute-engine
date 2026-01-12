// Unit tests for physics::black_holes::mod
use computational_engine::compute::physics::black_holes::mod::*;

use super::*;

    #[test]
    fn test_schwarzschild_creation() {
        let solar_mass = 1.989e30;
        let bh = BlackHoleConfig::schwarzschild(solar_mass);

        assert_eq!(bh.mass, solar_mass);
        assert_eq!(bh.spin, 0.0);
        assert_eq!(bh.bh_type, BlackHoleType::Schwarzschild);
    }

    #[test]
    fn test_schwarzschild_radius() {
        let solar_mass = 1.989e30;
        let bh = BlackHoleConfig::schwarzschild(solar_mass);
        let r_s = bh.schwarzschild_radius();

        // Solar mass black hole: r_s ≈ 3 km
        assert!((r_s - 2953.0).abs() < 100.0);
    }

    #[test]
    fn test_event_horizon() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_h = bh.event_horizon_radius();

        // Should equal Schwarzschild radius for non-rotating
        assert!((r_h - bh.schwarzschild_radius()).abs() < 1.0);
    }

    #[test]
    fn test_photon_sphere() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_ph = bh.photon_sphere_radius();
        let r_s = bh.schwarzschild_radius();

        // Photon sphere at 1.5 r_s
        assert!((r_ph - 1.5 * r_s).abs() < 10.0);
    }

    #[test]
    fn test_isco() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_isco = bh.isco_radius();
        let r_s = bh.schwarzschild_radius();

        // ISCO at 3 r_s (6 GM/c²)
        assert!((r_isco - 3.0 * r_s).abs() < 10.0);
    }

    #[test]
    fn test_hawking_temperature() {
        let solar_mass = 1.989e30;
        let bh = BlackHoleConfig::schwarzschild(solar_mass);
        let temp = bh.hawking_temperature();

        // Solar mass BH: T ~ 60 nanokelvin
        assert!(temp > 0.0);
        assert!(temp < 1e-6); // Very cold
    }

    #[test]
    fn test_kerr_creation() {
        let bh = BlackHoleConfig::kerr(1e30, 0.5e30);

        assert!(bh.spin > 0.0);
        assert!(bh.bh_type != BlackHoleType::Schwarzschild);
    }

    #[test]
    fn test_extremal_limit() {
        let mass = 1e30;
        let bh = BlackHoleConfig::kerr(mass, mass * 2.0); // Exceeds limit

        // Should be capped at extremal
        assert!(bh.spin <= mass);
    }

    #[test]
    fn test_surface_gravity() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let kappa = bh.surface_gravity();

        // Should be positive and finite
        assert!(kappa > 0.0);
        assert!(kappa.is_finite());
    }

    #[test]
    fn test_entropy() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let entropy = bh.entropy();

        // Entropy should be huge for macroscopic black hole
        assert!(entropy > 1e50);
    }

    #[test]
    fn test_evaporation_time() {
        let solar_mass = 1.989e30;
        let bh = BlackHoleConfig::schwarzschild(solar_mass);
        let t_evap = bh.evaporation_time();

        // Solar mass BH: ~10^67 years >> age of universe
        assert!(t_evap > 1e60); // seconds
    }
