// Unit tests for physics::wormholes::mod
use computational_engine::physics::wormholes::mod::*;

use super::*;

    #[test]
    fn test_wormhole_config_creation() {
        let config = WormholeConfig::new(100.0, 200.0);
        assert_eq!(config.throat_radius, 100.0);
        assert_eq!(config.length_scale, 200.0);
    }

    #[test]
    fn test_morris_thorne_config() {
        let config = WormholeConfig::morris_thorne(100.0);
        assert_eq!(config.throat_radius, 100.0);
        assert!(matches!(config.redshift_type, RedshiftFunction::Zero));
    }

    #[test]
    fn test_validation() {
        let config = WormholeConfig::new(100.0, 200.0);
        assert!(config.validate().is_ok());

        let bad_config = WormholeConfig::new(-10.0, 200.0);
        assert!(bad_config.validate().is_err());
    }

    #[test]
    fn test_schwarzschild_comparison() {
        let config = WormholeConfig::new(1000.0, 2000.0);
        let solar_mass = 1.989e30; // kg
        let r_s = config.schwarzschild_radius_equivalent(solar_mass);

        // Solar mass Schwarzschild radius is ~3 km
        assert!((r_s - 2953.0).abs() < 100.0);
    }
