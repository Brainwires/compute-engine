// Unit tests for physics::warp_drive::mod
use computational_engine::compute::physics::warp_drive::mod::*;

use super::*;

    #[test]
    fn test_warp_config_creation() {
        let config = WarpDriveConfig::new(1e6, 100.0, 10.0);
        assert_eq!(config.velocity, 1e6);
        assert_eq!(config.bubble_radius, 100.0);
        assert!(config.subluminal);
    }

    #[test]
    fn test_lorentz_factor_nonrelativistic() {
        let config = WarpDriveConfig::new(1e6, 100.0, 10.0); // ~0.003c
        let gamma = config.lorentz_factor();
        assert!((gamma - 1.0).abs() < 0.01); // Should be very close to 1
    }

    #[test]
    fn test_lorentz_factor_relativistic() {
        let config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);
        let gamma = config.lorentz_factor();
        assert!((gamma - 1.1547).abs() < 0.01); // γ(0.5c) ≈ 1.1547
    }

    #[test]
    fn test_subluminal_enforcement() {
        let config = WarpDriveConfig::subluminal(2.0 * C, 100.0, 10.0);
        assert!(config.velocity < C); // Should be capped below c
        assert!(config.subluminal);
    }

    #[test]
    fn test_validation() {
        let config = WarpDriveConfig::new(1e6, 100.0, 10.0);
        assert!(config.validate().is_ok());

        let bad_config = WarpDriveConfig::new(-1.0, 100.0, 10.0);
        assert!(bad_config.validate().is_err());
    }
