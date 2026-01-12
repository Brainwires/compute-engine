// Unit tests for physics::black_holes::hawking
use computational_engine::compute::physics::black_holes::hawking::*;

use super::*;
    use crate::compute::physics::black_holes::BlackHoleConfig;

    #[test]
    fn test_hawking_radiation() {
        let solar_mass = 1.989e30;
        let bh = BlackHoleConfig::schwarzschild(solar_mass);
        let radiation = hawking_radiation(&bh);

        // Solar mass BH: very cold
        assert!(radiation.temperature < 1e-6);
        assert!(radiation.luminosity > 0.0);
        assert!(radiation.evaporation_time > 1e60);
    }

    #[test]
    fn test_mass_loss_rate() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let dm_dt = mass_loss_rate(&bh);

        // Should be negative (losing mass)
        assert!(dm_dt < 0.0);
        assert!(dm_dt.is_finite());
    }
