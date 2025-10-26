// Unit tests for physics::black_holes::orbits
use computational_engine::physics::black_holes::orbits::*;

use super::*;
    use crate::physics::black_holes::BlackHoleConfig;

    #[test]
    fn test_circular_orbit() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_isco = bh.isco_radius();

        let orbit = circular_orbit(10.0 * r_isco, &bh);

        assert!(orbit.is_stable);
        assert!(orbit.angular_velocity > 0.0);
        assert!(orbit.orbital_period > 0.0);
    }

    #[test]
    fn test_orbital_velocity() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let v = orbital_velocity(bh.isco_radius(), &bh);

        assert!(v > 0.0);
        assert!(v < C);
    }
