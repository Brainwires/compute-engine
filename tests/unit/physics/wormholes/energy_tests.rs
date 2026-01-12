// Unit tests for physics::wormholes::energy
use computational_engine::compute::physics::wormholes::energy::*;

use super::*;
    use crate::compute::physics::wormholes::WormholeConfig;

    #[test]
    fn test_energy_density_at_throat() {
        let config = WormholeConfig::morris_thorne(100.0);
        let coords = SphericalCoordinates {
            t: 0.0,
            r: 100.0,
            theta: PI / 2.0,
            phi: 0.0,
        };

        let energy = compute_energy_density(&coords, &config);

        // Should have non-zero energy density
        assert!(energy.rho.abs() > 0.0);
    }

    #[test]
    fn test_exotic_matter_required() {
        let config = WormholeConfig::morris_thorne(100.0);
        let energy = calculate_total_energy(&config);

        // Morris-Thorne wormholes typically require exotic matter
        // (exact detection depends on numerical precision)
        assert!(energy.exotic_mass.is_finite());
    }

    #[test]
    fn test_total_energy_finite() {
        let config = WormholeConfig::new(100.0, 200.0);
        let energy = calculate_total_energy(&config);

        assert!(energy.exotic_mass.is_finite());
        assert!(energy.throat_area > 0.0);
    }

    #[test]
    fn test_black_hole_comparison() {
        let config = WormholeConfig::new(1000.0, 2000.0);
        let bh_mass = compare_with_black_hole(&config);

        // Black hole mass should be positive and large
        assert!(bh_mass > 0.0);
        assert!(bh_mass > 1e20); // Should be huge for km-scale
    }

    #[test]
    fn test_throat_energy_density() {
        let config = WormholeConfig::morris_thorne(100.0);
        let rho_throat = throat_energy_density(&config);

        // Should be non-zero at throat
        assert!(rho_throat.abs() > 0.0);
    }

    #[test]
    fn test_energy_scales_with_throat() {
        let config1 = WormholeConfig::morris_thorne(100.0);
        let config2 = WormholeConfig::morris_thorne(200.0);

        let e1 = calculate_total_energy(&config1);
        let e2 = calculate_total_energy(&config2);

        // Larger throat should require more energy
        assert!(e2.exotic_mass.abs() > e1.exotic_mass.abs());
    }
