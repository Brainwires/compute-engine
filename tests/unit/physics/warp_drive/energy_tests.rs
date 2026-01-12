// Unit tests for physics::warp_drive::energy
use computational_engine::compute::physics::warp_drive::energy::*;

use super::*;
    use crate::compute::physics::warp_drive::WarpDriveConfig;

    #[test]
    fn test_stress_energy_at_bubble_edge() {
        let config = WarpDriveConfig::new(1e7, 100.0, 10.0);
        let coords = Coordinates {
            t: 0.0,
            x: 100.0, // At bubble radius
            y: 0.0,
            z: 0.0,
        };

        let stress_energy = compute_stress_energy(&coords, &config);

        // Energy density should be non-zero at bubble wall
        assert!(stress_energy.energy_density.abs() > 0.0);
    }

    #[test]
    fn test_stress_energy_far_from_bubble() {
        let config = WarpDriveConfig::new(1e7, 100.0, 10.0);
        let coords = Coordinates {
            t: 0.0,
            x: 500.0, // Far from bubble
            y: 0.0,
            z: 0.0,
        };

        let stress_energy = compute_stress_energy(&coords, &config);

        // Energy density should be near zero far from bubble
        assert!(stress_energy.energy_density.abs() < 1e10);
    }

    #[test]
    fn test_total_energy_calculation() {
        let config = WarpDriveConfig::new(1e7, 100.0, 10.0);
        let requirements = calculate_total_energy(&config);

        // Should have finite total energy
        assert!(requirements.total_energy.is_finite());
        assert!(requirements.bubble_volume > 0.0);
    }

    #[test]
    fn test_superluminal_requires_exotic_matter() {
        let config = WarpDriveConfig::new(2.0 * C, 100.0, 10.0);
        let requirements = calculate_total_energy(&config);

        // Superluminal configurations typically require exotic matter
        assert!(requirements.requires_exotic_matter);
        assert_eq!(requirements.energy_type, EnergyType::Exotic);
    }

    #[test]
    fn test_subluminal_positive_energy() {
        let config = WarpDriveConfig::subluminal(0.1 * C, 100.0, 10.0);
        let energy = subluminal_positive_energy_estimate(&config);

        assert!(energy.is_ok());
        assert!(energy.unwrap() > 0.0); // Positive energy
    }

    #[test]
    fn test_alcubierre_estimate_negative() {
        let config = WarpDriveConfig::new(2.0 * C, 100.0, 10.0);
        let energy = alcubierre_energy_estimate(&config);

        // Classic Alcubierre requires negative energy
        assert!(energy < 0.0);
    }

    #[test]
    fn test_quantum_energy_bound() {
        let sigma = 1.0; // 1 meter wall thickness
        let bound = quantum_energy_bound(sigma);

        // Should be a small positive value
        assert!(bound > 0.0);
        assert!(bound < 1e-10); // Very small for meter-scale
    }

    #[test]
    fn test_energy_scales_with_velocity() {
        let config1 = WarpDriveConfig::new(1e7, 100.0, 10.0);
        let config2 = WarpDriveConfig::new(2e7, 100.0, 10.0);

        let e1 = alcubierre_energy_estimate(&config1).abs();
        let e2 = alcubierre_energy_estimate(&config2).abs();

        // Energy should scale with v²
        assert!(e2 > e1);
        assert!((e2 / e1 - 4.0).abs() < 0.5); // Should be ~4x for 2x velocity
    }
