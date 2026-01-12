// Unit tests for physics::relativity::mod
use computational_engine::compute::physics::relativity::mod::*;

use super::*;

    #[test]
    fn test_lorentz_factor() {
        // At 0.6c, γ = 1.25
        let gamma = lorentz_factor(0.6 * C).unwrap();
        assert!((gamma - 1.25).abs() < 0.01);
    }

    #[test]
    fn test_time_dilation() {
        let result = time_dilation(TimeDilationRequest {
            proper_time: 1.0,
            velocity: 0.6 * C,
        })
        .unwrap();

        assert!((result.gamma - 1.25).abs() < 0.01);
        assert!((result.dilated_time - 1.25).abs() < 0.01);
    }

    #[test]
    fn test_length_contraction() {
        let result = length_contraction(LengthContractionRequest {
            proper_length: 10.0,
            velocity: 0.6 * C,
        })
        .unwrap();

        assert!((result.contracted_length - 8.0).abs() < 0.1);
    }

    #[test]
    fn test_relativistic_energy() {
        let result = relativistic_energy(RelativisticEnergyRequest {
            mass: 1.0, // 1 kg
            velocity: 0.6 * C,
        })
        .unwrap();

        let expected_rest = C * C;
        assert!((result.rest_energy - expected_rest).abs() < 1e10);
        assert!(result.total_energy > result.rest_energy);
    }

    #[test]
    fn test_velocity_addition() {
        let result = velocity_addition(VelocityAdditionRequest {
            velocity1: 0.6 * C,
            velocity2: 0.6 * C,
        })
        .unwrap();

        // Should be less than c
        assert!(result.relativistic_sum < C);
        // Should be ~0.882c
        assert!((result.relativistic_sum / C - 0.882).abs() < 0.01);
    }

    #[test]
    fn test_schwarzschild_radius() {
        // Solar mass
        let solar_mass = 1.989e30; // kg
        let result = schwarzschild_metric(SchwarzschildRequest {
            mass: solar_mass,
            radius: 1e10, // 10,000 km
        })
        .unwrap();

        // Schwarzschild radius of sun ~2.95 km
        assert!((result.schwarzschild_radius - 2950.0).abs() < 100.0);
    }

    #[test]
    fn test_black_hole_properties() {
        let solar_mass = 1.989e30; // kg
        let result = black_hole_properties(BlackHoleRequest { mass: solar_mass }).unwrap();

        assert!(result.schwarzschild_radius > 0.0);
        assert!(result.hawking_temperature > 0.0);
        assert!(result.entropy > 0.0);
    }
