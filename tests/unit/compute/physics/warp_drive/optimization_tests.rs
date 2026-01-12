// Unit tests for physics::warp_drive::optimization
use computational_engine::compute::physics::warp_drive::optimization::*;

use super::*;
    use crate::compute::physics::warp_drive::C;

    #[test]
    fn test_optimization_minimize_energy() {
        let result = optimize_warp_drive(
            1e7,  // 10^7 m/s (~0.03c)
            100.0,
            OptimizationObjective::MinimizeEnergy,
            50,
        );

        assert!(result.iterations > 0);
        assert!(result.energy.total_energy.is_finite());
    }

    #[test]
    fn test_subluminal_optimization() {
        let result = optimize_warp_drive(
            0.1 * C,
            100.0,
            OptimizationObjective::PositiveEnergyOnly,
            50,
        );

        // For low velocities, should be able to achieve positive energy
        assert!(result.config.subluminal);
    }

    #[test]
    fn test_parameter_space_scan() {
        let scan = parameter_space_scan(
            1e7,
            50.0,
            150.0,
            5,
            5.0,
            15.0,
            5,
        );

        assert_eq!(scan.bubble_radius_range.len(), 5);
        assert_eq!(scan.wall_thickness_range.len(), 5);
        assert_eq!(scan.energy_grid.len(), 5);
        assert_eq!(scan.energy_grid[0].len(), 5);
    }

    #[test]
    fn test_find_subluminal_optimal() {
        let result = find_subluminal_optimal(0.05 * C, 100.0);

        // Should succeed for low velocities
        if let Ok(opt) = result {
            assert!(opt.config.velocity < C);
            assert_eq!(opt.config.subluminal, true);
        }
    }

    #[test]
    fn test_shape_function_comparison() {
        let comparison = compare_shape_functions(1e7, 100.0, 10.0);

        assert_eq!(comparison.results.len(), 4); // 4 shape functions

        // All should have finite energy
        for (_shape, energy) in &comparison.results {
            assert!(energy.total_energy.is_finite());
        }
    }

    #[test]
    fn test_evaluation_positive_energy_objective() {
        let mut energy = EnergyRequirements {
            total_energy: -1e30,
            solar_masses: 0.5,
            peak_energy_density: -1e10,
            bubble_volume: 1e6,
            requires_exotic_matter: true,
            energy_type: super::super::energy::EnergyType::Exotic,
        };

        let obj_val = evaluate_objective(&energy, OptimizationObjective::PositiveEnergyOnly);
        assert!(obj_val.is_infinite()); // Should heavily penalize exotic matter

        energy.requires_exotic_matter = false;
        let obj_val2 = evaluate_objective(&energy, OptimizationObjective::PositiveEnergyOnly);
        assert!(obj_val2.is_finite());
    }

    #[test]
    fn test_optimization_converges() {
        let result = optimize_warp_drive(
            1e7,
            100.0,
            OptimizationObjective::MinimizeEnergy,
            100,
        );

        // Should converge for reasonable parameters
        assert!(result.objective_value.is_finite());
        assert!(result.objective_value > 0.0);
    }
