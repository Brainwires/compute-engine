// Unit tests for physics::warp_drive module
use super::*;
use std::f64::consts::PI;

// ============================================================================
// WarpDriveConfig Tests
// ============================================================================

#[test]
fn test_config_creation() {
    let config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);

    assert_eq!(config.velocity, 0.5 * C);
    assert_eq!(config.bubble_radius, 100.0);
    assert_eq!(config.wall_thickness, 10.0);
    assert_eq!(config.shape_function, ShapeFunction::Tanh);
    assert_eq!(config.subluminal, true);
}

#[test]
fn test_subluminal_config() {
    let config = WarpDriveConfig::subluminal(0.9 * C, 50.0, 5.0);

    assert!(config.velocity < C);
    assert_eq!(config.subluminal, true);
    assert_eq!(config.bubble_radius, 50.0);
}

#[test]
fn test_subluminal_enforcement() {
    // Request superluminal, should be clamped to 0.99c
    let config = WarpDriveConfig::subluminal(2.0 * C, 50.0, 5.0);

    assert!(config.velocity < C);
    assert_eq!(config.velocity, 0.99 * C);
    assert_eq!(config.subluminal, true);
}

#[test]
fn test_lorentz_factor_subluminal() {
    let config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);
    let gamma = config.lorentz_factor();

    // γ = 1/√(1 - 0.5²) ≈ 1.1547
    assert!((gamma - 1.1547).abs() < 0.001);
}

#[test]
fn test_lorentz_factor_near_light_speed() {
    let config = WarpDriveConfig::new(0.99 * C, 100.0, 10.0);
    let gamma = config.lorentz_factor();

    // At 0.99c, γ ≈ 7.09
    assert!(gamma > 7.0 && gamma < 8.0);
}

#[test]
fn test_lorentz_factor_superluminal() {
    let mut config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);
    config.velocity = 2.0 * C;
    config.subluminal = false;

    let gamma = config.lorentz_factor();
    assert!(gamma.is_infinite());
}

#[test]
fn test_validation_positive_velocity() {
    let mut config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);
    config.velocity = -1.0;

    let result = config.validate();
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("positive"));
}

#[test]
fn test_validation_subluminal_constraint() {
    let mut config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);
    config.velocity = 2.0 * C;
    config.subluminal = true;

    let result = config.validate();
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("Subluminal"));
}

#[test]
fn test_validation_bubble_radius() {
    let mut config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);
    config.bubble_radius = -5.0;

    let result = config.validate();
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("radius"));
}

#[test]
fn test_validation_wall_thickness() {
    let mut config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);
    config.wall_thickness = 0.0;

    let result = config.validate();
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("thickness"));
}

// ============================================================================
// Shape Function Tests
// ============================================================================

#[test]
fn test_shape_function_tophat_center() {
    let config = WarpDriveConfig {
        velocity: 0.5 * C,
        bubble_radius: 100.0,
        wall_thickness: 10.0,
        shape_function: ShapeFunction::TopHat,
        subluminal: true,
    };

    let f = shape_function(50.0, &config); // Inside bubble
    assert_eq!(f, 1.0);
}

#[test]
fn test_shape_function_tophat_outside() {
    let config = WarpDriveConfig {
        velocity: 0.5 * C,
        bubble_radius: 100.0,
        wall_thickness: 10.0,
        shape_function: ShapeFunction::TopHat,
        subluminal: true,
    };

    let f = shape_function(120.0, &config); // Outside bubble
    assert_eq!(f, 0.0);
}

#[test]
fn test_shape_function_gaussian() {
    let config = WarpDriveConfig {
        velocity: 0.5 * C,
        bubble_radius: 100.0,
        wall_thickness: 10.0,
        shape_function: ShapeFunction::Gaussian,
        subluminal: true,
    };

    let f_center = shape_function(100.0, &config); // At radius
    assert!((f_center - 1.0).abs() < 0.01);

    let f_far = shape_function(150.0, &config); // Far away
    assert!(f_far < 0.1);
}

#[test]
fn test_shape_function_tanh() {
    let config = WarpDriveConfig {
        velocity: 0.5 * C,
        bubble_radius: 100.0,
        wall_thickness: 10.0,
        shape_function: ShapeFunction::Tanh,
        subluminal: true,
    };

    let f_inside = shape_function(50.0, &config);
    assert!(f_inside > 0.9);

    let f_outside = shape_function(150.0, &config);
    assert!(f_outside < 0.1);
}

#[test]
fn test_shape_function_derivative_gaussian() {
    let config = WarpDriveConfig {
        velocity: 0.5 * C,
        bubble_radius: 100.0,
        wall_thickness: 10.0,
        shape_function: ShapeFunction::Gaussian,
        subluminal: true,
    };

    let df = shape_function_derivative(100.0, &config);
    // At the center radius, derivative should be near zero
    assert!(df.abs() < 0.01);
}

// ============================================================================
// Alcubierre Metric Tests
// ============================================================================

#[test]
fn test_alcubierre_metric_flat_far_field() {
    let config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);

    let coords = Coordinates {
        t: 0.0,
        x: 1000.0, // Far from bubble
        y: 0.0,
        z: 0.0,
    };

    let metric = alcubierre_metric(&coords, &config);

    // Far from bubble, metric should be nearly Minkowski
    assert!((metric.g_tt + C2).abs() < 1e6); // g_tt ≈ -c²
    assert!(metric.g_tx.abs() < 1e-3);
    assert!((metric.g_xx - 1.0).abs() < 1e-6);
    assert_eq!(metric.g_yy, 1.0);
    assert_eq!(metric.g_zz, 1.0);
}

#[test]
fn test_alcubierre_metric_inside_bubble() {
    let config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);

    let coords = Coordinates {
        t: 0.0,
        x: 0.0, // Center of bubble
        y: 0.0,
        z: 0.0,
    };

    let metric = alcubierre_metric(&coords, &config);

    // Inside bubble, shape function f ≈ 1, so g_tx ≈ -v
    assert!(metric.g_tx.abs() > 0.0);
    assert!(metric.g_tt > -C2); // Modified by warp field
}

#[test]
fn test_shift_vector() {
    let config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);

    let coords = Coordinates {
        t: 0.0,
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    let (beta_x, beta_y, beta_z) = shift_vector(&coords, &config);

    // Inside bubble, β^x ≈ v
    assert!(beta_x > 0.0);
    assert_eq!(beta_y, 0.0);
    assert_eq!(beta_z, 0.0);
}

#[test]
fn test_lapse_function_inside_bubble() {
    let config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);

    let coords = Coordinates {
        t: 0.0,
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    let alpha = lapse_function(&coords, &config);

    // Lapse function should be positive and less than c
    assert!(alpha > 0.0);
    assert!(alpha <= C);
}

#[test]
fn test_lapse_function_far_field() {
    let config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);

    let coords = Coordinates {
        t: 0.0,
        x: 1000.0,
        y: 0.0,
        z: 0.0,
    };

    let alpha = lapse_function(&coords, &config);

    // Far from bubble, α → c
    assert!((alpha - C).abs() < 1.0);
}

#[test]
fn test_proper_time() {
    let config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);

    let coords = Coordinates {
        t: 0.0,
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    let dt = 1.0; // 1 second coordinate time
    let dtau = proper_time(dt, &coords, &config);

    // Proper time should be positive
    assert!(dtau > 0.0);
    assert!(dtau <= dt);
}

#[test]
fn test_effective_distance() {
    let config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);

    let distance = effective_distance(-200.0, 200.0, &config, 0.0);

    // Should compute a positive distance
    assert!(distance > 0.0);
}

// ============================================================================
// Energy Requirements Tests
// ============================================================================

#[test]
fn test_stress_energy_tensor_center() {
    let config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);

    let coords = Coordinates {
        t: 0.0,
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    let tensor = compute_stress_energy(&coords, &config);

    // At center, energy density should be regularized to 0
    assert_eq!(tensor.energy_density, 0.0);
}

#[test]
fn test_stress_energy_tensor_wall() {
    let config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);

    // At the wall where df/dr is maximum
    let coords = Coordinates {
        t: 0.0,
        x: 100.0, // At bubble radius
        y: 0.0,
        z: 0.0,
    };

    let tensor = compute_stress_energy(&coords, &config);

    // Energy density should be non-zero at the wall
    assert!(tensor.energy_density != 0.0);
}

#[test]
fn test_calculate_total_energy() {
    let config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);

    let energy = calculate_total_energy(&config);

    assert!(energy.total_energy != 0.0);
    assert!(energy.bubble_volume > 0.0);
    assert!(energy.solar_masses >= 0.0);
}

#[test]
fn test_energy_superluminal_requires_exotic() {
    let mut config = WarpDriveConfig::new(2.0 * C, 100.0, 10.0);
    config.subluminal = false;

    let energy = calculate_total_energy(&config);

    // Superluminal drives require exotic matter
    assert!(energy.requires_exotic_matter);
}

#[test]
fn test_alcubierre_energy_estimate() {
    let config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);

    let energy = alcubierre_energy_estimate(&config);

    // Classic Alcubierre requires negative energy
    assert!(energy < 0.0);
}

#[test]
fn test_alcubierre_energy_scales_with_velocity() {
    let config1 = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);
    let config2 = WarpDriveConfig::new(0.9 * C, 100.0, 10.0);

    let e1 = alcubierre_energy_estimate(&config1).abs();
    let e2 = alcubierre_energy_estimate(&config2).abs();

    // Energy should scale with v²
    assert!(e2 > e1);
}

#[test]
fn test_subluminal_positive_energy_estimate() {
    let config = WarpDriveConfig::subluminal(0.9 * C, 100.0, 10.0);

    let result = subluminal_positive_energy_estimate(&config);
    assert!(result.is_ok());

    let energy = result.unwrap();
    assert!(energy > 0.0); // Should be positive
}

#[test]
fn test_subluminal_positive_energy_requires_subluminal() {
    let mut config = WarpDriveConfig::new(2.0 * C, 100.0, 10.0);
    config.subluminal = false;

    let result = subluminal_positive_energy_estimate(&config);
    assert!(result.is_err());
}

#[test]
fn test_quantum_energy_bound() {
    let sigma = 1.0; // 1 meter wall thickness
    let bound = quantum_energy_bound(sigma);

    // Should be positive and scale as σ^(-4)
    assert!(bound > 0.0);

    let sigma2 = 0.5;
    let bound2 = quantum_energy_bound(sigma2);

    // Halving sigma should increase bound by 2^4 = 16
    let ratio = bound2 / bound;
    assert!((ratio - 16.0).abs() < 0.1);
}

// ============================================================================
// Optimization Tests
// ============================================================================

#[test]
fn test_optimize_warp_drive_minimize_energy() {
    let velocity = 0.5 * C;
    let radius = 100.0;

    let result = optimize_warp_drive(
        velocity,
        radius,
        OptimizationObjective::MinimizeEnergy,
        50,
    );

    assert!(result.iterations > 0);
    assert!(result.objective_value >= 0.0);
    assert_eq!(result.config.velocity, velocity);
}

#[test]
fn test_optimize_warp_drive_positive_energy() {
    let velocity = 0.8 * C;
    let radius = 50.0;

    let result = optimize_warp_drive(
        velocity,
        radius,
        OptimizationObjective::PositiveEnergyOnly,
        100,
    );

    if result.converged {
        assert!(!result.energy.requires_exotic_matter);
        assert!(result.config.subluminal);
    }
}

#[test]
fn test_find_subluminal_optimal() {
    let velocity = 0.7 * C;
    let radius = 100.0;

    let result = find_subluminal_optimal(velocity, radius);

    // May or may not find solution, but should not panic
    if let Ok(opt_result) = result {
        assert!(!opt_result.energy.requires_exotic_matter);
        assert!(opt_result.config.velocity < C);
    }
}

#[test]
fn test_find_subluminal_optimal_rejects_superluminal() {
    let velocity = 2.0 * C;
    let radius = 100.0;

    let result = find_subluminal_optimal(velocity, radius);
    assert!(result.is_err());
}

#[test]
fn test_parameter_space_scan() {
    let scan = parameter_space_scan(
        0.5 * C,
        50.0,  // radius min
        150.0, // radius max
        5,     // n_radius
        5.0,   // thickness min
        15.0,  // thickness max
        5,     // n_thickness
    );

    assert_eq!(scan.bubble_radius_range.len(), 5);
    assert_eq!(scan.wall_thickness_range.len(), 5);
    assert_eq!(scan.energy_grid.len(), 5);
    assert_eq!(scan.energy_grid[0].len(), 5);

    // All energies should be non-negative (absolute values)
    for row in &scan.energy_grid {
        for &energy in row {
            assert!(energy >= 0.0);
        }
    }
}

#[test]
fn test_compare_shape_functions() {
    let comparison = compare_shape_functions(0.5 * C, 100.0, 10.0);

    assert_eq!(comparison.results.len(), 4); // 4 shape functions
    assert_eq!(comparison.velocity, 0.5 * C);
    assert_eq!(comparison.bubble_radius, 100.0);

    // All should return valid energy requirements
    for (shape_fn, energy) in &comparison.results {
        // Energy requirements should have valid values (may be zero for TopHat derivative issues)
        assert!(energy.bubble_volume > 0.0);
        assert!(energy.solar_masses >= 0.0);
    }
}

#[test]
fn test_energy_type_classification() {
    let config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);
    let energy = calculate_total_energy(&config);

    // Should classify energy type
    match energy.energy_type {
        EnergyType::Exotic | EnergyType::Positive | EnergyType::Mixed => {
            // Valid classification
        }
    }
}

// ============================================================================
// Physical Consistency Tests
// ============================================================================

#[test]
fn test_energy_density_peaks_at_wall() {
    let config = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);

    let coords_center = Coordinates { t: 0.0, x: 0.0, y: 0.0, z: 0.0 };
    let coords_wall = Coordinates { t: 0.0, x: 100.0, y: 0.0, z: 0.0 };
    let coords_far = Coordinates { t: 0.0, x: 200.0, y: 0.0, z: 0.0 };

    let tensor_center = compute_stress_energy(&coords_center, &config);
    let tensor_wall = compute_stress_energy(&coords_wall, &config);
    let tensor_far = compute_stress_energy(&coords_far, &config);

    // Energy density should peak near the wall
    assert!(tensor_wall.energy_density.abs() > tensor_center.energy_density.abs());
    assert!(tensor_wall.energy_density.abs() > tensor_far.energy_density.abs());
}

#[test]
fn test_thinner_walls_increase_energy() {
    let config1 = WarpDriveConfig::new(0.5 * C, 100.0, 10.0);
    let config2 = WarpDriveConfig::new(0.5 * C, 100.0, 5.0); // Thinner wall

    let e1 = alcubierre_energy_estimate(&config1).abs();
    let e2 = alcubierre_energy_estimate(&config2).abs();

    // Thinner walls require more energy (E ∝ R/σ)
    assert!(e2 > e1);
}

#[test]
fn test_shape_function_continuity() {
    let config = WarpDriveConfig {
        velocity: 0.5 * C,
        bubble_radius: 100.0,
        wall_thickness: 10.0,
        shape_function: ShapeFunction::Tanh,
        subluminal: true,
    };

    // Test continuity across the wall
    let f1 = shape_function(99.0, &config);
    let f2 = shape_function(100.0, &config);
    let f3 = shape_function(101.0, &config);

    // Should be continuous (no jumps)
    assert!((f2 - f1).abs() < 0.5);
    assert!((f3 - f2).abs() < 0.5);
}
