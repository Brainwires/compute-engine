//! Warp Drive Optimization
//!
//! Find optimal warp field configurations to minimize energy requirements
//! while maintaining desired velocity and bubble properties.

use super::{WarpDriveConfig, ShapeFunction};
use super::energy::{calculate_total_energy, EnergyRequirements};
use serde::{Deserialize, Serialize};

/// Optimization objective
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum OptimizationObjective {
    /// Minimize total energy magnitude
    MinimizeEnergy,
    /// Minimize peak energy density
    MinimizePeakDensity,
    /// Achieve positive energy only (avoid exotic matter)
    PositiveEnergyOnly,
    /// Balance energy and field smoothness
    Balanced,
}

/// Optimization result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptimizationResult {
    /// Optimized configuration
    pub config: WarpDriveConfig,
    /// Energy requirements for optimized config
    pub energy: EnergyRequirements,
    /// Number of iterations performed
    pub iterations: usize,
    /// Final objective value
    pub objective_value: f64,
    /// Convergence achieved
    pub converged: bool,
}

/// Optimize warp drive configuration
pub fn optimize_warp_drive(
    initial_velocity: f64,
    target_bubble_radius: f64,
    objective: OptimizationObjective,
    max_iterations: usize,
) -> OptimizationResult {
    // Start with reasonable initial guess
    let mut best_config = WarpDriveConfig::new(initial_velocity, target_bubble_radius, target_bubble_radius * 0.1);
    let mut best_energy = calculate_total_energy(&best_config);
    let mut best_objective = evaluate_objective(&best_energy, objective);

    let mut converged = false;
    let mut iterations = 0;

    // Grid search over wall thickness and shape function
    let wall_thickness_range: Vec<f64> = (5..=20)
        .map(|i| target_bubble_radius * (i as f64) / 100.0)
        .collect();

    let shape_functions = vec![
        ShapeFunction::Gaussian,
        ShapeFunction::Tanh,
        ShapeFunction::Optimized,
    ];

    for shape_fn in &shape_functions {
        for &wall_thickness in &wall_thickness_range {
            iterations += 1;
            if iterations >= max_iterations {
                break;
            }

            let mut config = WarpDriveConfig::new(initial_velocity, target_bubble_radius, wall_thickness);
            config.shape_function = *shape_fn;

            // For positive energy objective, enforce subluminal
            if matches!(objective, OptimizationObjective::PositiveEnergyOnly) {
                config = WarpDriveConfig::subluminal(initial_velocity, target_bubble_radius, wall_thickness);
                config.shape_function = *shape_fn;
            }

            let energy = calculate_total_energy(&config);
            let objective_val = evaluate_objective(&energy, objective);

            // Check if this is better
            if objective_val < best_objective {
                // Additional check for positive energy objective
                if matches!(objective, OptimizationObjective::PositiveEnergyOnly) {
                    if !energy.requires_exotic_matter {
                        best_config = config;
                        best_energy = energy;
                        best_objective = objective_val;
                        converged = true;
                    }
                } else {
                    best_config = config;
                    best_energy = energy;
                    best_objective = objective_val;
                }
            }
        }
    }

    OptimizationResult {
        config: best_config,
        energy: best_energy,
        iterations,
        objective_value: best_objective,
        converged,
    }
}

/// Evaluate objective function
fn evaluate_objective(energy: &EnergyRequirements, objective: OptimizationObjective) -> f64 {
    match objective {
        OptimizationObjective::MinimizeEnergy => energy.total_energy.abs(),
        OptimizationObjective::MinimizePeakDensity => energy.peak_energy_density.abs(),
        OptimizationObjective::PositiveEnergyOnly => {
            if energy.requires_exotic_matter {
                f64::INFINITY // Penalize exotic matter heavily
            } else {
                energy.total_energy.abs()
            }
        }
        OptimizationObjective::Balanced => {
            // Balance between total energy and peak density
            energy.total_energy.abs() / 1e30 + energy.peak_energy_density.abs() / 1e10
        }
    }
}

/// Scan parameter space to understand energy landscape
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParameterScan {
    pub velocity: f64,
    pub bubble_radius_range: Vec<f64>,
    pub wall_thickness_range: Vec<f64>,
    pub energy_grid: Vec<Vec<f64>>, // [radius_idx][thickness_idx]
}

pub fn parameter_space_scan(
    velocity: f64,
    radius_min: f64,
    radius_max: f64,
    n_radius: usize,
    thickness_min: f64,
    thickness_max: f64,
    n_thickness: usize,
) -> ParameterScan {
    let mut bubble_radius_range = Vec::new();
    let mut wall_thickness_range = Vec::new();
    let mut energy_grid = Vec::new();

    let dr = (radius_max - radius_min) / (n_radius - 1) as f64;
    let dt = (thickness_max - thickness_min) / (n_thickness - 1) as f64;

    for i in 0..n_radius {
        let radius = radius_min + i as f64 * dr;
        bubble_radius_range.push(radius);

        let mut energy_row = Vec::new();

        for j in 0..n_thickness {
            let thickness = thickness_min + j as f64 * dt;

            if i == 0 {
                wall_thickness_range.push(thickness);
            }

            let config = WarpDriveConfig::new(velocity, radius, thickness);
            let energy = calculate_total_energy(&config);

            energy_row.push(energy.total_energy.abs());
        }

        energy_grid.push(energy_row);
    }

    ParameterScan {
        velocity,
        bubble_radius_range,
        wall_thickness_range,
        energy_grid,
    }
}

/// Find optimal subluminal configuration (positive energy)
pub fn find_subluminal_optimal(
    velocity: f64,
    target_radius: f64,
) -> Result<OptimizationResult, String> {
    if velocity >= super::C {
        return Err("Velocity must be subluminal (v < c)".to_string());
    }

    let result = optimize_warp_drive(
        velocity,
        target_radius,
        OptimizationObjective::PositiveEnergyOnly,
        100,
    );

    if result.energy.requires_exotic_matter {
        Err("Could not find positive-energy configuration".to_string())
    } else {
        Ok(result)
    }
}

/// Compare different shape functions for same parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ShapeFunctionComparison {
    pub velocity: f64,
    pub bubble_radius: f64,
    pub wall_thickness: f64,
    pub results: Vec<(ShapeFunction, EnergyRequirements)>,
}

pub fn compare_shape_functions(
    velocity: f64,
    bubble_radius: f64,
    wall_thickness: f64,
) -> ShapeFunctionComparison {
    let shape_functions = vec![
        ShapeFunction::TopHat,
        ShapeFunction::Gaussian,
        ShapeFunction::Tanh,
        ShapeFunction::Optimized,
    ];

    let mut results = Vec::new();

    for shape_fn in shape_functions {
        let mut config = WarpDriveConfig::new(velocity, bubble_radius, wall_thickness);
        config.shape_function = shape_fn;

        let energy = calculate_total_energy(&config);
        results.push((shape_fn, energy));
    }

    ShapeFunctionComparison {
        velocity,
        bubble_radius,
        wall_thickness,
        results,
    }
}

