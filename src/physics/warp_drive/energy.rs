//! Energy Requirements for Warp Drives
//!
//! Calculates the stress-energy tensor T_μν and total energy requirements
//! using Einstein field equations: G_μν = (8πG/c⁴) T_μν
//!
//! The classic Alcubierre drive requires **negative energy density** (exotic matter).
//! Recent work (Bobrick & Martire 2021) shows subluminal drives can use positive energy.

use super::{WarpDriveConfig, C, C2, C4, G};
use super::metric::{shape_function, shape_function_derivative, Coordinates};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Stress-energy tensor components T_μν
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StressEnergyTensor {
    /// Energy density ρ = -T⁰⁰ (negative for exotic matter)
    pub energy_density: f64,
    /// Pressure components
    pub pressure_x: f64,
    pub pressure_y: f64,
    pub pressure_z: f64,
    /// Indicates if energy density is negative (exotic matter required)
    pub requires_exotic_matter: bool,
}

/// Total energy requirements for the warp bubble
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnergyRequirements {
    /// Total energy (J) - can be negative for exotic matter
    pub total_energy: f64,
    /// Energy in solar masses
    pub solar_masses: f64,
    /// Peak energy density (J/m³)
    pub peak_energy_density: f64,
    /// Volume of warp bubble (m³)
    pub bubble_volume: f64,
    /// Whether configuration requires exotic (negative energy) matter
    pub requires_exotic_matter: bool,
    /// Energy type classification
    pub energy_type: EnergyType,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum EnergyType {
    /// Requires negative energy density (exotic matter)
    Exotic,
    /// Uses only positive energy (physically realizable)
    Positive,
    /// Mixed: some regions positive, some negative
    Mixed,
}

/// Compute stress-energy tensor at given coordinates
/// Using Einstein field equations: T_μν = (c⁴/8πG) G_μν
pub fn compute_stress_energy(coords: &Coordinates, config: &WarpDriveConfig) -> StressEnergyTensor {
    let vs = config.velocity;
    let xs = vs * coords.t;
    let rs = ((coords.x - xs).powi(2) + coords.y.powi(2) + coords.z.powi(2)).sqrt();

    let _f = shape_function(rs, config);
    let df_drs = shape_function_derivative(rs, config);

    // For Alcubierre metric, key components are:
    // Energy density: ρ ≈ -(c⁴/32πG) · vₛ² · (df/drₛ)² / rₛ²
    // This comes from Einstein tensor G_μν calculation

    let energy_density = if rs > 1e-6 {
        let factor = -(C4 / (32.0 * PI * G));
        let vs2 = vs * vs;
        let df2 = df_drs * df_drs;
        factor * vs2 * df2 / (rs * rs)
    } else {
        0.0 // Regularize at center
    };

    // Pressure components (simplified from full tensor analysis)
    // For detailed calculation, need full Ricci tensor computation
    let pressure_factor = energy_density / 3.0;

    StressEnergyTensor {
        energy_density,
        pressure_x: pressure_factor,
        pressure_y: pressure_factor,
        pressure_z: pressure_factor,
        requires_exotic_matter: energy_density < 0.0,
    }
}

/// Calculate total energy requirements by integrating over warp bubble volume
pub fn calculate_total_energy(config: &WarpDriveConfig) -> EnergyRequirements {
    let r_max = config.bubble_radius + 3.0 * config.wall_thickness;
    let n_r = 50;
    let n_theta = 20;
    let n_phi = 20;

    let dr = r_max / n_r as f64;
    let dtheta = PI / n_theta as f64;
    let dphi = 2.0 * PI / n_phi as f64;

    let mut total_energy = 0.0;
    let mut peak_density: f64 = 0.0;
    let mut has_exotic = false;
    let mut has_positive = false;

    // Integrate in spherical coordinates
    for i in 0..n_r {
        let r = (i as f64 + 0.5) * dr;

        for j in 0..n_theta {
            let theta = (j as f64 + 0.5) * dtheta;

            for _ in 0..n_phi {
                // Jacobian for spherical coordinates: r² sin(theta)
                let dv = r * r * theta.sin() * dr * dtheta * dphi;

                let coords = Coordinates {
                    t: 0.0,
                    x: r,
                    y: 0.0,
                    z: 0.0,
                };

                let stress_energy = compute_stress_energy(&coords, config);
                let rho = stress_energy.energy_density;

                // Total energy E = ∫ ρ dV
                total_energy += rho * dv;

                // Track peak density
                if rho.abs() > peak_density.abs() {
                    peak_density = rho;
                }

                // Track energy type
                if rho < 0.0 {
                    has_exotic = true;
                }
                if rho > 0.0 {
                    has_positive = true;
                }
            }
        }
    }

    // Convert to solar masses (1 solar mass = 1.989e30 kg)
    // E = mc² → m = E/c²
    let mass_kg = total_energy.abs() / C2;
    let solar_masses = mass_kg / 1.989e30;

    // Bubble volume (approximate sphere)
    let bubble_volume = (4.0 / 3.0) * PI * r_max.powi(3);

    let energy_type = if has_exotic && !has_positive {
        EnergyType::Exotic
    } else if has_positive && !has_exotic {
        EnergyType::Positive
    } else {
        EnergyType::Mixed
    };

    EnergyRequirements {
        total_energy,
        solar_masses,
        peak_energy_density: peak_density,
        bubble_volume,
        requires_exotic_matter: has_exotic,
        energy_type,
    }
}

/// Estimate energy for classic Alcubierre drive (analytic approximation)
/// E ≈ -c⁴/(32πG) · (vₛ/c)² · (R/σ) · 4πR²
pub fn alcubierre_energy_estimate(config: &WarpDriveConfig) -> f64 {
    let vs = config.velocity;
    let r = config.bubble_radius;
    let sigma = config.wall_thickness;

    let beta = vs / C;
    let surface_area = 4.0 * PI * r * r;

    // Negative energy required
    -(C4 / (32.0 * PI * G)) * beta * beta * (r / sigma) * surface_area
}

/// Subluminal positive-energy configuration analysis (Bobrick & Martire 2021)
/// For v < c and proper metric engineering, can achieve positive energy density
pub fn subluminal_positive_energy_estimate(config: &WarpDriveConfig) -> Result<f64, String> {
    if !config.subluminal {
        return Err("Configuration must be subluminal (v < c)".to_string());
    }

    let gamma = config.lorentz_factor();
    let r = config.bubble_radius;
    let sigma = config.wall_thickness;

    // For subluminal case with optimized metric:
    // E ≈ +mc²γ · (geometry factor)
    // Requires careful metric engineering to avoid exotic matter

    // Ship mass estimate (e.g., starship ~1e6 kg)
    let ship_mass = 1e6;

    // Positive energy proportional to relativistic mass
    let energy = ship_mass * C2 * gamma * (r / sigma);

    Ok(energy)
}

/// Quantum energy density lower bound (quantum inequality constraint)
/// Negative energy is limited by: |ρ| ≤ ћc/(σ⁴)
pub fn quantum_energy_bound(wall_thickness: f64) -> f64 {
    const HBAR: f64 = 1.054571817e-34; // Reduced Planck constant
    (HBAR * C) / wall_thickness.powi(4)
}

