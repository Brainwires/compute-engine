//! Wormhole Energy Requirements
//!
//! Calculate exotic matter requirements using Einstein field equations

use super::{WormholeConfig, G, C, C2, C4};
use super::metric::{shape_function, shape_function_derivative, redshift_function, SphericalCoordinates};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Energy density and stress at a point
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnergyDensity {
    /// Energy density ρ (J/m³)
    pub rho: f64,
    /// Radial pressure p_r
    pub pressure_radial: f64,
    /// Tangential pressure p_t
    pub pressure_tangential: f64,
    /// Total energy density (includes all components)
    pub total_density: f64,
    /// Is exotic matter (negative energy)?
    pub is_exotic: bool,
}

/// Total energy requirements
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WormholeEnergy {
    /// Total exotic matter mass (kg) - can be negative
    pub exotic_mass: f64,
    /// Total exotic matter in solar masses
    pub solar_masses: f64,
    /// Peak energy density (J/m³)
    pub peak_density: f64,
    /// Throat surface area (m²)
    pub throat_area: f64,
    /// Requires exotic matter
    pub requires_exotic: bool,
}

/// Compute energy density at given position
/// Using Einstein equations: T_μν = (c⁴/8πG) G_μν
pub fn compute_energy_density(coords: &SphericalCoordinates, config: &WormholeConfig) -> EnergyDensity {
    let r = coords.r;
    let b = shape_function(r, config);
    let db_dr = shape_function_derivative(r, config);
    let phi = redshift_function(r, config);

    // For Morris-Thorne wormhole, key energy-momentum components:
    // Energy density: ρ = -(c⁴/8πG) · b'/(r²)
    // Radial pressure: p_r = (c⁴/8πG) · b'/(r²)
    // Tangential pressure: p_t = (c⁴/8πG) · [b/(r³) - b'/(r²)]

    let factor = C4 / (8.0 * PI * G);

    let rho = -factor * db_dr / (r * r);
    let pressure_radial = factor * db_dr / (r * r);
    let pressure_tangential = factor * ((b / (r * r * r)) - (db_dr / (r * r)));

    let total_density = rho + 3.0 * pressure_radial.abs(); // Include pressure contributions

    EnergyDensity {
        rho,
        pressure_radial,
        pressure_tangential,
        total_density,
        is_exotic: rho < 0.0,
    }
}

/// Calculate total exotic matter requirement
pub fn calculate_total_energy(config: &WormholeConfig) -> WormholeEnergy {
    let r0 = config.throat_radius;
    let r_max = r0 * 5.0; // Integrate out to 5× throat radius

    let n_r = 100;
    let n_theta = 20;
    let n_phi = 20;

    let dr = (r_max - r0) / n_r as f64;
    let dtheta = PI / n_theta as f64;
    let dphi = 2.0 * PI / n_phi as f64;

    let mut total_mass = 0.0;
    let mut peak_density: f64 = 0.0;
    let mut has_exotic = false;

    // Integrate in spherical coordinates
    for i in 0..n_r {
        let r = r0 + (i as f64 + 0.5) * dr;

        for j in 0..n_theta {
            let theta = (j as f64 + 0.5) * dtheta;

            for _ in 0..n_phi {
                let coords = SphericalCoordinates {
                    t: 0.0,
                    r,
                    theta,
                    phi: 0.0,
                };

                let energy = compute_energy_density(&coords, config);

                // Volume element: r² sin(θ) dr dθ dφ
                let dv = r * r * theta.sin() * dr * dtheta * dphi;

                // Mass: M = ∫ (ρ/c²) dV
                total_mass += (energy.rho / C2) * dv;

                if energy.rho.abs() > peak_density.abs() {
                    peak_density = energy.rho;
                }

                if energy.is_exotic {
                    has_exotic = true;
                }
            }
        }
    }

    // Throat surface area
    let throat_area = 4.0 * PI * r0 * r0;

    // Convert to solar masses
    let solar_masses = total_mass / 1.989e30;

    WormholeEnergy {
        exotic_mass: total_mass,
        solar_masses,
        peak_density,
        throat_area,
        requires_exotic: has_exotic,
    }
}

/// Compare with Schwarzschild black hole of same throat size
/// A wormhole with throat radius r₀ requires less mass than
/// a black hole with Schwarzschild radius r₀
pub fn compare_with_black_hole(config: &WormholeConfig) -> f64 {
    let r0 = config.throat_radius;
    // Schwarzschild radius: r_s = 2GM/c²
    // Mass: M = r_s·c²/(2G)
    let black_hole_mass = r0 * C2 / (2.0 * G);
    black_hole_mass
}

/// Energy per unit throat area (surface energy density)
pub fn throat_energy_density(config: &WormholeConfig) -> f64 {
    let coords = SphericalCoordinates {
        t: 0.0,
        r: config.throat_radius,
        theta: PI / 2.0,
        phi: 0.0,
    };

    let energy = compute_energy_density(&coords, config);
    energy.rho
}

