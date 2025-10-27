//! Black Hole Physics
//!
//! Complete black hole analysis including Schwarzschild (non-rotating) and Kerr (rotating) solutions.
//!
//! **Schwarzschild Metric (non-rotating):**
//! ds² = -(1 - 2GM/rc²)c²dt² + dr²/(1 - 2GM/rc²) + r²(dθ² + sin²θ dφ²)
//!
//! **Kerr Metric (rotating):**
//! Complex metric with frame dragging, ergosphere, and multiple horizons
//!
//! **Key Features:**
//! - Event horizons and photon spheres
//! - Orbital mechanics (ISCO, photon orbits)
//! - Hawking radiation
//! - Penrose process (energy extraction)
//! - Gravitational time dilation
//! - Tidal forces (spaghettification)

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

pub mod schwarzschild;
pub mod kerr;
pub mod orbits;
pub mod hawking;

pub use schwarzschild::*;
pub use kerr::*;
pub use orbits::*;
pub use hawking::*;

use crate::physics::warp_drive::{C, G, C2};

/// Black hole configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlackHoleConfig {
    /// Mass (kg)
    pub mass: f64,
    /// Dimensionless spin parameter (0 ≤ spin ≤ 1), where spin = a/M_geom
    /// and M_geom = GM/c² is the geometric mass
    pub spin: f64,
    /// Black hole type
    pub bh_type: BlackHoleType,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum BlackHoleType {
    /// Non-rotating (Schwarzschild)
    Schwarzschild,
    /// Slowly rotating (|a| << M)
    SlowRotating,
    /// Rapidly rotating (a ~ M)
    RapidRotating,
    /// Extremal (a = M)
    Extremal,
}

impl BlackHoleConfig {
    /// Create Schwarzschild (non-rotating) black hole
    pub fn schwarzschild(mass: f64) -> Self {
        Self {
            mass,
            spin: 0.0,
            bh_type: BlackHoleType::Schwarzschild,
        }
    }

    /// Create Kerr (rotating) black hole
    /// spin: dimensionless spin parameter (0 ≤ spin ≤ 1)
    pub fn kerr(mass: f64, spin: f64) -> Self {
        let spin_normalized = spin.min(1.0).max(0.0); // Clamp to [0, 1]
        let bh_type = if spin_normalized < 0.1 {
            BlackHoleType::SlowRotating
        } else if (spin_normalized - 1.0).abs() < 0.01 {
            BlackHoleType::Extremal
        } else {
            BlackHoleType::RapidRotating
        };

        Self {
            mass,
            spin: spin_normalized,
            bh_type,
        }
    }

    /// Schwarzschild radius: r_s = 2GM/c²
    pub fn schwarzschild_radius(&self) -> f64 {
        2.0 * G * self.mass / C2
    }

    /// Event horizon radius (for Schwarzschild)
    pub fn event_horizon_radius(&self) -> f64 {
        if self.spin.abs() < 1e-10 {
            self.schwarzschild_radius()
        } else {
            // Kerr outer horizon: r+ = M + √(M² - a²) in geometric units
            // where a = spin * M (self.spin is dimensionless)
            let m_geom = G * self.mass / C2; // Geometric mass M = GM/c²
            let a_geom = self.spin * m_geom; // a = (spin parameter) * M
            let discriminant = m_geom * m_geom - a_geom * a_geom;
            m_geom + discriminant.max(0.0).sqrt()
        }
    }

    /// Photon sphere radius (where light can orbit)
    pub fn photon_sphere_radius(&self) -> f64 {
        if self.spin.abs() < 1e-10 {
            // Schwarzschild: r_ph = 3GM/c²
            1.5 * self.schwarzschild_radius()
        } else {
            // Kerr (approximate): depends on spin and orbit direction
            let r_s = self.schwarzschild_radius();
            r_s * (1.5 + 0.5 * (self.spin / self.mass))
        }
    }

    /// ISCO (Innermost Stable Circular Orbit) radius
    pub fn isco_radius(&self) -> f64 {
        if self.spin.abs() < 1e-10 {
            // Schwarzschild: r_isco = 6GM/c²
            3.0 * self.schwarzschild_radius()
        } else {
            // Kerr ISCO (approximate)
            let r_s = self.schwarzschild_radius();
            // self.spin is already dimensionless (0-1)
            // Prograde orbit: smaller ISCO for spinning BH
            r_s * (3.0 - 2.0 * self.spin)
        }
    }

    /// Hawking temperature: T = ℏc³/(8πGMk_B)
    pub fn hawking_temperature(&self) -> f64 {
        const HBAR: f64 = 1.054571817e-34; // Reduced Planck constant
        const K_B: f64 = 1.380649e-23; // Boltzmann constant

        (HBAR * C * C * C) / (8.0 * PI * G * self.mass * K_B)
    }

    /// Black hole lifetime from Hawking radiation
    /// t = 5120πG²M³/(ℏc⁴)
    pub fn evaporation_time(&self) -> f64 {
        const HBAR: f64 = 1.054571817e-34;
        let m3 = self.mass * self.mass * self.mass;
        (5120.0 * PI * G * G * m3) / (HBAR * C2 * C2)
    }

    /// Surface gravity at event horizon
    pub fn surface_gravity(&self) -> f64 {
        let r_h = self.event_horizon_radius();
        // κ = c²/(2r_h) for Schwarzschild
        C2 / (2.0 * r_h)
    }

    /// Surface area of event horizon
    pub fn horizon_area(&self) -> f64 {
        let r_h = self.event_horizon_radius();
        4.0 * PI * r_h * r_h
    }

    /// Bekenstein-Hawking entropy: S = (k_B c³ A)/(4ℏG)
    pub fn entropy(&self) -> f64 {
        const HBAR: f64 = 1.054571817e-34;
        const K_B: f64 = 1.380649e-23;

        let area = self.horizon_area();
        (K_B * C * C * C * area) / (4.0 * HBAR * G)
    }
}

#[cfg(test)]
#[path = "../../../tests/unit/physics/black_holes_tests.rs"]
mod tests;

