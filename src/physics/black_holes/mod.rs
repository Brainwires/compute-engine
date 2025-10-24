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
    /// Angular momentum parameter a = J/(Mc) (dimensionless, 0 ≤ a ≤ M)
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
    pub fn kerr(mass: f64, spin: f64) -> Self {
        let spin_normalized = spin.min(mass); // Can't exceed extremal
        let bh_type = if spin_normalized.abs() < mass * 0.1 {
            BlackHoleType::SlowRotating
        } else if (spin_normalized - mass).abs() < mass * 0.01 {
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
            let m_geom = G * self.mass / C2; // Geometric mass
            let a_geom = self.spin / (self.mass * C); // Geometric angular momentum
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
            let a_normalized = self.spin / self.mass;
            // Prograde orbit: smaller ISCO for spinning BH
            r_s * (3.0 - 2.0 * a_normalized)
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
mod tests {
    use super::*;

    #[test]
    fn test_schwarzschild_creation() {
        let solar_mass = 1.989e30;
        let bh = BlackHoleConfig::schwarzschild(solar_mass);

        assert_eq!(bh.mass, solar_mass);
        assert_eq!(bh.spin, 0.0);
        assert_eq!(bh.bh_type, BlackHoleType::Schwarzschild);
    }

    #[test]
    fn test_schwarzschild_radius() {
        let solar_mass = 1.989e30;
        let bh = BlackHoleConfig::schwarzschild(solar_mass);
        let r_s = bh.schwarzschild_radius();

        // Solar mass black hole: r_s ≈ 3 km
        assert!((r_s - 2953.0).abs() < 100.0);
    }

    #[test]
    fn test_event_horizon() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_h = bh.event_horizon_radius();

        // Should equal Schwarzschild radius for non-rotating
        assert!((r_h - bh.schwarzschild_radius()).abs() < 1.0);
    }

    #[test]
    fn test_photon_sphere() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_ph = bh.photon_sphere_radius();
        let r_s = bh.schwarzschild_radius();

        // Photon sphere at 1.5 r_s
        assert!((r_ph - 1.5 * r_s).abs() < 10.0);
    }

    #[test]
    fn test_isco() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_isco = bh.isco_radius();
        let r_s = bh.schwarzschild_radius();

        // ISCO at 3 r_s (6 GM/c²)
        assert!((r_isco - 3.0 * r_s).abs() < 10.0);
    }

    #[test]
    fn test_hawking_temperature() {
        let solar_mass = 1.989e30;
        let bh = BlackHoleConfig::schwarzschild(solar_mass);
        let temp = bh.hawking_temperature();

        // Solar mass BH: T ~ 60 nanokelvin
        assert!(temp > 0.0);
        assert!(temp < 1e-6); // Very cold
    }

    #[test]
    fn test_kerr_creation() {
        let bh = BlackHoleConfig::kerr(1e30, 0.5e30);

        assert!(bh.spin > 0.0);
        assert!(bh.bh_type != BlackHoleType::Schwarzschild);
    }

    #[test]
    fn test_extremal_limit() {
        let mass = 1e30;
        let bh = BlackHoleConfig::kerr(mass, mass * 2.0); // Exceeds limit

        // Should be capped at extremal
        assert!(bh.spin <= mass);
    }

    #[test]
    fn test_surface_gravity() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let kappa = bh.surface_gravity();

        // Should be positive and finite
        assert!(kappa > 0.0);
        assert!(kappa.is_finite());
    }

    #[test]
    fn test_entropy() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let entropy = bh.entropy();

        // Entropy should be huge for macroscopic black hole
        assert!(entropy > 1e50);
    }

    #[test]
    fn test_evaporation_time() {
        let solar_mass = 1.989e30;
        let bh = BlackHoleConfig::schwarzschild(solar_mass);
        let t_evap = bh.evaporation_time();

        // Solar mass BH: ~10^67 years >> age of universe
        assert!(t_evap > 1e60); // seconds
    }
}
