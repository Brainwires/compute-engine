//! Cosmology - Evolution and Structure of the Universe
//!
//! **Key Equations:**
//! - Friedmann equations (expansion of the universe)
//! - FLRW metric (spacetime geometry)
//! - Dark energy equation of state
//! - CMB temperature evolution
//! - Structure formation
//!
//! **References:**
//! - Planck Collaboration (2020) - Cosmological parameters
//! - Friedmann (1922) - Original cosmological solution

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

pub mod friedmann;
pub mod cmb;
pub mod dark_energy;

pub use friedmann::*;
pub use cmb::*;
pub use dark_energy::*;

/// Physical constants
pub const C: f64 = 299792458.0; // m/s
pub const G: f64 = 6.67430e-11; // m³/(kg·s²)
pub const H0: f64 = 67.4; // Hubble constant (km/s/Mpc) - Planck 2018
pub const K_B: f64 = 1.380649e-23; // Boltzmann constant

/// Cosmological parameters (Planck 2018)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CosmologyParams {
    pub omega_m: f64,      // Matter density parameter (0.315)
    pub omega_lambda: f64, // Dark energy density parameter (0.685)
    pub omega_b: f64,      // Baryon density parameter (0.049)
    pub omega_r: f64,      // Radiation density parameter (9.24e-5)
    pub h: f64,            // Hubble parameter H0/100 (0.674)
    pub t_cmb: f64,        // CMB temperature today (K) (2.7255)
}

impl Default for CosmologyParams {
    fn default() -> Self {
        Self {
            omega_m: 0.315,
            omega_lambda: 0.685,
            omega_b: 0.049,
            omega_r: 9.24e-5,
            h: 0.674,
            t_cmb: 2.7255,
        }
    }
}

impl CosmologyParams {
    /// Planck 2018 best-fit parameters
    pub fn planck_2018() -> Self {
        Self::default()
    }

    /// WMAP 9-year parameters
    pub fn wmap9() -> Self {
        Self {
            omega_m: 0.286,
            omega_lambda: 0.714,
            omega_b: 0.046,
            omega_r: 8.4e-5,
            h: 0.693,
            t_cmb: 2.7255,
        }
    }

    /// Total density parameter Ω_total = Ω_m + Ω_λ + Ω_r
    pub fn omega_total(&self) -> f64 {
        self.omega_m + self.omega_lambda + self.omega_r
    }

    /// Curvature parameter Ω_k = 1 - Ω_total
    pub fn omega_k(&self) -> f64 {
        1.0 - self.omega_total()
    }

    /// Hubble constant in SI units (1/s)
    pub fn hubble_constant_si(&self) -> f64 {
        self.h * 100.0 * 1000.0 / (3.086e22) // Convert km/s/Mpc to 1/s
    }

    /// Critical density today (kg/m³)
    pub fn critical_density(&self) -> f64 {
        let h_si = self.hubble_constant_si();
        3.0 * h_si * h_si / (8.0 * PI * G)
    }

    /// Age of the universe (seconds)
    pub fn age_of_universe(&self) -> f64 {
        // Approximate for flat ΛCDM
        let h_si = self.hubble_constant_si();
        (2.0 / 3.0) * (1.0 / h_si) * (1.0 / self.omega_m.sqrt())
    }
}

/// Redshift z = (a₀/a) - 1 where a is scale factor
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Redshift(pub f64);

impl Redshift {
    /// Convert redshift to scale factor: a = 1/(1+z)
    pub fn to_scale_factor(&self) -> f64 {
        1.0 / (1.0 + self.0)
    }

    /// Temperature at this redshift: T(z) = T₀(1+z)
    pub fn temperature(&self, params: &CosmologyParams) -> f64 {
        params.t_cmb * (1.0 + self.0)
    }

    /// Age of universe at this redshift (approximate)
    pub fn age(&self, params: &CosmologyParams) -> f64 {
        let h_si = params.hubble_constant_si();
        let om = params.omega_m;
        let ol = params.omega_lambda;
        let a = self.to_scale_factor();

        // Numerical integration would be more accurate
        // Simplified approximation
        (2.0 / (3.0 * h_si)) * (1.0 / (om * a.powf(-3.0) + ol).sqrt())
    }
}

impl From<f64> for Redshift {
    fn from(z: f64) -> Self {
        Redshift(z)
    }
}

#[cfg(test)]
#[path = "../../../../tests/unit/compute/physics/cosmology_tests.rs"]
mod tests;

