//! Traversable Wormholes
//!
//! Mathematical framework for Morris-Thorne traversable wormholes based on General Relativity.
//!
//! **Morris-Thorne Metric (1988):**
//! ds² = -e^(2Φ(r))c²dt² + dr²/(1-b(r)/r) + r²(dθ² + sin²θ dφ²)
//!
//! where:
//! - Φ(r) is the "redshift function" (gravitational potential)
//! - b(r) is the "shape function" defining the wormhole throat
//! - r is the proper radial coordinate
//! - b(r₀) = r₀ defines the throat radius
//!
//! **Key Physics:**
//! - Requires exotic matter (negative energy density) to hold throat open
//! - No event horizons (traversable in both directions)
//! - Tidal forces must be survivable
//! - Can connect distant regions of spacetime

use serde::{Deserialize, Serialize};

pub mod metric;
pub mod energy;
pub mod traversal;

pub use metric::*;
pub use energy::*;
pub use traversal::*;

use super::warp_drive::{C, G, C2, C4};

/// Wormhole configuration parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WormholeConfig {
    /// Throat radius (m) - minimum radius of the wormhole
    pub throat_radius: f64,
    /// Length parameter (m) - characteristic scale
    pub length_scale: f64,
    /// Redshift function type
    pub redshift_type: RedshiftFunction,
    /// Shape function type
    pub shape_type: ShapeFunction,
}

/// Redshift function Φ(r) types
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum RedshiftFunction {
    /// Zero redshift (simplest case): Φ(r) = 0
    Zero,
    /// Constant: Φ(r) = Φ₀
    Constant(f64),
    /// Exponential decay: Φ(r) = Φ₀ exp(-r/a)
    Exponential { phi0: f64, scale: f64 },
}

/// Shape function b(r) types
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum ShapeFunction {
    /// Simple polynomial: b(r) = r₀²/r
    Polynomial,
    /// Gaussian shell: b(r) = r₀ exp(-(r-r₀)²/σ²)
    Gaussian { sigma: f64 },
    /// Smooth cutoff
    SmoothCutoff,
}

impl WormholeConfig {
    /// Create a new wormhole configuration
    pub fn new(throat_radius: f64, length_scale: f64) -> Self {
        Self {
            throat_radius,
            length_scale,
            redshift_type: RedshiftFunction::Zero,
            shape_type: ShapeFunction::Polynomial,
        }
    }

    /// Create Morris-Thorne style configuration (zero redshift, polynomial shape)
    pub fn morris_thorne(throat_radius: f64) -> Self {
        Self {
            throat_radius,
            length_scale: throat_radius * 2.0,
            redshift_type: RedshiftFunction::Zero,
            shape_type: ShapeFunction::Polynomial,
        }
    }

    /// Validate configuration
    pub fn validate(&self) -> Result<(), String> {
        if self.throat_radius <= 0.0 {
            return Err("Throat radius must be positive".to_string());
        }
        if self.length_scale <= 0.0 {
            return Err("Length scale must be positive".to_string());
        }
        Ok(())
    }

    /// Get Schwarzschild radius for comparison (2GM/c²)
    /// If throat_radius < r_s, wormhole would collapse to black hole
    pub fn schwarzschild_radius_equivalent(&self, mass: f64) -> f64 {
        2.0 * G * mass / C2
    }
}

/// Physical constants from warp drive module
pub use super::warp_drive::C as SPEED_OF_LIGHT;
pub use super::warp_drive::G as GRAVITATIONAL_CONSTANT;

#[cfg(test)]
#[path = "../../../../tests/unit/physics/wormholes_tests.rs"]
mod tests;

