//! Alcubierre Warp Drive Physics
//!
//! Mathematical framework for warp drive geometries based on General Relativity.
//!
//! Implements:
//! - Classic Alcubierre metric (1994)
//! - Energy-momentum tensor calculations
//! - Energy requirement analysis
//! - Subluminal positive-energy configurations (Bobrick & Martire 2021)
//! - Warp field optimization
//!
//! **Mathematics based on Einstein field equations:**
//! G_μν = 8πG/c⁴ · T_μν
//!
//! where G_μν is the Einstein tensor and T_μν is the stress-energy tensor

use serde::{Deserialize, Serialize};

pub mod metric;
pub mod energy;
pub mod optimization;

pub use metric::*;
pub use energy::*;
pub use optimization::*;

/// Physical constants
pub const C: f64 = 299792458.0; // Speed of light (m/s)
pub const G: f64 = 6.67430e-11; // Gravitational constant (m³/kg/s²)
pub const C2: f64 = C * C;
pub const C4: f64 = C2 * C2;

/// Warp drive configuration parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WarpDriveConfig {
    /// Velocity of the warp bubble (m/s)
    pub velocity: f64,
    /// Radius of the warp bubble (m)
    pub bubble_radius: f64,
    /// Wall thickness parameter (m)
    pub wall_thickness: f64,
    /// Shape function type
    pub shape_function: ShapeFunction,
    /// Subluminal constraint (velocity < c)
    pub subluminal: bool,
}

/// Shape function types for the warp bubble
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum ShapeFunction {
    /// Classic Alcubierre top-hat function
    TopHat,
    /// Smooth Gaussian profile
    Gaussian,
    /// Hyperbolic tangent (smoother)
    Tanh,
    /// Optimized for minimal energy
    Optimized,
}

impl WarpDriveConfig {
    /// Create a new warp drive configuration
    pub fn new(velocity: f64, bubble_radius: f64, wall_thickness: f64) -> Self {
        Self {
            velocity,
            bubble_radius,
            wall_thickness,
            shape_function: ShapeFunction::Tanh,
            subluminal: velocity < C,
        }
    }

    /// Create a subluminal configuration (Bobrick & Martire style)
    pub fn subluminal(velocity: f64, bubble_radius: f64, wall_thickness: f64) -> Self {
        let v = velocity.min(C * 0.99); // Enforce subluminal
        Self {
            velocity: v,
            bubble_radius,
            wall_thickness,
            shape_function: ShapeFunction::Tanh,
            subluminal: true,
        }
    }

    /// Get Lorentz factor γ = 1/√(1 - v²/c²)
    pub fn lorentz_factor(&self) -> f64 {
        let beta = self.velocity / C;
        if beta >= 1.0 {
            return f64::INFINITY; // Superluminal case
        }
        1.0 / (1.0 - beta * beta).sqrt()
    }

    /// Validate configuration
    pub fn validate(&self) -> Result<(), String> {
        if self.velocity <= 0.0 {
            return Err("Velocity must be positive".to_string());
        }
        if self.subluminal && self.velocity >= C {
            return Err("Subluminal configuration requires v < c".to_string());
        }
        if self.bubble_radius <= 0.0 {
            return Err("Bubble radius must be positive".to_string());
        }
        if self.wall_thickness <= 0.0 {
            return Err("Wall thickness must be positive".to_string());
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "../../../../tests/unit/compute/physics/warp_drive_tests.rs"]
mod tests;

