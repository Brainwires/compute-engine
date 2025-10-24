//! Alcubierre Warp Drive Metric
//!
//! The Alcubierre metric describes a spacetime geometry where a "warp bubble"
//! contracts space in front and expands it behind, allowing apparent faster-than-light travel
//! while locally respecting special relativity.
//!
//! **Alcubierre Metric (in Cartesian coordinates):**
//!
//! ds² = -c²dt² + (dx - vₛf(rₛ)dt)² + dy² + dz²
//!
//! where:
//! - rₛ = √((x - xₛ(t))² + y² + z²) is the distance from bubble center
//! - xₛ(t) = vₛt is the bubble trajectory
//! - f(rₛ) is the "shape function" defining the warp field profile
//! - vₛ is the velocity of the bubble

use super::{WarpDriveConfig, ShapeFunction, C};
use serde::{Deserialize, Serialize};

/// Spacetime coordinates (t, x, y, z)
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Coordinates {
    pub t: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// Metric tensor components g_μν
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetricTensor {
    /// g_tt component
    pub g_tt: f64,
    /// g_tx component
    pub g_tx: f64,
    /// g_xx component
    pub g_xx: f64,
    /// g_yy component
    pub g_yy: f64,
    /// g_zz component
    pub g_zz: f64,
}

/// Shape function f(rₛ) determining warp field profile
pub fn shape_function(rs: f64, config: &WarpDriveConfig) -> f64 {
    let r = config.bubble_radius;
    let sigma = config.wall_thickness;

    match config.shape_function {
        ShapeFunction::TopHat => {
            // Classic Alcubierre step function
            if rs < r - sigma / 2.0 {
                1.0
            } else if rs > r + sigma / 2.0 {
                0.0
            } else {
                0.5
            }
        }
        ShapeFunction::Gaussian => {
            // Smooth Gaussian profile
            (-(rs - r).powi(2) / (2.0 * sigma * sigma)).exp()
        }
        ShapeFunction::Tanh => {
            // Hyperbolic tangent (continuously differentiable)
            0.5 * (1.0 - ((rs - r) / sigma).tanh())
        }
        ShapeFunction::Optimized => {
            // Smoother profile for reduced energy requirements
            if rs < r {
                1.0
            } else {
                let x = (rs - r) / sigma;
                (-x.powi(2) / (1.0 + x.powi(2))).exp()
            }
        }
    }
}

/// First derivative of shape function: df/drₛ
pub fn shape_function_derivative(rs: f64, config: &WarpDriveConfig) -> f64 {
    let r = config.bubble_radius;
    let sigma = config.wall_thickness;
    let eps = 1e-6;

    match config.shape_function {
        ShapeFunction::TopHat => {
            // Discontinuous - approximate with delta function
            0.0 // In practice, this needs regularization
        }
        ShapeFunction::Gaussian => {
            let arg = -(rs - r).powi(2) / (2.0 * sigma * sigma);
            -(rs - r) / (sigma * sigma) * arg.exp()
        }
        ShapeFunction::Tanh => {
            let arg = (rs - r) / sigma;
            -0.5 / sigma * (1.0 / arg.cosh()).powi(2)
        }
        ShapeFunction::Optimized => {
            // Numerical derivative for complex function
            let f_plus = shape_function(rs + eps, config);
            let f_minus = shape_function(rs - eps, config);
            (f_plus - f_minus) / (2.0 * eps)
        }
    }
}

/// Compute Alcubierre metric tensor at given coordinates
pub fn alcubierre_metric(coords: &Coordinates, config: &WarpDriveConfig) -> MetricTensor {
    let vs = config.velocity;

    // Bubble center position: xₛ = vₛ · t
    let xs = vs * coords.t;

    // Distance from bubble center
    let rs = ((coords.x - xs).powi(2) + coords.y.powi(2) + coords.z.powi(2)).sqrt();

    // Shape function value
    let f = shape_function(rs, config);

    // Metric components (Alcubierre 1994)
    // ds² = -c²dt² + (dx - vₛf dt)² + dy² + dz²
    //     = (-c² + vₛ²f²)dt² - 2vₛf dx dt + dx² + dy² + dz²

    let vs_f = vs * f;

    MetricTensor {
        g_tt: -(C * C) + vs_f * vs_f,
        g_tx: -vs_f,
        g_xx: 1.0,
        g_yy: 1.0,
        g_zz: 1.0,
    }
}

/// Compute shift vector β^i (representing frame dragging)
/// β^x = vₛ · f(rₛ), β^y = 0, β^z = 0
pub fn shift_vector(coords: &Coordinates, config: &WarpDriveConfig) -> (f64, f64, f64) {
    let vs = config.velocity;
    let xs = vs * coords.t;
    let rs = ((coords.x - xs).powi(2) + coords.y.powi(2) + coords.z.powi(2)).sqrt();
    let f = shape_function(rs, config);

    (vs * f, 0.0, 0.0)
}

/// Compute lapse function α (time dilation factor)
/// For Alcubierre metric: α = c√(1 - (vₛf/c)²)
pub fn lapse_function(coords: &Coordinates, config: &WarpDriveConfig) -> f64 {
    let vs = config.velocity;
    let xs = vs * coords.t;
    let rs = ((coords.x - xs).powi(2) + coords.y.powi(2) + coords.z.powi(2)).sqrt();
    let f = shape_function(rs, config);

    let beta_ratio = (vs * f) / C;
    C * (1.0 - beta_ratio * beta_ratio).abs().sqrt()
}

/// Calculate proper time τ experienced by observer at given position
/// dτ = α dt for stationary observer
pub fn proper_time(dt: f64, coords: &Coordinates, config: &WarpDriveConfig) -> f64 {
    let alpha = lapse_function(coords, config);
    (alpha / C) * dt
}

/// Calculate coordinate distance that appears contracted/expanded
pub fn effective_distance(start_x: f64, end_x: f64, config: &WarpDriveConfig, t: f64) -> f64 {
    // Integrate metric to get effective distance
    let n_steps = 100;
    let dx = (end_x - start_x) / n_steps as f64;
    let mut distance = 0.0;

    for i in 0..n_steps {
        let x = start_x + i as f64 * dx;
        let coords = Coordinates {
            t,
            x,
            y: 0.0,
            z: 0.0,
        };
        let metric = alcubierre_metric(&coords, config);
        // Spatial distance element: ds_spatial = √(g_xx) dx
        distance += metric.g_xx.sqrt() * dx.abs();
    }

    distance
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::warp_drive::WarpDriveConfig;

    #[test]
    fn test_shape_function_tanh() {
        let config = WarpDriveConfig::new(1e6, 100.0, 10.0);

        // At center: f should be ~1
        let f_center = shape_function(0.0, &config);
        assert!(f_center > 0.9);

        // At bubble edge: f should be ~0.5
        let f_edge = shape_function(100.0, &config);
        assert!((f_edge - 0.5).abs() < 0.1);

        // Far away: f should be ~0
        let f_far = shape_function(200.0, &config);
        assert!(f_far < 0.1);
    }

    #[test]
    fn test_shape_function_derivative() {
        let config = WarpDriveConfig::new(1e6, 100.0, 10.0);

        // Derivative should be negative (function decreasing from center to edge)
        let df = shape_function_derivative(100.0, &config);
        assert!(df < 0.0);
    }

    #[test]
    fn test_metric_flat_far_from_bubble() {
        let config = WarpDriveConfig::new(1e6, 100.0, 10.0);
        let coords = Coordinates {
            t: 0.0,
            x: 1000.0, // Far from bubble
            y: 0.0,
            z: 0.0,
        };

        let metric = alcubierre_metric(&coords, &config);

        // Should be approximately Minkowski metric far from bubble
        assert!((metric.g_tt + C * C).abs() < 1e6); // g_tt ≈ -c²
        assert!(metric.g_tx.abs() < 1e3);           // g_tx ≈ 0
        assert!((metric.g_xx - 1.0).abs() < 0.1);   // g_xx ≈ 1
    }

    #[test]
    fn test_lapse_function_positive() {
        let config = WarpDriveConfig::subluminal(0.1 * C, 100.0, 10.0);
        let coords = Coordinates {
            t: 0.0,
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };

        let alpha = lapse_function(&coords, &config);
        assert!(alpha > 0.0);
        assert!(alpha <= C);
    }

    #[test]
    fn test_shift_vector_inside_bubble() {
        let config = WarpDriveConfig::new(1e7, 100.0, 10.0);
        let coords = Coordinates {
            t: 0.0,
            x: 0.0, // Inside bubble
            y: 0.0,
            z: 0.0,
        };

        let (beta_x, beta_y, beta_z) = shift_vector(&coords, &config);

        assert!(beta_x.abs() > 0.0); // Should have x-component
        assert_eq!(beta_y, 0.0);      // No y-component
        assert_eq!(beta_z, 0.0);      // No z-component
    }

    #[test]
    fn test_proper_time_dilation() {
        let config = WarpDriveConfig::subluminal(0.5 * C, 100.0, 10.0);
        let coords = Coordinates {
            t: 0.0,
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };

        let dt = 1.0; // 1 second coordinate time
        let dtau = proper_time(dt, &coords, &config);

        // Proper time should be less than coordinate time (time dilation)
        assert!(dtau <= dt);
        assert!(dtau > 0.0);
    }

    #[test]
    fn test_effective_distance() {
        let config = WarpDriveConfig::new(1e7, 100.0, 10.0);

        // Measure effective distance across flat space (far from bubble)
        let dist = effective_distance(1000.0, 1100.0, &config, 0.0);

        // Should be approximately 100m in flat space
        assert!((dist - 100.0).abs() < 10.0);
    }
}
