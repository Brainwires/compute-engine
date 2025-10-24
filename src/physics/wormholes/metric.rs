//! Morris-Thorne Wormhole Metric
//!
//! The metric for a traversable wormhole in spherical coordinates:
//! ds² = -e^(2Φ(r))c²dt² + dr²/(1-b(r)/r) + r²(dθ² + sin²θ dφ²)

use super::{WormholeConfig, RedshiftFunction, ShapeFunction, C};
use serde::{Deserialize, Serialize};

/// Spherical coordinates (t, r, θ, φ)
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SphericalCoordinates {
    pub t: f64,
    pub r: f64,
    pub theta: f64,
    pub phi: f64,
}

/// Wormhole metric tensor components
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WormholeMetric {
    /// g_tt component
    pub g_tt: f64,
    /// g_rr component
    pub g_rr: f64,
    /// g_θθ component
    pub g_theta_theta: f64,
    /// g_φφ component
    pub g_phi_phi: f64,
}

/// Compute redshift function Φ(r)
pub fn redshift_function(r: f64, config: &WormholeConfig) -> f64 {
    match config.redshift_type {
        RedshiftFunction::Zero => 0.0,
        RedshiftFunction::Constant(phi0) => phi0,
        RedshiftFunction::Exponential { phi0, scale } => phi0 * (-r / scale).exp(),
    }
}

/// Compute shape function b(r)
/// Must satisfy: b(r₀) = r₀ (throat condition)
///              b'(r₀) < 1 (flaring-out condition)
pub fn shape_function(r: f64, config: &WormholeConfig) -> f64 {
    let r0 = config.throat_radius;

    match config.shape_type {
        ShapeFunction::Polynomial => {
            // b(r) = r₀²/r (simple Morris-Thorne)
            r0 * r0 / r
        }
        ShapeFunction::Gaussian { sigma } => {
            // Gaussian-like shape with throat at r₀
            let x = (r - r0) / sigma;
            r0 * (1.0 + x * x).powf(-0.5)
        }
        ShapeFunction::SmoothCutoff => {
            // Smooth transition using tanh
            let l = config.length_scale;
            r0 * (1.0 + (r - r0) / l).powf(-1.0)
        }
    }
}

/// Derivative of shape function: db/dr
pub fn shape_function_derivative(r: f64, config: &WormholeConfig) -> f64 {
    let eps = 1e-6;
    let b_plus = shape_function(r + eps, config);
    let b_minus = shape_function(r - eps, config);
    (b_plus - b_minus) / (2.0 * eps)
}

/// Compute Morris-Thorne metric at given coordinates
pub fn morris_thorne_metric(coords: &SphericalCoordinates, config: &WormholeConfig) -> WormholeMetric {
    let r = coords.r;
    let theta = coords.theta;

    let phi_r = redshift_function(r, config);
    let b_r = shape_function(r, config);

    // Metric components
    let g_tt = -(C * C) * (2.0 * phi_r).exp();
    let g_rr = 1.0 / (1.0 - b_r / r);
    let g_theta_theta = r * r;
    let g_phi_phi = r * r * theta.sin().powi(2);

    WormholeMetric {
        g_tt,
        g_rr,
        g_theta_theta,
        g_phi_phi,
    }
}

/// Check if flaring-out condition is satisfied: b'(r) < 1
pub fn satisfies_flaring_condition(r: f64, config: &WormholeConfig) -> bool {
    let db_dr = shape_function_derivative(r, config);
    db_dr < 1.0
}

/// Compute proper distance through wormhole throat
/// l = ∫ √(g_rr) dr from -r_max to +r_max
pub fn proper_distance(r_start: f64, r_end: f64, config: &WormholeConfig, n_steps: usize) -> f64 {
    let dr = (r_end - r_start) / n_steps as f64;
    let mut distance = 0.0;

    for i in 0..n_steps {
        let r = r_start + (i as f64 + 0.5) * dr;
        let coords = SphericalCoordinates {
            t: 0.0,
            r,
            theta: std::f64::consts::PI / 2.0,
            phi: 0.0,
        };
        let metric = morris_thorne_metric(&coords, config);
        distance += metric.g_rr.sqrt() * dr.abs();
    }

    distance
}

/// Compute embedding surface z(r) for visualization
/// The wormhole can be visualized as a surface of revolution
/// For Morris-Thorne: (dz/dr)² = (b/r)/(1 - b/r)
pub fn embedding_surface(r: f64, config: &WormholeConfig) -> f64 {
    let r0 = config.throat_radius;
    if r < r0 {
        return 0.0; // Not defined inside throat
    }

    let b = shape_function(r, config);
    let integrand = (b / r) / (1.0 - b / r);

    if integrand < 0.0 {
        0.0
    } else {
        // Numerical integration from r0 to r
        let n = 100;
        let dr = (r - r0) / n as f64;
        let mut z = 0.0;

        for i in 0..n {
            let r_i = r0 + (i as f64 + 0.5) * dr;
            let b_i = shape_function(r_i, config);
            let integrand_i = (b_i / r_i) / (1.0 - b_i / r_i);
            if integrand_i > 0.0 {
                z += integrand_i.sqrt() * dr;
            }
        }

        z
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::wormholes::WormholeConfig;

    #[test]
    fn test_redshift_function_zero() {
        let config = WormholeConfig::morris_thorne(100.0);
        let phi = redshift_function(50.0, &config);
        assert_eq!(phi, 0.0);
    }

    #[test]
    fn test_shape_function_at_throat() {
        let config = WormholeConfig::morris_thorne(100.0);
        let r0 = config.throat_radius;
        let b = shape_function(r0, &config);

        // At throat: b(r₀) = r₀
        assert!((b - r0).abs() < 1.0);
    }

    #[test]
    fn test_metric_computation() {
        let config = WormholeConfig::new(100.0, 200.0);
        let coords = SphericalCoordinates {
            t: 0.0,
            r: 150.0,
            theta: std::f64::consts::PI / 2.0,
            phi: 0.0,
        };

        let metric = morris_thorne_metric(&coords, &config);

        // Metric should have correct signs
        assert!(metric.g_tt < 0.0); // Timelike
        assert!(metric.g_rr > 0.0); // Spacelike
        assert!(metric.g_theta_theta > 0.0);
        assert!(metric.g_phi_phi > 0.0);
    }

    #[test]
    fn test_flaring_condition() {
        let config = WormholeConfig::morris_thorne(100.0);

        // Should satisfy flaring condition near throat
        let flares = satisfies_flaring_condition(100.0, &config);
        assert!(flares);
    }

    #[test]
    fn test_proper_distance_positive() {
        let config = WormholeConfig::new(100.0, 200.0);
        let distance = proper_distance(100.0, 200.0, &config, 100);

        // Should be positive and finite
        assert!(distance > 0.0);
        assert!(distance.is_finite());
    }

    #[test]
    fn test_embedding_surface() {
        let config = WormholeConfig::morris_thorne(100.0);

        // At throat, embedding z should be finite
        let z_throat = embedding_surface(100.0, &config);
        assert!(z_throat.is_finite());

        // Far from throat, z should also be finite
        let z_far = embedding_surface(200.0, &config);
        assert!(z_far.is_finite());
    }

    #[test]
    fn test_shape_function_decreases() {
        let config = WormholeConfig::morris_thorne(100.0);

        let b1 = shape_function(100.0, &config);
        let b2 = shape_function(200.0, &config);

        // Shape function should decrease with r for polynomial type
        assert!(b2 < b1);
    }
}
