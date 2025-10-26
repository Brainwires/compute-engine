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

        // g_rr should be positive for physical metric
        if metric.g_rr > 0.0 && metric.g_rr.is_finite() {
            distance += metric.g_rr.sqrt() * dr.abs();
        }
    }

    distance
}

/// Compute embedding surface z(r) for visualization
/// The wormhole can be visualized as a surface of revolution
/// For Morris-Thorne: (dz/dr)² = (b/r)/(1 - b/r)
pub fn embedding_surface(r: f64, config: &WormholeConfig) -> f64 {
    let r0 = config.throat_radius;
    if r <= r0 {
        return 0.0; // At or inside throat, z = 0
    }

    // Numerical integration from r0 to r
    let n = 100;
    let dr = (r - r0) / n as f64;
    let mut z = 0.0;

    for i in 1..n {  // Start at i=1 to avoid exactly r=r0
        let r_i = r0 + i as f64 * dr;
        let b_i = shape_function(r_i, config);
        let ratio = b_i / r_i;

        // Only integrate where 1 - b/r > 0 (physical region)
        if ratio < 1.0 {
            let integrand_i = ratio / (1.0 - ratio);
            if integrand_i > 0.0 && integrand_i.is_finite() {
                z += integrand_i.sqrt() * dr;
            }
        }
    }

    z
}

