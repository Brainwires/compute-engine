//! Kerr Black Hole (Rotating)
//!
//! Rotating black hole with angular momentum - much more complex than Schwarzschild

use super::BlackHoleConfig;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KerrMetric {
    pub g_tt: f64,
    pub g_rr: f64,
    pub g_theta_theta: f64,
    pub g_phi_phi: f64,
    pub g_t_phi: f64, // Frame dragging term
}

/// Ergosphere radius (stationary limit surface)
pub fn ergosphere_radius(theta: f64, config: &BlackHoleConfig) -> f64 {
    let m = config.mass;

    // config.spin is dimensionless (0-1), convert to geometric units
    let m_geom = super::G * m / super::C2; // M = GM/c²
    let a_geom = config.spin * m_geom; // a = (spin parameter) * M

    // r_ergo = M + √(M² - a²cos²θ) in geometric units
    let cos_theta = theta.cos();
    let discriminant = m_geom * m_geom - a_geom * a_geom * cos_theta * cos_theta;

    m_geom + discriminant.max(0.0).sqrt()
}

/// Frame dragging angular velocity (ZAMO frame)
pub fn frame_dragging_omega(r: f64, theta: f64, config: &BlackHoleConfig) -> f64 {
    let r_s = config.schwarzschild_radius();
    let a = config.spin;

    // Simplified: ω ~ 2Ma/(r³)
    (2.0 * r_s * a) / (r * r * r)
}

