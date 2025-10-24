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
    let r_s = config.schwarzschild_radius();
    let m = config.mass;

    // Convert spin to geometric units: a_geom = a/(Mc)
    let m_geom = super::G * m / super::C2;
    let a_geom = config.spin / (m * super::C);

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::black_holes::BlackHoleConfig;

    #[test]
    fn test_ergosphere() {
        let bh = BlackHoleConfig::kerr(1e30, 0.5e30);

        // At the pole (θ = 0), ergosphere should be larger than horizon
        let r_ergo_pole = ergosphere_radius(0.0, &bh);
        let r_h = bh.event_horizon_radius();

        // At equator (θ = π/2), they should be equal (cos²θ = 0)
        let r_ergo_equator = ergosphere_radius(std::f64::consts::PI / 2.0, &bh);

        eprintln!("Pole: {}, Horizon: {}, Equator: {}", r_ergo_pole, r_h, r_ergo_equator);

        // For a Kerr black hole, ergosphere >= horizon everywhere
        assert!(r_ergo_pole >= r_h - 1.0); // Allow small numerical error
        assert!((r_ergo_equator - r_h).abs() < 1.0); // Close to equal at equator
    }

    #[test]
    fn test_frame_dragging() {
        let bh = BlackHoleConfig::kerr(1e30, 0.5e30);
        let omega = frame_dragging_omega(10e3, std::f64::consts::PI / 2.0, &bh);
        assert!(omega.is_finite());
    }
}
