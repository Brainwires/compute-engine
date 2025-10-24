//! Schwarzschild Black Hole (Non-Rotating)
//!
//! The simplest black hole solution to Einstein's field equations.
//!
//! **Metric:**
//! ds² = -(1 - r_s/r)c²dt² + dr²/(1 - r_s/r) + r²dΩ²
//! where r_s = 2GM/c² is the Schwarzschild radius

use super::{BlackHoleConfig, C, C2};
use serde::{Deserialize, Serialize};

/// Schwarzschild metric components
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SchwarzschildMetric {
    pub g_tt: f64,
    pub g_rr: f64,
    pub g_theta_theta: f64,
    pub g_phi_phi: f64,
}

/// Calculate Schwarzschild metric at radius r
pub fn schwarzschild_metric(r: f64, config: &BlackHoleConfig) -> SchwarzschildMetric {
    let r_s = config.schwarzschild_radius();

    // Avoid singularity at r = r_s
    let factor = if r > r_s {
        1.0 - r_s / r
    } else {
        -1e-10 // Inside horizon
    };

    SchwarzschildMetric {
        g_tt: -factor * C2,
        g_rr: if factor.abs() > 1e-10 { 1.0 / factor } else { 1e10 },
        g_theta_theta: r * r,
        g_phi_phi: r * r, // At equator (θ = π/2)
    }
}

/// Gravitational time dilation factor
/// τ/t = √(1 - r_s/r) where τ is proper time, t is coordinate time
pub fn time_dilation_factor(r: f64, config: &BlackHoleConfig) -> f64 {
    let r_s = config.schwarzschild_radius();
    if r <= r_s {
        return 0.0; // Time stops at horizon
    }
    (1.0 - r_s / r).sqrt()
}

/// Gravitational redshift factor
/// λ_observed / λ_emitted = 1/√(1 - r_s/r)
pub fn redshift_factor(r: f64, config: &BlackHoleConfig) -> f64 {
    let factor = time_dilation_factor(r, config);
    if factor > 0.0 {
        1.0 / factor
    } else {
        f64::INFINITY
    }
}

/// Escape velocity at radius r
/// v_esc = c√(r_s/r)
pub fn escape_velocity(r: f64, config: &BlackHoleConfig) -> f64 {
    let r_s = config.schwarzschild_radius();
    if r <= r_s {
        return C; // At or inside horizon, need c to escape
    }
    C * (r_s / r).sqrt()
}

/// Tidal force (differential acceleration) across object of length L
/// Δa ≈ 2GM L / r³
pub fn tidal_acceleration(r: f64, length: f64, config: &BlackHoleConfig) -> f64 {
    2.0 * super::G * config.mass * length / (r * r * r)
}

/// Coordinate velocity of free-falling object (dropped from rest at infinity)
/// dr/dt = -c(r_s/r)√(1 - r_s/r)
pub fn freefall_velocity(r: f64, config: &BlackHoleConfig) -> f64 {
    let r_s = config.schwarzschild_radius();
    if r <= r_s {
        return C;
    }
    C * (r_s / r) * (1.0 - r_s / r).sqrt()
}

/// Proper time for radial freefall from r_0 to r
/// (approximate for r >> r_s)
pub fn freefall_time(r_start: f64, r_end: f64, config: &BlackHoleConfig) -> f64 {
    let r_s = config.schwarzschild_radius();

    // Simplified: integrate dt = dr / v(r)
    let n_steps = 100;
    let dr = (r_end - r_start) / n_steps as f64;
    let mut time = 0.0;

    for i in 0..n_steps {
        let r = r_start + (i as f64 + 0.5) * dr;
        if r > r_s {
            let v = freefall_velocity(r, config);
            if v.abs() > 1e-10 {
                time += dr.abs() / v;
            }
        }
    }

    time
}

/// Schwarzschild coordinates to Tortoise coordinates
/// r* = r + r_s ln|r/r_s - 1|
pub fn tortoise_coordinate(r: f64, config: &BlackHoleConfig) -> f64 {
    let r_s = config.schwarzschild_radius();
    if r <= r_s {
        return f64::NEG_INFINITY;
    }
    r + r_s * (r / r_s - 1.0).abs().ln()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::black_holes::BlackHoleConfig;

    #[test]
    fn test_metric_far_from_horizon() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_s = bh.schwarzschild_radius();

        let metric = schwarzschild_metric(10.0 * r_s, &bh);

        // Far from horizon, should approach Minkowski
        assert!(metric.g_tt < 0.0);
        assert!(metric.g_rr > 0.0);
        assert!((metric.g_tt / (-C2) - 0.9).abs() < 0.1);
    }

    #[test]
    fn test_time_dilation_at_horizon() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_s = bh.schwarzschild_radius();

        let factor = time_dilation_factor(r_s, &bh);

        // Time dilation infinite at horizon
        assert_eq!(factor, 0.0);
    }

    #[test]
    fn test_time_dilation_far_away() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_s = bh.schwarzschild_radius();

        let factor = time_dilation_factor(100.0 * r_s, &bh);

        // Should be close to 1 far from BH
        assert!((factor - 1.0).abs() < 0.1);
    }

    #[test]
    fn test_escape_velocity_at_horizon() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_s = bh.schwarzschild_radius();

        let v_esc = escape_velocity(r_s, &bh);

        // At horizon, escape velocity = c
        assert!((v_esc - C).abs() < 1e6);
    }

    #[test]
    fn test_escape_velocity_far_away() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_s = bh.schwarzschild_radius();

        let v_esc = escape_velocity(100.0 * r_s, &bh);

        // Far away, escape velocity should be small
        assert!(v_esc < C * 0.2);
    }

    #[test]
    fn test_tidal_acceleration() {
        let solar_mass = 1.989e30;
        let bh = BlackHoleConfig::schwarzschild(solar_mass);
        let r_s = bh.schwarzschild_radius();

        // Human height: 2m
        let tidal = tidal_acceleration(10.0 * r_s, 2.0, &bh);

        // Should be finite
        assert!(tidal.is_finite());
        assert!(tidal > 0.0);
    }

    #[test]
    fn test_freefall_velocity() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_s = bh.schwarzschild_radius();

        let v = freefall_velocity(3.0 * r_s, &bh);

        // Should be between 0 and c
        assert!(v > 0.0);
        assert!(v < C);
    }

    #[test]
    fn test_tortoise_coordinate() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_s = bh.schwarzschild_radius();

        let r_star = tortoise_coordinate(2.0 * r_s, &bh);

        // Should be finite and positive
        assert!(r_star.is_finite());
    }

    #[test]
    fn test_redshift_factor() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_s = bh.schwarzschild_radius();

        let z = redshift_factor(2.0 * r_s, &bh);

        // Near horizon, significant redshift
        assert!(z > 1.0);
    }
}
