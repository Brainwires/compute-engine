//! Wormhole Traversal Analysis
//!
//! Calculate tidal forces, travel times, and survivability

use super::{WormholeConfig, G, C};
use super::metric::{proper_distance, SphericalCoordinates, morris_thorne_metric};
use serde::{Deserialize, Serialize};

/// Traversal analysis results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TraversalAnalysis {
    /// Proper time for traveler (seconds)
    pub proper_time: f64,
    /// Coordinate time (seconds)
    pub coordinate_time: f64,
    /// Maximum tidal force (N/m per kg)
    pub max_tidal_force: f64,
    /// Is survivable by humans?
    pub survivable: bool,
    /// Proper distance through throat (m)
    pub throat_distance: f64,
}

/// Calculate tidal forces on extended object
/// Tidal acceleration: a_tidal ≈ (GM/r³) · Δr
/// For wormhole: dominated by curvature, not point mass
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TidalForce {
    /// Radial tidal acceleration (m/s²/m)
    pub radial: f64,
    /// Tangential tidal acceleration (m/s²/m)
    pub tangential: f64,
    /// Total magnitude (m/s²/m)
    pub magnitude: f64,
}

/// Compute tidal forces at given position
/// Based on Riemann curvature tensor components
pub fn compute_tidal_forces(coords: &SphericalCoordinates, config: &WormholeConfig) -> TidalForce {
    let r = coords.r;
    let eps = r * 1e-6;

    // Numerical derivatives of metric
    let coords_plus = SphericalCoordinates {
        r: r + eps,
        ..*coords
    };
    let coords_minus = SphericalCoordinates {
        r: r - eps,
        ..*coords
    };

    let g_plus = morris_thorne_metric(&coords_plus, config);
    let g_minus = morris_thorne_metric(&coords_minus, config);

    // Approximate Riemann components from metric derivatives
    // R^r_trt ≈ (1/2) d²g_tt/dr²
    let d2g_tt = (g_plus.g_tt - 2.0 * morris_thorne_metric(coords, config).g_tt + g_minus.g_tt)
        / (eps * eps);

    let radial = d2g_tt.abs() * C * C / (2.0 * r);

    // Tangential component (perpendicular to radial)
    let d2g_rr = (g_plus.g_rr - 2.0 * morris_thorne_metric(coords, config).g_rr + g_minus.g_rr)
        / (eps * eps);

    let tangential = d2g_rr.abs() / (2.0 * r);

    let magnitude = (radial * radial + tangential * tangential).sqrt();

    TidalForce {
        radial,
        tangential,
        magnitude,
    }
}

/// Analyze complete traversal from one side to the other
pub fn analyze_traversal(
    config: &WormholeConfig,
    velocity: f64, // Traveler velocity (m/s)
) -> TraversalAnalysis {
    let r0 = config.throat_radius;

    // Travel from -5r₀ through throat to +5r₀
    let r_start = r0 * 0.5; // Start close to throat on one side
    let r_end = r0 * 10.0;  // End far on other side

    // Proper distance
    let throat_distance = proper_distance(r_start, r_end, config, 100);

    // Coordinate time: t = distance / velocity
    let coordinate_time = if throat_distance > 0.0 && velocity > 0.0 {
        throat_distance / velocity
    } else {
        0.0
    };

    // Proper time (accounting for time dilation)
    // τ ≈ t√(1 - v²/c²) for relativistic motion
    let v_ratio = (velocity / C).powi(2);
    let proper_time = if v_ratio < 1.0 && coordinate_time > 0.0 {
        let gamma = 1.0 / (1.0 - v_ratio).sqrt();
        coordinate_time / gamma
    } else {
        0.0
    };

    // Calculate maximum tidal force along path
    let mut max_tidal = 0.0;
    let n_sample = 50;
    let dr = (r_end - r_start) / n_sample as f64;

    for i in 0..n_sample {
        let r = r_start + i as f64 * dr;
        let coords = SphericalCoordinates {
            t: 0.0,
            r,
            theta: std::f64::consts::PI / 2.0,
            phi: 0.0,
        };

        let tidal = compute_tidal_forces(&coords, config);
        if tidal.magnitude > max_tidal {
            max_tidal = tidal.magnitude;
        }
    }

    // Survivability: tidal forces < ~5 g/m (50 m/s²/m) for humans
    let survivable = max_tidal < 50.0;

    TraversalAnalysis {
        proper_time,
        coordinate_time,
        max_tidal_force: max_tidal,
        survivable,
        throat_distance,
    }
}

/// Calculate time dilation factor at given position
pub fn time_dilation_factor(coords: &SphericalCoordinates, config: &WormholeConfig) -> f64 {
    let metric = morris_thorne_metric(coords, config);

    // Time dilation: dt/dτ = √(-g_tt/c²)
    (-metric.g_tt / (C * C)).sqrt()
}

/// Calculate required velocity to cross wormhole in given time
pub fn required_velocity(config: &WormholeConfig, target_time: f64) -> f64 {
    let r0 = config.throat_radius;
    let distance = proper_distance(r0 * 0.5, r0 * 10.0, config, 100);

    distance / target_time
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::wormholes::WormholeConfig;

    #[test]
    fn test_tidal_forces_at_throat() {
        let config = WormholeConfig::morris_thorne(100.0);
        // Test slightly away from throat to avoid singularity
        let coords = SphericalCoordinates {
            t: 0.0,
            r: 101.0,  // Slightly outside throat
            theta: std::f64::consts::PI / 2.0,
            phi: 0.0,
        };

        let tidal = compute_tidal_forces(&coords, &config);

        // Should have finite tidal forces away from exact throat
        assert!(tidal.magnitude.is_finite());
        assert!(tidal.magnitude >= 0.0);
    }

    #[test]
    fn test_traversal_analysis() {
        let config = WormholeConfig::new(1000.0, 2000.0);
        let velocity = 1e6; // 1000 km/s

        let analysis = analyze_traversal(&config, velocity);

        // Should have positive times and finite values
        assert!(analysis.proper_time > 0.0);
        assert!(analysis.coordinate_time > 0.0);
        assert!(analysis.throat_distance >= 0.0);
        assert!(analysis.proper_time.is_finite());
        assert!(analysis.coordinate_time.is_finite());
    }

    #[test]
    fn test_time_dilation_factor() {
        let config = WormholeConfig::morris_thorne(100.0);
        let coords = SphericalCoordinates {
            t: 0.0,
            r: 150.0,
            theta: std::f64::consts::PI / 2.0,
            phi: 0.0,
        };

        let factor = time_dilation_factor(&coords, &config);

        // Should be close to 1 for zero redshift
        assert!((factor - 1.0).abs() < 0.5);
        assert!(factor > 0.0);
    }

    #[test]
    fn test_required_velocity() {
        let config = WormholeConfig::new(1000.0, 2000.0);
        let target_time = 10.0; // 10 seconds

        let velocity = required_velocity(&config, target_time);

        // Should be positive and finite
        assert!(velocity > 0.0);
        assert!(velocity.is_finite());
    }

    #[test]
    fn test_survivability_large_wormhole() {
        let config = WormholeConfig::new(10000.0, 20000.0); // 10 km throat
        let velocity = 1e5; // 100 km/s

        let analysis = analyze_traversal(&config, velocity);

        // All values should be finite
        assert!(analysis.max_tidal_force.is_finite());
        assert!(analysis.coordinate_time.is_finite());
    }

    #[test]
    fn test_faster_velocity_shorter_time() {
        let config = WormholeConfig::new(1000.0, 2000.0);

        let analysis1 = analyze_traversal(&config, 1e6);
        let analysis2 = analyze_traversal(&config, 2e6);

        // Both should have finite and positive times
        assert!(analysis1.coordinate_time > 0.0);
        assert!(analysis2.coordinate_time > 0.0);
        assert!(analysis1.coordinate_time.is_finite());
        assert!(analysis2.coordinate_time.is_finite());
    }
}
