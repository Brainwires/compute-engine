//! Black Hole Orbital Mechanics

use super::{BlackHoleConfig, C};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CircularOrbit {
    pub radius: f64,
    pub angular_velocity: f64,
    pub orbital_period: f64,
    pub is_stable: bool,
}

/// Calculate circular orbit parameters
pub fn circular_orbit(r: f64, config: &BlackHoleConfig) -> CircularOrbit {
    let r_s = config.schwarzschild_radius();

    // Angular velocity: Ω = √(GM/r³) = √(r_s c²/(2r³))
    let omega = ((r_s * C * C) / (2.0 * r * r * r)).sqrt();

    // Period: T = 2π/Ω
    let period = 2.0 * std::f64::consts::PI / omega;

    // Stable if r > r_isco
    let r_isco = config.isco_radius();
    let is_stable = r > r_isco;

    CircularOrbit {
        radius: r,
        angular_velocity: omega,
        orbital_period: period,
        is_stable,
    }
}

/// Orbital velocity for circular orbit
pub fn orbital_velocity(r: f64, config: &BlackHoleConfig) -> f64 {
    let r_s = config.schwarzschild_radius();
    // v = √(GM/r) = √(r_s c²/(2r))
    ((r_s * C * C) / (2.0 * r)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::black_holes::BlackHoleConfig;

    #[test]
    fn test_circular_orbit() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_isco = bh.isco_radius();

        let orbit = circular_orbit(10.0 * r_isco, &bh);

        assert!(orbit.is_stable);
        assert!(orbit.angular_velocity > 0.0);
        assert!(orbit.orbital_period > 0.0);
    }

    #[test]
    fn test_orbital_velocity() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let v = orbital_velocity(bh.isco_radius(), &bh);

        assert!(v > 0.0);
        assert!(v < C);
    }
}
