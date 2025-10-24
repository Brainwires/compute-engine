//! Chaotic Attractors

use super::Point3D;
use serde::{Deserialize, Serialize};

/// Lorenz attractor configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LorenzConfig {
    pub sigma: f64,  // Prandtl number (default: 10)
    pub rho: f64,    // Rayleigh number (default: 28)
    pub beta: f64,   // Geometric factor (default: 8/3)
}

impl Default for LorenzConfig {
    fn default() -> Self {
        Self {
            sigma: 10.0,
            rho: 28.0,
            beta: 8.0 / 3.0,
        }
    }
}

/// Lorenz attractor simulation result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AttractorTrajectory {
    pub points: Vec<Point3D>,
    pub times: Vec<f64>,
}

/// Simulate Lorenz attractor
/// dx/dt = σ(y - x)
/// dy/dt = x(ρ - z) - y
/// dz/dt = xy - βz
pub fn lorenz_attractor(
    initial: Point3D,
    config: &LorenzConfig,
    dt: f64,
    steps: usize,
) -> AttractorTrajectory {
    let mut points = vec![initial];
    let mut times = vec![0.0];
    let mut current = initial;

    for i in 1..=steps {
        // RK4 integration
        current = lorenz_rk4_step(current, config, dt);
        points.push(current);
        times.push(i as f64 * dt);
    }

    AttractorTrajectory { points, times }
}

fn lorenz_rk4_step(state: Point3D, config: &LorenzConfig, dt: f64) -> Point3D {
    let k1 = lorenz_derivative(state, config);
    let k2 = lorenz_derivative(
        Point3D::new(
            state.x + k1.x * dt / 2.0,
            state.y + k1.y * dt / 2.0,
            state.z + k1.z * dt / 2.0,
        ),
        config,
    );
    let k3 = lorenz_derivative(
        Point3D::new(
            state.x + k2.x * dt / 2.0,
            state.y + k2.y * dt / 2.0,
            state.z + k2.z * dt / 2.0,
        ),
        config,
    );
    let k4 = lorenz_derivative(
        Point3D::new(
            state.x + k3.x * dt,
            state.y + k3.y * dt,
            state.z + k3.z * dt,
        ),
        config,
    );

    Point3D::new(
        state.x + (k1.x + 2.0 * k2.x + 2.0 * k3.x + k4.x) * dt / 6.0,
        state.y + (k1.y + 2.0 * k2.y + 2.0 * k3.y + k4.y) * dt / 6.0,
        state.z + (k1.z + 2.0 * k2.z + 2.0 * k3.z + k4.z) * dt / 6.0,
    )
}

fn lorenz_derivative(state: Point3D, config: &LorenzConfig) -> Point3D {
    Point3D::new(
        config.sigma * (state.y - state.x),
        state.x * (config.rho - state.z) - state.y,
        state.x * state.y - config.beta * state.z,
    )
}

/// Rössler attractor
/// dx/dt = -y - z
/// dy/dt = x + ay
/// dz/dt = b + z(x - c)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RosslerConfig {
    pub a: f64,  // default: 0.2
    pub b: f64,  // default: 0.2
    pub c: f64,  // default: 5.7
}

impl Default for RosslerConfig {
    fn default() -> Self {
        Self { a: 0.2, b: 0.2, c: 5.7 }
    }
}

pub fn rossler_attractor(
    initial: Point3D,
    config: &RosslerConfig,
    dt: f64,
    steps: usize,
) -> AttractorTrajectory {
    let mut points = vec![initial];
    let mut times = vec![0.0];
    let mut current = initial;

    for i in 1..=steps {
        current = rossler_rk4_step(current, config, dt);
        points.push(current);
        times.push(i as f64 * dt);
    }

    AttractorTrajectory { points, times }
}

fn rossler_rk4_step(state: Point3D, config: &RosslerConfig, dt: f64) -> Point3D {
    let k1 = rossler_derivative(state, config);
    let k2 = rossler_derivative(
        Point3D::new(
            state.x + k1.x * dt / 2.0,
            state.y + k1.y * dt / 2.0,
            state.z + k1.z * dt / 2.0,
        ),
        config,
    );
    let k3 = rossler_derivative(
        Point3D::new(
            state.x + k2.x * dt / 2.0,
            state.y + k2.y * dt / 2.0,
            state.z + k2.z * dt / 2.0,
        ),
        config,
    );
    let k4 = rossler_derivative(
        Point3D::new(
            state.x + k3.x * dt,
            state.y + k3.y * dt,
            state.z + k3.z * dt,
        ),
        config,
    );

    Point3D::new(
        state.x + (k1.x + 2.0 * k2.x + 2.0 * k3.x + k4.x) * dt / 6.0,
        state.y + (k1.y + 2.0 * k2.y + 2.0 * k3.y + k4.y) * dt / 6.0,
        state.z + (k1.z + 2.0 * k2.z + 2.0 * k3.z + k4.z) * dt / 6.0,
    )
}

fn rossler_derivative(state: Point3D, config: &RosslerConfig) -> Point3D {
    Point3D::new(
        -state.y - state.z,
        state.x + config.a * state.y,
        config.b + state.z * (state.x - config.c),
    )
}

/// Logistic map bifurcation diagram
pub fn logistic_map_bifurcation(r_min: f64, r_max: f64, num_r: usize, iterations: usize, transient: usize) -> Vec<(f64, Vec<f64>)> {
    let mut results = vec![];
    let dr = (r_max - r_min) / (num_r as f64 - 1.0);

    for i in 0..num_r {
        let r = r_min + i as f64 * dr;
        let mut x = 0.5; // Initial condition

        // Discard transient
        for _ in 0..transient {
            x = r * x * (1.0 - x);
        }

        // Collect steady-state values
        let mut values = vec![];
        for _ in 0..iterations {
            x = r * x * (1.0 - x);
            values.push(x);
        }

        results.push((r, values));
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lorenz_attractor() {
        let initial = Point3D::new(1.0, 1.0, 1.0);
        let config = LorenzConfig::default();

        let trajectory = lorenz_attractor(initial, &config, 0.01, 100);

        assert_eq!(trajectory.points.len(), 101); // Initial + 100 steps
        assert_eq!(trajectory.times.len(), 101);

        // Check that trajectory doesn't blow up
        for point in &trajectory.points {
            assert!(point.x.is_finite());
            assert!(point.y.is_finite());
            assert!(point.z.is_finite());
        }
    }

    #[test]
    fn test_rossler_attractor() {
        let initial = Point3D::new(1.0, 1.0, 1.0);
        let config = RosslerConfig::default();

        let trajectory = rossler_attractor(initial, &config, 0.01, 100);

        assert_eq!(trajectory.points.len(), 101);

        // Check finite values
        for point in &trajectory.points {
            assert!(point.x.is_finite());
            assert!(point.y.is_finite());
            assert!(point.z.is_finite());
        }
    }

    #[test]
    fn test_logistic_map_fixed_point() {
        // For r < 1, logistic map should converge to 0
        let results = logistic_map_bifurcation(0.5, 0.9, 5, 100, 100);

        for (r, values) in &results {
            // All values should be close to 0
            let avg: f64 = values.iter().sum::<f64>() / values.len() as f64;
            assert!(avg < 0.1, "r={}, avg={}", r, avg);
        }
    }

    #[test]
    fn test_logistic_map_chaos() {
        // For r > 3.57, logistic map is chaotic
        let results = logistic_map_bifurcation(3.7, 3.9, 3, 100, 100);

        for (r, values) in &results {
            // Values should be spread out (chaotic)
            let max: f64 = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let min: f64 = values.iter().cloned().fold(f64::INFINITY, f64::min);
            let range = max - min;
            // At least some variation (not all identical)
            assert!(range > 0.01, "r={}, range={}", r, range);
        }
    }

    #[test]
    fn test_lorenz_butterfly_effect() {
        // Small perturbation should lead to divergence
        let initial1 = Point3D::new(1.0, 1.0, 1.0);
        let initial2 = Point3D::new(1.0, 1.0, 1.001);
        let config = LorenzConfig::default();

        let traj1 = lorenz_attractor(initial1, &config, 0.01, 500);
        let traj2 = lorenz_attractor(initial2, &config, 0.01, 500);

        // Distance should grow over time
        let dist_early = traj1.points[10].distance(&traj2.points[10]);
        let dist_late = traj1.points[500].distance(&traj2.points[500]);

        assert!(dist_late > dist_early);
    }
}
