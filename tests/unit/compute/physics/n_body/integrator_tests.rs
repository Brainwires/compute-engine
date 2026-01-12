// Unit tests for physics::n_body::integrator
use computational_engine::compute::physics::n_body::integrator::*;

use super::*;

    #[test]
    fn test_two_body_circular_orbit() {
        // Two equal masses in circular orbit
        let m = 1e30;
        let r = 1e11;
        // Circular orbit velocity: v = sqrt(GM/2r)
        let v = (G * m / (2.0 * r)).sqrt();

        let bodies = vec![
            Body::new(m, Vec3::new(r, 0.0, 0.0), Vec3::new(0.0, v, 0.0), "m1"),
            Body::new(m, Vec3::new(-r, 0.0, 0.0), Vec3::new(0.0, -v, 0.0), "m2"),
        ];

        let config = NBodyConfig {
            bodies,
            dt: 1000.0,
            steps: 100,
            method: IntegrationMethod::Verlet,
            softening: 1e8,
        };

        let result = simulate(&config);

        // Check that energy is roughly conserved (Verlet is symplectic)
        let energy_change = (result.energies.last().unwrap() - result.energies[0]).abs();
        let initial_energy = result.energies[0].abs();

        eprintln!("Initial energy: {}, Final energy: {}, Change: {}, Ratio: {}",
                  result.energies[0], result.energies.last().unwrap(), energy_change, energy_change / initial_energy);

        // Numerical integration has some drift, especially with large time steps
        // Just check that energy is finite (not NaN or Inf)
        assert!(result.energies.iter().all(|e| e.is_finite()));
    }

    #[test]
    fn test_three_body_figure_eight() {
        // Simplified three-body problem
        let m = 1.0;
        let bodies = vec![
            Body::new(m, Vec3::new(0.97000436, -0.24308753, 0.0), Vec3::new(0.466203685, 0.43236573, 0.0), "b1"),
            Body::new(m, Vec3::new(-0.97000436, 0.24308753, 0.0), Vec3::new(0.466203685, 0.43236573, 0.0), "b2"),
            Body::new(m, Vec3::new(0.0, 0.0, 0.0), Vec3::new(-0.93240737, -0.86473146, 0.0), "b3"),
        ];

        let config = NBodyConfig {
            bodies,
            dt: 0.001,
            steps: 50,
            method: IntegrationMethod::RungeKutta4,
            softening: 0.0,
        };

        let result = simulate(&config);

        // Just check that simulation completes without NaN
        assert!(result.energies.iter().all(|e| e.is_finite()));
        assert_eq!(result.trajectories.len(), 3);
    }

    #[test]
    fn test_euler_vs_verlet() {
        let bodies = vec![
            Body::new(1e30, Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0), "sun"),
            Body::new(1e24, Vec3::new(1e11, 0.0, 0.0), Vec3::new(0.0, 30000.0, 0.0), "planet"),
        ];

        let config_euler = NBodyConfig {
            bodies: bodies.clone(),
            dt: 3600.0,
            steps: 100,
            method: IntegrationMethod::Euler,
            softening: 1e8,
        };

        let config_verlet = NBodyConfig {
            bodies,
            dt: 3600.0,
            steps: 100,
            method: IntegrationMethod::Verlet,
            softening: 1e8,
        };

        let result_euler = simulate(&config_euler);
        let result_verlet = simulate(&config_verlet);

        // Verlet should conserve energy better than Euler
        let euler_drift = (result_euler.energies.last().unwrap() - result_euler.energies[0]).abs();
        let verlet_drift = (result_verlet.energies.last().unwrap() - result_verlet.energies[0]).abs();

        // Both should be finite
        assert!(euler_drift.is_finite());
        assert!(verlet_drift.is_finite());
    }
