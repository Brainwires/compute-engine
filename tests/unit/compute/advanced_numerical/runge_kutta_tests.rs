// Unit tests for specialized::advanced_numerical::runge_kutta
use computational_engine::compute::advanced_numerical::runge_kutta::*;

use super::*;

    // Test ODE: dy/dt = -2y, solution: y = y0 * exp(-2t)
    fn exponential_decay(_t: f64, y: &[f64]) -> Vec<f64> {
        vec![-2.0 * y[0]]
    }

    // Test ODE: dy/dt = y, solution: y = y0 * exp(t)
    fn exponential_growth(_t: f64, y: &[f64]) -> Vec<f64> {
        vec![y[0]]
    }

    #[test]
    fn test_rk4_exponential_decay() {
        let solver = RungeKuttaSolver::rk4();
        let solution = solver.solve_rk4(exponential_decay, 0.0, 1.0, vec![1.0], 0.1);

        // Check that solution decays
        assert!(solution.last().unwrap().1[0] < solution[0].1[0]);

        // Check approximate accuracy: y(1) ≈ exp(-2)
        let final_value = solution.last().unwrap().1[0];
        let expected = (-2.0_f64).exp();
        assert!((final_value - expected).abs() < 0.01);
    }

    #[test]
    fn test_rk4_exponential_growth() {
        let solver = RungeKuttaSolver::rk4();
        let solution = solver.solve_rk4(exponential_growth, 0.0, 1.0, vec![1.0], 0.1);

        // Check that solution grows
        assert!(solution.last().unwrap().1[0] > solution[0].1[0]);

        // Check approximate accuracy: y(1) ≈ e
        let final_value = solution.last().unwrap().1[0];
        let expected = 1.0_f64.exp();
        assert!((final_value - expected).abs() < 0.01);
    }

    #[test]
    fn test_rk45_adaptive() {
        let solver = RungeKuttaSolver::rk45(1e-6);
        let solution = solver.solve_rk45(exponential_decay, 0.0, 1.0, vec![1.0], 0.1);

        // Adaptive method should produce accurate results
        let final_value = solution.last().unwrap().1[0];
        let expected = (-2.0_f64).exp();
        assert!((final_value - expected).abs() < 1e-4);
    }

    #[test]
    fn test_rk4_system() {
        // Harmonic oscillator: d²x/dt² = -x
        // Convert to system: dy1/dt = y2, dy2/dt = -y1
        fn harmonic_oscillator(_t: f64, y: &[f64]) -> Vec<f64> {
            vec![y[1], -y[0]]
        }

        let solver = RungeKuttaSolver::rk4();
        let solution = solver.solve_rk4(harmonic_oscillator, 0.0, 2.0 * std::f64::consts::PI, vec![1.0, 0.0], 0.1);

        // After one period, should return close to initial position
        let final_pos = solution.last().unwrap().1[0];
        assert!((final_pos - 1.0).abs() < 0.1);
    }
