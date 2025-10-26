// Unit tests for specialized::control_theory::pid_controller
use computational_engine::specialized::control_theory::pid_controller::*;

use super::*;

    #[test]
    fn test_pid_creation() {
        let config = PIDConfig::default();
        let pid = PIDController::new(config);
        assert_eq!(pid.integral, 0.0);
    }

    #[test]
    fn test_proportional_only() {
        let config = PIDConfig {
            kp: 2.0,
            ki: 0.0,
            kd: 0.0,
            ..Default::default()
        };
        let mut pid = PIDController::new(config);

        let output = pid.update(10.0, 5.0, 0.1);
        // Error = 10 - 5 = 5, P = 2 * 5 = 10
        assert_eq!(output, 10.0);
    }

    #[test]
    fn test_integral_accumulation() {
        let config = PIDConfig {
            kp: 0.0,
            ki: 1.0,
            kd: 0.0,
            ..Default::default()
        };
        let mut pid = PIDController::new(config);

        // Step 1: error = 1, dt = 0.1 → integral = 0.1
        pid.update(1.0, 0.0, 0.1);
        assert!((pid.get_integral() - 0.1).abs() < 1e-10);

        // Step 2: error = 1, dt = 0.1 → integral = 0.2
        pid.update(1.0, 0.0, 0.1);
        assert!((pid.get_integral() - 0.2).abs() < 1e-10);
    }

    #[test]
    fn test_output_saturation() {
        let config = PIDConfig {
            kp: 10.0,
            ki: 0.0,
            kd: 0.0,
            output_limits: Some((-5.0, 5.0)),
            ..Default::default()
        };
        let mut pid = PIDController::new(config);

        // Large error → output should saturate
        let output = pid.update(100.0, 0.0, 0.1);
        assert_eq!(output, 5.0);
    }

    #[test]
    fn test_anti_windup() {
        let config = PIDConfig {
            kp: 0.0,
            ki: 10.0,
            kd: 0.0,
            output_limits: Some((-1.0, 1.0)),
            anti_windup: true,
            ..Default::default()
        };
        let mut pid = PIDController::new(config);

        // Large error with saturation - integral shouldn't wind up
        for _ in 0..10 {
            pid.update(10.0, 0.0, 0.1);
        }

        // Integral should be limited
        assert!(pid.get_integral() < 2.0);
    }

    #[test]
    fn test_reset() {
        let mut pid = PIDController::new(PIDConfig::default());

        pid.update(10.0, 0.0, 0.1);
        pid.reset();

        assert_eq!(pid.integral, 0.0);
        assert_eq!(pid.prev_error, 0.0);
    }

    #[test]
    fn test_ziegler_nichols() {
        let config = ziegler_nichols_tuning(2.0, 1.0);
        assert!((config.kp - 1.2).abs() < 1e-10);
        assert!((config.ki - 2.4).abs() < 1e-10);
        assert!((config.kd - 0.15).abs() < 1e-10);
    }
