//! Unit tests for pid_controller module

use super::{PIDController, PIDConfig, ziegler_nichols_tuning, cohen_coon_tuning};

#[test]
fn test_pid_creation() {
    let config = PIDConfig {
        kp: 1.0,
        ki: 0.1,
        kd: 0.05,
        ..Default::default()
    };
    let pid = PIDController::new(config);
    assert_eq!(pid.config.kp, 1.0);
}

#[test]
fn test_pid_update() {
    let config = PIDConfig {
        kp: 1.0,
        ki: 0.1,
        kd: 0.05,
        ..Default::default()
    };
    let mut pid = PIDController::new(config);

    let output = pid.update(10.0, 5.0, 0.1);
    assert!(output.is_finite());
}

#[test]
fn test_pid_reset() {
    let mut pid = PIDController::new(PIDConfig::default());
    pid.update(10.0, 5.0, 0.1);
    pid.reset();
    assert_eq!(pid.get_integral(), 0.0);
}

#[test]
fn test_ziegler_nichols_tuning() {
    let config = ziegler_nichols_tuning(2.0, 1.0);
    assert!(config.kp > 0.0);
    assert!(config.ki > 0.0);
    assert!(config.kd > 0.0);
}

#[test]
fn test_cohen_coon_tuning() {
    let config = cohen_coon_tuning(1.5, 2.0, 0.5);
    assert!(config.kp.is_finite());
    assert!(config.ki.is_finite());
    assert!(config.kd.is_finite());
}
