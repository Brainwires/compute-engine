//! Control Theory Module
//!
//! Provides control systems analysis and design tools:
//! - State-space representation
//! - Transfer functions
//! - PID controllers
//! - LQR (Linear Quadratic Regulator)
//! - Kalman filtering
//! - Controllability and observability analysis

pub mod state_space;
pub mod transfer_function;
pub mod pid_controller;
pub mod lqr;
pub mod kalman_filter;
pub mod analysis;

pub use state_space::*;
pub use transfer_function::*;
pub use pid_controller::*;
pub use lqr::*;
pub use kalman_filter::*;
pub use analysis::*;

use serde::{Deserialize, Serialize};

/// Control system types
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum SystemType {
    /// Continuous-time system
    Continuous,
    /// Discrete-time system
    Discrete,
}

/// Controller types
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum ControllerType {
    /// Proportional controller
    P,
    /// Proportional-Integral controller
    PI,
    /// Proportional-Derivative controller
    PD,
    /// Proportional-Integral-Derivative controller
    PID,
}

/// System response characteristics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResponseCharacteristics {
    /// Rise time (10% to 90%)
    pub rise_time: Option<f64>,
    /// Settling time (within 2% of final value)
    pub settling_time: Option<f64>,
    /// Peak time
    pub peak_time: Option<f64>,
    /// Overshoot (percentage)
    pub overshoot: Option<f64>,
    /// Steady-state error
    pub steady_state_error: Option<f64>,
}

/// Stability result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StabilityResult {
    /// Whether the system is stable
    pub is_stable: bool,
    /// Eigenvalues (poles)
    pub eigenvalues: Vec<num_complex::Complex64>,
    /// Stability margin (distance to instability)
    pub stability_margin: Option<f64>,
}
