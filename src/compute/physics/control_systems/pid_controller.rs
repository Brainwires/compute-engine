//! PID Controller Implementation
//!
//! Proportional-Integral-Derivative controller:
//! u(t) = Kp*e(t) + Ki*∫e(τ)dτ + Kd*de(t)/dt
//!
//! Features:
//! - Anti-windup for integral term
//! - Derivative filtering
//! - Setpoint weighting

use serde::{Deserialize, Serialize};

/// PID controller configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PIDConfig {
    /// Proportional gain
    pub kp: f64,
    /// Integral gain
    pub ki: f64,
    /// Derivative gain
    pub kd: f64,
    /// Setpoint weighting for proportional term (0 to 1)
    pub setpoint_weight: f64,
    /// Output limits (min, max)
    pub output_limits: Option<(f64, f64)>,
    /// Anti-windup enabled
    pub anti_windup: bool,
    /// Derivative filter time constant
    pub derivative_filter_tau: Option<f64>,
}

impl Default for PIDConfig {
    fn default() -> Self {
        Self {
            kp: 1.0,
            ki: 0.0,
            kd: 0.0,
            setpoint_weight: 1.0,
            output_limits: None,
            anti_windup: true,
            derivative_filter_tau: None,
        }
    }
}

/// PID controller state
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PIDController {
    /// Configuration
    pub config: PIDConfig,
    /// Integral term accumulator
    integral: f64,
    /// Previous error (for derivative)
    prev_error: f64,
    /// Previous output (for anti-windup)
    prev_output: f64,
    /// Filtered derivative
    filtered_derivative: f64,
}

impl PIDController {
    /// Create a new PID controller
    ///
    /// # Example
    /// ```
    /// use computational_engine::compute::physics::control_systems::{PIDController, PIDConfig};
    ///
    /// let config = PIDConfig {
    ///     kp: 1.0,
    ///     ki: 0.1,
    ///     kd: 0.05,
    ///     ..Default::default()
    /// };
    /// let pid = PIDController::new(config);
    /// ```
    pub fn new(config: PIDConfig) -> Self {
        Self {
            config,
            integral: 0.0,
            prev_error: 0.0,
            prev_output: 0.0,
            filtered_derivative: 0.0,
        }
    }

    /// Reset controller state
    pub fn reset(&mut self) {
        self.integral = 0.0;
        self.prev_error = 0.0;
        self.prev_output = 0.0;
        self.filtered_derivative = 0.0;
    }

    /// Compute control output
    ///
    /// # Arguments
    /// * `setpoint` - Desired value
    /// * `measurement` - Current measured value
    /// * `dt` - Time step
    ///
    /// # Returns
    /// Control output u(t)
    pub fn update(&mut self, setpoint: f64, measurement: f64, dt: f64) -> f64 {
        // Error
        let error = setpoint - measurement;

        // Proportional term with setpoint weighting
        let proportional = self.config.kp
            * (self.config.setpoint_weight * setpoint - measurement);

        // Integral term
        self.integral += error * dt;

        // Anti-windup: don't integrate if output is saturated
        if self.config.anti_windup {
            if let Some((min, max)) = self.config.output_limits {
                if (self.prev_output >= max && error > 0.0)
                    || (self.prev_output <= min && error < 0.0)
                {
                    self.integral -= error * dt; // Undo integration
                }
            }
        }

        let integral = self.config.ki * self.integral;

        // Derivative term
        let derivative_raw = if dt > 0.0 {
            (error - self.prev_error) / dt
        } else {
            0.0
        };

        // Apply derivative filtering if configured
        let derivative = if let Some(tau) = self.config.derivative_filter_tau {
            // First-order low-pass filter
            let alpha = dt / (tau + dt);
            self.filtered_derivative =
                alpha * derivative_raw + (1.0 - alpha) * self.filtered_derivative;
            self.config.kd * self.filtered_derivative
        } else {
            self.config.kd * derivative_raw
        };

        // Compute output
        let mut output = proportional + integral + derivative;

        // Apply output limits
        if let Some((min, max)) = self.config.output_limits {
            output = output.clamp(min, max);
        }

        // Store state for next iteration
        self.prev_error = error;
        self.prev_output = output;

        output
    }

    /// Get current integral term (for monitoring)
    pub fn get_integral(&self) -> f64 {
        self.integral
    }

    /// Manually set integral term (for bumpless transfer)
    pub fn set_integral(&mut self, value: f64) {
        self.integral = value;
    }
}

/// Auto-tune PID using Ziegler-Nichols method
///
/// Based on ultimate gain (Ku) and ultimate period (Pu)
pub fn ziegler_nichols_tuning(ku: f64, pu: f64) -> PIDConfig {
    PIDConfig {
        kp: 0.6 * ku,
        ki: 1.2 * ku / pu,
        kd: 0.075 * ku * pu,
        ..Default::default()
    }
}

/// Auto-tune PID using Cohen-Coon method
///
/// Based on process parameters from step response
pub fn cohen_coon_tuning(k: f64, tau: f64, theta: f64) -> PIDConfig {
    let r = theta / tau;

    PIDConfig {
        kp: (1.0 / k) * (1.35 / r + 0.25),
        ki: (1.0 / k) * (1.35 / r + 0.25) * (tau / (2.5 + 2.0 * r)) * (1.0 + 0.6 * r),
        kd: (1.0 / k) * (1.35 / r + 0.25) * (0.37 * tau / (1.0 + 0.2 * r)),
        ..Default::default()
    }
}

// Test module
#[cfg(test)]
#[path = "../../../../tests/unit/compute/control/control_theory_pid_controller_tests.rs"]
mod tests;

