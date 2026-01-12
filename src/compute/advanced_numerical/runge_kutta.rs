//! Runge-Kutta Methods
//!
//! High-order numerical methods for solving ordinary differential equations.

use serde::{Deserialize, Serialize};
use super::ODESystem;

/// Runge-Kutta solver for ODEs
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RungeKuttaSolver {
    /// Order of the method (4 for RK4, 5 for RK45)
    pub order: usize,
    /// Adaptive step size control
    pub adaptive: bool,
    /// Tolerance for adaptive methods
    pub tolerance: f64,
}

impl RungeKuttaSolver {
    /// Create a new RK4 solver
    pub fn rk4() -> Self {
        Self {
            order: 4,
            adaptive: false,
            tolerance: 1e-6,
        }
    }

    /// Create a new adaptive RK45 solver
    pub fn rk45(tolerance: f64) -> Self {
        Self {
            order: 5,
            adaptive: true,
            tolerance,
        }
    }

    /// Solve ODE system using RK4
    pub fn solve_rk4(
        &self,
        f: ODESystem,
        t0: f64,
        tf: f64,
        y0: Vec<f64>,
        dt: f64,
    ) -> Vec<(f64, Vec<f64>)> {
        let mut solution = Vec::new();
        let mut t = t0;
        let mut y = y0;

        solution.push((t, y.clone()));

        while t < tf {
            let h = dt.min(tf - t);
            y = self.rk4_step(f, t, &y, h);
            t += h;
            solution.push((t, y.clone()));
        }

        solution
    }

    /// Single RK4 step
    fn rk4_step(&self, f: ODESystem, t: f64, y: &[f64], h: f64) -> Vec<f64> {
        let n = y.len();

        // k1 = f(t, y)
        let k1 = f(t, y);

        // k2 = f(t + h/2, y + h*k1/2)
        let y2: Vec<f64> = y.iter().zip(&k1).map(|(&yi, &k1i)| yi + h * k1i / 2.0).collect();
        let k2 = f(t + h / 2.0, &y2);

        // k3 = f(t + h/2, y + h*k2/2)
        let y3: Vec<f64> = y.iter().zip(&k2).map(|(&yi, &k2i)| yi + h * k2i / 2.0).collect();
        let k3 = f(t + h / 2.0, &y3);

        // k4 = f(t + h, y + h*k3)
        let y4: Vec<f64> = y.iter().zip(&k3).map(|(&yi, &k3i)| yi + h * k3i).collect();
        let k4 = f(t + h, &y4);

        // y_next = y + h/6 * (k1 + 2*k2 + 2*k3 + k4)
        let mut y_next = vec![0.0; n];
        for i in 0..n {
            y_next[i] = y[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }

        y_next
    }

    /// Solve ODE system using adaptive RK45 (Runge-Kutta-Fehlberg)
    pub fn solve_rk45(
        &self,
        f: ODESystem,
        t0: f64,
        tf: f64,
        y0: Vec<f64>,
        h_initial: f64,
    ) -> Vec<(f64, Vec<f64>)> {
        let mut solution = Vec::new();
        let mut t = t0;
        let mut y = y0;
        let mut h = h_initial;

        solution.push((t, y.clone()));

        while t < tf {
            // Limit step size
            h = h.min(tf - t);

            // Compute both 4th and 5th order solutions
            let (y_4th, y_5th) = self.rk45_step(f, t, &y, h);

            // Estimate error
            let error: f64 = y_4th
                .iter()
                .zip(&y_5th)
                .map(|(&y4, &y5)| (y4 - y5).abs())
                .fold(0.0, f64::max);

            // Check if step is acceptable
            if error < self.tolerance || h < 1e-12 {
                // Accept step
                t += h;
                y = y_5th;
                solution.push((t, y.clone()));

                // Increase step size if error is small
                if error < self.tolerance / 10.0 {
                    h *= 2.0;
                }
            } else {
                // Reject step and decrease step size
                h *= 0.5;
            }
        }

        solution
    }

    /// Single RK45 step (Runge-Kutta-Fehlberg)
    /// Returns both 4th and 5th order solutions
    fn rk45_step(&self, f: ODESystem, t: f64, y: &[f64], h: f64) -> (Vec<f64>, Vec<f64>) {
        let n = y.len();

        // Butcher tableau coefficients for RKF45
        let k1 = f(t, y);

        let y2: Vec<f64> = y
            .iter()
            .zip(&k1)
            .map(|(&yi, &k1i)| yi + h * k1i / 4.0)
            .collect();
        let k2 = f(t + h / 4.0, &y2);

        let y3: Vec<f64> = y
            .iter()
            .enumerate()
            .map(|(i, &yi)| yi + h * (3.0 * k1[i] + 9.0 * k2[i]) / 32.0)
            .collect();
        let k3 = f(t + 3.0 * h / 8.0, &y3);

        let y4: Vec<f64> = y
            .iter()
            .enumerate()
            .map(|(i, &yi)| {
                yi + h * (1932.0 * k1[i] - 7200.0 * k2[i] + 7296.0 * k3[i]) / 2197.0
            })
            .collect();
        let k4 = f(t + 12.0 * h / 13.0, &y4);

        let y5: Vec<f64> = y
            .iter()
            .enumerate()
            .map(|(i, &yi)| {
                yi + h
                    * (439.0 * k1[i] / 216.0 - 8.0 * k2[i] + 3680.0 * k3[i] / 513.0
                        - 845.0 * k4[i] / 4104.0)
            })
            .collect();
        let k5 = f(t + h, &y5);

        let y6: Vec<f64> = y
            .iter()
            .enumerate()
            .map(|(i, &yi)| {
                yi + h
                    * (-8.0 * k1[i] / 27.0 + 2.0 * k2[i] - 3544.0 * k3[i] / 2565.0
                        + 1859.0 * k4[i] / 4104.0
                        - 11.0 * k5[i] / 40.0)
            })
            .collect();
        let k6 = f(t + h / 2.0, &y6);

        // 4th order solution
        let mut y_4th = vec![0.0; n];
        for i in 0..n {
            y_4th[i] = y[i]
                + h * (25.0 * k1[i] / 216.0 + 1408.0 * k3[i] / 2565.0
                    + 2197.0 * k4[i] / 4104.0
                    - k5[i] / 5.0);
        }

        // 5th order solution
        let mut y_5th = vec![0.0; n];
        for i in 0..n {
            y_5th[i] = y[i]
                + h * (16.0 * k1[i] / 135.0 + 6656.0 * k3[i] / 12825.0
                    + 28561.0 * k4[i] / 56430.0
                    - 9.0 * k5[i] / 50.0
                    + 2.0 * k6[i] / 55.0);
        }

        (y_4th, y_5th)
    }
}

