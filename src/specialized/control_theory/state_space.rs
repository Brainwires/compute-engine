//! State-Space Representation
//!
//! State-space models of the form:
//! - ẋ = Ax + Bu (continuous-time)
//! - x[k+1] = Ax[k] + Bu[k] (discrete-time)
//! - y = Cx + Du
//!
//! Where:
//! - x: state vector
//! - u: input vector
//! - y: output vector
//! - A: state matrix
//! - B: input matrix
//! - C: output matrix
//! - D: feedthrough matrix

use nalgebra::{DMatrix, DVector};
use serde::{Deserialize, Serialize};

/// State-space system representation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StateSpace {
    /// State matrix A (n × n)
    pub a: Vec<Vec<f64>>,
    /// Input matrix B (n × m)
    pub b: Vec<Vec<f64>>,
    /// Output matrix C (p × n)
    pub c: Vec<Vec<f64>>,
    /// Feedthrough matrix D (p × m)
    pub d: Vec<Vec<f64>>,
    /// Current state
    pub state: Vec<f64>,
}

impl StateSpace {
    /// Create a new state-space system
    ///
    /// # Example
    /// ```
    /// use computational_engine::specialized::control_theory::StateSpace;
    ///
    /// // Simple integrator: ẋ = u, y = x
    /// let sys = StateSpace::new(
    ///     vec![vec![0.0]],  // A
    ///     vec![vec![1.0]],  // B
    ///     vec![vec![1.0]],  // C
    ///     vec![vec![0.0]],  // D
    /// ).unwrap();
    /// ```
    pub fn new(
        a: Vec<Vec<f64>>,
        b: Vec<Vec<f64>>,
        c: Vec<Vec<f64>>,
        d: Vec<Vec<f64>>,
    ) -> Result<Self, String> {
        // Validate dimensions
        let n = a.len(); // Number of states
        if n == 0 {
            return Err("Matrix A cannot be empty".to_string());
        }

        // Check A is square
        for row in &a {
            if row.len() != n {
                return Err("Matrix A must be square".to_string());
            }
        }

        // Check B dimensions
        let m = b[0].len(); // Number of inputs
        if b.len() != n {
            return Err(format!("Matrix B must have {} rows", n));
        }

        // Check C dimensions
        let p = c.len(); // Number of outputs
        for row in &c {
            if row.len() != n {
                return Err(format!("Matrix C must have {} columns", n));
            }
        }

        // Check D dimensions
        if d.len() != p {
            return Err(format!("Matrix D must have {} rows", p));
        }
        for row in &d {
            if row.len() != m {
                return Err(format!("Matrix D must have {} columns", m));
            }
        }

        Ok(Self {
            a,
            b,
            c,
            d,
            state: vec![0.0; n],
        })
    }

    /// Get number of states
    pub fn num_states(&self) -> usize {
        self.a.len()
    }

    /// Get number of inputs
    pub fn num_inputs(&self) -> usize {
        self.b[0].len()
    }

    /// Get number of outputs
    pub fn num_outputs(&self) -> usize {
        self.c.len()
    }

    /// Set the current state
    pub fn set_state(&mut self, state: Vec<f64>) -> Result<(), String> {
        if state.len() != self.num_states() {
            return Err(format!(
                "State vector must have {} elements",
                self.num_states()
            ));
        }
        self.state = state;
        Ok(())
    }

    /// Compute output y = Cx + Du
    pub fn output(&self, input: &[f64]) -> Result<Vec<f64>, String> {
        if input.len() != self.num_inputs() {
            return Err(format!(
                "Input vector must have {} elements",
                self.num_inputs()
            ));
        }

        let n = self.num_states();
        let p = self.num_outputs();
        let m = self.num_inputs();

        let mut y = vec![0.0; p];

        // y = Cx
        for i in 0..p {
            for j in 0..n {
                y[i] += self.c[i][j] * self.state[j];
            }
        }

        // y += Du
        for i in 0..p {
            for j in 0..m {
                y[i] += self.d[i][j] * input[j];
            }
        }

        Ok(y)
    }

    /// Simulate one step: x_next = Ax + Bu (discrete-time)
    pub fn step(&mut self, input: &[f64], dt: f64) -> Result<Vec<f64>, String> {
        if input.len() != self.num_inputs() {
            return Err(format!(
                "Input vector must have {} elements",
                self.num_inputs()
            ));
        }

        let n = self.num_states();
        let m = self.num_inputs();

        // Discretize if continuous (using forward Euler)
        // x[k+1] = x[k] + dt * (Ax[k] + Bu[k])
        let mut x_dot = vec![0.0; n];

        // x_dot = Ax
        for i in 0..n {
            for j in 0..n {
                x_dot[i] += self.a[i][j] * self.state[j];
            }
        }

        // x_dot += Bu
        for i in 0..n {
            for j in 0..m {
                x_dot[i] += self.b[i][j] * input[j];
            }
        }

        // Update state: x = x + dt * x_dot
        for i in 0..n {
            self.state[i] += dt * x_dot[i];
        }

        // Return output
        self.output(input)
    }

    /// Simulate multiple steps
    pub fn simulate(
        &mut self,
        inputs: &[Vec<f64>],
        dt: f64,
    ) -> Result<Vec<Vec<f64>>, String> {
        let mut outputs = Vec::new();

        for input in inputs {
            let output = self.step(input, dt)?;
            outputs.push(output);
        }

        Ok(outputs)
    }

    /// Convert to nalgebra matrices for advanced operations
    pub fn to_matrices(&self) -> (DMatrix<f64>, DMatrix<f64>, DMatrix<f64>, DMatrix<f64>) {
        let n = self.num_states();
        let m = self.num_inputs();
        let p = self.num_outputs();

        let a_flat: Vec<f64> = self.a.iter().flatten().copied().collect();
        let b_flat: Vec<f64> = self.b.iter().flatten().copied().collect();
        let c_flat: Vec<f64> = self.c.iter().flatten().copied().collect();
        let d_flat: Vec<f64> = self.d.iter().flatten().copied().collect();

        (
            DMatrix::from_row_slice(n, n, &a_flat),
            DMatrix::from_row_slice(n, m, &b_flat),
            DMatrix::from_row_slice(p, n, &c_flat),
            DMatrix::from_row_slice(p, m, &d_flat),
        )
    }
}

/// Create a SISO (Single-Input Single-Output) state-space system
pub fn siso_state_space(
    a: Vec<Vec<f64>>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: f64,
) -> Result<StateSpace, String> {
    let n = a.len();

    let b_matrix = b.into_iter().map(|val| vec![val]).collect();
    let c_matrix = vec![c];
    let d_matrix = vec![vec![d]];

    StateSpace::new(a, b_matrix, c_matrix, d_matrix)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_state_space_creation() {
        // Simple integrator
        let sys = StateSpace::new(
            vec![vec![0.0]],
            vec![vec![1.0]],
            vec![vec![1.0]],
            vec![vec![0.0]],
        );

        assert!(sys.is_ok());
        let sys = sys.unwrap();
        assert_eq!(sys.num_states(), 1);
        assert_eq!(sys.num_inputs(), 1);
        assert_eq!(sys.num_outputs(), 1);
    }

    #[test]
    fn test_output_computation() {
        // y = 2x (C = [2], D = [0])
        let sys = StateSpace::new(
            vec![vec![0.0]],
            vec![vec![1.0]],
            vec![vec![2.0]],
            vec![vec![0.0]],
        )
        .unwrap();

        let mut sys_mut = sys.clone();
        sys_mut.set_state(vec![5.0]).unwrap();

        let output = sys_mut.output(&[0.0]).unwrap();
        assert_eq!(output[0], 10.0);
    }

    #[test]
    fn test_integrator_simulation() {
        // Integrator: ẋ = u
        let mut sys = StateSpace::new(
            vec![vec![0.0]],
            vec![vec![1.0]],
            vec![vec![1.0]],
            vec![vec![0.0]],
        )
        .unwrap();

        let dt = 0.1;
        let input = vec![1.0]; // Constant input

        // Step 10 times
        for _ in 0..10 {
            sys.step(&input, dt).unwrap();
        }

        // After 10 steps with dt=0.1 and u=1, x should be ~1.0
        assert!((sys.state[0] - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_siso_creation() {
        let sys = siso_state_space(vec![vec![-1.0]], vec![1.0], vec![1.0], 0.0);

        assert!(sys.is_ok());
        let sys = sys.unwrap();
        assert_eq!(sys.num_states(), 1);
        assert_eq!(sys.num_inputs(), 1);
        assert_eq!(sys.num_outputs(), 1);
    }

    #[test]
    fn test_dimension_validation() {
        // Mismatched B dimensions
        let sys = StateSpace::new(
            vec![vec![0.0, 1.0], vec![-1.0, 0.0]],
            vec![vec![1.0]],        // Should be 2×1, but only 1×1
            vec![vec![1.0, 0.0]],
            vec![vec![0.0]],
        );

        assert!(sys.is_err());
    }
}
