//! Linear Programming Module
//!
//! Solve linear optimization problems using the Simplex algorithm and related methods.
//!
//! # Linear Programming Problem
//!
//! Maximize (or minimize): c^T · x
//! Subject to:
//!   A · x ≤ b  (inequality constraints)
//!   A_eq · x = b_eq  (equality constraints)
//!   x ≥ 0  (non-negativity)
//!
//! Where:
//! - x: decision variables
//! - c: objective coefficients
//! - A, b: inequality constraint matrix and bounds
//! - A_eq, b_eq: equality constraint matrix and bounds

pub mod simplex;
pub mod dual;

pub use simplex::*;
pub use dual::*;

use serde::{Deserialize, Serialize};

/// Linear programming problem specification
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LinearProgram {
    /// Objective function coefficients (to maximize c^T · x)
    pub objective: Vec<f64>,

    /// Inequality constraint matrix A (A · x ≤ b)
    pub inequality_constraints: Option<Vec<Vec<f64>>>,

    /// Inequality constraint bounds b
    pub inequality_bounds: Option<Vec<f64>>,

    /// Equality constraint matrix A_eq (A_eq · x = b_eq)
    pub equality_constraints: Option<Vec<Vec<f64>>>,

    /// Equality constraint bounds b_eq
    pub equality_bounds: Option<Vec<f64>>,

    /// Variable bounds (lower, upper). None means unbounded
    pub variable_bounds: Option<Vec<(Option<f64>, Option<f64>)>>,

    /// Whether to maximize (true) or minimize (false)
    pub maximize: bool,
}

/// Solution to a linear program
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LinearProgramSolution {
    /// Optimal solution vector
    pub solution: Vec<f64>,

    /// Optimal objective value
    pub objective_value: f64,

    /// Status of the solution
    pub status: SolutionStatus,

    /// Number of iterations performed
    pub iterations: usize,

    /// Dual variable values (shadow prices)
    pub dual_variables: Option<Vec<f64>>,

    /// Slack variable values
    pub slack_variables: Option<Vec<f64>>,
}

/// Status of the linear programming solution
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SolutionStatus {
    /// Optimal solution found
    Optimal,

    /// Problem is infeasible (no solution exists)
    Infeasible,

    /// Problem is unbounded (objective can increase indefinitely)
    Unbounded,

    /// Maximum iterations reached without convergence
    MaxIterations,

    /// Numerical error occurred
    NumericalError,
}

impl LinearProgram {
    /// Create a new linear program in standard form
    pub fn new(objective: Vec<f64>, maximize: bool) -> Self {
        Self {
            objective,
            inequality_constraints: None,
            inequality_bounds: None,
            equality_constraints: None,
            equality_bounds: None,
            variable_bounds: None,
            maximize,
        }
    }

    /// Add inequality constraints (A · x ≤ b)
    pub fn with_inequality_constraints(mut self, a: Vec<Vec<f64>>, b: Vec<f64>) -> Self {
        self.inequality_constraints = Some(a);
        self.inequality_bounds = Some(b);
        self
    }

    /// Add equality constraints (A_eq · x = b_eq)
    pub fn with_equality_constraints(mut self, a_eq: Vec<Vec<f64>>, b_eq: Vec<f64>) -> Self {
        self.equality_constraints = Some(a_eq);
        self.equality_bounds = Some(b_eq);
        self
    }

    /// Set variable bounds
    pub fn with_variable_bounds(mut self, bounds: Vec<(Option<f64>, Option<f64>)>) -> Self {
        self.variable_bounds = Some(bounds);
        self
    }

    /// Get number of decision variables
    pub fn num_variables(&self) -> usize {
        self.objective.len()
    }

    /// Get number of inequality constraints
    pub fn num_inequality_constraints(&self) -> usize {
        self.inequality_constraints.as_ref().map_or(0, |c| c.len())
    }

    /// Get number of equality constraints
    pub fn num_equality_constraints(&self) -> usize {
        self.equality_constraints.as_ref().map_or(0, |c| c.len())
    }
}


// Test module
#[cfg(test)]
#[path = "../../../../tests/unit/solve/specialized/linear_programming_tests.rs"]
mod tests;
