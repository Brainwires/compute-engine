//! Advanced Numerical Methods
//!
//! Sophisticated numerical techniques for solving differential equations,
//! boundary value problems, and partial differential equations.

use serde::{Deserialize, Serialize};

pub mod runge_kutta;
pub mod finite_element;
pub mod spectral;
pub mod multigrid;

pub use runge_kutta::*;
pub use finite_element::*;
pub use spectral::*;
pub use multigrid::*;

/// Type alias for ODE system: f(t, y) -> dy/dt
pub type ODESystem = fn(f64, &[f64]) -> Vec<f64>;

/// Type alias for PDE: represents spatial discretization
pub type PDEField = Vec<f64>;

/// Boundary condition types
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum BoundaryCondition {
    /// Dirichlet: u(boundary) = value
    Dirichlet(f64),
    /// Neumann: du/dn(boundary) = value
    Neumann(f64),
    /// Robin: a*u + b*du/dn = value
    Robin { a: f64, b: f64, value: f64 },
    /// Periodic boundary conditions
    Periodic,
}

/// Mesh types for numerical discretization
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MeshType {
    /// Uniform grid spacing
    Uniform,
    /// Non-uniform adaptive mesh
    Adaptive,
    /// Chebyshev nodes for spectral methods
    Chebyshev,
}
