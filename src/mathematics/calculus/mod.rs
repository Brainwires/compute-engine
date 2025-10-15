//! Calculus Module
//!
//! Comprehensive calculus operations including:
//! - Fractional calculus (non-integer order derivatives and integrals)
//! - Special functions (Riemann zeta, elliptic integrals, hypergeometric, etc.)
//! - Variational calculus (Euler-Lagrange equations)
//! - Stochastic calculus (Itô and Stratonovich integrals, SDEs)
//! - Symbolic integration

pub mod fractional;
pub mod special_functions;
pub mod variational;
pub mod stochastic;
pub mod symbolic_integration;

// Re-export main functions for convenience
pub use fractional::*;
pub use special_functions::*;
pub use variational::*;
pub use stochastic::*;
pub use symbolic_integration::*;
