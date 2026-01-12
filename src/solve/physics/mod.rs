//! Physics-Specific Solvers
//!
//! This module contains solvers for various physics domains:
//! - Einstein field equations (general relativity)
//! - Electromagnetic equations
//! - Fluid dynamics
//! - Chemical equations

pub mod chemical;
pub mod einstein;
pub mod electromagnetic;
pub mod fluid;

pub use chemical::solve_chemical;
pub use einstein::solve_einstein;
pub use electromagnetic::solve_electromagnetic;
pub use fluid::solve_fluid;
