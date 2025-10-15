//! Fluid dynamics module
//!
//! Provides computational fluid dynamics capabilities:
//! - Navier-Stokes equation solvers
//! - Analytical flow solutions
//! - Boundary condition handling
//! - Cavity flow simulations
//! - Flow field analysis

pub mod analytical;
pub mod boundary_conditions;
pub mod cavity_flow;
pub mod flow_analysis;
pub mod grid;
pub mod navier_stokes;

// Re-export main types and functions
pub use analytical::*;
pub use boundary_conditions::*;
pub use cavity_flow::*;
pub use flow_analysis::*;
pub use grid::*;
pub use navier_stokes::*;
