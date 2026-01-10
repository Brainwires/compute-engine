//! Fluid dynamics module
//!
//! Provides computational fluid dynamics capabilities:
//! - Navier-Stokes equation solvers (2D classical, 1D/2D quantum)
//! - Quantum Navier-Stokes with Bohm potential corrections
//! - Analytical flow solutions
//! - Boundary condition handling
//! - Cavity flow simulations
//! - Flow field analysis
//! - Spectral analysis (energy spectrum E(k))

pub mod analytical;
pub mod boundary_conditions;
pub mod cavity_flow;
pub mod flow_analysis;
pub mod grid;
pub mod navier_stokes;
pub mod navier_stokes_3d;
pub mod quantum_navier_stokes_1d;
pub mod quantum_navier_stokes_2d;
pub mod spectral_analysis;

// Re-export main types and functions
pub use analytical::*;
pub use boundary_conditions::*;
pub use cavity_flow::*;
pub use flow_analysis::*;
pub use grid::*;
pub use navier_stokes::*;
pub use navier_stokes_3d::*;
pub use quantum_navier_stokes_1d::*;
pub use quantum_navier_stokes_2d::*;
pub use spectral_analysis::*;
