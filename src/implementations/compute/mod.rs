//! Compute operation modules
//!
//! This module organizes compute operations into domain-specific submodules
//! for better code organization and maintainability.

pub mod datetime;
pub mod engineering;
pub mod geophysics;
pub mod tensor;

// Core mathematical modules
pub mod graph;
pub mod information;
pub mod matrix;
pub mod number_theory;
pub mod scientific;

// Physics modules
pub mod control;
pub mod nuclear;
pub mod physics;
pub mod quantum;
pub mod relativity;
pub mod statistical_physics;

// Re-export the compute functions for convenience
pub use datetime::compute_datetime;
pub use engineering::compute_engineering;
pub use geophysics::compute_geophysics;
pub use tensor::compute_tensor;

// Re-export core mathematical functions
pub use graph::compute_graph;
pub use information::compute_information;
pub use matrix::{compute_matrix_decomp, compute_matrix_op};
pub use number_theory::compute_number_theory;
pub use scientific::{compute_biology, compute_chemistry, compute_optics, compute_thermodynamics};

// Re-export physics compute functions
pub use control::compute_control_systems;
pub use nuclear::compute_nuclear_physics;
pub use physics::{compute_em, compute_fourier_series, compute_physics};
pub use quantum::compute_quantum_mechanics;
pub use relativity::compute_relativity;
pub use statistical_physics::compute_statistical_physics;
