//! Quantum physics module
//!
//! Provides quantum physics simulation capabilities:
//! - Quantum state evolution and measurement
//! - Particle physics simulations
//! - Symbolic quantum mathematics
//! - GPU acceleration support (when enabled)
//!
//! Note: Tensor computations have been moved to mathematics::tensor_calculus::quantum_tensors

pub mod particle_physics;
pub mod quantum_simulation;
pub mod symbolic_mathematics;

// Re-export main types and functions
pub use particle_physics::*;
pub use quantum_simulation::*;
pub use symbolic_mathematics::*;
