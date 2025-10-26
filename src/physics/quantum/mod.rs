//! Quantum physics module
//!
//! Provides quantum physics simulation capabilities:
//! - Quantum state evolution and measurement
//! - Particle physics simulations
//! - Quantum Field Theory (QFT): propagators, Feynman rules, cross sections
//! - Particle decays and lifetimes
//! - Symbolic quantum mathematics
//! - GPU acceleration support (when enabled)
//!
//! For non-relativistic quantum mechanics, see physics::quantum_mechanics
//! Note: Tensor computations have been moved to mathematics::tensor_calculus::quantum_tensors

pub mod cross_sections;
pub mod feynman_rules;
pub mod particle_decays;
pub mod particle_physics;
pub mod propagators;
pub mod quantum_simulation;
pub mod symbolic_mathematics;

// Re-export main types and functions
pub use cross_sections::*;
pub use feynman_rules::*;
pub use particle_decays::*;
pub use particle_physics::*;
pub use propagators::*;
pub use quantum_simulation::*;
pub use symbolic_mathematics::*;
