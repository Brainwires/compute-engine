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

/// Fine structure constant (shared constant for all quantum submodules)
pub const ALPHA_EM: f64 = 1.0 / 137.036;

// Re-export main types and functions
pub use cross_sections::*;
pub use feynman_rules::*;
pub use particle_decays::*;
pub use particle_physics::*;
// Note: propagators re-exported selectively to avoid ALPHA_EM conflict
pub use propagators::{
    coulomb_potential, fermion_propagator_trace, photon_propagator, scalar_propagator,
    vector_boson_propagator, yukawa_potential,
};
pub use quantum_simulation::*;
pub use symbolic_mathematics::*;
