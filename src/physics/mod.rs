//! Physics Domain
//!
//! Physical simulation and modeling modules

pub mod control_systems;
pub mod electromagnetism;
pub mod fluid_dynamics;
pub mod nuclear_physics;
pub mod quantum;
pub mod quantum_mechanics;
pub mod relativity;
pub mod statistical_physics;

// Re-export quantum_physics as quantum for backwards compatibility
pub use quantum as quantum_physics;
