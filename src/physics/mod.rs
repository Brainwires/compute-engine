//! Physics Domain
//!
//! Physical simulation and modeling modules

pub mod fluid_dynamics;
pub mod quantum;
pub mod electromagnetism;
pub mod relativity;
pub mod statistical_physics;
pub mod quantum_mechanics;
pub mod control_systems;
pub mod nuclear_physics;

// Re-export quantum_physics as quantum for backwards compatibility
pub use quantum as quantum_physics;
