//! Physics Domain
//!
//! Physical simulation and modeling modules

pub mod black_holes;
pub mod control_systems;
pub mod cosmology;
pub mod electromagnetism;
pub mod fluid_dynamics;
pub mod gravitational_waves;
pub mod n_body;
pub mod nuclear_physics;
pub mod plasma;
pub mod quantum;
pub mod quantum_mechanics;
pub mod relativity;
pub mod statistical_physics;
pub mod warp_drive;
pub mod wormholes;

// Re-export quantum_physics as quantum for backwards compatibility
pub use quantum as quantum_physics;
