//! Physics computation operations
//!
//! This module handles general physics computations including electromagnetism
//! and Fourier series calculations.
//!
//! ## Submodules
//! - `black_holes` - Black hole physics (Schwarzschild, Kerr)
//! - `control_systems` - Control system theory
//! - `cosmology` - Cosmological models and parameters
//! - `electromagnetism` - Electromagnetic field calculations
//! - `gravitational_waves` - Gravitational wave physics
//! - `n_body` - N-body simulations
//! - `nuclear_physics` - Nuclear physics calculations
//! - `plasma` - Plasma physics
//! - `quantum` - Quantum field theory
//! - `quantum_mechanics` - Quantum mechanics
//! - `relativity` - Special and general relativity
//! - `statistical_physics` - Statistical mechanics
//! - `warp_drive` - Alcubierre warp drive physics
//! - `wormholes` - Traversable wormhole physics

use crate::engine::*;

// Physics submodules
pub mod black_holes;
pub mod control_systems;
pub mod cosmology;
pub mod electromagnetism;
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

// Re-export quantum as quantum_physics for backwards compatibility
pub use quantum as quantum_physics;

/// Compute electromagnetism operations
///
/// Handles electromagnetic field calculations including:
/// - Poynting vector: Energy flux in electromagnetic fields
/// - Skin effect: Current distribution in conductors at high frequencies
pub fn compute_em(op: &EMComputation, input: &ComputeInput) -> ToolResult<ComputeOutput> {
    use electromagnetism::*;

    let result_json = match op {
        EMComputation::PoyntingVector => {
            let req: PoyntingRequest = serde_json::from_value(
                serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
            )
            .map_err(|e| format!("Failed to parse Poynting vector request: {}", e))?;

            let result =
                poynting_vector(req).map_err(|e| format!("Poynting vector error: {}", e))?;
            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
        EMComputation::SkinEffect => {
            let req: SkinEffectRequest = serde_json::from_value(
                serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
            )
            .map_err(|e| format!("Failed to parse skin effect request: {}", e))?;

            let result = skin_effect(req).map_err(|e| format!("Skin effect error: {}", e))?;
            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
    };

    Ok(ComputeOutput {
        result: result_json,
        additional: None,
        metadata: None,
    })
}

/// Compute Fourier series
///
/// Calculates the Fourier series representation of a periodic function,
/// decomposing it into sine and cosine components.
pub fn compute_fourier_series(input: &ComputeInput) -> ToolResult<ComputeOutput> {
    use super::transforms::*;

    let req: FourierSeriesRequest = serde_json::from_value(
        serde_json::to_value(&input.parameters)
            .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
    )
    .map_err(|e| format!("Failed to parse Fourier series request: {}", e))?;

    let result = compute_fourier_series(req).map_err(|e| format!("Fourier series error: {}", e))?;

    let result_json =
        serde_json::to_value(result).map_err(|e| format!("Failed to serialize result: {}", e))?;

    Ok(ComputeOutput {
        result: result_json,
        additional: None,
        metadata: None,
    })
}

/// Compute physics operations dispatcher
///
/// Routes physics computation requests to the appropriate specialized handler.
pub fn compute_physics(op: &PhysicsOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
    match op {
        PhysicsOp::Relativity(relativity_op) => {
            super::relativity::compute_relativity(relativity_op, input)
        }
        PhysicsOp::StatisticalPhysics(stat_phys_op) => {
            super::statistical_physics::compute_statistical_physics(stat_phys_op, input)
        }
        PhysicsOp::QuantumMechanics(qm_op) => {
            super::quantum::compute_quantum_mechanics(qm_op, input)
        }
        PhysicsOp::ControlSystems(cs_op) => super::control::compute_control_systems(cs_op, input),
        PhysicsOp::NuclearPhysics(np_op) => super::nuclear::compute_nuclear_physics(np_op, input),
    }
}
