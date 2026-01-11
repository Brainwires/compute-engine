//! Physics computation operations
//!
//! This module handles general physics computations including electromagnetism
//! and Fourier series calculations.

use crate::engine::*;

/// Compute electromagnetism operations
///
/// Handles electromagnetic field calculations including:
/// - Poynting vector: Energy flux in electromagnetic fields
/// - Skin effect: Current distribution in conductors at high frequencies
pub fn compute_em(op: &EMComputation, input: &ComputeInput) -> ToolResult<ComputeOutput> {
    use crate::physics::electromagnetism::*;

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
    use crate::tools::signal_processing::*;

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
