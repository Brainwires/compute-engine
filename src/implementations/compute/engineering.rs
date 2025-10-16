//! Engineering computation operations
//!
//! Handles engineering calculations across multiple disciplines including acoustics,
//! materials science, fluid mechanics, and control theory.

use crate::engine::*;

/// Compute engineering operations
pub fn compute_engineering(op: &EngineeringOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
    use crate::engineering::*;

    // Engineering operations are organized by discipline, infer from operation
    let (discipline, params) = match op {
        EngineeringOp::SoundPressureLevel
        | EngineeringOp::DopplerEffect
        | EngineeringOp::ReverberationTime => (
            EngineeringDiscipline::Acoustics,
            serde_json::to_value(&input.parameters)
                .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
        ),
        EngineeringOp::Stress | EngineeringOp::Strain | EngineeringOp::FractureMechanics => (
            EngineeringDiscipline::Materials,
            serde_json::to_value(&input.parameters)
                .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
        ),
        EngineeringOp::Bernoulli | EngineeringOp::Poiseuille | EngineeringOp::Drag => (
            EngineeringDiscipline::FluidMechanics,
            serde_json::to_value(&input.parameters)
                .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
        ),
        EngineeringOp::PidController | EngineeringOp::FirstOrderResponse => (
            EngineeringDiscipline::ControlTheory,
            serde_json::to_value(&input.parameters)
                .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
        ),
    };

    let eng_params: EngineeringParams = serde_json::from_value(params)
        .map_err(|e| format!("Failed to parse engineering parameters: {}", e))?;

    let eng_input = EngineeringInput {
        discipline,
        parameters: eng_params,
    };

    let result = calculate_engineering(eng_input)
        .map_err(|e| format!("Engineering calculation error: {}", e))?;

    Ok(ComputeOutput {
        result: serde_json::json!({
            "value": result.value,
            "unit": result.unit,
            "formula_used": result.formula_used,
            "classification": result.classification,
            "interpretation": result.interpretation,
            "additional": result.additional
        }),
        additional: None,
        metadata: None,
    })
}
