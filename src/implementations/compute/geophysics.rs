//! Geophysics computation operations
//!
//! Handles geophysics calculations including seismology, atmospheric science,
//! radiometric dating, and planetary science.

use crate::engine::*;

/// Compute geophysics operations
pub fn compute_geophysics(op: &GeophysicsOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
    use crate::geophysics::*;

    let category = match op {
        GeophysicsOp::Seismology => GeophysicsCategory::Seismology,
        GeophysicsOp::Atmosphere => GeophysicsCategory::Atmosphere,
        GeophysicsOp::RadiometricDating => GeophysicsCategory::RadiometricDating,
        GeophysicsOp::PlanetaryScience => GeophysicsCategory::PlanetaryScience,
    };

    let params: GeophysicsParams = serde_json::from_value(
        serde_json::to_value(&input.parameters)
            .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
    )
    .map_err(|e| format!("Failed to parse geophysics parameters: {}", e))?;

    let geo_input = GeophysicsInput {
        category,
        parameters: params,
    };

    let result = calculate_geophysics(geo_input)
        .map_err(|e| format!("Geophysics calculation error: {}", e))?;

    Ok(ComputeOutput {
        result: serde_json::json!({
            "value": result.value,
            "unit": result.unit,
            "formula_used": result.formula_used,
            "uncertainty": result.uncertainty,
            "interpretation": result.interpretation,
            "additional_data": result.additional_data
        }),
        additional: None,
        metadata: None,
    })
}
