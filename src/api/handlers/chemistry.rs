use crate::api::types::ComputationResponse;
use crate::chemistry::{ChemistryInput, calculate_chemistry};
/**
 * Chemistry API Handler - Updated to use new chemistry module
 */
use serde_json::{Value, json};

pub fn handle(request: &crate::api::types::ComputationRequest) -> ComputationResponse {
    // Convert HashMap to serde_json::Map
    let params_map: serde_json::Map<String, Value> =
        request.parameters.clone().into_iter().collect();

    // Parse the input from the request parameters
    let input: ChemistryInput = match serde_json::from_value(serde_json::Value::Object(params_map))
    {
        Ok(input) => input,
        Err(e) => {
            return ComputationResponse::error(
                request.module.clone(),
                request.operation.clone(),
                format!("Invalid chemistry request: {}", e),
            );
        }
    };

    // Execute calculation
    match calculate_chemistry(input) {
        Ok(result) => ComputationResponse::success(
            request.module.clone(),
            request.operation.clone(),
            json!({
                "value": result.value,
                "unit": result.unit,
                "formula_used": result.formula_used,
                "interpretation": result.interpretation
            }),
        ),
        Err(e) => ComputationResponse::error(request.module.clone(), request.operation.clone(), e),
    }
}

// Legacy function for backwards compatibility
pub fn handle_chemistry_request(request: serde_json::Map<String, Value>) -> Value {
    // Parse the input
    let input: ChemistryInput = match serde_json::from_value(serde_json::Value::Object(request)) {
        Ok(input) => input,
        Err(e) => {
            return json!({
                "error": format!("Invalid chemistry request: {}", e),
                "success": false
            });
        }
    };

    // Execute calculation
    match calculate_chemistry(input) {
        Ok(result) => json!({
            "success": true,
            "value": result.value,
            "unit": result.unit,
            "formula": result.formula_used,
            "interpretation": result.interpretation
        }),
        Err(e) => json!({
            "error": e,
            "success": false
        }),
    }
}
