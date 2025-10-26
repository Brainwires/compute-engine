use crate::api::types::ComputationResponse;
use crate::biology::{BiologyInput, calculate_biology};
/**
 * Biology API Handler
 */
use serde_json::{Value, json};

pub fn handle(request: &crate::api::types::ComputationRequest) -> ComputationResponse {
    // Convert HashMap to serde_json::Map
    let params_map: serde_json::Map<String, Value> =
        request.parameters.clone().into_iter().collect();

    // Parse the input from the request parameters
    let input: BiologyInput = match serde_json::from_value(serde_json::Value::Object(params_map)) {
        Ok(input) => input,
        Err(e) => {
            return ComputationResponse::error(
                request.module.clone(),
                request.operation.clone(),
                format!("Invalid biology request: {}", e),
            );
        }
    };

    // Execute calculation
    match calculate_biology(input) {
        Ok(result) => ComputationResponse::success(
            request.module.clone(),
            request.operation.clone(),
            json!({
                "value": result.value,
                "unit": result.unit,
                "formula_used": result.formula_used,
                "interpretation": result.interpretation,
                "additional_data": result.additional_data
            }),
        ),
        Err(e) => ComputationResponse::error(request.module.clone(), request.operation.clone(), e),
    }
}

#[cfg(test)]
#[path = "../../../tests/unit/api/handlers/biology_handler_tests.rs"]
mod tests;
