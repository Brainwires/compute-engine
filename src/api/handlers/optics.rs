use crate::api::types::ComputationResponse;
use crate::optics::{OpticsInput, calculate_optics};
/**
 * Optics API Handler
 */
use serde_json::{Value, json};

pub fn handle(request: &crate::api::types::ComputationRequest) -> ComputationResponse {
    // Convert HashMap to serde_json::Map
    let params_map: serde_json::Map<String, Value> =
        request.parameters.clone().into_iter().collect();

    // Parse the input from the request parameters
    let input: OpticsInput = match serde_json::from_value(serde_json::Value::Object(params_map)) {
        Ok(input) => input,
        Err(e) => {
            return ComputationResponse::error(
                request.module.clone(),
                request.operation.clone(),
                format!("Invalid optics request: {}", e),
            );
        }
    };

    // Execute calculation
    match calculate_optics(input) {
        Ok(result) => ComputationResponse::success(
            request.module.clone(),
            request.operation.clone(),
            json!({
                "primary_value": result.primary_value,
                "unit": result.unit,
                "formula_used": result.formula_used,
                "secondary_values": result.secondary_values,
                "interpretation": result.interpretation
            }),
        ),
        Err(e) => ComputationResponse::error(request.module.clone(), request.operation.clone(), e),
    }
}
