/**
 * Thermodynamics API Handler
 */

use serde_json::{json, Value};
use crate::thermodynamics::{ThermodynamicsInput, calculate_thermodynamics};
use crate::api::types::ComputationResponse;

pub fn handle(request: &crate::api::types::ComputationRequest) -> ComputationResponse {
    // Convert HashMap to serde_json::Map
    let params_map: serde_json::Map<String, Value> = request.parameters.clone().into_iter().collect();

    // Parse the input from the request parameters
    let input: ThermodynamicsInput = match serde_json::from_value(serde_json::Value::Object(params_map)) {
        Ok(input) => input,
        Err(e) => {
            return ComputationResponse::error(
                request.module.clone(),
                request.operation.clone(),
                format!("Invalid thermodynamics request: {}", e)
            );
        }
    };

    // Execute calculation
    match calculate_thermodynamics(input) {
        Ok(result) => ComputationResponse::success(
            request.module.clone(),
            request.operation.clone(),
            json!({
                "value": result.value,
                "unit": result.unit,
                "formula_used": result.formula_used,
                "interpretation": result.interpretation,
                "additional_info": result.additional_info
            })
        ),
        Err(e) => ComputationResponse::error(request.module.clone(), request.operation.clone(), e),
    }
}
