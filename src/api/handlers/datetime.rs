use crate::api::types::ComputationResponse;
use crate::datetime::{DateTimeInput, calculate_datetime};
/**
 * Date/Time API Handler
 */
use serde_json::{Value, json};

pub fn handle(request: &crate::api::types::ComputationRequest) -> ComputationResponse {
    // Convert HashMap to serde_json::Map
    let params_map: serde_json::Map<String, Value> =
        request.parameters.clone().into_iter().collect();

    // Parse the input from the request parameters
    let input: DateTimeInput = match serde_json::from_value(serde_json::Value::Object(params_map)) {
        Ok(input) => input,
        Err(e) => {
            return ComputationResponse::error(
                request.module.clone(),
                request.operation.clone(),
                format!("Invalid datetime request: {}", e),
            );
        }
    };

    // Execute calculation
    match calculate_datetime(input) {
        Ok(result) => ComputationResponse::success(
            request.module.clone(),
            request.operation.clone(),
            json!({
                "value": result.value,
                "numeric_value": result.numeric_value,
                "unit": result.unit,
                "operation_used": result.operation_used,
                "interpretation": result.interpretation,
                "additional_info": result.additional_info
            }),
        ),
        Err(e) => ComputationResponse::error(request.module.clone(), request.operation.clone(), e),
    }
}

#[cfg(test)]
#[path = "../../../tests/unit/api/handlers/datetime_handler_tests.rs"]
mod tests;
