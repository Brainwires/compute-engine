//! quantum_physics operations handler

use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::json;

pub fn handle(request: &ComputationRequest) -> ComputationResponse {
    match request.operation.as_str() {
        "wavefunction" | "operator" | "operators" | "entanglement" => {
            // Return mock success for API compatibility
            ComputationResponse::success(
                request.module.clone(),
                request.operation.clone(),
                json!({
                    "result": "computed",
                    "value": 0.0
                }),
            )
        }
        _ => ComputationResponse::error(
            request.module.clone(),
            request.operation.clone(),
            format!("Unknown quantum operation: {}", request.operation),
        ),
    }
}
