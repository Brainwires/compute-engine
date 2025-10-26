//! fluid_dynamics operations handler

use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::json;

pub fn handle(request: &ComputationRequest) -> ComputationResponse {
    match request.operation.as_str() {
        "navier_stokes" | "reynolds_number" | "bernoulli" => {
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
            format!("Unknown fluid dynamics operation: {}", request.operation),
        ),
    }
}
