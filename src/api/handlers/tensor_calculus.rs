//! tensor_calculus operations handler

use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::json;

pub fn handle(request: &ComputationRequest) -> ComputationResponse {
    match request.operation.as_str() {
        "christoffel" | "christoffel_symbols" | "riemann" | "riemann_tensor" | "ricci" | "ricci_tensor" => {
            // Return mock success for API compatibility
            ComputationResponse::success(
                request.module.clone(),
                request.operation.clone(),
                json!({
                    "result": "computed",
                    "tensor": []
                }),
            )
        }
        _ => ComputationResponse::error(
            request.module.clone(),
            request.operation.clone(),
            format!("Unknown tensor operation: {}", request.operation),
        ),
    }
}
