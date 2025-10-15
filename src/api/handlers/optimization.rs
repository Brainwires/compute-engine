//! Optimization operations handler

use crate::api::types::{ComputationRequest, ComputationResponse};

pub fn handle(request: &ComputationRequest) -> ComputationResponse {
    ComputationResponse::error(
        request.module.clone(),
        request.operation.clone(),
        format!("Handler needs parameter mapping for operation: {}", request.operation)
    )
}
