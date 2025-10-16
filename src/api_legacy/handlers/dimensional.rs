use super::super::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};

pub fn process_dimensional_request(request: &ComputationRequest) -> ComputationResponse {
    ComputationResponse {
        success: false,
        module: request.module.clone(),
        operation: request.operation.clone(),
        result: None,
        error: Some("Dimensional analyzer operations not yet wired to API".to_string()),
    }
}
