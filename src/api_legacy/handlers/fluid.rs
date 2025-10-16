use super::super::{ComputationRequest, ComputationResponse};

pub fn process_fluid_request(request: &ComputationRequest) -> ComputationResponse {
    ComputationResponse {
        success: false,
        module: request.module.clone(),
        operation: request.operation.clone(),
        result: None,
        error: Some("Fluid dynamics operations not yet implemented in API".to_string()),
    }
}
