use super::super::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};

pub fn process_optimization_request(request: &ComputationRequest) -> ComputationResponse {
    use crate::optimization::*;

    let result = match request.operation.as_str() {
        "curve_fitting" => {
            let req: CurveFitRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse {
                    success: false,
                    module: request.module.clone(),
                    operation: request.operation.clone(),
                    result: None,
                    error: Some(format!("Invalid parameters: {}", e)),
                },
            };
            match curve_fitting(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        _ => Err(format!("Unknown optimization operation: {}", request.operation)),
    };

    match result {
        Ok(value) => ComputationResponse {
            success: true,
            module: request.module.clone(),
            operation: request.operation.clone(),
            result: Some(value),
            error: None,
        },
        Err(e) => ComputationResponse {
            success: false,
            module: request.module.clone(),
            operation: request.operation.clone(),
            result: None,
            error: Some(e),
        },
    }
}
