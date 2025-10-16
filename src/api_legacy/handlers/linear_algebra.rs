use super::super::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};

pub fn process_linear_algebra_request(request: &ComputationRequest) -> ComputationResponse {
    use crate::linear_algebra::*;

    let result = match request.operation.as_str() {
        "compute_svd" => {
            let req: MatrixInput = match serde_json::from_value(serde_json::Value::Object(
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
            match compute_svd(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "compute_eigendecomposition" => {
            let req: MatrixInput = match serde_json::from_value(serde_json::Value::Object(
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
            match compute_eigendecomposition(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "compute_pca" => {
            let req: MatrixInput = match serde_json::from_value(serde_json::Value::Object(
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
            let n_components = request.parameters.get("n_components")
                .and_then(|v| v.as_u64())
                .map(|n| n as usize);
            match compute_pca(req, n_components) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "compute_matrix_rank" => {
            let req: MatrixInput = match serde_json::from_value(serde_json::Value::Object(
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
            match compute_matrix_rank(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "compute_pseudoinverse" => {
            let req: MatrixInput = match serde_json::from_value(serde_json::Value::Object(
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
            match compute_pseudoinverse(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "matrix_operations" => {
            let req: MatrixInput = match serde_json::from_value(serde_json::Value::Object(
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
            match matrix_operations(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        _ => Err(format!("Unknown linear_algebra operation: {}", request.operation)),
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
