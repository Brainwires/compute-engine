//! Linear algebra operations handler

use crate::api::types::{ComputationRequest, ComputationResponse};
use crate::mathematics::linear_algebra::*;
use serde_json::json;

pub fn handle(request: &ComputationRequest) -> ComputationResponse {
    let result = match request.operation.as_str() {
        "compute_svd" => {
            let req: MatrixInput = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect(),
            )) {
                Ok(r) => r,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid parameters: {}", e),
                    );
                }
            };
            match compute_svd(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        }
        "compute_eigendecomposition" | "eigendecomposition" => {
            let req: MatrixInput = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect(),
            )) {
                Ok(r) => r,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid parameters: {}", e),
                    );
                }
            };
            match compute_eigendecomposition(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        }
        "compute_pca" | "pca" => {
            let req: MatrixInput = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect(),
            )) {
                Ok(r) => r,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid parameters: {}", e),
                    );
                }
            };
            let n_components = request
                .parameters
                .get("n_components")
                .and_then(|v| v.as_u64())
                .map(|n| n as usize);
            match compute_pca(req, n_components) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        }
        "compute_matrix_rank" | "matrix_rank" => {
            let req: MatrixInput = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect(),
            )) {
                Ok(r) => r,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid parameters: {}", e),
                    );
                }
            };
            match compute_matrix_rank(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        }
        "compute_pseudoinverse" | "pseudoinverse" => {
            let req: MatrixInput = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect(),
            )) {
                Ok(r) => r,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid parameters: {}", e),
                    );
                }
            };
            match compute_pseudoinverse(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        }
        "matrix_operations" | "matrix_ops" => {
            let req: MatrixInput = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect(),
            )) {
                Ok(r) => r,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid parameters: {}", e),
                    );
                }
            };
            match matrix_operations(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        }
        "matrix_norm" | "norm" => {
            let req: MatrixNormRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect(),
            )) {
                Ok(r) => r,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid parameters: {}", e),
                    );
                }
            };
            match matrix_norm(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        }
        "matrix_power" | "power" => {
            let req: MatrixPowerRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect(),
            )) {
                Ok(r) => r,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid parameters: {}", e),
                    );
                }
            };
            match matrix_power(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        }
        "matrix_exp" | "exp" => {
            let req: MatrixExpRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect(),
            )) {
                Ok(r) => r,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid parameters: {}", e),
                    );
                }
            };
            match matrix_exp(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        }
        "cholesky" | "cholesky_decomposition" => {
            let req: CholeskyRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect(),
            )) {
                Ok(r) => r,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid parameters: {}", e),
                    );
                }
            };
            match cholesky_decomposition(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        }
        _ => {
            return ComputationResponse::error(
                request.module.clone(),
                request.operation.clone(),
                format!("Unknown operation: {}", request.operation),
            );
        }
    };

    match result {
        Ok(result_value) => ComputationResponse::success(
            request.module.clone(),
            request.operation.clone(),
            result_value,
        ),
        Err(error_msg) => {
            ComputationResponse::error(request.module.clone(), request.operation.clone(), error_msg)
        }
    }
}

// Test module
#[cfg(test)]
#[path = "../../../tests/unit/api/handlers/linear_algebra_handler_tests.rs"]
mod tests;
