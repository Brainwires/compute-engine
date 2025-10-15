//! Numerical Methods request handler

use crate::api::types::{ComputationRequest, ComputationResponse};
use crate::numerical_methods::*;
use serde_json::json;

pub fn handle(request: &ComputationRequest) -> ComputationResponse {
    let result = match request.operation.as_str() {
        "solve_ode" => {
            let req: ODESolverRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match solve_ode(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "find_root" => {
            let req: RootFindingRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match find_root(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "integrate" => {
            let req: IntegrationRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match integrate(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "interpolate" => {
            let req: InterpolationRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match interpolate(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "solve_linear_system" => {
            let req: LinearSystemRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match solve_linear_system(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "differentiate" => {
            let req: DifferentiationRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match differentiate(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "solve_pde" => {
            let req: PDESolverRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match solve_pde(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        _ => Err(format!("Unknown numerical_methods operation: {}", request.operation)),
    };

    match result {
        Ok(value) => ComputationResponse::success(
            request.module.clone(),
            request.operation.clone(),
            value,
        ),
        Err(e) => ComputationResponse::error(
            request.module.clone(),
            request.operation.clone(),
            e,
        ),
    }
}
