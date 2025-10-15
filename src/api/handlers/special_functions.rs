//! Special Functions request handler

use crate::api::types::{ComputationRequest, ComputationResponse};
use crate::special_functions::*;
use serde_json::json;

pub fn handle(request: &ComputationRequest) -> ComputationResponse {
    let result = match request.operation.as_str() {
        "bessel" => {
            let req: BesselRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match bessel_function(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "gamma" => {
            let req: GammaRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match gamma_function(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "erf" | "error_function" => {
            let req: ErrorFunctionRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match error_function(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "elliptic" => {
            let req: EllipticIntegralRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match elliptic_integral(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "polynomial" | "orthogonal_polynomial" => {
            let req: PolynomialRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match orthogonal_polynomial(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "airy" => {
            let req: AiryRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match airy_function(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        _ => Err(format!("Unknown special_functions operation: {}", request.operation)),
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
