//! Dimensional analysis operations handler

use crate::api::types::{ComputationRequest, ComputationResponse};
use crate::tools::dimensional_analysis::*;
use serde_json::json;
use std::collections::HashMap;

pub fn handle(request: &ComputationRequest) -> ComputationResponse {
    let result = match request.operation.as_str() {
        "check_dimensions" | "analyze" => {
            let equation = match request.parameters.get("equation") {
                Some(v) => match v.as_str() {
                    Some(s) => s.to_string(),
                    None => return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        "Invalid equation parameter".to_string()
                    ),
                },
                None => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    "Missing equation parameter".to_string()
                ),
            };

            let units: HashMap<String, String> = match request.parameters.get("units") {
                Some(v) => match serde_json::from_value(v.clone()) {
                    Ok(u) => u,
                    Err(e) => return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid units: {}", e)
                    ),
                },
                None => HashMap::new(),
            };

            let target = request.parameters.get("target_dimension")
                .and_then(|v| v.as_str())
                .map(|s| s.to_string());

            match dimensional_analysis(equation, units, target) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e.to_string()),
            }
        },
        _ => {
            return ComputationResponse::error(
                request.module.clone(),
                request.operation.clone(),
                format!("Unknown operation: {}", request.operation)
            );
        }
    };

    match result {
        Ok(result_value) => ComputationResponse::success(
            request.module.clone(),
            request.operation.clone(),
            result_value
        ),
        Err(error_msg) => ComputationResponse::error(
            request.module.clone(),
            request.operation.clone(),
            error_msg
        ),
    }
}
