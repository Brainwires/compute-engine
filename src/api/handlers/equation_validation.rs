//! Equation validation operations handler

use crate::api::types::{ComputationRequest, ComputationResponse};
use crate::tools::equation_validation::*;
use serde_json::json;
use std::collections::HashMap;

pub fn handle(request: &ComputationRequest) -> ComputationResponse {
    let result = match request.operation.as_str() {
        "validate" => {
            let equation = match request.parameters.get("equation") {
                Some(v) => match v.as_str() {
                    Some(s) => s.to_string(),
                    None => {
                        return ComputationResponse::error(
                            request.module.clone(),
                            request.operation.clone(),
                            "Invalid equation parameter".to_string(),
                        );
                    }
                },
                None => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        "Missing equation parameter".to_string(),
                    );
                }
            };

            let domain = match request.parameters.get("domain") {
                Some(v) => match v.as_str() {
                    Some(s) => s.to_string(),
                    None => "general".to_string(),
                },
                None => "general".to_string(),
            };

            let units: Option<HashMap<String, String>> = request
                .parameters
                .get("units")
                .and_then(|v| serde_json::from_value(v.clone()).ok());

            let conservation_laws: Option<Vec<String>> = request
                .parameters
                .get("conservation_laws")
                .and_then(|v| serde_json::from_value(v.clone()).ok());

            let symmetries: Option<Vec<String>> = request
                .parameters
                .get("symmetries")
                .and_then(|v| serde_json::from_value(v.clone()).ok());

            match validate_equation(equation, domain, units, conservation_laws, symmetries) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e.to_string()),
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
