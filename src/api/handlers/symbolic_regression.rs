//! symbolic_regression operations handler

use crate::api::types::{ComputationRequest, ComputationResponse};
use crate::mathematics::symbolic_regression::*;
use serde_json::json;
use std::collections::HashMap;

pub fn handle(request: &ComputationRequest) -> ComputationResponse {
    let result = match request.operation.as_str() {
        "fit" | "discover" => {
            // Extract x and y arrays from parameters
            let x_data: Vec<f64> = match request.parameters.get("x") {
                Some(v) => match serde_json::from_value(v.clone()) {
                    Ok(data) => data,
                    Err(e) => {
                        return ComputationResponse::error(
                            request.module.clone(),
                            request.operation.clone(),
                            format!("Invalid x data: {}", e),
                        );
                    }
                },
                None => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        "Missing x parameter".to_string(),
                    );
                }
            };

            let y_data: Vec<f64> = match request.parameters.get("y") {
                Some(v) => match serde_json::from_value(v.clone()) {
                    Ok(data) => data,
                    Err(e) => {
                        return ComputationResponse::error(
                            request.module.clone(),
                            request.operation.clone(),
                            format!("Invalid y data: {}", e),
                        );
                    }
                },
                None => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        "Missing y parameter".to_string(),
                    );
                }
            };

            // Convert to InputData format (single variable)
            let data_points: Vec<Vec<f64>> = x_data.iter().map(|&x| vec![x]).collect();
            let input_data = InputData {
                variables: vec!["x".to_string()],
                data_points,
                target_values: y_data,
            };

            let max_complexity = request
                .parameters
                .get("max_complexity")
                .and_then(|v| v.as_u64())
                .unwrap_or(5) as u32;

            let domain = request
                .parameters
                .get("domain")
                .and_then(|v| v.as_str())
                .unwrap_or("general")
                .to_string();

            let units: Option<HashMap<String, String>> = request
                .parameters
                .get("units")
                .and_then(|v| serde_json::from_value(v.clone()).ok());

            match discover_equations(input_data, domain, max_complexity, units) {
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

#[cfg(test)]
#[path = "../../../tests/unit/api/handlers/symbolic_regression_handler_tests.rs"]
mod tests;
