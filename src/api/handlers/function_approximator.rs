//! Function approximator handler

use crate::api::types::{ComputationRequest, ComputationResponse};
use crate::tools::function_approximator::*;
use serde_json::json;

pub fn handle(request: &ComputationRequest) -> ComputationResponse {
    let result = match request.operation.as_str() {
        "approximate" | "discover_function" => {
            // Extract inputs
            let inputs: Vec<f64> = match request.parameters.get("inputs") {
                Some(v) => match serde_json::from_value(v.clone()) {
                    Ok(i) => i,
                    Err(e) => {
                        return ComputationResponse::error(
                            request.module.clone(),
                            request.operation.clone(),
                            format!("Invalid inputs: {}", e),
                        )
                    }
                },
                None => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        "Missing inputs parameter".to_string(),
                    )
                }
            };

            // Extract outputs
            let outputs: Vec<f64> = match request.parameters.get("outputs") {
                Some(v) => match serde_json::from_value(v.clone()) {
                    Ok(o) => o,
                    Err(e) => {
                        return ComputationResponse::error(
                            request.module.clone(),
                            request.operation.clone(),
                            format!("Invalid outputs: {}", e),
                        )
                    }
                },
                None => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        "Missing outputs parameter".to_string(),
                    )
                }
            };

            // Build config
            let mut config = FunctionApproximatorConfig::default();

            if let Some(v) = request.parameters.get("population_size") {
                if let Some(size) = v.as_u64() {
                    config.population_size = size as usize;
                }
            }

            if let Some(v) = request.parameters.get("generations") {
                if let Some(gens) = v.as_u64() {
                    config.generations = gens as usize;
                }
            }

            if let Some(v) = request.parameters.get("complexity") {
                if let Some(c) = v.as_f64() {
                    config.complexity_penalty = c;
                }
            }

            if let Some(v) = request.parameters.get("death_probability") {
                if let Some(d) = v.as_f64() {
                    config.death_probability = d;
                }
            }

            if let Some(v) = request.parameters.get("diversity") {
                if let Some(d) = v.as_u64() {
                    config.diversity_factor = d as usize;
                }
            }

            // Create approximator and run
            let approximator = match FunctionApproximator::new(inputs.clone(), outputs.clone(), config) {
                Ok(a) => a,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        e,
                    )
                }
            };

            match approximator.approximate() {
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
        Err(error_msg) => ComputationResponse::error(
            request.module.clone(),
            request.operation.clone(),
            error_msg,
        ),
    }
}
