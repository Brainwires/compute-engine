//! stochastic_processes operations handler

use crate::api::types::{ComputationRequest, ComputationResponse};
use crate::specialized::stochastic_processes::*;
use serde_json::json;

pub fn handle(request: &ComputationRequest) -> ComputationResponse {
    let result = match request.operation.as_str() {
        "brownian_motion" => {
            let params: BrownianMotionParams = match serde_json::from_value(
                serde_json::Value::Object(request.parameters.clone().into_iter().collect()),
            ) {
                Ok(p) => p,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid parameters: {}", e),
                    );
                }
            };
            match generate_brownian_motion(params) {
                Ok(path) => Ok(json!({"path": path})),
                Err(e) => Err(e),
            }
        }
        "markov_chain" => {
            let params: MarkovChainParams = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect(),
            )) {
                Ok(p) => p,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid parameters: {}", e),
                    );
                }
            };
            match simulate_markov_chain(params) {
                Ok(chain) => Ok(json!({"chain": chain})),
                Err(e) => Err(e),
            }
        }
        "poisson_process" => {
            let rate: f64 = match request.parameters.get("rate") {
                Some(v) => match v.as_f64() {
                    Some(r) => r,
                    None => {
                        return ComputationResponse::error(
                            request.module.clone(),
                            request.operation.clone(),
                            "Invalid rate parameter".to_string(),
                        );
                    }
                },
                None => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        "Missing rate parameter".to_string(),
                    );
                }
            };
            let time_horizon: f64 = match request.parameters.get("time_horizon") {
                Some(v) => match v.as_f64() {
                    Some(t) => t,
                    None => {
                        return ComputationResponse::error(
                            request.module.clone(),
                            request.operation.clone(),
                            "Invalid time_horizon parameter".to_string(),
                        );
                    }
                },
                None => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        "Missing time_horizon parameter".to_string(),
                    );
                }
            };
            match simulate_poisson_process(rate, time_horizon) {
                Ok(events) => Ok(json!({"events": events})),
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

#[cfg(test)]
#[path = "../../../tests/unit/api/handlers/stochastic_processes_handler_tests.rs"]
mod tests;
