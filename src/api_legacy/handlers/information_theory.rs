use super::super::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};

pub fn process_information_theory_request(request: &ComputationRequest) -> ComputationResponse {
    use crate::information_theory::*;

    let result = match request.operation.as_str() {
        "shannon_entropy" | "entropy" => {
            let req: EntropyRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match shannon_entropy(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "mutual_information" => {
            let req: MutualInfoRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match mutual_information(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "channel_capacity" => {
            let req: ChannelCapacityRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match channel_capacity(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "huffman_coding" => {
            let req: HuffmanRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match huffman_coding(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "kolmogorov_complexity" => {
            let req: KolmogorovRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match kolmogorov_complexity(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "conditional_entropy" => {
            let req: ConditionalEntropyRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match conditional_entropy(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "relative_entropy" | "kl_divergence" => {
            let req: RelativeEntropyRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match relative_entropy(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        _ => Err(format!("Unknown information_theory operation: {}", request.operation)),
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
