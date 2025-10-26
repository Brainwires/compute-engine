//! Graph Theory request handler

use crate::api::types::{ComputationRequest, ComputationResponse};
use crate::graph_theory::*;
use serde_json::json;

pub fn handle(request: &ComputationRequest) -> ComputationResponse {
    let result = match request.operation.as_str() {
        "shortest_path" => {
            let req: ShortestPathRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match shortest_path(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        }
        "minimum_spanning_tree" => {
            let req: MSTRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match minimum_spanning_tree(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        }
        "connected_components" => {
            let req: ConnectedComponentsRequest = match serde_json::from_value(
                serde_json::Value::Object(request.parameters.clone().into_iter().collect()),
            ) {
                Ok(r) => r,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid parameters: {}", e),
                    );
                }
            };
            match connected_components(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        }
        "graph_properties" => {
            let req: GraphPropertiesRequest = match serde_json::from_value(
                serde_json::Value::Object(request.parameters.clone().into_iter().collect()),
            ) {
                Ok(r) => r,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid parameters: {}", e),
                    );
                }
            };
            match graph_properties(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        }
        "topological_sort" => {
            let req: TopologicalSortRequest = match serde_json::from_value(
                serde_json::Value::Object(request.parameters.clone().into_iter().collect()),
            ) {
                Ok(r) => r,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid parameters: {}", e),
                    );
                }
            };
            match topological_sort(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        }
        _ => Err(format!(
            "Unknown graph_theory operation: {}",
            request.operation
        )),
    };

    match result {
        Ok(value) => {
            ComputationResponse::success(request.module.clone(), request.operation.clone(), value)
        }
        Err(e) => ComputationResponse::error(request.module.clone(), request.operation.clone(), e),
    }
}

#[cfg(test)]
#[path = "../../../tests/unit/api/handlers/graph_theory_handler_tests.rs"]
mod tests;
