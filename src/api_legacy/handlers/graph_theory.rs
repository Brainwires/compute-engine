use super::super::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};

pub fn process_graph_theory_request(request: &ComputationRequest) -> ComputationResponse {
    use crate::graph_theory::*;

    let result = match request.operation.as_str() {
        "shortest_path" => {
            let req: ShortestPathRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match shortest_path(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "minimum_spanning_tree" => {
            let req: MSTRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match minimum_spanning_tree(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "connected_components" => {
            let req: ConnectedComponentsRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match connected_components(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "graph_properties" => {
            let req: GraphPropertiesRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match graph_properties(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "topological_sort" => {
            let req: TopologicalSortRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match topological_sort(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        _ => Err(format!("Unknown graph_theory operation: {}", request.operation)),
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
