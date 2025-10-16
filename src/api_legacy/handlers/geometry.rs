use super::super::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};

pub fn process_geometry_request(request: &ComputationRequest) -> ComputationResponse {
    use crate::computational_geometry::*;

    let result = match request.operation.as_str() {
        "convex_hull" => {
            let req: ConvexHullRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match convex_hull(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "delaunay_triangulation" => {
            let req: DelaunayRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match delaunay_triangulation(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "voronoi_diagram" => {
            let req: VoronoiRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match voronoi_diagram(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "polygon_area" => {
            let req: PolygonAreaRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match polygon_area(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "point_in_polygon" => {
            let req: PointInPolygonRequest = match serde_json::from_value(serde_json::Value::Object(
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
            match point_in_polygon(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        _ => Err(format!("Unknown computational_geometry operation: {}", request.operation)),
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
