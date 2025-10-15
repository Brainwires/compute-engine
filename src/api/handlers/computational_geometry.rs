//! Computational geometry operations handler

use crate::api::types::{ComputationRequest, ComputationResponse};
use crate::tools::computational_geometry::*;
use serde_json::json;

pub fn handle(request: &ComputationRequest) -> ComputationResponse {
    let result = match request.operation.as_str() {
        "convex_hull" => {
            let req: ConvexHullRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e)
                ),
            };
            match convex_hull(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "delaunay_triangulation" | "delaunay" => {
            let req: DelaunayRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e)
                ),
            };
            match delaunay_triangulation(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "voronoi_diagram" | "voronoi" => {
            let req: VoronoiRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e)
                ),
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
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e)
                ),
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
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e)
                ),
            };
            match point_in_polygon(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
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
