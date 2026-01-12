//! Computational geometry module
//!
//! Provides geometric algorithms and spatial computations:
//! - Convex hull algorithms (2D and 3D)
//! - Voronoi diagrams
//! - Triangulation
//! - Spatial indexing (KD-trees)

mod advanced;
pub mod spatial_3d;

pub use advanced::*;
pub use spatial_3d::*;

use crate::engine::*;

/// Compute geometry operations
pub fn compute_geometry(op: &GeometryOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
    let result_json = match op {
        GeometryOp::ConvexHull => {
            let req: ConvexHullRequest = serde_json::from_value(
                serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
            )
            .map_err(|e| format!("Failed to parse convex hull request: {}", e))?;

            let result = convex_hull(req).map_err(|e| format!("Convex hull error: {}", e))?;
            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
        GeometryOp::Delaunay => {
            let req: DelaunayRequest = serde_json::from_value(
                serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
            )
            .map_err(|e| format!("Failed to parse Delaunay request: {}", e))?;

            let result = delaunay_triangulation(req)
                .map_err(|e| format!("Delaunay triangulation error: {}", e))?;
            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
        GeometryOp::Voronoi => {
            let req: VoronoiRequest = serde_json::from_value(
                serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
            )
            .map_err(|e| format!("Failed to parse Voronoi request: {}", e))?;

            let result =
                voronoi_diagram(req).map_err(|e| format!("Voronoi diagram error: {}", e))?;
            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
        GeometryOp::PolygonArea => {
            let req: PolygonAreaRequest = serde_json::from_value(
                serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
            )
            .map_err(|e| format!("Failed to parse polygon area request: {}", e))?;

            let result = polygon_area(req).map_err(|e| format!("Polygon area error: {}", e))?;
            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
        GeometryOp::PointInPolygon => {
            let req: PointInPolygonRequest = serde_json::from_value(
                serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
            )
            .map_err(|e| format!("Failed to parse point in polygon request: {}", e))?;

            let result =
                point_in_polygon(req).map_err(|e| format!("Point in polygon error: {}", e))?;
            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
    };

    Ok(ComputeOutput {
        result: result_json,
        additional: None,
        metadata: None,
    })
}
