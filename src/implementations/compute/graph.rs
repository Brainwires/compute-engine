//! Graph theory operations
//!
//! Provides shortest path, minimum spanning tree, topological sort, and related operations

use crate::engine::*;
use crate::specialized::graph_theory::*;

/// Compute graph theory operations
pub fn compute_graph(op: &GraphOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
    let result_json = match op {
        GraphOp::ShortestPath => {
            let req: ShortestPathRequest = serde_json::from_value(
                serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
            )
            .map_err(|e| format!("Failed to parse shortest path request: {}", e))?;

            let result = shortest_path(req).map_err(|e| format!("Shortest path error: {}", e))?;
            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
        GraphOp::MinimumSpanningTree => {
            let req: MSTRequest = serde_json::from_value(
                serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
            )
            .map_err(|e| format!("Failed to parse MST request: {}", e))?;

            let result = minimum_spanning_tree(req).map_err(|e| format!("MST error: {}", e))?;
            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
        GraphOp::TopologicalSort => {
            let req: TopologicalSortRequest = serde_json::from_value(
                serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
            )
            .map_err(|e| format!("Failed to parse topological sort request: {}", e))?;

            let result =
                topological_sort(req).map_err(|e| format!("Topological sort error: {}", e))?;
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
