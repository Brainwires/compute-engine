//! Unit tests for graph_theory API handler
//!
//! Tests all handler operations including:
//! - Shortest path (Dijkstra's algorithm)
//! - Minimum spanning tree (Kruskal's algorithm)
//! - Connected components
//! - Graph properties (density, degree distribution, connectivity)
//! - Topological sort (Kahn's algorithm)

use super::*; // Import from parent module (graph_theory handler)
use crate::api::types::ComputationRequest;
use serde_json::json;
use std::collections::HashMap;

fn create_request(operation: &str, parameters: serde_json::Value) -> ComputationRequest {
    let params_map: HashMap<String, serde_json::Value> = serde_json::from_value(parameters)
        .expect("Failed to convert parameters to HashMap");

    ComputationRequest {
        module: "graph_theory".to_string(),
        operation: operation.to_string(),
        parameters: params_map,
    }
}

// ============================================================================
// Shortest Path Tests
// ============================================================================

#[test]
fn test_shortest_path_simple_graph() {
    let params = json!({
        "graph": {
            "nodes": ["A", "B", "C", "D"],
            "edges": [
                {"from": "A", "to": "B", "weight": 1.0},
                {"from": "B", "to": "C", "weight": 2.0},
                {"from": "A", "to": "D", "weight": 5.0},
                {"from": "D", "to": "C", "weight": 1.0}
            ]
        },
        "source": "A",
        "target": "C"
    });
    let request = create_request("shortest_path", params);
    let response = handle(&request);

    assert!(response.success, "Shortest path should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("path").is_some(), "Should have path");
    assert!(result.get("distance").is_some(), "Should have distance");
    assert!(result.get("found").is_some(), "Should have found flag");

    let distance = result.get("distance").unwrap().as_f64().unwrap();
    assert_eq!(distance, 3.0, "Distance A->B->C should be 3.0");
}

#[test]
fn test_shortest_path_no_path_exists() {
    let params = json!({
        "graph": {
            "nodes": ["A", "B", "C", "D"],
            "edges": [
                {"from": "A", "to": "B", "weight": 1.0},
                {"from": "C", "to": "D", "weight": 1.0}
            ]
        },
        "source": "A",
        "target": "D"
    });
    let request = create_request("shortest_path", params);
    let response = handle(&request);

    assert!(response.success, "Should succeed even when no path exists");
    let result = response.result.unwrap();

    let found = result.get("found").unwrap().as_bool().unwrap();
    assert!(!found, "Should indicate path not found");
}

#[test]
fn test_shortest_path_with_unweighted_edges() {
    let params = json!({
        "graph": {
            "nodes": ["A", "B", "C"],
            "edges": [
                {"from": "A", "to": "B"},
                {"from": "B", "to": "C"}
            ]
        },
        "source": "A",
        "target": "C"
    });
    let request = create_request("shortest_path", params);
    let response = handle(&request);

    assert!(response.success, "Should handle unweighted edges (default weight 1.0)");
    let result = response.result.unwrap();

    let distance = result.get("distance").unwrap().as_f64().unwrap();
    assert_eq!(distance, 2.0, "Unweighted edges should default to 1.0");
}

#[test]
fn test_shortest_path_invalid_params() {
    let params = json!({
        "invalid_field": "test"
    });
    let request = create_request("shortest_path", params);
    let response = handle(&request);

    assert!(!response.success, "Should fail with invalid parameters");
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Invalid parameters"));
}

// ============================================================================
// Minimum Spanning Tree Tests
// ============================================================================

#[test]
fn test_minimum_spanning_tree_simple() {
    let params = json!({
        "graph": {
            "nodes": ["A", "B", "C", "D"],
            "edges": [
                {"from": "A", "to": "B", "weight": 1.0},
                {"from": "B", "to": "C", "weight": 2.0},
                {"from": "C", "to": "D", "weight": 1.0},
                {"from": "A", "to": "D", "weight": 5.0},
                {"from": "B", "to": "D", "weight": 3.0}
            ]
        }
    });
    let request = create_request("minimum_spanning_tree", params);
    let response = handle(&request);

    assert!(response.success, "MST computation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("edges").is_some(), "Should have edges");
    assert!(result.get("total_weight").is_some(), "Should have total weight");

    let edges = result.get("edges").unwrap().as_array().unwrap();
    assert_eq!(edges.len(), 3, "MST of 4 nodes should have 3 edges");

    let total_weight = result.get("total_weight").unwrap().as_f64().unwrap();
    assert_eq!(total_weight, 4.0, "MST weight should be 1+2+1=4");
}

#[test]
fn test_minimum_spanning_tree_single_node() {
    let params = json!({
        "graph": {
            "nodes": ["A"],
            "edges": []
        }
    });
    let request = create_request("minimum_spanning_tree", params);
    let response = handle(&request);

    assert!(response.success, "Should handle single node graph");
    let result = response.result.unwrap();

    let edges = result.get("edges").unwrap().as_array().unwrap();
    assert_eq!(edges.len(), 0, "Single node MST has no edges");
}

#[test]
fn test_minimum_spanning_tree_disconnected_graph() {
    let params = json!({
        "graph": {
            "nodes": ["A", "B", "C", "D"],
            "edges": [
                {"from": "A", "to": "B", "weight": 1.0},
                {"from": "C", "to": "D", "weight": 2.0}
            ]
        }
    });
    let request = create_request("minimum_spanning_tree", params);
    let response = handle(&request);

    assert!(response.success, "Should handle disconnected graph");
    let result = response.result.unwrap();

    let edges = result.get("edges").unwrap().as_array().unwrap();
    assert_eq!(edges.len(), 2, "Disconnected graph MST includes edges from all components");
}

#[test]
fn test_minimum_spanning_tree_invalid_params() {
    let params = json!({
        "wrong_param": "value"
    });
    let request = create_request("minimum_spanning_tree", params);
    let response = handle(&request);

    assert!(!response.success, "Should fail with invalid parameters");
    assert!(response.error.is_some());
}

// ============================================================================
// Connected Components Tests
// ============================================================================

#[test]
fn test_connected_components_single_component() {
    let params = json!({
        "graph": {
            "nodes": ["A", "B", "C", "D"],
            "edges": [
                {"from": "A", "to": "B"},
                {"from": "B", "to": "C"},
                {"from": "C", "to": "D"}
            ]
        }
    });
    let request = create_request("connected_components", params);
    let response = handle(&request);

    assert!(response.success, "Connected components should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("components").is_some());
    assert!(result.get("count").is_some());

    let count = result.get("count").unwrap().as_u64().unwrap();
    assert_eq!(count, 1, "Fully connected graph has 1 component");
}

#[test]
fn test_connected_components_multiple_components() {
    let params = json!({
        "graph": {
            "nodes": ["A", "B", "C", "D", "E", "F"],
            "edges": [
                {"from": "A", "to": "B"},
                {"from": "C", "to": "D"},
                {"from": "E", "to": "F"}
            ]
        }
    });
    let request = create_request("connected_components", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();

    let count = result.get("count").unwrap().as_u64().unwrap();
    assert_eq!(count, 3, "Should have 3 separate components");

    let components = result.get("components").unwrap().as_array().unwrap();
    assert_eq!(components.len(), 3);
}

#[test]
fn test_connected_components_isolated_nodes() {
    let params = json!({
        "graph": {
            "nodes": ["A", "B", "C"],
            "edges": []
        }
    });
    let request = create_request("connected_components", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();

    let count = result.get("count").unwrap().as_u64().unwrap();
    assert_eq!(count, 3, "Each isolated node is its own component");
}

#[test]
fn test_connected_components_invalid_params() {
    let params = json!({
        "not_a_graph": "invalid"
    });
    let request = create_request("connected_components", params);
    let response = handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
}

// ============================================================================
// Graph Properties Tests
// ============================================================================

#[test]
fn test_graph_properties_undirected() {
    let params = json!({
        "graph": {
            "nodes": ["A", "B", "C"],
            "edges": [
                {"from": "A", "to": "B"},
                {"from": "B", "to": "C"}
            ]
        },
        "directed": false
    });
    let request = create_request("graph_properties", params);
    let response = handle(&request);

    assert!(response.success, "Graph properties should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("node_count").is_some());
    assert!(result.get("edge_count").is_some());
    assert!(result.get("density").is_some());
    assert!(result.get("degree_distribution").is_some());
    assert!(result.get("average_degree").is_some());
    assert!(result.get("is_connected").is_some());

    let node_count = result.get("node_count").unwrap().as_u64().unwrap();
    assert_eq!(node_count, 3);

    let edge_count = result.get("edge_count").unwrap().as_u64().unwrap();
    assert_eq!(edge_count, 2);

    let is_connected = result.get("is_connected").unwrap().as_bool().unwrap();
    assert!(is_connected, "Graph should be connected");
}

#[test]
fn test_graph_properties_directed() {
    let params = json!({
        "graph": {
            "nodes": ["A", "B", "C", "D"],
            "edges": [
                {"from": "A", "to": "B"},
                {"from": "B", "to": "C"},
                {"from": "C", "to": "D"}
            ]
        },
        "directed": true
    });
    let request = create_request("graph_properties", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();

    let density = result.get("density").unwrap().as_f64().unwrap();
    // For directed graph: 3 edges / (4 * 3) = 0.25
    assert!((density - 0.25).abs() < 1e-10, "Directed graph density should be 0.25");
}

#[test]
fn test_graph_properties_disconnected() {
    let params = json!({
        "graph": {
            "nodes": ["A", "B", "C", "D"],
            "edges": [
                {"from": "A", "to": "B"},
                {"from": "C", "to": "D"}
            ]
        },
        "directed": false
    });
    let request = create_request("graph_properties", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();

    let is_connected = result.get("is_connected").unwrap().as_bool().unwrap();
    assert!(!is_connected, "Graph should not be connected");
}

#[test]
fn test_graph_properties_invalid_params() {
    let params = json!({
        "graph": {"nodes": [], "edges": []}
    });
    let request = create_request("graph_properties", params);
    let response = handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
}

// ============================================================================
// Topological Sort Tests
// ============================================================================

#[test]
fn test_topological_sort_dag() {
    let params = json!({
        "graph": {
            "nodes": ["A", "B", "C", "D"],
            "edges": [
                {"from": "A", "to": "B"},
                {"from": "A", "to": "C"},
                {"from": "B", "to": "D"},
                {"from": "C", "to": "D"}
            ]
        }
    });
    let request = create_request("topological_sort", params);
    let response = handle(&request);

    assert!(response.success, "Topological sort should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("order").is_some());
    assert!(result.get("has_cycle").is_some());

    let has_cycle = result.get("has_cycle").unwrap().as_bool().unwrap();
    assert!(!has_cycle, "DAG should not have cycle");

    let order = result.get("order").unwrap().as_array().unwrap();
    assert_eq!(order.len(), 4, "All nodes should be in topological order");

    // Verify A comes before B and C, and both B and C come before D
    let a_pos = order.iter().position(|v| v.as_str().unwrap() == "A").unwrap();
    let b_pos = order.iter().position(|v| v.as_str().unwrap() == "B").unwrap();
    let c_pos = order.iter().position(|v| v.as_str().unwrap() == "C").unwrap();
    let d_pos = order.iter().position(|v| v.as_str().unwrap() == "D").unwrap();

    assert!(a_pos < b_pos, "A should come before B");
    assert!(a_pos < c_pos, "A should come before C");
    assert!(b_pos < d_pos, "B should come before D");
    assert!(c_pos < d_pos, "C should come before D");
}

#[test]
fn test_topological_sort_with_cycle() {
    let params = json!({
        "graph": {
            "nodes": ["A", "B", "C"],
            "edges": [
                {"from": "A", "to": "B"},
                {"from": "B", "to": "C"},
                {"from": "C", "to": "A"}
            ]
        }
    });
    let request = create_request("topological_sort", params);
    let response = handle(&request);

    assert!(response.success, "Should succeed but detect cycle");
    let result = response.result.unwrap();

    let has_cycle = result.get("has_cycle").unwrap().as_bool().unwrap();
    assert!(has_cycle, "Should detect cycle");

    let order = result.get("order").unwrap().as_array().unwrap();
    assert!(order.len() < 3, "Should not include all nodes when cycle exists");
}

#[test]
fn test_topological_sort_linear_chain() {
    let params = json!({
        "graph": {
            "nodes": ["1", "2", "3", "4", "5"],
            "edges": [
                {"from": "1", "to": "2"},
                {"from": "2", "to": "3"},
                {"from": "3", "to": "4"},
                {"from": "4", "to": "5"}
            ]
        }
    });
    let request = create_request("topological_sort", params);
    let response = handle(&request);

    assert!(response.success);
    let result = response.result.unwrap();

    let has_cycle = result.get("has_cycle").unwrap().as_bool().unwrap();
    assert!(!has_cycle);

    let order = result.get("order").unwrap().as_array().unwrap();
    assert_eq!(order.len(), 5);
    // Should be in order: 1, 2, 3, 4, 5
    assert_eq!(order[0].as_str().unwrap(), "1");
    assert_eq!(order[4].as_str().unwrap(), "5");
}

#[test]
fn test_topological_sort_invalid_params() {
    let params = json!({
        "invalid": "data"
    });
    let request = create_request("topological_sort", params);
    let response = handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
}

// ============================================================================
// Unknown Operation Test
// ============================================================================

#[test]
fn test_unknown_operation() {
    let params = json!({
        "graph": {
            "nodes": ["A", "B"],
            "edges": []
        }
    });
    let request = create_request("unknown_operation", params);
    let response = handle(&request);

    assert!(!response.success, "Unknown operation should fail");
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Unknown graph_theory operation"));
}
