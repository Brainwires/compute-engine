//! Unit tests for compute::graph (graph theory)

use crate::compute::graph::*;

#[test]
fn test_shortest_path_simple() {
    let request = ShortestPathRequest {
        graph: Graph {
            nodes: vec!["0".to_string(), "1".to_string(), "2".to_string()],
            edges: vec![
                Edge { from: "0".to_string(), to: "1".to_string(), weight: Some(1.0) },
                Edge { from: "1".to_string(), to: "2".to_string(), weight: Some(1.0) },
            ],
        },
        source: "0".to_string(),
        target: "2".to_string(),
    };

    let result = shortest_path(request).unwrap();
    assert!(result.found);
    assert_eq!(result.distance, 2.0);
    assert_eq!(result.path, vec!["0", "1", "2"]);
}

#[test]
fn test_shortest_path_same_node() {
    let request = ShortestPathRequest {
        graph: Graph {
            nodes: vec!["0".to_string()],
            edges: vec![],
        },
        source: "0".to_string(),
        target: "0".to_string(),
    };

    let result = shortest_path(request).unwrap();
    assert!(result.found);
    assert_eq!(result.distance, 0.0);
    assert_eq!(result.path, vec!["0"]);
}

#[test]
fn test_shortest_path_unreachable() {
    let request = ShortestPathRequest {
        graph: Graph {
            nodes: vec!["0".to_string(), "1".to_string()],
            edges: vec![],
        },
        source: "0".to_string(),
        target: "1".to_string(),
    };

    let result = shortest_path(request).unwrap();
    assert!(!result.found);
    assert!(result.distance.is_infinite());
}

#[test]
fn test_mst_basic() {
    let request = MSTRequest {
        graph: Graph {
            nodes: vec!["0".to_string(), "1".to_string(), "2".to_string()],
            edges: vec![
                Edge { from: "0".to_string(), to: "1".to_string(), weight: Some(2.0) },
                Edge { from: "1".to_string(), to: "2".to_string(), weight: Some(3.0) },
                Edge { from: "0".to_string(), to: "2".to_string(), weight: Some(4.0) },
            ],
        },
    };

    let result = minimum_spanning_tree(request).unwrap();
    assert_eq!(result.total_weight, 5.0); // 2 + 3
    assert_eq!(result.edges.len(), 2);
}

#[test]
fn test_mst_single_node() {
    let request = MSTRequest {
        graph: Graph {
            nodes: vec!["0".to_string()],
            edges: vec![],
        },
    };

    let result = minimum_spanning_tree(request).unwrap();
    assert_eq!(result.total_weight, 0.0);
    assert_eq!(result.edges.len(), 0);
}

#[test]
fn test_connected_components_single() {
    let request = ConnectedComponentsRequest {
        graph: Graph {
            nodes: vec!["0".to_string(), "1".to_string(), "2".to_string()],
            edges: vec![
                Edge { from: "0".to_string(), to: "1".to_string(), weight: None },
                Edge { from: "1".to_string(), to: "2".to_string(), weight: None },
            ],
        },
    };

    let result = connected_components(request).unwrap();
    assert_eq!(result.count, 1);
    assert_eq!(result.components[0].size, 3);
}

#[test]
fn test_connected_components_multiple() {
    let request = ConnectedComponentsRequest {
        graph: Graph {
            nodes: vec!["0".to_string(), "1".to_string(), "2".to_string(), "3".to_string()],
            edges: vec![
                Edge { from: "0".to_string(), to: "1".to_string(), weight: None },
                Edge { from: "2".to_string(), to: "3".to_string(), weight: None },
            ],
        },
    };

    let result = connected_components(request).unwrap();
    assert_eq!(result.count, 2);
}

#[test]
fn test_graph_properties_basic() {
    let request = GraphPropertiesRequest {
        graph: Graph {
            nodes: vec!["0".to_string(), "1".to_string(), "2".to_string()],
            edges: vec![
                Edge { from: "0".to_string(), to: "1".to_string(), weight: None },
                Edge { from: "1".to_string(), to: "2".to_string(), weight: None },
            ],
        },
        directed: false,
    };

    let result = graph_properties(request).unwrap();
    assert_eq!(result.node_count, 3);
    assert_eq!(result.edge_count, 2);
    assert!(result.is_connected);
}

#[test]
fn test_topological_sort_basic() {
    let request = TopologicalSortRequest {
        graph: Graph {
            nodes: vec!["0".to_string(), "1".to_string(), "2".to_string()],
            edges: vec![
                Edge { from: "0".to_string(), to: "1".to_string(), weight: None },
                Edge { from: "1".to_string(), to: "2".to_string(), weight: None },
            ],
        },
    };

    let result = topological_sort(request).unwrap();
    assert_eq!(result.order, vec!["0", "1", "2"]);
    assert!(!result.has_cycle);
}

#[test]
fn test_topological_sort_empty() {
    let request = TopologicalSortRequest {
        graph: Graph {
            nodes: vec!["0".to_string(), "1".to_string(), "2".to_string()],
            edges: vec![],
        },
    };

    let result = topological_sort(request).unwrap();
    assert_eq!(result.order.len(), 3);
    assert!(!result.has_cycle);
}

#[test]
fn test_topological_sort_cycle() {
    let request = TopologicalSortRequest {
        graph: Graph {
            nodes: vec!["0".to_string(), "1".to_string(), "2".to_string()],
            edges: vec![
                Edge { from: "0".to_string(), to: "1".to_string(), weight: None },
                Edge { from: "1".to_string(), to: "2".to_string(), weight: None },
                Edge { from: "2".to_string(), to: "0".to_string(), weight: None },
            ],
        },
    };

    let result = topological_sort(request).unwrap();
    assert!(result.has_cycle);
}
