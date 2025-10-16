use computational_engine::graph_theory::*;

#[test]
fn test_shortest_path() {
    // Simple directed graph: A -> B -> D
    let edges = vec![
        Edge {
            from: "A".to_string(),
            to: "B".to_string(),
            weight: Some(1.0),
        },
        Edge {
            from: "B".to_string(),
            to: "D".to_string(),
            weight: Some(1.0),
        },
    ];

    let graph = Graph {
        nodes: vec!["A".to_string(), "B".to_string(), "D".to_string()],
        edges,
    };

    let result = shortest_path(ShortestPathRequest {
        graph,
        source: "A".to_string(),
        target: "D".to_string(),
    })
    .unwrap();

    println!(
        "Result: found={}, distance={}, path={:?}",
        result.found, result.distance, result.path
    );

    assert!(result.found, "Path should be found");
    assert_eq!(result.distance, 2.0);
    assert_eq!(result.path, vec!["A", "B", "D"]);
}

#[test]
fn test_shortest_path_weighted() {
    let edges = vec![
        Edge {
            from: "A".to_string(),
            to: "B".to_string(),
            weight: Some(7.0),
        },
        Edge {
            from: "A".to_string(),
            to: "C".to_string(),
            weight: Some(9.0),
        },
        Edge {
            from: "A".to_string(),
            to: "F".to_string(),
            weight: Some(14.0),
        },
        Edge {
            from: "B".to_string(),
            to: "C".to_string(),
            weight: Some(10.0),
        },
        Edge {
            from: "B".to_string(),
            to: "D".to_string(),
            weight: Some(15.0),
        },
        Edge {
            from: "C".to_string(),
            to: "D".to_string(),
            weight: Some(11.0),
        },
        Edge {
            from: "C".to_string(),
            to: "F".to_string(),
            weight: Some(2.0),
        },
        Edge {
            from: "D".to_string(),
            to: "E".to_string(),
            weight: Some(6.0),
        },
        Edge {
            from: "E".to_string(),
            to: "F".to_string(),
            weight: Some(9.0),
        },
    ];

    let graph = Graph {
        nodes: vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
            "E".to_string(),
            "F".to_string(),
        ],
        edges,
    };

    let result = shortest_path(ShortestPathRequest {
        graph,
        source: "A".to_string(),
        target: "E".to_string(),
    })
    .unwrap();

    assert!(result.found);
    assert!(result.distance < 30.0); // Should find efficient path
}

#[test]
fn test_shortest_path_no_path() {
    let edges = vec![
        Edge {
            from: "A".to_string(),
            to: "B".to_string(),
            weight: Some(1.0),
        },
        Edge {
            from: "C".to_string(),
            to: "D".to_string(),
            weight: Some(1.0),
        },
    ];

    let graph = Graph {
        nodes: vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ],
        edges,
    };

    let result = shortest_path(ShortestPathRequest {
        graph,
        source: "A".to_string(),
        target: "D".to_string(),
    })
    .unwrap();

    assert!(!result.found);
    assert!(result.distance.is_infinite());
}

#[test]
fn test_shortest_path_self_loop() {
    let edges = vec![
        Edge {
            from: "A".to_string(),
            to: "A".to_string(),
            weight: Some(5.0),
        },
        Edge {
            from: "A".to_string(),
            to: "B".to_string(),
            weight: Some(3.0),
        },
    ];

    let graph = Graph {
        nodes: vec!["A".to_string(), "B".to_string()],
        edges,
    };

    let result = shortest_path(ShortestPathRequest {
        graph,
        source: "A".to_string(),
        target: "B".to_string(),
    })
    .unwrap();

    assert!(result.found);
    assert_eq!(result.distance, 3.0);
}

#[test]
fn test_minimum_spanning_tree() {
    let edges = vec![
        Edge {
            from: "A".to_string(),
            to: "B".to_string(),
            weight: Some(4.0),
        },
        Edge {
            from: "A".to_string(),
            to: "C".to_string(),
            weight: Some(2.0),
        },
        Edge {
            from: "B".to_string(),
            to: "C".to_string(),
            weight: Some(1.0),
        },
        Edge {
            from: "B".to_string(),
            to: "D".to_string(),
            weight: Some(5.0),
        },
        Edge {
            from: "C".to_string(),
            to: "D".to_string(),
            weight: Some(8.0),
        },
    ];

    let graph = Graph {
        nodes: vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ],
        edges,
    };

    let result = minimum_spanning_tree(MSTRequest { graph }).unwrap();

    // MST should have n-1 edges for n vertices
    assert_eq!(result.edges.len(), 3);
    assert!(result.total_weight > 0.0);
}

#[test]
fn test_minimum_spanning_tree_complete_graph() {
    let edges = vec![
        Edge {
            from: "A".to_string(),
            to: "B".to_string(),
            weight: Some(1.0),
        },
        Edge {
            from: "A".to_string(),
            to: "C".to_string(),
            weight: Some(3.0),
        },
        Edge {
            from: "B".to_string(),
            to: "C".to_string(),
            weight: Some(2.0),
        },
    ];

    let graph = Graph {
        nodes: vec!["A".to_string(), "B".to_string(), "C".to_string()],
        edges,
    };

    let result = minimum_spanning_tree(MSTRequest { graph }).unwrap();

    assert_eq!(result.edges.len(), 2);
    assert_eq!(result.total_weight, 3.0); // Should pick edges with weights 1 and 2
}

#[test]
fn test_minimum_spanning_tree_single_node() {
    let edges = vec![];
    let graph = Graph {
        nodes: vec!["A".to_string()],
        edges,
    };

    let result = minimum_spanning_tree(MSTRequest { graph }).unwrap();

    assert_eq!(result.edges.len(), 0);
    assert_eq!(result.total_weight, 0.0);
}

#[test]
fn test_connected_components() {
    let edges = vec![
        Edge {
            from: "A".to_string(),
            to: "B".to_string(),
            weight: None,
        },
        Edge {
            from: "C".to_string(),
            to: "D".to_string(),
            weight: None,
        },
    ];

    let graph = Graph {
        nodes: vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ],
        edges,
    };

    let result = connected_components(ConnectedComponentsRequest { graph }).unwrap();

    // Should have 2 components
    assert_eq!(result.count, 2);
}

#[test]
fn test_connected_components_single() {
    let edges = vec![
        Edge {
            from: "A".to_string(),
            to: "B".to_string(),
            weight: None,
        },
        Edge {
            from: "B".to_string(),
            to: "C".to_string(),
            weight: None,
        },
        Edge {
            from: "C".to_string(),
            to: "D".to_string(),
            weight: None,
        },
    ];

    let graph = Graph {
        nodes: vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ],
        edges,
    };

    let result = connected_components(ConnectedComponentsRequest { graph }).unwrap();

    assert_eq!(result.count, 1);
    assert_eq!(result.components[0].size, 4);
}

#[test]
fn test_connected_components_isolated_nodes() {
    let edges = vec![];
    let graph = Graph {
        nodes: vec!["A".to_string(), "B".to_string(), "C".to_string()],
        edges,
    };

    let result = connected_components(ConnectedComponentsRequest { graph }).unwrap();

    assert_eq!(result.count, 3);
    for component in result.components {
        assert_eq!(component.size, 1);
    }
}

#[test]
fn test_topological_sort() {
    let edges = vec![
        Edge {
            from: "A".to_string(),
            to: "B".to_string(),
            weight: None,
        },
        Edge {
            from: "A".to_string(),
            to: "C".to_string(),
            weight: None,
        },
        Edge {
            from: "B".to_string(),
            to: "D".to_string(),
            weight: None,
        },
        Edge {
            from: "C".to_string(),
            to: "D".to_string(),
            weight: None,
        },
    ];

    let graph = Graph {
        nodes: vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ],
        edges,
    };

    let result = topological_sort(TopologicalSortRequest { graph }).unwrap();

    assert!(result.order.len() >= 4);
    // A should come before B and C
    let a_pos = result.order.iter().position(|x| x == "A").unwrap();
    let b_pos = result.order.iter().position(|x| x == "B").unwrap();
    assert!(a_pos < b_pos);
}

#[test]
fn test_topological_sort_linear() {
    let edges = vec![
        Edge {
            from: "A".to_string(),
            to: "B".to_string(),
            weight: None,
        },
        Edge {
            from: "B".to_string(),
            to: "C".to_string(),
            weight: None,
        },
        Edge {
            from: "C".to_string(),
            to: "D".to_string(),
            weight: None,
        },
    ];

    let graph = Graph {
        nodes: vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ],
        edges,
    };

    let result = topological_sort(TopologicalSortRequest { graph }).unwrap();

    assert!(!result.has_cycle);
    assert_eq!(result.order, vec!["A", "B", "C", "D"]);
}

#[test]
fn test_topological_sort_cycle_detection() {
    let edges = vec![
        Edge {
            from: "A".to_string(),
            to: "B".to_string(),
            weight: None,
        },
        Edge {
            from: "B".to_string(),
            to: "C".to_string(),
            weight: None,
        },
        Edge {
            from: "C".to_string(),
            to: "A".to_string(),
            weight: None,
        },
    ];

    let graph = Graph {
        nodes: vec!["A".to_string(), "B".to_string(), "C".to_string()],
        edges,
    };

    let result = topological_sort(TopologicalSortRequest { graph }).unwrap();

    assert!(result.has_cycle);
}

#[test]
fn test_graph_properties() {
    let edges = vec![
        Edge {
            from: "A".to_string(),
            to: "B".to_string(),
            weight: None,
        },
        Edge {
            from: "B".to_string(),
            to: "C".to_string(),
            weight: None,
        },
        Edge {
            from: "C".to_string(),
            to: "A".to_string(),
            weight: None,
        },
    ];

    let graph = Graph {
        nodes: vec!["A".to_string(), "B".to_string(), "C".to_string()],
        edges,
    };

    let result = graph_properties(GraphPropertiesRequest {
        graph,
        directed: false,
    })
    .unwrap();

    assert_eq!(result.node_count, 3);
    assert_eq!(result.edge_count, 3);
    assert!(result.density > 0.0);
}

#[test]
fn test_graph_properties_directed() {
    let edges = vec![
        Edge {
            from: "A".to_string(),
            to: "B".to_string(),
            weight: None,
        },
        Edge {
            from: "B".to_string(),
            to: "C".to_string(),
            weight: None,
        },
    ];

    let graph = Graph {
        nodes: vec!["A".to_string(), "B".to_string(), "C".to_string()],
        edges,
    };

    let result = graph_properties(GraphPropertiesRequest {
        graph,
        directed: true,
    })
    .unwrap();

    assert_eq!(result.node_count, 3);
    assert_eq!(result.edge_count, 2);
    assert_eq!(result.degree_distribution.get("A").unwrap(), &1);
    assert_eq!(result.degree_distribution.get("B").unwrap(), &1);
    assert_eq!(result.degree_distribution.get("C").unwrap(), &0);
}

#[test]
fn test_graph_properties_complete_graph() {
    let edges = vec![
        Edge {
            from: "A".to_string(),
            to: "B".to_string(),
            weight: None,
        },
        Edge {
            from: "A".to_string(),
            to: "C".to_string(),
            weight: None,
        },
        Edge {
            from: "B".to_string(),
            to: "C".to_string(),
            weight: None,
        },
    ];

    let graph = Graph {
        nodes: vec!["A".to_string(), "B".to_string(), "C".to_string()],
        edges,
    };

    let result = graph_properties(GraphPropertiesRequest {
        graph,
        directed: false,
    })
    .unwrap();

    assert_eq!(result.density, 1.0); // Complete graph has maximum density
    assert!(result.is_connected);
}

#[test]
fn test_graph_properties_disconnected() {
    let edges = vec![Edge {
        from: "A".to_string(),
        to: "B".to_string(),
        weight: None,
    }];

    let graph = Graph {
        nodes: vec!["A".to_string(), "B".to_string(), "C".to_string()],
        edges,
    };

    let result = graph_properties(GraphPropertiesRequest {
        graph,
        directed: false,
    })
    .unwrap();

    assert!(!result.is_connected);
}
