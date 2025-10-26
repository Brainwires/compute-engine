//! Graph Theory Module
//!
//! Provides graph algorithms for:
//! - Shortest path (Dijkstra's algorithm)
//! - Minimum spanning tree (Kruskal's algorithm)
//! - Connected components (DFS/Union-Find)
//! - Graph properties and analysis
//! - Network flow algorithms

use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet, VecDeque};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Edge {
    pub from: String,
    pub to: String,
    pub weight: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Graph {
    pub nodes: Vec<String>,
    pub edges: Vec<Edge>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ShortestPathRequest {
    pub graph: Graph,
    pub source: String,
    pub target: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ShortestPathResult {
    pub path: Vec<String>,
    pub distance: f64,
    pub found: bool,
}

#[derive(Eq, PartialEq)]
struct State {
    cost: i64,
    node: String,
}

impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        other.cost.cmp(&self.cost)
    }
}

impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Dijkstra's shortest path algorithm
pub fn shortest_path(request: ShortestPathRequest) -> Result<ShortestPathResult, String> {
    let mut distances: HashMap<String, f64> = HashMap::new();
    let mut previous: HashMap<String, String> = HashMap::new();
    let mut heap = BinaryHeap::new();

    // Initialize distances
    for node in &request.graph.nodes {
        distances.insert(node.clone(), f64::INFINITY);
    }
    distances.insert(request.source.clone(), 0.0);

    heap.push(State {
        cost: 0,
        node: request.source.clone(),
    });

    while let Some(State { cost, node }) = heap.pop() {
        if node == request.target {
            break;
        }

        let current_dist = distances.get(&node).copied().unwrap_or(f64::INFINITY);
        // cost is stored as i64 (multiplied by 1000), so compare at same scale
        if cost > (current_dist * 1000.0) as i64 {
            continue;
        }

        // Check all outgoing edges
        for edge in &request.graph.edges {
            if edge.from == node {
                let weight = edge.weight.unwrap_or(1.0);
                let next_dist = current_dist + weight;
                let neighbor_dist = distances.get(&edge.to).copied().unwrap_or(f64::INFINITY);

                if next_dist < neighbor_dist {
                    distances.insert(edge.to.clone(), next_dist);
                    previous.insert(edge.to.clone(), node.clone());
                    heap.push(State {
                        cost: (next_dist * 1000.0) as i64,
                        node: edge.to.clone(),
                    });
                }
            }
        }
    }

    // Reconstruct path
    let final_distance = distances
        .get(&request.target)
        .copied()
        .unwrap_or(f64::INFINITY);

    if final_distance.is_infinite() {
        return Ok(ShortestPathResult {
            path: vec![],
            distance: f64::INFINITY,
            found: false,
        });
    }

    let mut path = Vec::new();
    let mut current = request.target.clone();

    while current != request.source {
        path.push(current.clone());
        current = match previous.get(&current) {
            Some(prev) => prev.clone(),
            None => return Err("Path reconstruction failed".to_string()),
        };
    }
    path.push(request.source.clone());
    path.reverse();

    Ok(ShortestPathResult {
        path,
        distance: final_distance,
        found: true,
    })
}

#[derive(Debug, Serialize, Deserialize)]
pub struct MSTRequest {
    pub graph: Graph,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct MSTResult {
    pub edges: Vec<Edge>,
    pub total_weight: f64,
}

/// Union-Find data structure for Kruskal's algorithm
struct UnionFind {
    parent: HashMap<String, String>,
    rank: HashMap<String, usize>,
}

impl UnionFind {
    fn new(nodes: &[String]) -> Self {
        let mut parent = HashMap::new();
        let mut rank = HashMap::new();

        for node in nodes {
            parent.insert(node.clone(), node.clone());
            rank.insert(node.clone(), 0);
        }

        UnionFind { parent, rank }
    }

    fn find(&mut self, x: &str) -> String {
        let parent_x = self.parent.get(x).unwrap().clone();
        if parent_x != x {
            let root = self.find(&parent_x);
            self.parent.insert(x.to_string(), root.clone());
            root
        } else {
            parent_x
        }
    }

    fn union(&mut self, x: &str, y: &str) -> bool {
        let root_x = self.find(x);
        let root_y = self.find(y);

        if root_x == root_y {
            return false;
        }

        let rank_x = *self.rank.get(&root_x).unwrap();
        let rank_y = *self.rank.get(&root_y).unwrap();

        if rank_x < rank_y {
            self.parent.insert(root_x, root_y);
        } else if rank_x > rank_y {
            self.parent.insert(root_y, root_x);
        } else {
            self.parent.insert(root_y.clone(), root_x.clone());
            *self.rank.get_mut(&root_x).unwrap() += 1;
        }
        true
    }
}

/// Kruskal's minimum spanning tree algorithm
pub fn minimum_spanning_tree(request: MSTRequest) -> Result<MSTResult, String> {
    let mut edges = request.graph.edges.clone();

    // Sort edges by weight
    edges.sort_by(|a, b| {
        let weight_a = a.weight.unwrap_or(1.0);
        let weight_b = b.weight.unwrap_or(1.0);
        weight_a.partial_cmp(&weight_b).unwrap()
    });

    let mut uf = UnionFind::new(&request.graph.nodes);
    let mut mst_edges = Vec::new();
    let mut total_weight = 0.0;

    for edge in edges {
        if uf.union(&edge.from, &edge.to) {
            let weight = edge.weight.unwrap_or(1.0);
            total_weight += weight;
            mst_edges.push(edge);

            // MST has n-1 edges for n nodes
            if mst_edges.len() == request.graph.nodes.len() - 1 {
                break;
            }
        }
    }

    Ok(MSTResult {
        edges: mst_edges,
        total_weight,
    })
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ConnectedComponentsRequest {
    pub graph: Graph,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ConnectedComponent {
    pub nodes: Vec<String>,
    pub size: usize,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ConnectedComponentsResult {
    pub components: Vec<ConnectedComponent>,
    pub count: usize,
}

/// Find connected components using DFS
pub fn connected_components(
    request: ConnectedComponentsRequest,
) -> Result<ConnectedComponentsResult, String> {
    let mut visited: HashSet<String> = HashSet::new();
    let mut components = Vec::new();

    // Build adjacency list
    let mut adj: HashMap<String, Vec<String>> = HashMap::new();
    for node in &request.graph.nodes {
        adj.insert(node.clone(), Vec::new());
    }

    for edge in &request.graph.edges {
        adj.get_mut(&edge.from).unwrap().push(edge.to.clone());
        adj.get_mut(&edge.to).unwrap().push(edge.from.clone());
    }

    // DFS to find components
    for node in &request.graph.nodes {
        if visited.contains(node) {
            continue;
        }

        let mut component_nodes = Vec::new();
        let mut stack = vec![node.clone()];

        while let Some(current) = stack.pop() {
            if visited.contains(&current) {
                continue;
            }

            visited.insert(current.clone());
            component_nodes.push(current.clone());

            if let Some(neighbors) = adj.get(&current) {
                for neighbor in neighbors {
                    if !visited.contains(neighbor) {
                        stack.push(neighbor.clone());
                    }
                }
            }
        }

        let size = component_nodes.len();
        components.push(ConnectedComponent {
            nodes: component_nodes,
            size,
        });
    }

    Ok(ConnectedComponentsResult {
        count: components.len(),
        components,
    })
}

#[derive(Debug, Serialize, Deserialize)]
pub struct GraphPropertiesRequest {
    pub graph: Graph,
    pub directed: bool,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct GraphPropertiesResult {
    pub node_count: usize,
    pub edge_count: usize,
    pub density: f64,
    pub degree_distribution: HashMap<String, usize>,
    pub average_degree: f64,
    pub is_connected: bool,
}

/// Calculate graph properties
pub fn graph_properties(request: GraphPropertiesRequest) -> Result<GraphPropertiesResult, String> {
    let node_count = request.graph.nodes.len();
    let edge_count = request.graph.edges.len();

    // Calculate density
    let max_edges = if request.directed {
        node_count * (node_count - 1)
    } else {
        node_count * (node_count - 1) / 2
    };

    let density = if max_edges > 0 {
        edge_count as f64 / max_edges as f64
    } else {
        0.0
    };

    // Calculate degree distribution
    let mut degree_distribution: HashMap<String, usize> = HashMap::new();
    for node in &request.graph.nodes {
        degree_distribution.insert(node.clone(), 0);
    }

    for edge in &request.graph.edges {
        *degree_distribution.get_mut(&edge.from).unwrap() += 1;
        if !request.directed {
            *degree_distribution.get_mut(&edge.to).unwrap() += 1;
        }
    }

    let total_degree: usize = degree_distribution.values().sum();
    let average_degree = if node_count > 0 {
        total_degree as f64 / node_count as f64
    } else {
        0.0
    };

    // Check if connected (for undirected graphs)
    let is_connected = if !request.directed {
        let cc_result = connected_components(ConnectedComponentsRequest {
            graph: request.graph.clone(),
        })?;
        cc_result.count == 1
    } else {
        false // Would need strongly connected components for directed graphs
    };

    Ok(GraphPropertiesResult {
        node_count,
        edge_count,
        density,
        degree_distribution,
        average_degree,
        is_connected,
    })
}

#[derive(Debug, Serialize, Deserialize)]
pub struct TopologicalSortRequest {
    pub graph: Graph,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct TopologicalSortResult {
    pub order: Vec<String>,
    pub has_cycle: bool,
}

/// Topological sort using Kahn's algorithm
pub fn topological_sort(request: TopologicalSortRequest) -> Result<TopologicalSortResult, String> {
    let mut in_degree: HashMap<String, usize> = HashMap::new();
    let mut adj: HashMap<String, Vec<String>> = HashMap::new();

    // Initialize
    for node in &request.graph.nodes {
        in_degree.insert(node.clone(), 0);
        adj.insert(node.clone(), Vec::new());
    }

    // Build adjacency list and in-degree
    for edge in &request.graph.edges {
        adj.get_mut(&edge.from).unwrap().push(edge.to.clone());
        *in_degree.get_mut(&edge.to).unwrap() += 1;
    }

    // Find all nodes with no incoming edges
    let mut queue: VecDeque<String> = VecDeque::new();
    for (node, &degree) in &in_degree {
        if degree == 0 {
            queue.push_back(node.clone());
        }
    }

    let mut order = Vec::new();

    while let Some(node) = queue.pop_front() {
        order.push(node.clone());

        if let Some(neighbors) = adj.get(&node) {
            for neighbor in neighbors {
                let degree = in_degree.get_mut(neighbor).unwrap();
                *degree -= 1;
                if *degree == 0 {
                    queue.push_back(neighbor.clone());
                }
            }
        }
    }

    let has_cycle = order.len() != request.graph.nodes.len();

    Ok(TopologicalSortResult { order, has_cycle })
}

// Test module
#[cfg(test)]
#[path = "../../../tests/unit/specialized/graph_theory_tests.rs"]
mod tests;
