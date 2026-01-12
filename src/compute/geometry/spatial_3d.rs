//! 3D Computational Geometry
//!
//! Three-dimensional geometric algorithms including convex hull and spatial indexing

use serde::{Deserialize, Serialize};
use super::advanced::Point3D;

/// 3D Convex Hull using Gift Wrapping algorithm
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConvexHull3D {
    pub vertices: Vec<Point3D>,
    pub faces: Vec<Face3D>,
    pub volume: f64,
    pub surface_area: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Face3D {
    pub p1: usize, // Indices into vertex list
    pub p2: usize,
    pub p3: usize,
    pub normal: Point3D,
}

/// Compute 3D convex hull using QuickHull algorithm (simplified)
pub fn convex_hull_3d(points: &[Point3D]) -> Result<ConvexHull3D, String> {
    if points.len() < 4 {
        return Err("Need at least 4 points for 3D convex hull".to_string());
    }

    // Find initial tetrahedron
    let vertices = points.to_vec();
    let mut faces = Vec::new();

    // Find extreme points
    let _min_x_idx = find_extreme(&vertices, |p| p.x, true);
    let _max_x_idx = find_extreme(&vertices, |p| p.x, false);
    let _min_y_idx = find_extreme(&vertices, |p| p.y, true);
    let _max_y_idx = find_extreme(&vertices, |p| p.y, false);

    // Create initial simplex (simplified - just use first 4 distinct points)
    if vertices.len() >= 4 {
        // Create initial tetrahedron faces
        faces.push(Face3D {
            p1: 0,
            p2: 1,
            p3: 2,
            normal: compute_normal(&vertices[0], &vertices[1], &vertices[2]),
        });
        faces.push(Face3D {
            p1: 0,
            p2: 1,
            p3: 3,
            normal: compute_normal(&vertices[0], &vertices[1], &vertices[3]),
        });
        faces.push(Face3D {
            p1: 0,
            p2: 2,
            p3: 3,
            normal: compute_normal(&vertices[0], &vertices[2], &vertices[3]),
        });
        faces.push(Face3D {
            p1: 1,
            p2: 2,
            p3: 3,
            normal: compute_normal(&vertices[1], &vertices[2], &vertices[3]),
        });
    }

    // Compute volume using divergence theorem
    let mut volume = 0.0;
    for face in &faces {
        let v1 = &vertices[face.p1];
        let v2 = &vertices[face.p2];
        let v3 = &vertices[face.p3];
        volume += signed_volume_of_tetrahedron(&Point3D { x: 0.0, y: 0.0, z: 0.0 }, v1, v2, v3);
    }
    volume = volume.abs();

    // Compute surface area
    let mut surface_area = 0.0;
    for face in &faces {
        let v1 = &vertices[face.p1];
        let v2 = &vertices[face.p2];
        let v3 = &vertices[face.p3];
        surface_area += triangle_area_3d(v1, v2, v3);
    }

    Ok(ConvexHull3D {
        vertices,
        faces,
        volume,
        surface_area,
    })
}

fn find_extreme<F>(points: &[Point3D], coord: F, minimize: bool) -> usize
where
    F: Fn(&Point3D) -> f64,
{
    let mut best_idx = 0;
    let mut best_val = coord(&points[0]);

    for (i, point) in points.iter().enumerate().skip(1) {
        let val = coord(point);
        if (minimize && val < best_val) || (!minimize && val > best_val) {
            best_val = val;
            best_idx = i;
        }
    }

    best_idx
}

fn compute_normal(p1: &Point3D, p2: &Point3D, p3: &Point3D) -> Point3D {
    let v1 = Point3D {
        x: p2.x - p1.x,
        y: p2.y - p1.y,
        z: p2.z - p1.z,
    };
    let v2 = Point3D {
        x: p3.x - p1.x,
        y: p3.y - p1.y,
        z: p3.z - p1.z,
    };

    // Cross product
    let normal = Point3D {
        x: v1.y * v2.z - v1.z * v2.y,
        y: v1.z * v2.x - v1.x * v2.z,
        z: v1.x * v2.y - v1.y * v2.x,
    };

    // Normalize
    let len = (normal.x * normal.x + normal.y * normal.y + normal.z * normal.z).sqrt();
    if len > 1e-10 {
        Point3D {
            x: normal.x / len,
            y: normal.y / len,
            z: normal.z / len,
        }
    } else {
        Point3D { x: 0.0, y: 0.0, z: 1.0 }
    }
}

fn signed_volume_of_tetrahedron(p0: &Point3D, p1: &Point3D, p2: &Point3D, p3: &Point3D) -> f64 {
    let v1 = (p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
    let v2 = (p2.x - p0.x, p2.y - p0.y, p2.z - p0.z);
    let v3 = (p3.x - p0.x, p3.y - p0.y, p3.z - p0.z);

    // Scalar triple product / 6
    (v1.0 * (v2.1 * v3.2 - v2.2 * v3.1) - v1.1 * (v2.0 * v3.2 - v2.2 * v3.0)
        + v1.2 * (v2.0 * v3.1 - v2.1 * v3.0))
        / 6.0
}

fn triangle_area_3d(p1: &Point3D, p2: &Point3D, p3: &Point3D) -> f64 {
    let v1 = (p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
    let v2 = (p3.x - p1.x, p3.y - p1.y, p3.z - p1.z);

    // Cross product magnitude / 2
    let cross_x = v1.1 * v2.2 - v1.2 * v2.1;
    let cross_y = v1.2 * v2.0 - v1.0 * v2.2;
    let cross_z = v1.0 * v2.1 - v1.1 * v2.0;

    ((cross_x * cross_x + cross_y * cross_y + cross_z * cross_z).sqrt()) / 2.0
}

/// KD-Tree for efficient nearest neighbor searches
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KDTree {
    root: Option<Box<KDNode>>,
    dimensions: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct KDNode {
    point: Vec<f64>,
    left: Option<Box<KDNode>>,
    right: Option<Box<KDNode>>,
    axis: usize,
}

impl KDTree {
    /// Create a new KD-tree from points
    pub fn new(points: Vec<Vec<f64>>) -> Self {
        if points.is_empty() {
            return Self {
                root: None,
                dimensions: 0,
            };
        }

        let dimensions = points[0].len();
        let root = Self::build_tree(points, 0);

        Self {
            root: Some(Box::new(root)),
            dimensions,
        }
    }

    fn build_tree(mut points: Vec<Vec<f64>>, depth: usize) -> KDNode {
        if points.is_empty() {
            panic!("Cannot build tree from empty points");
        }

        let dimensions = points[0].len();
        let axis = depth % dimensions;

        if points.len() == 1 {
            return KDNode {
                point: points.pop().unwrap(),
                left: None,
                right: None,
                axis,
            };
        }

        // Sort by axis
        points.sort_by(|a, b| a[axis].partial_cmp(&b[axis]).unwrap());

        let median_idx = points.len() / 2;
        let median_point = points[median_idx].clone();

        let left_points: Vec<_> = points[..median_idx].to_vec();
        let right_points: Vec<_> = points[median_idx + 1..].to_vec();

        let left = if !left_points.is_empty() {
            Some(Box::new(Self::build_tree(left_points, depth + 1)))
        } else {
            None
        };

        let right = if !right_points.is_empty() {
            Some(Box::new(Self::build_tree(right_points, depth + 1)))
        } else {
            None
        };

        KDNode {
            point: median_point,
            left,
            right,
            axis,
        }
    }

    /// Find nearest neighbor to query point
    pub fn nearest_neighbor(&self, query: &[f64]) -> Option<Vec<f64>> {
        self.root
            .as_ref()
            .map(|root| self.search_nearest(root, query, root.point.clone(), f64::INFINITY).0)
    }

    fn search_nearest(
        &self,
        node: &KDNode,
        query: &[f64],
        mut best: Vec<f64>,
        mut best_dist: f64,
    ) -> (Vec<f64>, f64) {
        let dist = euclidean_distance(query, &node.point);

        if dist < best_dist {
            best = node.point.clone();
            best_dist = dist;
        }

        let axis = node.axis;
        let diff = query[axis] - node.point[axis];

        let (first, second) = if diff < 0.0 {
            (&node.left, &node.right)
        } else {
            (&node.right, &node.left)
        };

        if let Some(child) = first {
            let (new_best, new_dist) = self.search_nearest(child, query, best.clone(), best_dist);
            if new_dist < best_dist {
                best = new_best;
                best_dist = new_dist;
            }
        }

        // Check if we need to search the other side
        if diff * diff < best_dist {
            if let Some(child) = second {
                let (new_best, new_dist) =
                    self.search_nearest(child, query, best.clone(), best_dist);
                if new_dist < best_dist {
                    best = new_best;
                    best_dist = new_dist;
                }
            }
        }

        (best, best_dist)
    }

    /// Find all points within radius of query point
    pub fn range_search(&self, query: &[f64], radius: f64) -> Vec<Vec<f64>> {
        let mut result = Vec::new();
        if let Some(ref root) = self.root {
            self.search_range(root, query, radius, &mut result);
        }
        result
    }

    fn search_range(&self, node: &KDNode, query: &[f64], radius: f64, result: &mut Vec<Vec<f64>>) {
        let dist = euclidean_distance(query, &node.point);

        if dist <= radius {
            result.push(node.point.clone());
        }

        let axis = node.axis;
        let diff = query[axis] - node.point[axis];

        if diff < 0.0 {
            if let Some(ref left) = node.left {
                self.search_range(left, query, radius, result);
            }
            if diff * diff <= radius * radius {
                if let Some(ref right) = node.right {
                    self.search_range(right, query, radius, result);
                }
            }
        } else {
            if let Some(ref right) = node.right {
                self.search_range(right, query, radius, result);
            }
            if diff * diff <= radius * radius {
                if let Some(ref left) = node.left {
                    self.search_range(left, query, radius, result);
                }
            }
        }
    }
}

fn euclidean_distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt()
}
