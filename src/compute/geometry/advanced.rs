//! Advanced Computational Geometry Algorithms
//!
//! Implements advanced geometric algorithms including:
//! - Convex hull computation (Graham scan)
//! - Delaunay triangulation
//! - Voronoi diagrams
//! - Polygon operations
//! - Spatial indexing

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Point2D {
    pub x: f64,
    pub y: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Triangle {
    pub p1: Point2D,
    pub p2: Point2D,
    pub p3: Point2D,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VoronoiCell {
    pub site: Point2D,
    pub vertices: Vec<Point2D>,
}

// Request/Response types
#[derive(Debug, Deserialize)]
pub struct ConvexHullRequest {
    pub points: Vec<Point2D>,
}

#[derive(Debug, Serialize)]
pub struct ConvexHullResult {
    pub hull: Vec<Point2D>,
    pub area: f64,
    pub perimeter: f64,
}

#[derive(Debug, Deserialize)]
pub struct DelaunayRequest {
    pub points: Vec<Point2D>,
}

#[derive(Debug, Serialize)]
pub struct DelaunayResult {
    pub triangles: Vec<Triangle>,
    pub num_triangles: usize,
}

#[derive(Debug, Deserialize)]
pub struct VoronoiRequest {
    pub points: Vec<Point2D>,
    pub bounds: Option<(f64, f64, f64, f64)>, // (min_x, min_y, max_x, max_y)
}

#[derive(Debug, Serialize)]
pub struct VoronoiResult {
    pub cells: Vec<VoronoiCell>,
    pub num_cells: usize,
}

#[derive(Debug, Deserialize)]
pub struct PolygonAreaRequest {
    pub vertices: Vec<Point2D>,
}

#[derive(Debug, Serialize)]
pub struct PolygonAreaResult {
    pub area: f64,
    pub perimeter: f64,
    pub centroid: Point2D,
    pub is_convex: bool,
}

#[derive(Debug, Deserialize)]
pub struct PointInPolygonRequest {
    pub point: Point2D,
    pub polygon: Vec<Point2D>,
}

#[derive(Debug, Serialize)]
pub struct PointInPolygonResult {
    pub inside: bool,
    pub distance_to_boundary: f64,
}

// Helper functions
fn cross_product(o: &Point2D, a: &Point2D, b: &Point2D) -> f64 {
    (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)
}

fn distance(a: &Point2D, b: &Point2D) -> f64 {
    ((a.x - b.x).powi(2) + (a.y - b.y).powi(2)).sqrt()
}

fn polar_angle(p: &Point2D, reference: &Point2D) -> f64 {
    (p.y - reference.y).atan2(p.x - reference.x)
}

/// Compute convex hull using Graham scan algorithm
pub fn convex_hull(request: ConvexHullRequest) -> Result<ConvexHullResult, String> {
    let mut points = request.points;

    if points.len() < 3 {
        return Err("Need at least 3 points for convex hull".to_string());
    }

    // Find lowest point (break ties by x-coordinate)
    let start_idx = points
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            a.y.partial_cmp(&b.y)
                .unwrap()
                .then_with(|| a.x.partial_cmp(&b.x).unwrap())
        })
        .map(|(i, _)| i)
        .unwrap();

    let start = points.swap_remove(start_idx);

    // Sort by polar angle
    points.sort_by(|a, b| {
        polar_angle(a, &start)
            .partial_cmp(&polar_angle(b, &start))
            .unwrap()
            .then_with(|| {
                distance(a, &start)
                    .partial_cmp(&distance(b, &start))
                    .unwrap()
            })
    });

    let mut hull = vec![start.clone()];

    for point in points {
        while hull.len() > 1 {
            let len = hull.len();
            if cross_product(&hull[len - 2], &hull[len - 1], &point) <= 0.0 {
                hull.pop();
            } else {
                break;
            }
        }
        hull.push(point);
    }

    // Calculate area using shoelace formula
    let mut area = 0.0;
    let n = hull.len();
    for i in 0..n {
        let j = (i + 1) % n;
        area += hull[i].x * hull[j].y;
        area -= hull[j].x * hull[i].y;
    }
    area = area.abs() / 2.0;

    // Calculate perimeter
    let mut perimeter = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        perimeter += distance(&hull[i], &hull[j]);
    }

    Ok(ConvexHullResult {
        hull,
        area,
        perimeter,
    })
}

/// Bowyer-Watson algorithm for Delaunay triangulation
pub fn delaunay_triangulation(request: DelaunayRequest) -> Result<DelaunayResult, String> {
    let points = request.points;

    if points.len() < 3 {
        return Err("Need at least 3 points for triangulation".to_string());
    }

    // Find bounding box
    let min_x = points.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let max_x = points.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let min_y = points.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let max_y = points.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);

    let dx = max_x - min_x;
    let dy = max_y - min_y;
    let delta_max = dx.max(dy);
    let mid_x = (min_x + max_x) / 2.0;
    let mid_y = (min_y + max_y) / 2.0;

    // Create super-triangle
    let p1 = Point2D {
        x: mid_x - 20.0 * delta_max,
        y: mid_y - delta_max,
    };
    let p2 = Point2D {
        x: mid_x,
        y: mid_y + 20.0 * delta_max,
    };
    let p3 = Point2D {
        x: mid_x + 20.0 * delta_max,
        y: mid_y - delta_max,
    };

    let mut triangles = vec![Triangle { p1, p2, p3 }];

    // Add points one at a time
    for point in &points {
        let mut bad_triangles = Vec::new();

        for (i, tri) in triangles.iter().enumerate() {
            if point_in_circumcircle(point, tri) {
                bad_triangles.push(i);
            }
        }

        let mut polygon = Vec::new();

        for &idx in &bad_triangles {
            let tri = &triangles[idx];
            polygon.push((tri.p1.clone(), tri.p2.clone()));
            polygon.push((tri.p2.clone(), tri.p3.clone()));
            polygon.push((tri.p3.clone(), tri.p1.clone()));
        }

        // Remove bad triangles (in reverse to maintain indices)
        for &idx in bad_triangles.iter().rev() {
            triangles.remove(idx);
        }

        // Remove duplicate edges
        let mut unique_edges = Vec::new();
        for (i, edge) in polygon.iter().enumerate() {
            let mut is_duplicate = false;
            for (j, other) in polygon.iter().enumerate() {
                if i != j
                    && edge.0.x == other.1.x
                    && edge.0.y == other.1.y
                    && edge.1.x == other.0.x
                    && edge.1.y == other.0.y
                {
                    is_duplicate = true;
                    break;
                }
            }
            if !is_duplicate {
                unique_edges.push(edge.clone());
            }
        }

        // Create new triangles
        for (p1, p2) in unique_edges {
            triangles.push(Triangle {
                p1,
                p2,
                p3: point.clone(),
            });
        }
    }

    // Remove triangles that share vertices with super-triangle
    let super_points = vec![
        Point2D {
            x: mid_x - 20.0 * delta_max,
            y: mid_y - delta_max,
        },
        Point2D {
            x: mid_x,
            y: mid_y + 20.0 * delta_max,
        },
        Point2D {
            x: mid_x + 20.0 * delta_max,
            y: mid_y - delta_max,
        },
    ];

    triangles.retain(|tri| {
        !super_points.iter().any(|sp| {
            (sp.x == tri.p1.x && sp.y == tri.p1.y)
                || (sp.x == tri.p2.x && sp.y == tri.p2.y)
                || (sp.x == tri.p3.x && sp.y == tri.p3.y)
        })
    });

    Ok(DelaunayResult {
        num_triangles: triangles.len(),
        triangles,
    })
}

fn point_in_circumcircle(point: &Point2D, triangle: &Triangle) -> bool {
    let ax = triangle.p1.x - point.x;
    let ay = triangle.p1.y - point.y;
    let bx = triangle.p2.x - point.x;
    let by = triangle.p2.y - point.y;
    let cx = triangle.p3.x - point.x;
    let cy = triangle.p3.y - point.y;

    let det = (ax * ax + ay * ay) * (bx * cy - cx * by) - (bx * bx + by * by) * (ax * cy - cx * ay)
        + (cx * cx + cy * cy) * (ax * by - bx * ay);

    det > 0.0
}

/// Compute Voronoi diagram from Delaunay triangulation
pub fn voronoi_diagram(request: VoronoiRequest) -> Result<VoronoiResult, String> {
    let points = request.points.clone();

    if points.len() < 3 {
        return Err("Need at least 3 points for Voronoi diagram".to_string());
    }

    // First compute Delaunay triangulation
    let delaunay = delaunay_triangulation(DelaunayRequest {
        points: points.clone(),
    })?;

    // Compute circumcenters of triangles (these become Voronoi vertices)
    let mut cells: Vec<VoronoiCell> = points
        .iter()
        .map(|p| VoronoiCell {
            site: p.clone(),
            vertices: Vec::new(),
        })
        .collect();

    for triangle in &delaunay.triangles {
        let circumcenter = compute_circumcenter(triangle)?;

        // Add circumcenter to cells of triangle vertices
        // (Simplified - in real implementation, need proper edge tracking)
        for (i, site) in points.iter().enumerate() {
            if points_equal(site, &triangle.p1)
                || points_equal(site, &triangle.p2)
                || points_equal(site, &triangle.p3)
            {
                cells[i].vertices.push(circumcenter.clone());
            }
        }
    }

    Ok(VoronoiResult {
        num_cells: cells.len(),
        cells,
    })
}

fn compute_circumcenter(triangle: &Triangle) -> Result<Point2D, String> {
    let ax = triangle.p1.x;
    let ay = triangle.p1.y;
    let bx = triangle.p2.x;
    let by = triangle.p2.y;
    let cx = triangle.p3.x;
    let cy = triangle.p3.y;

    let d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));

    if d.abs() < 1e-10 {
        return Err("Degenerate triangle".to_string());
    }

    let ux = ((ax * ax + ay * ay) * (by - cy)
        + (bx * bx + by * by) * (cy - ay)
        + (cx * cx + cy * cy) * (ay - by))
        / d;

    let uy = ((ax * ax + ay * ay) * (cx - bx)
        + (bx * bx + by * by) * (ax - cx)
        + (cx * cx + cy * cy) * (bx - ax))
        / d;

    Ok(Point2D { x: ux, y: uy })
}

fn points_equal(a: &Point2D, b: &Point2D) -> bool {
    (a.x - b.x).abs() < 1e-10 && (a.y - b.y).abs() < 1e-10
}

/// Calculate polygon properties
pub fn polygon_area(request: PolygonAreaRequest) -> Result<PolygonAreaResult, String> {
    let vertices = request.vertices;

    if vertices.len() < 3 {
        return Err("Polygon must have at least 3 vertices".to_string());
    }

    let n = vertices.len();

    // Shoelace formula for area
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        area += vertices[i].x * vertices[j].y;
        area -= vertices[j].x * vertices[i].y;
    }
    area = area.abs() / 2.0;

    // Perimeter
    let mut perimeter = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        perimeter += distance(&vertices[i], &vertices[j]);
    }

    // Centroid
    let mut cx = 0.0;
    let mut cy = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        let factor = vertices[i].x * vertices[j].y - vertices[j].x * vertices[i].y;
        cx += (vertices[i].x + vertices[j].x) * factor;
        cy += (vertices[i].y + vertices[j].y) * factor;
    }
    let area_factor = 6.0 * if area > 0.0 { area } else { -area };
    cx /= area_factor;
    cy /= area_factor;

    // Check if convex
    let mut is_convex = true;
    for i in 0..n {
        let a = &vertices[i];
        let b = &vertices[(i + 1) % n];
        let c = &vertices[(i + 2) % n];
        if cross_product(a, b, c) < 0.0 {
            is_convex = false;
            break;
        }
    }

    Ok(PolygonAreaResult {
        area,
        perimeter,
        centroid: Point2D { x: cx, y: cy },
        is_convex,
    })
}

/// Point-in-polygon test using ray casting algorithm
pub fn point_in_polygon(request: PointInPolygonRequest) -> Result<PointInPolygonResult, String> {
    let point = request.point;
    let polygon = request.polygon;

    if polygon.len() < 3 {
        return Err("Polygon must have at least 3 vertices".to_string());
    }

    let mut inside = false;
    let n = polygon.len();

    for i in 0..n {
        let j = (i + 1) % n;
        let xi = polygon[i].x;
        let yi = polygon[i].y;
        let xj = polygon[j].x;
        let yj = polygon[j].y;

        let intersect = ((yi > point.y) != (yj > point.y))
            && (point.x < (xj - xi) * (point.y - yi) / (yj - yi) + xi);

        if intersect {
            inside = !inside;
        }
    }

    // Calculate distance to boundary
    let mut min_dist = f64::INFINITY;
    for i in 0..n {
        let j = (i + 1) % n;
        let dist = point_to_segment_distance(&point, &polygon[i], &polygon[j]);
        min_dist = min_dist.min(dist);
    }

    Ok(PointInPolygonResult {
        inside,
        distance_to_boundary: min_dist,
    })
}

fn point_to_segment_distance(p: &Point2D, a: &Point2D, b: &Point2D) -> f64 {
    let px = b.x - a.x;
    let py = b.y - a.y;
    let norm = px * px + py * py;

    if norm == 0.0 {
        return distance(p, a);
    }

    let u = ((p.x - a.x) * px + (p.y - a.y) * py) / norm;

    if u < 0.0 {
        distance(p, a)
    } else if u > 1.0 {
        distance(p, b)
    } else {
        let closest = Point2D {
            x: a.x + u * px,
            y: a.y + u * py,
        };
        distance(p, &closest)
    }
}
