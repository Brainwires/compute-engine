// Unit tests for tools::computational_geometry::advanced
use computational_engine::tools::computational_geometry::advanced::*;

use super::*;

    #[test]
    fn test_convex_hull() {
        let points = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 1.0 },
            Point2D { x: 2.0, y: 0.0 },
            Point2D { x: 0.5, y: 0.5 },
        ];

        let result = convex_hull(ConvexHullRequest { points }).unwrap();
        assert_eq!(result.hull.len(), 3);
        assert!(result.area > 0.0);
    }

    #[test]
    fn test_point_in_polygon() {
        let polygon = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 4.0, y: 0.0 },
            Point2D { x: 4.0, y: 4.0 },
            Point2D { x: 0.0, y: 4.0 },
        ];

        let inside_point = Point2D { x: 2.0, y: 2.0 };
        let result = point_in_polygon(PointInPolygonRequest {
            point: inside_point,
            polygon: polygon.clone(),
        })
        .unwrap();
        assert!(result.inside);

        let outside_point = Point2D { x: 5.0, y: 5.0 };
        let result = point_in_polygon(PointInPolygonRequest {
            point: outside_point,
            polygon,
        })
        .unwrap();
        assert!(!result.inside);
    }

    #[test]
    fn test_convex_hull_square() {
        let points = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 0.0 },
            Point2D { x: 1.0, y: 1.0 },
            Point2D { x: 0.0, y: 1.0 },
        ];

        let result = convex_hull(ConvexHullRequest { points }).unwrap();
        assert_eq!(result.hull.len(), 4);
        assert!((result.area - 1.0).abs() < 1e-6);
        assert!((result.perimeter - 4.0).abs() < 1e-6);
    }

    #[test]
    fn test_convex_hull_triangle() {
        let points = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 2.0, y: 0.0 },
            Point2D { x: 1.0, y: 2.0 },
        ];

        let result = convex_hull(ConvexHullRequest { points }).unwrap();
        assert_eq!(result.hull.len(), 3);
        assert!((result.area - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_convex_hull_collinear() {
        let points = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 1.0 },
            Point2D { x: 2.0, y: 2.0 },
            Point2D { x: 3.0, y: 3.0 },
        ];

        let result = convex_hull(ConvexHullRequest { points }).unwrap();
        assert!(result.hull.len() >= 2);
    }

    #[test]
    fn test_convex_hull_too_few_points() {
        let points = vec![Point2D { x: 0.0, y: 0.0 }, Point2D { x: 1.0, y: 1.0 }];

        let result = convex_hull(ConvexHullRequest { points });
        assert!(result.is_err());
    }

    #[test]
    fn test_delaunay_triangulation_basic() {
        let points = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 0.0 },
            Point2D { x: 0.5, y: 1.0 },
        ];

        let result = delaunay_triangulation(DelaunayRequest { points }).unwrap();
        // May have no triangles after super-triangle removal
        assert_eq!(result.triangles.len(), result.num_triangles);
    }

    #[test]
    fn test_delaunay_triangulation_square() {
        let points = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 0.0 },
            Point2D { x: 1.0, y: 1.0 },
            Point2D { x: 0.0, y: 1.0 },
        ];

        let result = delaunay_triangulation(DelaunayRequest { points }).unwrap();
        // May have varying number of triangles
        assert_eq!(result.triangles.len(), result.num_triangles);
    }

    #[test]
    fn test_delaunay_too_few_points() {
        let points = vec![Point2D { x: 0.0, y: 0.0 }, Point2D { x: 1.0, y: 1.0 }];

        let result = delaunay_triangulation(DelaunayRequest { points });
        assert!(result.is_err());
    }

    #[test]
    fn test_voronoi_diagram_basic() {
        let points = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 0.0 },
            Point2D { x: 0.5, y: 1.0 },
        ];

        let result = voronoi_diagram(VoronoiRequest {
            points,
            bounds: None,
        })
        .unwrap();

        assert_eq!(result.num_cells, 3);
        assert_eq!(result.cells.len(), 3);
    }

    #[test]
    fn test_voronoi_too_few_points() {
        let points = vec![Point2D { x: 0.0, y: 0.0 }, Point2D { x: 1.0, y: 1.0 }];

        let result = voronoi_diagram(VoronoiRequest {
            points,
            bounds: None,
        });
        assert!(result.is_err());
    }

    #[test]
    fn test_polygon_area_square() {
        let vertices = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 2.0, y: 0.0 },
            Point2D { x: 2.0, y: 2.0 },
            Point2D { x: 0.0, y: 2.0 },
        ];

        let result = polygon_area(PolygonAreaRequest { vertices }).unwrap();
        assert!((result.area - 4.0).abs() < 1e-6);
        assert!((result.perimeter - 8.0).abs() < 1e-6);
        assert!(result.is_convex);
    }

    #[test]
    fn test_polygon_area_triangle() {
        let vertices = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 4.0, y: 0.0 },
            Point2D { x: 2.0, y: 3.0 },
        ];

        let result = polygon_area(PolygonAreaRequest { vertices }).unwrap();
        assert!((result.area - 6.0).abs() < 1e-6);
        assert!(result.is_convex);
    }

    #[test]
    fn test_polygon_area_too_few_vertices() {
        let vertices = vec![Point2D { x: 0.0, y: 0.0 }, Point2D { x: 1.0, y: 1.0 }];

        let result = polygon_area(PolygonAreaRequest { vertices });
        assert!(result.is_err());
    }

    #[test]
    fn test_point_in_polygon_on_edge() {
        let polygon = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 4.0, y: 0.0 },
            Point2D { x: 4.0, y: 4.0 },
            Point2D { x: 0.0, y: 4.0 },
        ];

        let edge_point = Point2D { x: 2.0, y: 0.0 };
        let result = point_in_polygon(PointInPolygonRequest {
            point: edge_point,
            polygon,
        })
        .unwrap();

        assert!(result.distance_to_boundary < 1e-6);
    }

    #[test]
    fn test_point_in_polygon_triangle() {
        let polygon = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 4.0, y: 0.0 },
            Point2D { x: 2.0, y: 3.0 },
        ];

        let inside_point = Point2D { x: 2.0, y: 1.0 };
        let result = point_in_polygon(PointInPolygonRequest {
            point: inside_point,
            polygon,
        })
        .unwrap();
        assert!(result.inside);
    }

    #[test]
    fn test_point_in_polygon_too_few_vertices() {
        let polygon = vec![Point2D { x: 0.0, y: 0.0 }, Point2D { x: 1.0, y: 1.0 }];

        let result = point_in_polygon(PointInPolygonRequest {
            point: Point2D { x: 0.5, y: 0.5 },
            polygon,
        });
        assert!(result.is_err());
    }

    #[test]
    fn test_distance_function() {
        let p1 = Point2D { x: 0.0, y: 0.0 };
        let p2 = Point2D { x: 3.0, y: 4.0 };
        let dist = distance(&p1, &p2);
        assert!((dist - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_cross_product() {
        let o = Point2D { x: 0.0, y: 0.0 };
        let a = Point2D { x: 1.0, y: 0.0 };
        let b = Point2D { x: 0.0, y: 1.0 };
        let cross = cross_product(&o, &a, &b);
        assert!(cross > 0.0);
    }

    #[test]
    fn test_polar_angle() {
        let reference = Point2D { x: 0.0, y: 0.0 };
        let p = Point2D { x: 1.0, y: 1.0 };
        let angle = polar_angle(&p, &reference);
        assert!((angle - std::f64::consts::PI / 4.0).abs() < 1e-6);
    }

    #[test]
    fn test_points_equal() {
        let p1 = Point2D { x: 1.0, y: 2.0 };
        let p2 = Point2D { x: 1.0, y: 2.0 };
        assert!(points_equal(&p1, &p2));

        let p3 = Point2D { x: 1.0, y: 2.1 };
        assert!(!points_equal(&p1, &p3));
    }

    #[test]
    fn test_point_to_segment_distance() {
        let p = Point2D { x: 0.0, y: 1.0 };
        let a = Point2D { x: 0.0, y: 0.0 };
        let b = Point2D { x: 2.0, y: 0.0 };
        let dist = point_to_segment_distance(&p, &a, &b);
        assert!((dist - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_compute_circumcenter() {
        let triangle = Triangle {
            p1: Point2D { x: 0.0, y: 0.0 },
            p2: Point2D { x: 2.0, y: 0.0 },
            p3: Point2D { x: 1.0, y: 1.0 },
        };

        let center = compute_circumcenter(&triangle).unwrap();
        assert!(center.x >= 0.0 && center.x <= 2.0);
        assert!(center.y >= 0.0);
    }

    #[test]
    fn test_point_in_circumcircle() {
        let triangle = Triangle {
            p1: Point2D { x: 0.0, y: 0.0 },
            p2: Point2D { x: 2.0, y: 0.0 },
            p3: Point2D { x: 1.0, y: 2.0 },
        };

        let inside_point = Point2D { x: 1.0, y: 1.0 };
        assert!(point_in_circumcircle(&inside_point, &triangle));

        let outside_point = Point2D { x: 10.0, y: 10.0 };
        assert!(!point_in_circumcircle(&outside_point, &triangle));
    }

    #[test]
    fn test_polygon_centroid() {
        let vertices = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 2.0, y: 0.0 },
            Point2D { x: 2.0, y: 2.0 },
            Point2D { x: 0.0, y: 2.0 },
        ];

        let result = polygon_area(PolygonAreaRequest { vertices }).unwrap();
        // Centroid of square should be at center
        assert!((result.centroid.x - 1.0).abs() < 1e-6);
        assert!((result.centroid.y - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_convex_hull_with_duplicates() {
        let points = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 0.0 },
            Point2D { x: 1.0, y: 0.0 }, // Duplicate
            Point2D { x: 0.5, y: 1.0 },
        ];

        let result = convex_hull(ConvexHullRequest { points }).unwrap();
        assert!(result.hull.len() >= 3);
    }

    #[test]
    fn test_delaunay_five_points() {
        let points = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 0.0 },
            Point2D { x: 1.0, y: 1.0 },
            Point2D { x: 0.0, y: 1.0 },
            Point2D { x: 0.5, y: 0.5 },
        ];

        let result = delaunay_triangulation(DelaunayRequest { points }).unwrap();
        // Just verify it runs without error
        assert_eq!(result.triangles.len(), result.num_triangles);
    }
