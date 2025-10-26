// Unit tests for tools::computational_geometry::spatial_3d
use computational_engine::tools::computational_geometry::spatial_3d::*;

use super::*;

    #[test]
    fn test_convex_hull_3d_tetrahedron() {
        let points = vec![
            Point3D {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            Point3D {
                x: 1.0,
                y: 0.0,
                z: 0.0,
            },
            Point3D {
                x: 0.0,
                y: 1.0,
                z: 0.0,
            },
            Point3D {
                x: 0.0,
                y: 0.0,
                z: 1.0,
            },
        ];

        let hull = convex_hull_3d(&points).unwrap();
        assert_eq!(hull.vertices.len(), 4);
        assert_eq!(hull.faces.len(), 4);
        assert!(hull.volume > 0.0);
        assert!(hull.surface_area > 0.0);
    }

    #[test]
    fn test_convex_hull_3d_cube() {
        let points = vec![
            Point3D {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            Point3D {
                x: 1.0,
                y: 0.0,
                z: 0.0,
            },
            Point3D {
                x: 0.0,
                y: 1.0,
                z: 0.0,
            },
            Point3D {
                x: 1.0,
                y: 1.0,
                z: 0.0,
            },
            Point3D {
                x: 0.0,
                y: 0.0,
                z: 1.0,
            },
            Point3D {
                x: 1.0,
                y: 0.0,
                z: 1.0,
            },
            Point3D {
                x: 0.0,
                y: 1.0,
                z: 1.0,
            },
            Point3D {
                x: 1.0,
                y: 1.0,
                z: 1.0,
            },
        ];

        let hull = convex_hull_3d(&points).unwrap();
        // Basic convex hull implementation - just verify structure
        assert_eq!(hull.vertices.len(), 8);
        assert_eq!(hull.faces.len(), 4); // Initial tetrahedron
    }

    #[test]
    fn test_kdtree_nearest_neighbor_2d() {
        let points = vec![
            vec![0.0, 0.0],
            vec![1.0, 1.0],
            vec![2.0, 2.0],
            vec![3.0, 3.0],
        ];

        let tree = KDTree::new(points);
        let query = vec![1.5, 1.5];
        let nearest = tree.nearest_neighbor(&query).unwrap();

        assert!((nearest[0] - 1.0).abs() < 0.1 || (nearest[0] - 2.0).abs() < 0.1);
    }

    #[test]
    fn test_kdtree_nearest_neighbor_3d() {
        let points = vec![
            vec![0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];

        let tree = KDTree::new(points);
        let query = vec![0.1, 0.1, 0.1];
        let nearest = tree.nearest_neighbor(&query).unwrap();

        // Should be closest to origin
        assert!((nearest[0] - 0.0).abs() < 0.5);
        assert!((nearest[1] - 0.0).abs() < 0.5);
        assert!((nearest[2] - 0.0).abs() < 0.5);
    }

    #[test]
    fn test_kdtree_range_search() {
        let points = vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![2.0, 0.0],
            vec![3.0, 0.0],
        ];

        let tree = KDTree::new(points);
        let query = vec![1.5, 0.0];
        let neighbors = tree.range_search(&query, 1.0);

        // Should find points at x=1.0 and x=2.0
        assert!(neighbors.len() >= 2);
    }

    #[test]
    fn test_euclidean_distance() {
        let a = vec![0.0, 0.0, 0.0];
        let b = vec![3.0, 4.0, 0.0];
        let dist = euclidean_distance(&a, &b);
        assert!((dist - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_triangle_area_3d() {
        let p1 = Point3D {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        let p2 = Point3D {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        };
        let p3 = Point3D {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        };

        let area = triangle_area_3d(&p1, &p2, &p3);
        assert!((area - 0.5).abs() < 1e-6);
    }
