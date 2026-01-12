// Unit tests for physics::fluid_dynamics::grid
use computational_engine::simulate::fluids::grid::*;

use super::*;

    #[test]
    fn test_grid_creation() {
        let grid = Grid2D::uniform(11, 11, 1.0, 1.0);
        assert_eq!(grid.nx, 11);
        assert_eq!(grid.ny, 11);
        assert!((grid.dx - 0.1).abs() < 1e-10);
        assert!((grid.dy - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_coordinates() {
        let grid = Grid2D::uniform(11, 11, 1.0, 1.0);
        assert!((grid.x_coord(0) - 0.0).abs() < 1e-10);
        assert!((grid.x_coord(10) - 1.0).abs() < 1e-10);
        assert!((grid.y_coord(5) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_boundary_detection() {
        let grid = Grid2D::uniform(5, 5, 1.0, 1.0);
        assert!(grid.is_boundary(0, 2));
        assert!(grid.is_boundary(4, 2));
        assert!(grid.is_boundary(2, 0));
        assert!(grid.is_boundary(2, 4));
        assert!(grid.is_interior(2, 2));
        assert!(!grid.is_interior(0, 0));
    }
