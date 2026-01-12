// Unit tests for physics::n_body::barnes_hut
use computational_engine::compute::physics::n_body::barnes_hut::*;

use super::*;

    #[test]
    fn test_octree_build() {
        let bodies = vec![
            Body::new(1.0, Vec3::new(0.0, 0.0, 0.0), Vec3::zero(), "b1"),
            Body::new(1.0, Vec3::new(1.0, 0.0, 0.0), Vec3::zero(), "b2"),
            Body::new(1.0, Vec3::new(0.0, 1.0, 0.0), Vec3::zero(), "b3"),
        ];

        let tree = build_octree(&bodies);

        assert_eq!(tree.total_mass, 3.0);
        assert!(tree.center_of_mass.x > 0.0);
    }

    #[test]
    fn test_octree_force() {
        let bodies = vec![
            Body::new(1e30, Vec3::new(0.0, 0.0, 0.0), Vec3::zero(), "sun"),
            Body::new(1e24, Vec3::new(1e11, 0.0, 0.0), Vec3::zero(), "earth"),
        ];

        let tree = build_octree(&bodies);
        let test_body = Body::new(1.0, Vec3::new(5e10, 0.0, 0.0), Vec3::zero(), "probe");

        let force = tree.calculate_force(&test_body, 0.5, 0.0);

        // Force should point away from origin (towards sun's mass concentration)
        assert!(force.mag() > 0.0);
    }
