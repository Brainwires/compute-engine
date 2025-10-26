// Unit tests for mathematics::chaos::attractors
use computational_engine::mathematics::chaos::attractors::*;

use super::*;

    #[test]
    fn test_lorenz_attractor() {
        let initial = Point3D::new(1.0, 1.0, 1.0);
        let config = LorenzConfig::default();

        let trajectory = lorenz_attractor(initial, &config, 0.01, 100);

        assert_eq!(trajectory.points.len(), 101); // Initial + 100 steps
        assert_eq!(trajectory.times.len(), 101);

        // Check that trajectory doesn't blow up
        for point in &trajectory.points {
            assert!(point.x.is_finite());
            assert!(point.y.is_finite());
            assert!(point.z.is_finite());
        }
    }

    #[test]
    fn test_rossler_attractor() {
        let initial = Point3D::new(1.0, 1.0, 1.0);
        let config = RosslerConfig::default();

        let trajectory = rossler_attractor(initial, &config, 0.01, 100);

        assert_eq!(trajectory.points.len(), 101);

        // Check finite values
        for point in &trajectory.points {
            assert!(point.x.is_finite());
            assert!(point.y.is_finite());
            assert!(point.z.is_finite());
        }
    }

    #[test]
    fn test_logistic_map_fixed_point() {
        // For r < 1, logistic map should converge to 0
        let results = logistic_map_bifurcation(0.5, 0.9, 5, 100, 100);

        for (r, values) in &results {
            // All values should be close to 0
            let avg: f64 = values.iter().sum::<f64>() / values.len() as f64;
            assert!(avg < 0.1, "r={}, avg={}", r, avg);
        }
    }

    #[test]
    fn test_logistic_map_chaos() {
        // For r > 3.57, logistic map is chaotic
        let results = logistic_map_bifurcation(3.7, 3.9, 3, 100, 100);

        for (r, values) in &results {
            // Values should be spread out (chaotic)
            let max: f64 = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let min: f64 = values.iter().cloned().fold(f64::INFINITY, f64::min);
            let range = max - min;
            // At least some variation (not all identical)
            assert!(range > 0.01, "r={}, range={}", r, range);
        }
    }

    #[test]
    fn test_lorenz_butterfly_effect() {
        // Small perturbation should lead to divergence
        let initial1 = Point3D::new(1.0, 1.0, 1.0);
        let initial2 = Point3D::new(1.0, 1.0, 1.001);
        let config = LorenzConfig::default();

        let traj1 = lorenz_attractor(initial1, &config, 0.01, 500);
        let traj2 = lorenz_attractor(initial2, &config, 0.01, 500);

        // Distance should grow over time
        let dist_early = traj1.points[10].distance(&traj2.points[10]);
        let dist_late = traj1.points[500].distance(&traj2.points[500]);

        assert!(dist_late > dist_early);
    }
