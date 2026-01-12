// Unit tests for physics::warp_drive::metric
use computational_engine::compute::physics::warp_drive::metric::*;

use super::*;
    use crate::compute::physics::warp_drive::WarpDriveConfig;

    #[test]
    fn test_shape_function_tanh() {
        let config = WarpDriveConfig::new(1e6, 100.0, 10.0);

        // At center: f should be ~1
        let f_center = shape_function(0.0, &config);
        assert!(f_center > 0.9);

        // At bubble edge: f should be ~0.5
        let f_edge = shape_function(100.0, &config);
        assert!((f_edge - 0.5).abs() < 0.1);

        // Far away: f should be ~0
        let f_far = shape_function(200.0, &config);
        assert!(f_far < 0.1);
    }

    #[test]
    fn test_shape_function_derivative() {
        let config = WarpDriveConfig::new(1e6, 100.0, 10.0);

        // Derivative should be negative (function decreasing from center to edge)
        let df = shape_function_derivative(100.0, &config);
        assert!(df < 0.0);
    }

    #[test]
    fn test_metric_flat_far_from_bubble() {
        let config = WarpDriveConfig::new(1e6, 100.0, 10.0);
        let coords = Coordinates {
            t: 0.0,
            x: 1000.0, // Far from bubble
            y: 0.0,
            z: 0.0,
        };

        let metric = alcubierre_metric(&coords, &config);

        // Should be approximately Minkowski metric far from bubble
        assert!((metric.g_tt + C * C).abs() < 1e6); // g_tt ≈ -c²
        assert!(metric.g_tx.abs() < 1e3);           // g_tx ≈ 0
        assert!((metric.g_xx - 1.0).abs() < 0.1);   // g_xx ≈ 1
    }

    #[test]
    fn test_lapse_function_positive() {
        let config = WarpDriveConfig::subluminal(0.1 * C, 100.0, 10.0);
        let coords = Coordinates {
            t: 0.0,
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };

        let alpha = lapse_function(&coords, &config);
        assert!(alpha > 0.0);
        assert!(alpha <= C);
    }

    #[test]
    fn test_shift_vector_inside_bubble() {
        let config = WarpDriveConfig::new(1e7, 100.0, 10.0);
        let coords = Coordinates {
            t: 0.0,
            x: 0.0, // Inside bubble
            y: 0.0,
            z: 0.0,
        };

        let (beta_x, beta_y, beta_z) = shift_vector(&coords, &config);

        assert!(beta_x.abs() > 0.0); // Should have x-component
        assert_eq!(beta_y, 0.0);      // No y-component
        assert_eq!(beta_z, 0.0);      // No z-component
    }

    #[test]
    fn test_proper_time_dilation() {
        let config = WarpDriveConfig::subluminal(0.5 * C, 100.0, 10.0);
        let coords = Coordinates {
            t: 0.0,
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };

        let dt = 1.0; // 1 second coordinate time
        let dtau = proper_time(dt, &coords, &config);

        // Proper time should be less than coordinate time (time dilation)
        assert!(dtau <= dt);
        assert!(dtau > 0.0);
    }

    #[test]
    fn test_effective_distance() {
        let config = WarpDriveConfig::new(1e7, 100.0, 10.0);

        // Measure effective distance across flat space (far from bubble)
        let dist = effective_distance(1000.0, 1100.0, &config, 0.0);

        // Should be approximately 100m in flat space
        assert!((dist - 100.0).abs() < 10.0);
    }
