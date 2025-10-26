// Unit tests for physics::wormholes::metric
use computational_engine::physics::wormholes::metric::*;

use super::*;
    use crate::physics::wormholes::WormholeConfig;

    #[test]
    fn test_redshift_function_zero() {
        let config = WormholeConfig::morris_thorne(100.0);
        let phi = redshift_function(50.0, &config);
        assert_eq!(phi, 0.0);
    }

    #[test]
    fn test_shape_function_at_throat() {
        let config = WormholeConfig::morris_thorne(100.0);
        let r0 = config.throat_radius;
        let b = shape_function(r0, &config);

        // At throat: b(r₀) = r₀
        assert!((b - r0).abs() < 1.0);
    }

    #[test]
    fn test_metric_computation() {
        let config = WormholeConfig::new(100.0, 200.0);
        let coords = SphericalCoordinates {
            t: 0.0,
            r: 150.0,
            theta: std::f64::consts::PI / 2.0,
            phi: 0.0,
        };

        let metric = morris_thorne_metric(&coords, &config);

        // Metric should have correct signs
        assert!(metric.g_tt < 0.0); // Timelike
        assert!(metric.g_rr > 0.0); // Spacelike
        assert!(metric.g_theta_theta > 0.0);
        assert!(metric.g_phi_phi > 0.0);
    }

    #[test]
    fn test_flaring_condition() {
        let config = WormholeConfig::morris_thorne(100.0);

        // Should satisfy flaring condition near throat
        let flares = satisfies_flaring_condition(100.0, &config);
        assert!(flares);
    }

    #[test]
    fn test_proper_distance_positive() {
        let config = WormholeConfig::new(100.0, 200.0);
        let distance = proper_distance(100.0, 200.0, &config, 100);

        // Should be positive and finite
        assert!(distance > 0.0);
        assert!(distance.is_finite());
    }

    #[test]
    fn test_embedding_surface() {
        let config = WormholeConfig::morris_thorne(100.0);

        // At throat, embedding z should be finite
        let z_throat = embedding_surface(100.0, &config);
        assert!(z_throat.is_finite());

        // Far from throat, z should also be finite
        let z_far = embedding_surface(200.0, &config);
        assert!(z_far.is_finite());
    }

    #[test]
    fn test_shape_function_decreases() {
        let config = WormholeConfig::morris_thorne(100.0);

        let b1 = shape_function(100.0, &config);
        let b2 = shape_function(200.0, &config);

        // Shape function should decrease with r for polynomial type
        assert!(b2 < b1);
    }
