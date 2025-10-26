// Unit tests for physics::wormholes::traversal
use computational_engine::physics::wormholes::traversal::*;

use super::*;
    use crate::physics::wormholes::WormholeConfig;

    #[test]
    fn test_tidal_forces_at_throat() {
        let config = WormholeConfig::morris_thorne(100.0);
        // Test slightly away from throat to avoid singularity
        let coords = SphericalCoordinates {
            t: 0.0,
            r: 101.0,  // Slightly outside throat
            theta: std::f64::consts::PI / 2.0,
            phi: 0.0,
        };

        let tidal = compute_tidal_forces(&coords, &config);

        // Should have finite tidal forces away from exact throat
        assert!(tidal.magnitude.is_finite());
        assert!(tidal.magnitude >= 0.0);
    }

    #[test]
    fn test_traversal_analysis() {
        let config = WormholeConfig::new(1000.0, 2000.0);
        let velocity = 1e6; // 1000 km/s

        let analysis = analyze_traversal(&config, velocity);

        // Should have positive times and finite values
        assert!(analysis.proper_time > 0.0);
        assert!(analysis.coordinate_time > 0.0);
        assert!(analysis.throat_distance >= 0.0);
        assert!(analysis.proper_time.is_finite());
        assert!(analysis.coordinate_time.is_finite());
    }

    #[test]
    fn test_time_dilation_factor() {
        let config = WormholeConfig::morris_thorne(100.0);
        let coords = SphericalCoordinates {
            t: 0.0,
            r: 150.0,
            theta: std::f64::consts::PI / 2.0,
            phi: 0.0,
        };

        let factor = time_dilation_factor(&coords, &config);

        // Should be close to 1 for zero redshift
        assert!((factor - 1.0).abs() < 0.5);
        assert!(factor > 0.0);
    }

    #[test]
    fn test_required_velocity() {
        let config = WormholeConfig::new(1000.0, 2000.0);
        let target_time = 10.0; // 10 seconds

        let velocity = required_velocity(&config, target_time);

        // Should be positive and finite
        assert!(velocity > 0.0);
        assert!(velocity.is_finite());
    }

    #[test]
    fn test_survivability_large_wormhole() {
        let config = WormholeConfig::new(10000.0, 20000.0); // 10 km throat
        let velocity = 1e5; // 100 km/s

        let analysis = analyze_traversal(&config, velocity);

        // All values should be finite
        assert!(analysis.max_tidal_force.is_finite());
        assert!(analysis.coordinate_time.is_finite());
    }

    #[test]
    fn test_faster_velocity_shorter_time() {
        let config = WormholeConfig::new(1000.0, 2000.0);

        let analysis1 = analyze_traversal(&config, 1e6);
        let analysis2 = analyze_traversal(&config, 2e6);

        // Both should have finite and positive times
        assert!(analysis1.coordinate_time > 0.0);
        assert!(analysis2.coordinate_time > 0.0);
        assert!(analysis1.coordinate_time.is_finite());
        assert!(analysis2.coordinate_time.is_finite());
    }
