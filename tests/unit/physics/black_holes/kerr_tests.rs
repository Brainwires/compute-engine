// Unit tests for physics::black_holes::kerr
use computational_engine::physics::black_holes::kerr::*;

use super::*;
    use crate::physics::black_holes::BlackHoleConfig;

    #[test]
    fn test_ergosphere() {
        let bh = BlackHoleConfig::kerr(1e30, 0.5e30);

        // At the pole (θ = 0), ergosphere should be larger than horizon
        let r_ergo_pole = ergosphere_radius(0.0, &bh);
        let r_h = bh.event_horizon_radius();

        // At equator (θ = π/2), they should be equal (cos²θ = 0)
        let r_ergo_equator = ergosphere_radius(std::f64::consts::PI / 2.0, &bh);

        eprintln!("Pole: {}, Horizon: {}, Equator: {}", r_ergo_pole, r_h, r_ergo_equator);

        // For a Kerr black hole, ergosphere >= horizon everywhere
        assert!(r_ergo_pole >= r_h - 1.0); // Allow small numerical error
        assert!((r_ergo_equator - r_h).abs() < 1.0); // Close to equal at equator
    }

    #[test]
    fn test_frame_dragging() {
        let bh = BlackHoleConfig::kerr(1e30, 0.5e30);
        let omega = frame_dragging_omega(10e3, std::f64::consts::PI / 2.0, &bh);
        assert!(omega.is_finite());
    }
