// Unit tests for physics::black_holes::schwarzschild
use computational_engine::physics::black_holes::schwarzschild::*;

use super::*;
    use crate::physics::black_holes::BlackHoleConfig;

    #[test]
    fn test_metric_far_from_horizon() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_s = bh.schwarzschild_radius();

        let metric = schwarzschild_metric(10.0 * r_s, &bh);

        // Far from horizon, should approach Minkowski
        assert!(metric.g_tt < 0.0);
        assert!(metric.g_rr > 0.0);
        assert!((metric.g_tt / (-C2) - 0.9).abs() < 0.1);
    }

    #[test]
    fn test_time_dilation_at_horizon() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_s = bh.schwarzschild_radius();

        let factor = time_dilation_factor(r_s, &bh);

        // Time dilation infinite at horizon
        assert_eq!(factor, 0.0);
    }

    #[test]
    fn test_time_dilation_far_away() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_s = bh.schwarzschild_radius();

        let factor = time_dilation_factor(100.0 * r_s, &bh);

        // Should be close to 1 far from BH
        assert!((factor - 1.0).abs() < 0.1);
    }

    #[test]
    fn test_escape_velocity_at_horizon() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_s = bh.schwarzschild_radius();

        let v_esc = escape_velocity(r_s, &bh);

        // At horizon, escape velocity = c
        assert!((v_esc - C).abs() < 1e6);
    }

    #[test]
    fn test_escape_velocity_far_away() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_s = bh.schwarzschild_radius();

        let v_esc = escape_velocity(100.0 * r_s, &bh);

        // Far away, escape velocity should be small
        assert!(v_esc < C * 0.2);
    }

    #[test]
    fn test_tidal_acceleration() {
        let solar_mass = 1.989e30;
        let bh = BlackHoleConfig::schwarzschild(solar_mass);
        let r_s = bh.schwarzschild_radius();

        // Human height: 2m
        let tidal = tidal_acceleration(10.0 * r_s, 2.0, &bh);

        // Should be finite
        assert!(tidal.is_finite());
        assert!(tidal > 0.0);
    }

    #[test]
    fn test_freefall_velocity() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_s = bh.schwarzschild_radius();

        let v = freefall_velocity(3.0 * r_s, &bh);

        // Should be between 0 and c
        assert!(v > 0.0);
        assert!(v < C);
    }

    #[test]
    fn test_tortoise_coordinate() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_s = bh.schwarzschild_radius();

        let r_star = tortoise_coordinate(2.0 * r_s, &bh);

        // Should be finite and positive
        assert!(r_star.is_finite());
    }

    #[test]
    fn test_redshift_factor() {
        let bh = BlackHoleConfig::schwarzschild(1e30);
        let r_s = bh.schwarzschild_radius();

        let z = redshift_factor(2.0 * r_s, &bh);

        // Near horizon, significant redshift
        assert!(z > 1.0);
    }
