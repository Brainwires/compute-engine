// Unit tests for electrical::impedance
use computational_engine::electrical::impedance::*;

use super::*;

    #[test]
    fn test_rectangular_polar_conversion() {
        let z = Complex64::new(3.0, 4.0);
        let (mag, phase) = rectangular_to_polar(z);
        assert!((mag - 5.0).abs() < 0.001);
        assert!((phase - 0.927).abs() < 0.01); // atan(4/3) ≈ 0.927 rad

        let z2 = polar_to_rectangular(5.0, 0.927);
        assert!((z2.re - 3.0).abs() < 0.01);
        assert!((z2.im - 4.0).abs() < 0.01);
    }

    #[test]
    fn test_series_parallel() {
        let z1 = Complex64::new(10.0, 5.0);
        let z2 = Complex64::new(20.0, -10.0);

        let z_series = series_impedances(&[z1, z2]);
        assert_eq!(z_series.re, 30.0);
        assert_eq!(z_series.im, -5.0);

        let z_parallel = parallel_impedances(&[
            Complex64::new(10.0, 0.0),
            Complex64::new(10.0, 0.0),
        ]);
        assert!((z_parallel.re - 5.0).abs() < 0.01);
    }

    #[test]
    fn test_reflection_coefficient() {
        let z_load = Complex64::new(100.0, 0.0);
        let z_0 = Complex64::new(50.0, 0.0);
        let gamma = reflection_coefficient(z_load, z_0);

        assert!((gamma.re - 0.333).abs() < 0.01);
        assert!(gamma.im.abs() < 0.01);

        let vswr = vswr_from_reflection(gamma);
        assert!((vswr - 2.0).abs() < 0.01);
    }

    #[test]
    fn test_quarter_wave_transformer() {
        let z_0 = quarter_wave_transformer_impedance(50.0, 75.0);
        assert!((z_0 - 61.237).abs() < 0.01);
    }
