// Unit tests for materials_science::xrd_analysis
use computational_engine::materials_science::xrd_analysis::*;

use super::*;

    #[test]
    fn test_two_theta_calculation() {
        // Silicon (111): d = 3.135 Å with Cu Kα
        let two_theta = calculate_two_theta(3.135, wavelengths::CU_KA).unwrap();
        assert!((two_theta - 28.44).abs() < 0.1); // Should be ~28.44°
    }

    #[test]
    fn test_d_spacing_from_two_theta() {
        // Reverse calculation
        let d = calculate_d_spacing(28.44, wavelengths::CU_KA);
        assert!((d - 3.135).abs() < 0.01);
    }

    #[test]
    fn test_scherrer_equation() {
        // Typical values: λ = 1.54 Å, FWHM = 0.5°, θ = 15°
        let wavelength = 1.54;
        let fwhm_rad = 0.5_f64.to_radians();
        let theta_rad = 15.0_f64.to_radians();

        let size = scherrer_crystallite_size(wavelength, fwhm_rad, theta_rad, 0.9);
        assert!(size > 100.0 && size < 200.0); // Should be ~150 Å
    }

    #[test]
    fn test_lattice_parameter_cubic() {
        // Silicon: FCC, a = 5.431 Å
        // (111) reflection: d = a/√3
        let d = 5.431 / 3.0_f64.sqrt();
        let a = lattice_parameter_from_d_cubic(d, 1, 1, 1);
        assert!((a - 5.431).abs() < 0.001);
    }

    #[test]
    fn test_generate_cubic_pattern() {
        // Silicon FCC with Cu Kα
        let pattern = generate_cubic_pattern(
            5.431,
            BravaisLattice::FaceCenteredCubic,
            wavelengths::CU_KA,
            90.0,
        );

        // Should have several peaks
        assert!(pattern.len() > 5);

        // First peak should be (111) around 28°
        assert!((pattern[0].two_theta - 28.44).abs() < 1.0);
        assert_eq!(pattern[0].hkl, (1, 1, 1));
    }

    #[test]
    fn test_multiplicity() {
        assert_eq!(calculate_multiplicity((1, 0, 0)), 6); // {100}
        assert_eq!(calculate_multiplicity((1, 1, 0)), 24); // {110}
        assert_eq!(calculate_multiplicity((1, 1, 1)), 8); // {111}
        assert_eq!(calculate_multiplicity((2, 1, 0)), 12); // {210} - has zero
    }

    #[test]
    fn test_instrumental_correction() {
        let observed = 0.5; // degrees
        let instrumental = 0.1; // degrees
        let corrected = correct_instrumental_broadening(observed, instrumental);

        assert!((corrected - 0.4899).abs() < 0.001); // √(0.25 - 0.01) ≈ 0.4899
    }

    #[test]
    fn test_lorentz_polarization() {
        // At 2θ = 30°
        let lp = lorentz_polarization_factor(30.0_f64.to_radians());
        assert!(lp > 25.0 && lp < 30.0); // ~27 for 2θ=30°
    }

    #[test]
    fn test_degree_of_crystallinity() {
        let xc = degree_of_crystallinity(80.0, 20.0);
        assert!((xc - 0.8).abs() < 0.001); // 80% crystalline
    }
