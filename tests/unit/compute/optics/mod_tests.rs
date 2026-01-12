// Unit tests for optics::mod
use computational_engine::optics::mod::*;

use super::*;

    // Thin Lens Tests (4 tests)
    #[test]
    fn test_thin_lens_converging() {
        let params = OpticsParams {
            focal_length: Some(0.1),    // 10 cm converging lens
            object_distance: Some(0.2), // 20 cm
            ..Default::default()
        };

        let result = calculate_thin_lens(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // 1/di = 1/f - 1/do = 1/0.1 - 1/0.2 = 10 - 5 = 5 → di = 0.2 m
        assert!((secondary["image_distance"] - 0.2).abs() < 0.01);
        assert!(result.interpretation.contains("real"));
    }

    #[test]
    fn test_thin_lens_diverging() {
        let params = OpticsParams {
            focal_length: Some(-0.15),  // 15 cm diverging lens
            object_distance: Some(0.3), // 30 cm
            ..Default::default()
        };

        let result = calculate_thin_lens(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // 1/di = 1/f - 1/do = 1/-0.15 - 1/0.3 = -6.67 - 3.33 = -10 → di = -0.1 m (virtual)
        assert!(secondary["image_distance"] < 0.0);
        assert!(result.interpretation.contains("virtual"));
    }

    #[test]
    fn test_thin_lens_magnification() {
        let params = OpticsParams {
            focal_length: Some(0.05),   // 5 cm
            object_distance: Some(0.1), // 10 cm (2f)
            ..Default::default()
        };

        let result = calculate_thin_lens(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // At 2f, image forms at 2f with magnification = -1
        assert!((secondary["image_distance"] - 0.1).abs() < 0.01);
        assert!((secondary["magnification"].abs() - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_thin_lens_virtual_image() {
        let params = OpticsParams {
            focal_length: Some(0.2),    // 20 cm
            object_distance: Some(0.1), // 10 cm (inside focal length)
            ..Default::default()
        };

        let result = calculate_thin_lens(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // Inside f gives virtual, upright, enlarged image
        assert!(secondary["image_distance"] < 0.0);
        assert!(secondary["magnification"] > 0.0); // Upright
        assert!(secondary["magnification"].abs() > 1.0); // Enlarged
    }

    // Snell's Law Tests (4 tests)
    #[test]
    fn test_snells_law_air_to_glass() {
        let params = OpticsParams {
            n1: Some(1.0), // air
            n2: Some(1.5), // glass
            theta1: Some(30.0),
            ..Default::default()
        };

        let result = calculate_snells_law(&params).unwrap();
        // Should bend toward normal (smaller angle)
        assert!(result.primary_value < 30.0);
        // n1·sin(30°) = 1.5·sin(θ2) → sin(θ2) = 0.5/1.5 = 0.333 → θ2 ≈ 19.47°
        assert!((result.primary_value - 19.47).abs() < 0.5);
    }

    #[test]
    fn test_snells_law_total_internal_reflection() {
        let params = OpticsParams {
            n1: Some(1.5),      // glass
            n2: Some(1.0),      // air
            theta1: Some(50.0), // Greater than critical angle
            ..Default::default()
        };

        let result = calculate_snells_law(&params).unwrap();
        // Should report total internal reflection
        assert!(result.interpretation.contains("Total internal reflection"));
        // Critical angle = arcsin(1.0/1.5) ≈ 41.8°
        assert!((result.primary_value - 41.8).abs() < 1.0);
    }

    #[test]
    fn test_snells_law_reverse_calculation() {
        let params = OpticsParams {
            n1: Some(1.0),
            n2: Some(1.33),     // water
            theta2: Some(22.0), // Refracted angle given
            ..Default::default()
        };

        let result = calculate_snells_law(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // Calculate incident angle from refracted angle
        // sin(θ1) = 1.33·sin(22°) / 1.0 ≈ 0.498 → θ1 ≈ 29.9°
        assert!(secondary["incident_angle"] > 22.0);
        assert!((secondary["incident_angle"] - 29.9).abs() < 1.0);
    }

    #[test]
    fn test_snells_law_critical_angle() {
        let params = OpticsParams {
            n1: Some(1.5),
            n2: Some(1.0),
            theta1: Some(30.0), // Below critical angle
            ..Default::default()
        };

        let result = calculate_snells_law(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // Critical angle should be calculated
        assert!(secondary.contains_key("critical_angle"));
        // Critical = arcsin(1.0/1.5) ≈ 41.8°
        assert!((secondary["critical_angle"] - 41.8).abs() < 1.0);
    }

    // Diffraction Grating Tests (3 tests)
    #[test]
    fn test_diffraction_first_order() {
        let params = OpticsParams {
            grating_spacing: Some(2e-6), // 2 μm (500 lines/mm)
            wavelength: Some(600e-9),    // 600 nm (red light)
            order: Some(1),
            ..Default::default()
        };

        let result = calculate_diffraction(&params).unwrap();
        // sin(θ) = m·λ/d = 1·600e-9/2e-6 = 0.3 → θ ≈ 17.46°
        assert!((result.primary_value - 17.46).abs() < 1.0);
        assert_eq!(result.unit, "degrees");
    }

    #[test]
    fn test_diffraction_higher_order() {
        let params = OpticsParams {
            grating_spacing: Some(1e-6), // 1 μm
            wavelength: Some(500e-9),    // 500 nm
            order: Some(2),              // Second order
            ..Default::default()
        };

        let result = calculate_diffraction(&params).unwrap();
        // sin(θ) = 2·500e-9/1e-6 = 1.0 → θ = 90° (grazing angle)
        assert!((result.primary_value - 90.0).abs() < 1.0);
        assert!(result.interpretation.contains("Order 2"));
    }

    #[test]
    fn test_diffraction_max_order() {
        let params = OpticsParams {
            grating_spacing: Some(3e-6), // 3 μm
            wavelength: Some(600e-9),    // 600 nm
            order: Some(1),
            ..Default::default()
        };

        let result = calculate_diffraction(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // m_max = d/λ = 3e-6/600e-9 = 5
        assert_eq!(secondary["max_order"] as i32, 5);
    }

    // Fresnel Equations Tests (3 tests)
    #[test]
    fn test_fresnel_s_polarization() {
        let params = OpticsParams {
            n1: Some(1.0), // air
            n2: Some(1.5), // glass
            theta1: Some(45.0),
            polarization: Some("s".to_string()),
            ..Default::default()
        };

        let result = calculate_fresnel(&params).unwrap();
        // Should have reflection and transmission coefficients
        let secondary = result.secondary_values.unwrap();
        assert!(secondary.contains_key("reflection_coefficient"));
        assert!(secondary.contains_key("transmission_coefficient"));
        // Sum should be close to 1 (energy conservation)
        let sum = secondary["reflection_coefficient"] + secondary["transmission_coefficient"];
        assert!((sum - 1.0).abs() < 0.1);
    }

    #[test]
    fn test_fresnel_p_polarization() {
        let params = OpticsParams {
            n1: Some(1.0),
            n2: Some(1.5),
            theta1: Some(30.0),
            polarization: Some("p".to_string()),
            ..Default::default()
        };

        let result = calculate_fresnel(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // p-polarization should include Brewster's angle
        assert!(secondary.contains_key("brewsters_angle"));
        // Brewster = arctan(n2/n1) = arctan(1.5) ≈ 56.3°
        assert!((secondary["brewsters_angle"] - 56.3).abs() < 1.0);
    }

    #[test]
    fn test_fresnel_brewster_angle() {
        let params = OpticsParams {
            n1: Some(1.0),
            n2: Some(1.5),
            theta1: Some(56.3), // At Brewster's angle
            polarization: Some("p".to_string()),
            ..Default::default()
        };

        let result = calculate_fresnel(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // At Brewster's angle, p-polarized reflection should be minimal
        assert!(secondary["reflection_coefficient"] < 0.01);
        // Transmission should be nearly 100%
        assert!(secondary["transmission_coefficient"] > 0.9);
    }
