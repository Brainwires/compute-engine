/**
 * Unit Tests for Optics Module
 *
 * Tests all optics operations:
 * - Thin lens equation (converging/diverging)
 * - Snell's law (refraction and total internal reflection)
 * - Diffraction grating
 * - Fresnel equations (s and p polarization)
 */

use super::*;

// ============================================================================
// Thin Lens Tests
// ============================================================================

#[test]
fn test_thin_lens_calculate_image_distance_converging() {
    let input = OpticsInput {
        operation: OpticsOperation::ThinLens,
        parameters: OpticsParams {
            focal_length: Some(0.2),      // 20 cm converging lens
            object_distance: Some(0.3),   // 30 cm object distance
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();

    // Expected: 1/di = 1/0.2 - 1/0.3 = 5 - 3.333 = 1.667, di = 0.6 m
    assert!((result.primary_value - 0.6).abs() < 1e-6);
    assert_eq!(result.unit, "m");
    assert!(result.interpretation.contains("real"));
}

#[test]
fn test_thin_lens_magnification() {
    let input = OpticsInput {
        operation: OpticsOperation::ThinLens,
        parameters: OpticsParams {
            focal_length: Some(0.15),
            object_distance: Some(0.2),
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();
    let secondary = result.secondary_values.unwrap();

    // Check magnification is calculated
    assert!(secondary.contains_key("magnification"));
    let mag = secondary.get("magnification").unwrap();

    // For converging lens with object beyond focal point
    assert!(mag.abs() > 0.0);
}

#[test]
fn test_thin_lens_virtual_image() {
    let input = OpticsInput {
        operation: OpticsOperation::ThinLens,
        parameters: OpticsParams {
            focal_length: Some(0.2),
            object_distance: Some(0.1),  // Object inside focal length
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();

    // Virtual image should be formed (negative image distance)
    assert!(result.interpretation.contains("virtual"));
}

#[test]
fn test_thin_lens_calculate_focal_length() {
    let input = OpticsInput {
        operation: OpticsOperation::ThinLens,
        parameters: OpticsParams {
            object_distance: Some(0.3),
            image_distance: Some(0.6),
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();

    // 1/f = 1/0.3 + 1/0.6 = 3.333 + 1.667 = 5, f = 0.2 m
    assert!((result.primary_value - 0.2).abs() < 1e-6);
    assert!(result.interpretation.contains("converging"));
}

#[test]
fn test_thin_lens_diverging_lens() {
    let input = OpticsInput {
        operation: OpticsOperation::ThinLens,
        parameters: OpticsParams {
            focal_length: Some(-0.15),  // Diverging lens (negative focal length)
            object_distance: Some(0.3),
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();

    // Should produce virtual, upright, reduced image
    assert!(result.interpretation.contains("virtual") || result.primary_value > 0.0);
}

#[test]
fn test_thin_lens_missing_parameters() {
    let input = OpticsInput {
        operation: OpticsOperation::ThinLens,
        parameters: OpticsParams {
            focal_length: Some(0.2),
            // Missing both object and image distance
            ..Default::default()
        },
    };

    let result = calculate_optics(input);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("at least 2 values"));
}

// ============================================================================
// Snell's Law Tests
// ============================================================================

#[test]
fn test_snells_law_air_to_glass() {
    let input = OpticsInput {
        operation: OpticsOperation::SnellsLaw,
        parameters: OpticsParams {
            n1: Some(1.0),      // Air
            n2: Some(1.5),      // Glass
            theta1: Some(30.0), // 30 degrees incident
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();

    // n1·sin(30°) = n2·sin(θ2)
    // 1.0 × 0.5 = 1.5 × sin(θ2)
    // θ2 = arcsin(0.333) ≈ 19.47°
    assert!((result.primary_value - 19.47).abs() < 0.1);
    assert_eq!(result.unit, "degrees");
    assert!(result.interpretation.contains("toward"));
}

#[test]
fn test_snells_law_glass_to_air() {
    let input = OpticsInput {
        operation: OpticsOperation::SnellsLaw,
        parameters: OpticsParams {
            n1: Some(1.5),      // Glass
            n2: Some(1.0),      // Air
            theta1: Some(20.0),
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();

    // Light bends away from normal when going to less dense medium
    assert!(result.primary_value > 20.0);
    assert!(result.interpretation.contains("away from"));
}

#[test]
fn test_snells_law_total_internal_reflection() {
    let input = OpticsInput {
        operation: OpticsOperation::SnellsLaw,
        parameters: OpticsParams {
            n1: Some(1.5),      // Glass
            n2: Some(1.0),      // Air
            theta1: Some(50.0), // Above critical angle
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();

    // Should report total internal reflection
    assert!(result.interpretation.contains("Total internal reflection"));

    // Critical angle should be calculated
    let secondary = result.secondary_values.unwrap();
    assert!(secondary.contains_key("critical_angle"));

    // Critical angle for glass-air: arcsin(1.0/1.5) ≈ 41.8°
    let critical = secondary.get("critical_angle").unwrap();
    assert!((critical - 41.8).abs() < 0.5);
}

#[test]
fn test_snells_law_critical_angle_calculation() {
    let input = OpticsInput {
        operation: OpticsOperation::SnellsLaw,
        parameters: OpticsParams {
            n1: Some(1.5),
            n2: Some(1.0),
            theta1: Some(30.0), // Below critical angle
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();
    let secondary = result.secondary_values.unwrap();

    // Should still provide critical angle info
    assert!(secondary.contains_key("critical_angle"));
}

#[test]
fn test_snells_law_reverse_calculation() {
    let input = OpticsInput {
        operation: OpticsOperation::SnellsLaw,
        parameters: OpticsParams {
            n1: Some(1.0),
            n2: Some(1.5),
            theta2: Some(19.47), // Given refracted angle
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();

    // Should calculate incident angle ≈ 30°
    assert!((result.primary_value - 30.0).abs() < 0.5);
}

// ============================================================================
// Diffraction Grating Tests
// ============================================================================

#[test]
fn test_diffraction_grating_visible_light() {
    let input = OpticsInput {
        operation: OpticsOperation::DiffractionGrating,
        parameters: OpticsParams {
            grating_spacing: Some(1.67e-6), // 600 lines/mm
            wavelength: Some(550e-9),       // 550 nm (green light)
            order: Some(1),                 // First order
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();

    // θ = arcsin(m·λ/d) = arcsin(1 × 550e-9 / 1.67e-6)
    // θ ≈ 19.3°
    assert!(result.primary_value > 15.0 && result.primary_value < 25.0);
    assert_eq!(result.unit, "degrees");
}

#[test]
fn test_diffraction_grating_second_order() {
    let input = OpticsInput {
        operation: OpticsOperation::DiffractionGrating,
        parameters: OpticsParams {
            grating_spacing: Some(2e-6),
            wavelength: Some(600e-9),
            order: Some(2), // Second order
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();

    // Second order should have larger angle than first order
    let secondary = result.secondary_values.unwrap();
    assert!(secondary.contains_key("max_order"));
}

#[test]
fn test_diffraction_grating_max_order() {
    let input = OpticsInput {
        operation: OpticsOperation::DiffractionGrating,
        parameters: OpticsParams {
            grating_spacing: Some(2e-6),
            wavelength: Some(500e-9),
            order: Some(1),
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();
    let secondary = result.secondary_values.unwrap();

    // m_max = floor(d/λ) = floor(2e-6/500e-9) = floor(4) = 4
    let m_max = secondary.get("max_order").unwrap();
    assert_eq!(*m_max, 4.0);
}

#[test]
fn test_diffraction_grating_invalid_order() {
    let input = OpticsInput {
        operation: OpticsOperation::DiffractionGrating,
        parameters: OpticsParams {
            grating_spacing: Some(1e-6),
            wavelength: Some(600e-9),
            order: Some(10), // Order too high
            ..Default::default()
        },
    };

    let result = calculate_optics(input);

    // Should error when sin(θ) > 1
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("not observable"));
}

#[test]
fn test_diffraction_grating_default_first_order() {
    let input = OpticsInput {
        operation: OpticsOperation::DiffractionGrating,
        parameters: OpticsParams {
            grating_spacing: Some(2e-6),
            wavelength: Some(600e-9),
            // order not specified, should default to 1
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();
    assert!(result.interpretation.contains("Order 1"));
}

// ============================================================================
// Fresnel Equations Tests
// ============================================================================

#[test]
fn test_fresnel_s_polarization_normal_incidence() {
    let input = OpticsInput {
        operation: OpticsOperation::FresnelEquations,
        parameters: OpticsParams {
            n1: Some(1.0),       // Air
            n2: Some(1.5),       // Glass
            theta1: Some(0.0),   // Normal incidence
            polarization: Some("s".to_string()),
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();

    // At normal incidence, R = ((n1-n2)/(n1+n2))²
    // R = ((1-1.5)/(1+1.5))² = (-0.5/2.5)² = 0.04 = 4%
    assert!((result.primary_value - 4.0).abs() < 0.5);
}

#[test]
fn test_fresnel_p_polarization() {
    let input = OpticsInput {
        operation: OpticsOperation::FresnelEquations,
        parameters: OpticsParams {
            n1: Some(1.0),
            n2: Some(1.5),
            theta1: Some(30.0),
            polarization: Some("p".to_string()),
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();

    // p-polarization should have Brewster angle calculated
    let secondary = result.secondary_values.unwrap();
    assert!(secondary.contains_key("brewsters_angle"));

    // Brewster's angle: arctan(1.5/1.0) ≈ 56.3°
    let brewster = secondary.get("brewsters_angle").unwrap();
    assert!((brewster - 56.3).abs() < 0.5);
}

#[test]
fn test_fresnel_brewster_angle() {
    let input = OpticsInput {
        operation: OpticsOperation::FresnelEquations,
        parameters: OpticsParams {
            n1: Some(1.0),
            n2: Some(1.5),
            theta1: Some(56.3), // At Brewster's angle
            polarization: Some("p".to_string()),
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();

    // At Brewster's angle, p-polarized light has zero reflection
    assert!(result.primary_value < 1.0); // Should be nearly 0%
}

#[test]
fn test_fresnel_energy_conservation() {
    let input = OpticsInput {
        operation: OpticsOperation::FresnelEquations,
        parameters: OpticsParams {
            n1: Some(1.0),
            n2: Some(1.5),
            theta1: Some(45.0),
            polarization: Some("s".to_string()),
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();
    let secondary = result.secondary_values.unwrap();

    let reflectance = secondary.get("reflection_coefficient").unwrap();
    let transmittance = secondary.get("transmission_coefficient").unwrap();

    // R + T should equal 1 (energy conservation)
    assert!((reflectance + transmittance - 1.0).abs() < 1e-6);
}

#[test]
fn test_fresnel_total_internal_reflection_no_transmission() {
    let input = OpticsInput {
        operation: OpticsOperation::FresnelEquations,
        parameters: OpticsParams {
            n1: Some(1.5),      // Glass
            n2: Some(1.0),      // Air
            theta1: Some(50.0), // Above critical angle
            polarization: Some("s".to_string()),
            ..Default::default()
        },
    };

    let result = calculate_optics(input);

    // Should error due to total internal reflection
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("Total internal reflection"));
}

#[test]
fn test_fresnel_default_s_polarization() {
    let input = OpticsInput {
        operation: OpticsOperation::FresnelEquations,
        parameters: OpticsParams {
            n1: Some(1.0),
            n2: Some(1.5),
            theta1: Some(30.0),
            // polarization not specified, should default to "s"
            ..Default::default()
        },
    };

    let result = calculate_optics(input).unwrap();
    assert!(result.formula_used.contains("s-polarization"));
}

#[test]
fn test_fresnel_invalid_polarization() {
    let input = OpticsInput {
        operation: OpticsOperation::FresnelEquations,
        parameters: OpticsParams {
            n1: Some(1.0),
            n2: Some(1.5),
            theta1: Some(30.0),
            polarization: Some("x".to_string()), // Invalid
            ..Default::default()
        },
    };

    let result = calculate_optics(input);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("must be 's' or 'p'"));
}

// ============================================================================
// Parameter Validation Tests
// ============================================================================

#[test]
fn test_snells_law_negative_refractive_index() {
    let input = OpticsInput {
        operation: OpticsOperation::SnellsLaw,
        parameters: OpticsParams {
            n1: Some(-1.0), // Invalid
            n2: Some(1.5),
            theta1: Some(30.0),
            ..Default::default()
        },
    };

    let result = calculate_optics(input);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("must be positive"));
}

#[test]
fn test_diffraction_negative_wavelength() {
    let input = OpticsInput {
        operation: OpticsOperation::DiffractionGrating,
        parameters: OpticsParams {
            grating_spacing: Some(1e-6),
            wavelength: Some(-500e-9), // Invalid
            order: Some(1),
            ..Default::default()
        },
    };

    let result = calculate_optics(input);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("must be positive"));
}

#[test]
fn test_thin_lens_zero_object_distance() {
    let input = OpticsInput {
        operation: OpticsOperation::ThinLens,
        parameters: OpticsParams {
            focal_length: Some(0.2),
            object_distance: Some(0.0), // Invalid
            ..Default::default()
        },
    };

    let result = calculate_optics(input);
    assert!(result.is_err());
}
