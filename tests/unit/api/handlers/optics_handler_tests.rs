//! Unit tests for optics API handler
//!
//! Tests all optics operations including:
//! - Thin lens equation (focal length, image distance, object distance)
//! - Snell's law (refraction, total internal reflection, critical angle)
//! - Diffraction grating (angle calculation, maximum order)
//! - Fresnel equations (reflection/transmission coefficients, polarization)
//! - Error handling for invalid inputs

use crate::api::handlers::optics;
use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};
use std::collections::HashMap;

/// Helper function to create an optics request
fn create_request(operation: &str, parameters: HashMap<String, Value>) -> ComputationRequest {
    ComputationRequest {
        module: "optics".to_string(),
        operation: operation.to_string(),
        parameters,
    }
}

/// Helper function to extract result value from response
fn get_result(response: &ComputationResponse) -> &Value {
    assert!(response.success, "Expected success response");
    response.result.as_ref().expect("Expected result value")
}

// ============================================================================
// Thin Lens Tests
// ============================================================================

#[test]
fn test_thin_lens_calculate_image_distance() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("thin_lens"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "focal_length": 0.1,     // 10 cm = 0.1 m
            "object_distance": 0.3,  // 30 cm = 0.3 m
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // 1/di = 1/f - 1/do = 1/0.1 - 1/0.3 = 10 - 3.33 = 6.67
    // di = 0.15 m
    let primary_value = result["primary_value"].as_f64().unwrap();
    assert!((primary_value - 0.15).abs() < 0.01);

    assert_eq!(result["unit"].as_str().unwrap(), "m");
    assert!(result["interpretation"].as_str().unwrap().contains("image"));
}

#[test]
fn test_thin_lens_calculate_focal_length() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("thin_lens"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "object_distance": 0.2,  // 20 cm
            "image_distance": 0.1,   // 10 cm
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // 1/f = 1/do + 1/di = 1/0.2 + 1/0.1 = 5 + 10 = 15
    // f = 1/15 ≈ 0.0667 m
    let primary_value = result["primary_value"].as_f64().unwrap();
    assert!((primary_value - 0.0667).abs() < 0.001);

    assert!(result["interpretation"].as_str().unwrap().contains("lens"));
}

#[test]
fn test_thin_lens_calculate_object_distance() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("thin_lens"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "focal_length": 0.05,    // 5 cm
            "image_distance": 0.15,  // 15 cm
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // 1/do = 1/f - 1/di = 1/0.05 - 1/0.15 = 20 - 6.67 = 13.33
    // do = 0.075 m
    let primary_value = result["primary_value"].as_f64().unwrap();
    assert!((primary_value - 0.075).abs() < 0.001);
}

#[test]
fn test_thin_lens_magnification() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("thin_lens"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "focal_length": 0.1,
            "object_distance": 0.15,
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Check secondary values contain magnification
    let secondary = result["secondary_values"].as_object().unwrap();
    let magnification = secondary["magnification"].as_f64().unwrap();
    assert!(magnification.abs() > 0.0);
}

#[test]
fn test_thin_lens_insufficient_parameters() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("thin_lens"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "focal_length": 0.1,
            // Missing both object_distance and image_distance
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("at least 2 values"));
}

// ============================================================================
// Snell's Law Tests
// ============================================================================

#[test]
fn test_snells_law_calculate_refracted_angle() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("snells_law"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "n1": 1.0,    // Air
            "n2": 1.5,    // Glass
            "theta1": 30.0,  // degrees
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // sin(θ2) = (n1/n2) * sin(θ1) = (1/1.5) * sin(30°) = 0.667 * 0.5 = 0.333
    // θ2 = arcsin(0.333) ≈ 19.47°
    let theta2 = result["primary_value"].as_f64().unwrap();
    assert!((theta2 - 19.47).abs() < 0.1);

    assert_eq!(result["unit"].as_str().unwrap(), "degrees");
    assert!(result["interpretation"].as_str().unwrap().contains("refracted"));
}

#[test]
fn test_snells_law_calculate_incident_angle() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("snells_law"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "n1": 1.0,
            "n2": 1.5,
            "theta2": 20.0,
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // sin(θ1) = (n2/n1) * sin(θ2) = 1.5 * sin(20°)
    let theta1 = result["primary_value"].as_f64().unwrap();
    assert!(theta1 > 20.0); // Should be larger due to n2 > n1

    assert!(result["interpretation"].as_str().unwrap().contains("incident"));
}

#[test]
fn test_snells_law_total_internal_reflection() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("snells_law"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "n1": 1.5,    // Glass
            "n2": 1.0,    // Air
            "theta1": 60.0,  // Large angle from high to low index
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Critical angle = arcsin(n2/n1) = arcsin(1/1.5) ≈ 41.8°
    // Since θ1 = 60° > critical angle, total internal reflection occurs
    assert!(result["interpretation"].as_str().unwrap().contains("Total internal reflection"));

    let secondary = result["secondary_values"].as_object().unwrap();
    let critical_angle = secondary["critical_angle"].as_f64().unwrap();
    assert!((critical_angle - 41.8).abs() < 0.5);
}

#[test]
fn test_snells_law_critical_angle() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("snells_law"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "n1": 1.5,
            "n2": 1.0,
            "theta1": 30.0,  // Below critical angle
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should include critical angle in secondary values
    let secondary = result["secondary_values"].as_object().unwrap();
    assert!(secondary.contains_key("critical_angle"));
}

#[test]
fn test_snells_law_missing_refractive_indices() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("snells_law"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "theta1": 30.0,
            // Missing n1 and n2
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Refractive index"));
}

// ============================================================================
// Diffraction Grating Tests
// ============================================================================

#[test]
fn test_diffraction_grating_first_order() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("diffraction_grating"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "grating_spacing": 1.67e-6,  // 1.67 μm (600 lines/mm)
            "wavelength": 500e-9,        // 500 nm (green light)
            "order": 1,
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // sin(θ) = m·λ/d = 1 * 500e-9 / 1.67e-6 = 0.299
    // θ = arcsin(0.299) ≈ 17.4°
    let angle = result["primary_value"].as_f64().unwrap();
    assert!((angle - 17.4).abs() < 1.0);

    assert_eq!(result["unit"].as_str().unwrap(), "degrees");
    assert!(result["interpretation"].as_str().unwrap().contains("Order 1"));
}

#[test]
fn test_diffraction_grating_second_order() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("diffraction_grating"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "grating_spacing": 2e-6,
            "wavelength": 600e-9,
            "order": 2,
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // sin(θ) = 2 * 600e-9 / 2e-6 = 0.6
    let angle = result["primary_value"].as_f64().unwrap();
    assert!(angle > 0.0);
    assert!(result["interpretation"].as_str().unwrap().contains("Order 2"));
}

#[test]
fn test_diffraction_grating_max_order() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("diffraction_grating"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "grating_spacing": 1.67e-6,
            "wavelength": 500e-9,
            "order": 1,
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Check that max_order is calculated
    let secondary = result["secondary_values"].as_object().unwrap();
    let max_order = secondary["max_order"].as_f64().unwrap();
    assert!(max_order >= 1.0);
}

#[test]
fn test_diffraction_grating_order_not_observable() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("diffraction_grating"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "grating_spacing": 500e-9,   // Small spacing
            "wavelength": 600e-9,        // Large wavelength
            "order": 3,                  // High order
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    // Should fail because sin(θ) > 1
    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("not observable"));
}

#[test]
fn test_diffraction_grating_missing_wavelength() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("diffraction_grating"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "grating_spacing": 1.67e-6,
            // Missing wavelength
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Wavelength"));
}

// ============================================================================
// Fresnel Equations Tests
// ============================================================================

#[test]
fn test_fresnel_s_polarization() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("fresnel_equations"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "n1": 1.0,
            "n2": 1.5,
            "theta1": 30.0,
            "polarization": "s",
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Primary value is reflectance percentage
    let reflectance_percent = result["primary_value"].as_f64().unwrap();
    assert!(reflectance_percent >= 0.0);
    assert!(reflectance_percent <= 100.0);

    assert_eq!(result["unit"].as_str().unwrap(), "% reflected");
    assert!(result["formula_used"].as_str().unwrap().contains("s-polarization"));

    // Check secondary values
    let secondary = result["secondary_values"].as_object().unwrap();
    assert!(secondary.contains_key("reflection_coefficient"));
    assert!(secondary.contains_key("transmission_coefficient"));
}

#[test]
fn test_fresnel_p_polarization() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("fresnel_equations"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "n1": 1.0,
            "n2": 1.5,
            "theta1": 30.0,
            "polarization": "p",
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    assert!(result["formula_used"].as_str().unwrap().contains("p-polarization"));

    // p-polarization should include Brewster's angle
    let secondary = result["secondary_values"].as_object().unwrap();
    assert!(secondary.contains_key("brewsters_angle"));

    let brewster = secondary["brewsters_angle"].as_f64().unwrap();
    // Brewster's angle = arctan(n2/n1) = arctan(1.5) ≈ 56.3°
    assert!((brewster - 56.3).abs() < 1.0);
}

#[test]
fn test_fresnel_conservation_of_energy() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("fresnel_equations"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "n1": 1.0,
            "n2": 1.5,
            "theta1": 30.0,
            "polarization": "s",
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // R + T should equal 1 (energy conservation)
    let secondary = result["secondary_values"].as_object().unwrap();
    let r = secondary["reflection_coefficient"].as_f64().unwrap();
    let t = secondary["transmission_coefficient"].as_f64().unwrap();

    assert!((r + t - 1.0).abs() < 0.01, "Energy not conserved: R={}, T={}", r, t);
}

#[test]
fn test_fresnel_total_internal_reflection() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("fresnel_equations"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "n1": 1.5,    // Glass
            "n2": 1.0,    // Air
            "theta1": 60.0,  // Above critical angle
            "polarization": "s",
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    // Should fail with total internal reflection message
    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Total internal reflection"));
}

#[test]
fn test_fresnel_invalid_polarization() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("fresnel_equations"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "n1": 1.0,
            "n2": 1.5,
            "theta1": 30.0,
            "polarization": "invalid",
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Polarization must be"));
}

#[test]
fn test_fresnel_missing_incident_angle() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("fresnel_equations"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "n1": 1.0,
            "n2": 1.5,
            // Missing theta1
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Incident angle"));
}

// ============================================================================
// Error Handling Tests
// ============================================================================

#[test]
fn test_invalid_json_structure() {
    let mut params = HashMap::new();
    params.insert("invalid_field".to_string(), json!("value"));

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Invalid optics request"));
}

#[test]
fn test_response_structure() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("thin_lens"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "focal_length": 0.1,
            "object_distance": 0.3,
        }),
    );

    let request = create_request("optics", params);
    let response = optics::handle(&request);

    // Verify response structure
    assert!(response.success);
    assert_eq!(response.module, "optics");
    assert_eq!(response.operation, "optics");
    assert!(response.result.is_some());
    assert!(response.error.is_none());

    let result = get_result(&response);
    assert!(result["primary_value"].is_f64());
    assert!(result["unit"].is_string());
    assert!(result["formula_used"].is_string());
    assert!(result["interpretation"].is_string());
}
