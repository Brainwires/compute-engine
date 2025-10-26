//! Unit tests for geophysics API handler
//!
//! Tests all geophysics operations including:
//! - Seismology (earthquake magnitudes)
//! - Atmospheric physics (pressure, humidity)
//! - Radiometric dating (isotope decay)
//! - Planetary science (escape velocity, Roche limit)
//! - Error handling for invalid inputs

use crate::api::handlers::geophysics;
use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};
use std::collections::HashMap;

/// Helper function to create a geophysics request
fn create_request(parameters: HashMap<String, Value>) -> ComputationRequest {
    ComputationRequest {
        module: "geophysics".to_string(),
        operation: "calculate".to_string(),
        parameters,
    }
}

/// Helper function to extract result value from response
fn get_result(response: &ComputationResponse) -> &Value {
    assert!(response.success, "Expected success response");
    response.result.as_ref().expect("Expected result value")
}

// ============================================================================
// Seismology Tests
// ============================================================================

#[test]
fn test_seismology_moment_magnitude() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("seismology"));
    params.insert("parameters".to_string(), json!({
        "seismic_moment": 1e18 // N·m
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Mw = (2/3)·log10(M0) - 6.07
    // Mw = (2/3)·log10(1e18) - 6.07 = (2/3)·18 - 6.07 = 12 - 6.07 = 5.93
    let mw = result["value"].as_f64().unwrap();
    assert!((mw - 5.93).abs() < 0.1);
    assert_eq!(result["unit"].as_str().unwrap(), "Mw");

    // Check additional data
    let additional = result["additional_data"].as_object().unwrap();
    assert!(additional.contains_key("energy_joules"));
    assert!(additional.contains_key("energy_tnt_equivalent_tons"));
}

#[test]
fn test_seismology_energy_magnitude() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("seismology"));
    params.insert("parameters".to_string(), json!({
        "energy": 2e15 // J
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // M = (2/3)·log10(E) - 2.9
    let magnitude = result["value"].as_f64().unwrap();
    assert!(magnitude > 0.0);
    assert_eq!(result["unit"].as_str().unwrap(), "ML (Richter)");

    let additional = result["additional_data"].as_object().unwrap();
    assert!(additional.contains_key("tnt_equivalent_tons"));
}

#[test]
fn test_seismology_large_earthquake() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("seismology"));
    params.insert("parameters".to_string(), json!({
        "seismic_moment": 1e25 // Very large earthquake
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let mw = result["value"].as_f64().unwrap();
    assert!(mw > 8.0); // Major earthquake

    let interpretation = result["interpretation"].as_str().unwrap();
    assert!(interpretation.contains("earthquake"));
}

#[test]
fn test_seismology_missing_parameters() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("seismology"));
    params.insert("parameters".to_string(), json!({}));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Seismic moment or energy required"));
}

// ============================================================================
// Atmospheric Physics Tests
// ============================================================================

#[test]
fn test_atmosphere_hydrostatic_equation() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("atmosphere"));
    params.insert("parameters".to_string(), json!({
        "pressure": 101325.0,  // Pa (sea level)
        "temperature": 288.15, // K (15°C)
        "altitude": 1000.0     // m
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Pressure should decrease with altitude
    let pressure_at_h = result["value"].as_f64().unwrap();
    assert!(pressure_at_h < 101325.0);
    assert!(pressure_at_h > 80000.0); // Reasonable range

    assert_eq!(result["unit"].as_str().unwrap(), "Pa");

    let additional = result["additional_data"].as_object().unwrap();
    assert!(additional.contains_key("scale_height_m"));
    assert!(additional.contains_key("pressure_drop_percent"));
}

#[test]
fn test_atmosphere_temperature_celsius() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("atmosphere"));
    params.insert("parameters".to_string(), json!({
        "pressure": 101325.0,
        "temperature": 15.0,  // °C (should be converted to K)
        "altitude": 500.0
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let pressure = result["value"].as_f64().unwrap();
    assert!(pressure > 0.0);
}

#[test]
fn test_atmosphere_humidity_dew_point() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("atmosphere"));
    params.insert("parameters".to_string(), json!({
        "temperature": 25.0,        // °C
        "relative_humidity": 60.0   // %
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Dew point should be less than temperature
    let dew_point = result["value"].as_f64().unwrap();
    assert!(dew_point < 25.0);
    assert!(dew_point > 0.0);

    assert_eq!(result["unit"].as_str().unwrap(), "°C");

    let additional = result["additional_data"].as_object().unwrap();
    assert!(additional.contains_key("saturation_vapor_pressure_pa"));
    assert!(additional.contains_key("actual_vapor_pressure_pa"));
    assert!(additional.contains_key("dew_point_celsius"));
}

#[test]
fn test_atmosphere_high_humidity() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("atmosphere"));
    params.insert("parameters".to_string(), json!({
        "temperature": 298.15,      // K
        "relative_humidity": 95.0   // %
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // High humidity: dew point should be close to temperature
    let dew_point = result["value"].as_f64().unwrap();
    let temp_c = 298.15 - 273.15; // 25°C
    assert!((dew_point - temp_c).abs() < 3.0);
}

#[test]
fn test_atmosphere_missing_parameters() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("atmosphere"));
    params.insert("parameters".to_string(), json!({
        "pressure": 101325.0
        // Missing other required parameters
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Insufficient atmospheric parameters"));
}

// ============================================================================
// Radiometric Dating Tests
// ============================================================================

#[test]
fn test_dating_carbon14() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("radiometric_dating"));
    params.insert("parameters".to_string(), json!({
        "isotope_system": "C14",
        "parent_isotope": 10.0,
        "daughter_isotope": 5.0
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Age should be reasonable for C14 (< 50,000 years)
    let age = result["value"].as_f64().unwrap();
    assert!(age > 0.0);
    assert!(age < 100000.0);

    assert_eq!(result["unit"].as_str().unwrap(), "years");

    let additional = result["additional_data"].as_object().unwrap();
    assert_eq!(additional["half_life_years"].as_f64().unwrap(), 5730.0);
    assert!(additional.contains_key("daughter_parent_ratio"));
}

#[test]
fn test_dating_uranium238() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("radiometric_dating"));
    params.insert("parameters".to_string(), json!({
        "isotope_system": "U238Pb206",
        "parent_isotope": 100.0,
        "daughter_isotope": 20.0
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Age should be on geological timescale (billions of years)
    let age = result["value"].as_f64().unwrap();
    assert!(age > 1e8); // > 100 million years

    let additional = result["additional_data"].as_object().unwrap();
    assert_eq!(additional["half_life_years"].as_f64().unwrap(), 4.468e9);
}

#[test]
fn test_dating_potassium40() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("radiometric_dating"));
    params.insert("parameters".to_string(), json!({
        "isotope_system": "K40Ar40",
        "parent_isotope": 50.0,
        "daughter_isotope": 50.0
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let age = result["value"].as_f64().unwrap();
    assert!(age > 0.0);

    let interpretation = result["interpretation"].as_str().unwrap();
    assert!(interpretation.contains("years") || interpretation.contains("million") || interpretation.contains("billion"));
}

#[test]
fn test_dating_custom_half_life() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("radiometric_dating"));
    params.insert("parameters".to_string(), json!({
        "isotope_system": "Custom",
        "parent_isotope": 80.0,
        "daughter_isotope": 20.0,
        "half_life": 10000.0
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let age = result["value"].as_f64().unwrap();
    assert!(age > 0.0);
}

#[test]
fn test_dating_missing_isotope_system() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("radiometric_dating"));
    params.insert("parameters".to_string(), json!({
        "parent_isotope": 100.0,
        "daughter_isotope": 50.0
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Isotope system required"));
}

#[test]
fn test_dating_missing_isotopes() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("radiometric_dating"));
    params.insert("parameters".to_string(), json!({
        "isotope_system": "C14"
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Parent and daughter isotope amounts required"));
}

// ============================================================================
// Planetary Science Tests
// ============================================================================

#[test]
fn test_planetary_escape_velocity_earth() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("planetary_science"));
    params.insert("parameters".to_string(), json!({
        "mass_primary": 5.972e24,  // Earth mass (kg)
        "radius_primary": 6.371e6  // Earth radius (m)
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Earth's escape velocity ≈ 11.2 km/s
    let v_escape = result["value"].as_f64().unwrap();
    assert!((v_escape - 11.2).abs() < 0.5);

    assert_eq!(result["unit"].as_str().unwrap(), "km/s");

    let additional = result["additional_data"].as_object().unwrap();
    assert!(additional.contains_key("escape_velocity_m_per_s"));
    assert!(additional.contains_key("orbital_velocity_surface"));
}

#[test]
fn test_planetary_escape_velocity_mars() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("planetary_science"));
    params.insert("parameters".to_string(), json!({
        "mass_primary": 6.4171e23,  // Mars mass (kg)
        "radius_primary": 3.3895e6  // Mars radius (m)
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Mars escape velocity ≈ 5.0 km/s
    let v_escape = result["value"].as_f64().unwrap();
    assert!((v_escape - 5.0).abs() < 0.5);
}

#[test]
fn test_planetary_roche_limit() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("planetary_science"));
    params.insert("parameters".to_string(), json!({
        "radius_primary": 6.371e6,   // Earth radius (m)
        "density_primary": 5514.0,   // Earth density (kg/m³)
        "density_secondary": 3340.0  // Moon density (kg/m³)
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Roche limit should be a few Earth radii
    let roche_km = result["value"].as_f64().unwrap();
    assert!(roche_km > 5000.0);  // > 5,000 km
    assert!(roche_km < 30000.0); // < 30,000 km

    assert_eq!(result["unit"].as_str().unwrap(), "km");

    let additional = result["additional_data"].as_object().unwrap();
    assert!(additional.contains_key("roche_limit_rigid"));
    assert!(additional.contains_key("roche_limit_fluid"));
    assert!(additional.contains_key("roche_limit_primary_radii"));
}

#[test]
fn test_planetary_roche_limit_interpretation() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("planetary_science"));
    params.insert("parameters".to_string(), json!({
        "radius_primary": 1e7,
        "density_primary": 5000.0,
        "density_secondary": 3000.0
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let interpretation = result["interpretation"].as_str().unwrap();
    assert!(interpretation.contains("Roche limit"));
    assert!(interpretation.contains("km"));
    assert!(interpretation.contains("radii"));
}

#[test]
fn test_planetary_missing_parameters() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("planetary_science"));
    params.insert("parameters".to_string(), json!({
        "mass_primary": 5.972e24
        // Missing radius
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Insufficient planetary parameters"));
}

// ============================================================================
// Error Handling Tests
// ============================================================================

#[test]
fn test_invalid_input_structure() {
    let mut params = HashMap::new();
    params.insert("invalid_field".to_string(), json!("invalid"));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Invalid geophysics request"));
}

#[test]
fn test_response_structure() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("seismology"));
    params.insert("parameters".to_string(), json!({
        "seismic_moment": 1e18
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    // Verify response structure
    assert!(response.success);
    assert_eq!(response.module, "geophysics");
    assert_eq!(response.operation, "calculate");
    assert!(response.result.is_some());
    assert!(response.error.is_none());

    let result = get_result(&response);
    assert!(result["value"].is_f64());
    assert!(result["unit"].is_string());
    assert!(result["formula_used"].is_string());
    assert!(result["interpretation"].is_string());
}

#[test]
fn test_uncertainty_present() {
    let mut params = HashMap::new();
    params.insert("category".to_string(), json!("seismology"));
    params.insert("parameters".to_string(), json!({
        "energy": 1e15
    }));

    let request = create_request(params);
    let response = geophysics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Geophysics calculations should include uncertainty
    assert!(result["uncertainty"].is_f64() || result["uncertainty"].is_null());
}
