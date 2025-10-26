//! Unit tests for electromagnetism API handler
//!
//! Tests coverage for src/api/handlers/electromagnetism.rs
//!
//! This test file covers:
//! - maxwell_equations (gauss_electric, gauss_magnetic, faraday, ampere)
//! - em_wave (vacuum, air, dielectric mediums)
//! - antenna_analysis (dipole, monopole, patch)
//! - transmission_line (coax, microstrip, stripline)
//! - waveguide (rectangular, circular, coaxial, dielectric)
//! - scattering (sphere, cylinder, plane)
//! - poynting_vector
//! - skin_effect
//! - Error handling for invalid inputs

use crate::api::handlers::electromagnetism::handle;
use crate::api::types::ComputationRequest;
use serde_json::{json, Value};
use std::collections::HashMap;

// Helper function to convert json! to HashMap
fn to_params(value: Value) -> HashMap<String, Value> {
    match value {
        Value::Object(map) => map.into_iter().collect(),
        _ => HashMap::new(),
    }
}

// ============================================================================
// MAXWELL'S EQUATIONS TESTS
// ============================================================================

#[test]
fn test_maxwell_gauss_electric() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "maxwell".to_string(),
        parameters: to_params(json!({
            "equation": "gauss_electric",
            "electric_field": [1.0, 2.0, 3.0],
            "position": [0.0, 0.0, 0.0]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Gauss electric law should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("divergence").is_some());
    assert!(result.get("field_strength").is_some());
    assert_eq!(result.get("equation").unwrap(), "Gauss's Law (Electric)");
}

#[test]
fn test_maxwell_gauss_magnetic() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "maxwell_equations".to_string(),
        parameters: to_params(json!({
            "equation": "gauss_magnetic",
            "magnetic_field": [0.5, 1.5, 2.5],
            "position": [0.0, 0.0, 0.0]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Gauss magnetic law should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("divergence").is_some());
    assert_eq!(result.get("equation").unwrap(), "Gauss's Law (Magnetic)");
}

#[test]
fn test_maxwell_faraday() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "maxwell".to_string(),
        parameters: to_params(json!({
            "equation": "faraday",
            "electric_field": [1.0, 0.0, -1.0],
            "position": [0.0, 0.0, 0.0]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Faraday's law should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("curl").is_some());
    assert_eq!(result.get("equation").unwrap(), "Faraday's Law");
}

#[test]
fn test_maxwell_ampere() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "maxwell_equations".to_string(),
        parameters: to_params(json!({
            "equation": "ampere",
            "magnetic_field": [2.0, 1.0, 0.0],
            "position": [0.0, 0.0, 0.0]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Ampere-Maxwell law should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("curl").is_some());
    assert_eq!(result.get("equation").unwrap(), "Ampere-Maxwell Law");
}

#[test]
fn test_maxwell_missing_field() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "maxwell".to_string(),
        parameters: to_params(json!({
            "equation": "gauss_electric",
            "position": [0.0, 0.0, 0.0]
            // Missing electric_field
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Should fail with missing electric field");
    assert!(response.error.is_some());
}

#[test]
fn test_maxwell_unknown_equation() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "maxwell".to_string(),
        parameters: to_params(json!({
            "equation": "unknown_equation",
            "electric_field": [1.0, 2.0, 3.0],
            "position": [0.0, 0.0, 0.0]
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Should fail with unknown equation");
    assert!(response.error.is_some());
}

// ============================================================================
// ELECTROMAGNETIC WAVE TESTS
// ============================================================================

#[test]
fn test_em_wave_vacuum() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "em_wave".to_string(),
        parameters: to_params(json!({
            "frequency": 2.4e9,  // 2.4 GHz (WiFi frequency)
            "medium": "vacuum"
        })),
    };

    let response = handle(&request);
    assert!(response.success, "EM wave in vacuum should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("wavelength").is_some());
    assert!(result.get("phase_velocity").is_some());
    assert!(result.get("impedance").is_some());
}

#[test]
fn test_em_wave_air() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "wave".to_string(),
        parameters: to_params(json!({
            "frequency": 1.0e9,  // 1 GHz
            "medium": "air"
        })),
    };

    let response = handle(&request);
    assert!(response.success, "EM wave in air should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_em_wave_dielectric() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "em_wave".to_string(),
        parameters: to_params(json!({
            "frequency": 5.0e9,  // 5 GHz
            "medium": "dielectric",
            "permittivity": 2.0e-11,
            "permeability": 1.257e-6
        })),
    };

    let response = handle(&request);
    assert!(response.success, "EM wave in dielectric should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_em_wave_unknown_medium() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "wave".to_string(),
        parameters: to_params(json!({
            "frequency": 1.0e9,
            "medium": "unknown_medium"
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Should fail with unknown medium");
    assert!(response.error.is_some());
}

// ============================================================================
// ANTENNA ANALYSIS TESTS
// ============================================================================

#[test]
fn test_antenna_dipole() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "antenna".to_string(),
        parameters: to_params(json!({
            "antenna_type": "dipole",
            "frequency": 2.4e9,  // 2.4 GHz
            "distance": 10.0,    // 10 meters
            "theta": 1.57,       // 90 degrees in radians
            "power": 1.0         // 1 Watt
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Dipole antenna analysis should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("gain").is_some());
    assert!(result.get("directivity").is_some());
    assert!(result.get("radiation_resistance").is_some());
    assert!(result.get("efficiency").is_some());
    assert!(result.get("beam_width").is_some());
}

#[test]
fn test_antenna_monopole() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "antenna".to_string(),
        parameters: to_params(json!({
            "antenna_type": "monopole",
            "frequency": 900e6,  // 900 MHz
            "distance": 5.0,
            "theta": 0.0,
            "power": 2.0
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Monopole antenna analysis should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_antenna_patch() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "antenna".to_string(),
        parameters: to_params(json!({
            "antenna_type": "patch",
            "frequency": 5.8e9,  // 5.8 GHz
            "distance": 15.0,
            "theta": 0.785,      // 45 degrees
            "power": 0.5
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Patch antenna analysis should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_antenna_unknown_type() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "antenna".to_string(),
        parameters: to_params(json!({
            "antenna_type": "unknown_type",
            "frequency": 1.0e9,
            "distance": 10.0,
            "theta": 0.0,
            "power": 1.0
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Should fail with unknown antenna type");
    assert!(response.error.is_some());
}

// ============================================================================
// TRANSMISSION LINE TESTS
// ============================================================================

#[test]
fn test_transmission_line_matched() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "transmission_line".to_string(),
        parameters: to_params(json!({
            "line_type": "coax",
            "frequency": 1.0e9,
            "length": 1.0,
            "z0": 50.0,
            "load_impedance": 50.0  // Matched load
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Matched transmission line should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("vswr").is_some());
    assert!(result.get("reflection_coefficient").is_some());
    assert!(result.get("return_loss").is_some());
}

#[test]
fn test_transmission_line_mismatched() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "transmission_line".to_string(),
        parameters: to_params(json!({
            "line_type": "microstrip",
            "frequency": 2.4e9,
            "length": 0.5,
            "z0": 50.0,
            "load_impedance": 75.0,  // Mismatched load
            "attenuation": 0.1
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Mismatched transmission line should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_transmission_line_default_load() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "transmission_line".to_string(),
        parameters: to_params(json!({
            "line_type": "stripline",
            "frequency": 5.0e9,
            "length": 0.25,
            "z0": 50.0
            // No load_impedance specified (should default to z0)
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Transmission line with default load should succeed");
    assert!(response.result.is_some());
}

// ============================================================================
// WAVEGUIDE TESTS
// ============================================================================

#[test]
fn test_waveguide_rectangular() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "waveguide".to_string(),
        parameters: to_params(json!({
            "guide_type": "rectangular",
            "frequency": 10.0e9,  // 10 GHz (well above cutoff)
            "width": 0.02286,     // WR-90 standard
            "height": 0.01016,
            "mode": "TE10"
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Rectangular waveguide should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("cutoff_frequency").is_some());
    assert!(result.get("guide_wavelength").is_some());
    assert!(result.get("phase_velocity").is_some());
    assert!(result.get("group_velocity").is_some());
}

#[test]
fn test_waveguide_circular() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "waveguide".to_string(),
        parameters: to_params(json!({
            "guide_type": "circular",
            "frequency": 12.0e9,
            "width": 0.015,  // radius
            "mode": "TE11"
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Circular waveguide should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_waveguide_coaxial() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "waveguide".to_string(),
        parameters: to_params(json!({
            "guide_type": "coaxial",
            "frequency": 5.0e9,
            "width": 0.001,   // inner radius
            "height": 0.003,  // outer radius
            "mode": "TEM"
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Coaxial waveguide should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_waveguide_dielectric() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "waveguide".to_string(),
        parameters: to_params(json!({
            "guide_type": "dielectric",
            "frequency": 30.0e9,  // 30 GHz (well above cutoff)
            "width": 0.005,       // slab thickness (5mm)
            "height": 1.5,        // refractive index of core
            "mode": "TE0"
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Dielectric waveguide should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_waveguide_below_cutoff() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "waveguide".to_string(),
        parameters: to_params(json!({
            "guide_type": "rectangular",
            "frequency": 1.0e9,  // Too low for typical waveguide
            "width": 0.02286,
            "height": 0.01016,
            "mode": "TE10"
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Should fail below cutoff frequency");
    assert!(response.error.is_some());
}

#[test]
fn test_waveguide_unknown_type() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "waveguide".to_string(),
        parameters: to_params(json!({
            "guide_type": "unknown_type",
            "frequency": 10.0e9,
            "width": 0.02,
            "mode": "TE10"
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Should fail with unknown waveguide type");
    assert!(response.error.is_some());
}

// ============================================================================
// SCATTERING TESTS
// ============================================================================

#[test]
fn test_scattering_rayleigh_regime() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "scattering".to_string(),
        parameters: to_params(json!({
            "scatterer_type": "sphere",
            "radius": 0.001,      // 1 mm
            "wavelength": 0.1,    // 10 cm (ka << 1)
            "incident_angle": 0.0,
            "polarization": "TE"
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Rayleigh scattering should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("radar_cross_section").is_some());
    assert!(result.get("scattering_pattern").is_some());
}

#[test]
fn test_scattering_geometric_regime() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "scattering".to_string(),
        parameters: to_params(json!({
            "scatterer_type": "sphere",
            "radius": 1.0,        // 1 m
            "wavelength": 0.01,   // 1 cm (ka >> 1)
            "incident_angle": 1.57,
            "polarization": "TM"
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Geometric optics scattering should succeed");
    assert!(response.result.is_some());
}

// ============================================================================
// POYNTING VECTOR TESTS
// ============================================================================

#[test]
fn test_poynting_vector() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "poynting".to_string(),
        parameters: to_params(json!({
            "electric_field": [1.0, 0.0, 0.0],
            "magnetic_field": [0.0, 1.0, 0.0]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Poynting vector calculation should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("poynting_vector").is_some());
    assert!(result.get("power_density").is_some());
    assert!(result.get("intensity").is_some());
}

#[test]
fn test_poynting_vector_alias() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "poynting_vector".to_string(),
        parameters: to_params(json!({
            "electric_field": [2.0, 1.0, 0.0],
            "magnetic_field": [0.0, 0.0, 3.0]
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Poynting vector (alias) should succeed");
    assert!(response.result.is_some());
}

#[test]
fn test_poynting_vector_invalid_dimensions() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "poynting".to_string(),
        parameters: to_params(json!({
            "electric_field": [1.0, 0.0],  // Only 2D
            "magnetic_field": [0.0, 1.0]
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Should fail with 2D vectors");
    assert!(response.error.is_some());
}

// ============================================================================
// SKIN EFFECT TESTS
// ============================================================================

#[test]
fn test_skin_effect_copper() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "skin_effect".to_string(),
        parameters: to_params(json!({
            "frequency": 1.0e6,        // 1 MHz
            "conductivity": 5.96e7     // Copper conductivity (S/m)
        })),
    };

    let response = handle(&request);
    assert!(response.success, "Skin effect for copper should succeed");
    assert!(response.result.is_some());

    let result = response.result.unwrap();
    assert!(result.get("skin_depth").is_some());
    assert!(result.get("surface_resistance").is_some());
}

#[test]
fn test_skin_effect_high_frequency() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "skin_effect".to_string(),
        parameters: to_params(json!({
            "frequency": 1.0e9,        // 1 GHz
            "conductivity": 3.5e7,     // Aluminum conductivity
            "permeability": 1.257e-6   // Free space permeability
        })),
    };

    let response = handle(&request);
    assert!(response.success, "High frequency skin effect should succeed");
    assert!(response.result.is_some());
}

// ============================================================================
// ERROR HANDLING TESTS
// ============================================================================

#[test]
fn test_unknown_operation() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "unknown_operation".to_string(),
        parameters: to_params(json!({})),
    };

    let response = handle(&request);
    assert!(!response.success, "Should fail with unknown operation");
    assert!(response.error.is_some());
    assert!(response.error.unwrap().contains("Unknown electromagnetism operation"));
}

#[test]
fn test_invalid_json_parameters() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "em_wave".to_string(),
        parameters: to_params(json!({
            "frequency": "not_a_number",  // Invalid type
            "medium": "vacuum"
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Should fail with invalid parameter types");
    assert!(response.error.is_some());
}

#[test]
fn test_missing_required_parameters() {
    let request = ComputationRequest {
        module: "electromagnetism".to_string(),
        operation: "antenna".to_string(),
        parameters: to_params(json!({
            "antenna_type": "dipole"
            // Missing frequency, distance, theta, power
        })),
    };

    let response = handle(&request);
    assert!(!response.success, "Should fail with missing required parameters");
    assert!(response.error.is_some());
}
