//! Tests for old API handlers that have 0% coverage
//!
//! This test file covers the following API handlers:
//! - biology.rs
//! - datetime.rs
//! - engineering.rs
//! - fluid_dynamics.rs
//! - geophysics.rs
//! - optics.rs
//! - thermodynamics.rs

use computational_engine::api::types::ComputationRequest;
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
// BIOLOGY API HANDLER TESTS
// ============================================================================

#[test]
fn test_biology_michaelis_menten() {
    let request = ComputationRequest {
        module: "biology".to_string(),
        operation: "enzyme_kinetics".to_string(),
        parameters: to_params(json!({
            "operation": "michaelis_menten",
            "parameters": {
                "vmax": 100.0,
                "km": 0.5,
                "substrate_concentration": 2.0
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Michaelis-Menten calculation should succeed");
}

#[test]
fn test_biology_hardy_weinberg() {
    let request = ComputationRequest {
        module: "biology".to_string(),
        operation: "population_genetics".to_string(),
        parameters: to_params(json!({
            "operation": "hardy_weinberg",
            "parameters": {
                "allele_frequency_p": 0.6
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Hardy-Weinberg calculation should succeed");
}

#[test]
fn test_biology_invalid_input() {
    let request = ComputationRequest {
        module: "biology".to_string(),
        operation: "enzyme_kinetics".to_string(),
        parameters: to_params(json!({
            "operation": "michaelis_menten",
            "parameters": {
                // Missing required parameters
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(!response.success, "Should fail with missing parameters");
}

// ============================================================================
// DATETIME API HANDLER TESTS
// ============================================================================

#[test]
fn test_datetime_add_interval() {
    let request = ComputationRequest {
        module: "datetime".to_string(),
        operation: "date_arithmetic".to_string(),
        parameters: to_params(json!({
            "operation": "add_interval",
            "parameters": {
                "date1": "2024-10-26",
                "days": 7
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "DateTime add interval should succeed");
}

#[test]
fn test_datetime_is_leap_year() {
    let request = ComputationRequest {
        module: "datetime".to_string(),
        operation: "calendar_info".to_string(),
        parameters: to_params(json!({
            "operation": "is_leap_year",
            "parameters": {
                "year": 2024
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Leap year check should succeed");
}

#[test]
fn test_datetime_business_days() {
    let request = ComputationRequest {
        module: "datetime".to_string(),
        operation: "date_difference".to_string(),
        parameters: to_params(json!({
            "operation": "business_days",
            "parameters": {
                "date1": "2024-10-21",
                "date2": "2024-10-25"
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Business days calculation should succeed");
}

// ============================================================================
// ENGINEERING API HANDLER TESTS
// ============================================================================

#[test]
fn test_engineering_sound_pressure_level() {
    let request = ComputationRequest {
        module: "engineering".to_string(),
        operation: "acoustics".to_string(),
        parameters: to_params(json!({
            "discipline": "acoustics",
            "parameters": {
                "pressure_rms": 1.0
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Sound pressure level should succeed");
}

#[test]
fn test_engineering_stress() {
    let request = ComputationRequest {
        module: "engineering".to_string(),
        operation: "materials".to_string(),
        parameters: to_params(json!({
            "discipline": "materials",
            "parameters": {
                "youngs_modulus": 200e9, // 200 GPa (steel)
                "strain": 0.001 // 0.1% strain
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Stress calculation should succeed");
}

#[test]
fn test_engineering_fluid_mechanics() {
    let request = ComputationRequest {
        module: "engineering".to_string(),
        operation: "fluid_mechanics".to_string(),
        parameters: to_params(json!({
            "discipline": "fluid_mechanics",
            "parameters": {
                "density": 1.225,
                "velocity": 30.0,
                "drag_coefficient": 0.3,
                "cross_sectional_area": 2.0
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Fluid mechanics calculation should succeed");
}

// ============================================================================
// FLUID DYNAMICS API HANDLER TESTS
// ============================================================================

#[test]
fn test_fluid_navier_stokes() {
    let request = ComputationRequest {
        module: "fluid_dynamics".to_string(),
        operation: "navier_stokes".to_string(),
        parameters: to_params(json!({})),
    };

    let response = computational_engine::api::process_request(&request);
    // Should return mock success
    assert!(response.success, "Navier-Stokes should return mock success");
}

#[test]
fn test_fluid_reynolds_number() {
    let request = ComputationRequest {
        module: "fluid_dynamics".to_string(),
        operation: "reynolds_number".to_string(),
        parameters: to_params(json!({})),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Reynolds number should return mock success");
}

#[test]
fn test_fluid_bernoulli() {
    let request = ComputationRequest {
        module: "fluid_dynamics".to_string(),
        operation: "bernoulli".to_string(),
        parameters: to_params(json!({})),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Bernoulli should return mock success");
}

// ============================================================================
// GEOPHYSICS API HANDLER TESTS
// ============================================================================

#[test]
fn test_geophysics_seismology() {
    let request = ComputationRequest {
        module: "geophysics".to_string(),
        operation: "seismology".to_string(),
        parameters: to_params(json!({
            "category": "seismology",
            "parameters": {
                "amplitude": 1000.0,
                "distance": 100.0,
                "seismic_moment": 1e20 // Added: seismic moment required
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Seismology calculation should succeed");
}

#[test]
fn test_geophysics_atmosphere() {
    let request = ComputationRequest {
        module: "geophysics".to_string(),
        operation: "atmosphere".to_string(),
        parameters: to_params(json!({
            "category": "atmosphere",
            "parameters": {
                "altitude": 10000.0,
                "pressure": 101325.0,
                "temperature": 288.15
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Atmosphere calculation should succeed");
}

#[test]
fn test_geophysics_radiometric_dating() {
    let request = ComputationRequest {
        module: "geophysics".to_string(),
        operation: "radiometric_dating".to_string(),
        parameters: to_params(json!({
            "category": "radiometric_dating",
            "parameters": {
                "parent_isotope": 1000.0,
                "daughter_isotope": 500.0,
                "half_life": 5730.0,
                "isotope_system": "C14"
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Radiometric dating should succeed");
}

// ============================================================================
// OPTICS API HANDLER TESTS
// ============================================================================

#[test]
fn test_optics_thin_lens() {
    let request = ComputationRequest {
        module: "optics".to_string(),
        operation: "thin_lens".to_string(),
        parameters: to_params(json!({
            "operation": "thin_lens",
            "parameters": {
                "focal_length": 0.5,
                "object_distance": 1.0
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Thin lens calculation should succeed");
}

#[test]
fn test_optics_snells_law() {
    let request = ComputationRequest {
        module: "optics".to_string(),
        operation: "snells_law".to_string(),
        parameters: to_params(json!({
            "operation": "snells_law",
            "parameters": {
                "n1": 1.0,
                "n2": 1.5,
                "theta1": 30.0 // Changed from angle1 to theta1
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Snell's law calculation should succeed");
}

#[test]
fn test_optics_diffraction_grating() {
    let request = ComputationRequest {
        module: "optics".to_string(),
        operation: "diffraction_grating".to_string(),
        parameters: to_params(json!({
            "operation": "diffraction_grating",
            "parameters": {
                "wavelength": 5.5e-7,
                "grating_spacing": 1.0e-6,
                "order": 1 // Changed from 1.0 to 1 (integer)
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Diffraction grating should succeed");
}

// ============================================================================
// THERMODYNAMICS API HANDLER TESTS
// ============================================================================

#[test]
fn test_thermodynamics_conduction() {
    let request = ComputationRequest {
        module: "thermodynamics".to_string(),
        operation: "conduction".to_string(),
        parameters: to_params(json!({
            "operation": "conduction",
            "parameters": {
                "thermal_conductivity": 0.6,
                "area": 1.0,
                "temp_hot": 320.0,
                "temp_cold": 300.0,
                "thickness": 0.1
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Conduction calculation should succeed");
}

#[test]
fn test_thermodynamics_convection() {
    let request = ComputationRequest {
        module: "thermodynamics".to_string(),
        operation: "convection".to_string(),
        parameters: to_params(json!({
            "operation": "convection",
            "parameters": {
                "heat_transfer_coefficient": 10.0,
                "area": 1.0,
                "surface_temp": 100.0,
                "fluid_temp": 20.0
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Convection calculation should succeed");
}

#[test]
fn test_thermodynamics_radiation() {
    let request = ComputationRequest {
        module: "thermodynamics".to_string(),
        operation: "radiation".to_string(),
        parameters: to_params(json!({
            "operation": "radiation",
            "parameters": {
                "emissivity": 0.9,
                "area": 1.0,
                "surface_temp_1": 400.0,
                "surface_temp_2": 300.0
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Radiation calculation should succeed");
}

#[test]
fn test_thermodynamics_entropy() {
    let request = ComputationRequest {
        module: "thermodynamics".to_string(),
        operation: "entropy".to_string(),
        parameters: to_params(json!({
            "operation": "entropy",
            "parameters": {
                "heat_transfer": 1000.0, // Changed from heat to heat_transfer
                "temperature": 300.0
            }
        })),
    };

    let response = computational_engine::api::process_request(&request);
    assert!(response.success, "Entropy calculation should succeed");
}
