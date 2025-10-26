//! Unit tests for engineering API handler
//!
//! Tests all engineering operations including:
//! - Acoustics (SPL, Doppler effect, reverberation time)
//! - Materials science (stress-strain, fracture mechanics)
//! - Fluid mechanics (Bernoulli, Poiseuille, drag force)
//! - Control theory (PID control calculations)
//! - Error handling for invalid inputs

use crate::api::handlers::engineering;
use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};
use std::collections::HashMap;

/// Helper function to create an engineering request
fn create_request(operation: &str, parameters: HashMap<String, Value>) -> ComputationRequest {
    ComputationRequest {
        module: "engineering".to_string(),
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
// Acoustics Tests - Sound Pressure Level (SPL)
// ============================================================================

#[test]
fn test_acoustics_spl_whisper_quiet() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("acoustics"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "pressure_rms": 0.0002 // 0.2 mPa, should be ~20 dB SPL
        }),
    );

    let request = create_request("acoustics", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let spl = result["value"].as_f64().unwrap();
    assert!((spl - 20.0).abs() < 5.0); // ~20 dB SPL

    let classification = result["classification"].as_str().unwrap();
    assert!(classification.contains("Whisper") || classification.contains("Quiet"));
}

#[test]
fn test_acoustics_spl_moderate() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("acoustics"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "pressure_rms": 0.02 // 20 mPa, should be ~60 dB SPL
        }),
    );

    let request = create_request("acoustics", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let spl = result["value"].as_f64().unwrap();
    assert!((spl - 60.0).abs() < 5.0); // ~60 dB SPL
    assert_eq!(result["unit"].as_str().unwrap(), "dB SPL");
}

#[test]
fn test_acoustics_spl_dangerous() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("acoustics"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "pressure_rms": 20.0 // High pressure, dangerous levels
        }),
    );

    let request = create_request("acoustics", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let spl = result["value"].as_f64().unwrap();
    assert!(spl > 100.0); // Should be > 100 dB

    let classification = result["classification"].as_str().unwrap();
    assert!(classification.contains("damage") || classification.contains("Dangerous") || classification.contains("Pain"));
}

// ============================================================================
// Acoustics Tests - Doppler Effect
// ============================================================================

#[test]
fn test_acoustics_doppler_approaching() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("acoustics"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "frequency": 1000.0,           // 1 kHz
            "sound_speed": 343.0,           // m/s at 20°C
            "velocity_source": 30.0,        // Source approaching at 30 m/s
            "velocity_observer": 0.0        // Stationary observer
        }),
    );

    let request = create_request("acoustics", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let f_observed = result["value"].as_f64().unwrap();

    // Approaching source should have higher frequency
    assert!(f_observed > 1000.0);

    // f' = f * (v + v_obs) / (v - v_src) = 1000 * 343 / (343 - 30) ≈ 1096 Hz
    assert!((f_observed - 1096.0).abs() < 5.0);

    assert_eq!(result["unit"].as_str().unwrap(), "Hz");
}

#[test]
fn test_acoustics_doppler_receding() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("acoustics"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "frequency": 1000.0,
            "sound_speed": 343.0,
            "velocity_source": -20.0,       // Source receding (negative velocity)
            "velocity_observer": 0.0
        }),
    );

    let request = create_request("acoustics", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let f_observed = result["value"].as_f64().unwrap();

    // Receding source should have lower frequency
    assert!(f_observed < 1000.0);
}

// ============================================================================
// Acoustics Tests - Reverberation Time
// ============================================================================

#[test]
fn test_acoustics_reverberation_dry_room() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("acoustics"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "room_volume": 100.0,              // 100 m³
            "absorption_coefficient": 0.5      // High absorption (dry)
        }),
    );

    let request = create_request("acoustics", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let rt60 = result["value"].as_f64().unwrap();

    // High absorption should give short RT60
    assert!(rt60 < 1.0);
    assert_eq!(result["unit"].as_str().unwrap(), "seconds");

    let classification = result["classification"].as_str().unwrap();
    assert!(classification.contains("dry") || classification.contains("Dry"));
}

#[test]
fn test_acoustics_reverberation_cathedral() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("acoustics"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "room_volume": 10000.0,            // Large space
            "absorption_coefficient": 0.05     // Low absorption (reflective)
        }),
    );

    let request = create_request("acoustics", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let rt60 = result["value"].as_f64().unwrap();

    // Low absorption and large volume should give long RT60
    assert!(rt60 > 2.0);
}

// ============================================================================
// Materials Science Tests - Hooke's Law (Stress-Strain)
// ============================================================================

#[test]
fn test_materials_stress_from_strain() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("materials"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "youngs_modulus": 200e9,  // 200 GPa (steel)
            "strain": 0.001           // 0.1% strain
        }),
    );

    let request = create_request("materials", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // σ = E·ε = 200e9 * 0.001 = 200e6 Pa = 200 MPa
    let stress_mpa = result["value"].as_f64().unwrap();
    assert!((stress_mpa - 200.0).abs() < 5.0);
    assert_eq!(result["unit"].as_str().unwrap(), "MPa");

    let classification = result["classification"].as_str().unwrap();
    assert!(classification.contains("Elastic"));
}

#[test]
fn test_materials_strain_from_stress() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("materials"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "youngs_modulus": 200e9,  // 200 GPa
            "stress": 100e6           // 100 MPa
        }),
    );

    let request = create_request("materials", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // ε = σ/E = 100e6 / 200e9 = 0.0005 = 0.05%
    let strain_percent = result["value"].as_f64().unwrap();
    assert!((strain_percent - 0.05).abs() < 0.01);
    assert_eq!(result["unit"].as_str().unwrap(), "% strain");
}

#[test]
fn test_materials_plastic_deformation() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("materials"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "youngs_modulus": 200e9,
            "strain": 0.005,              // High strain
            "yield_strength": 250e6       // 250 MPa yield strength
        }),
    );

    let request = create_request("materials", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // σ = 200e9 * 0.005 = 1000 MPa, exceeds yield strength
    let classification = result["classification"].as_str().unwrap();
    assert!(classification.contains("Plastic") || classification.contains("deformation"));
}

// ============================================================================
// Materials Science Tests - Fracture Mechanics
// ============================================================================

#[test]
fn test_materials_fracture_mechanics_safe() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("materials"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "stress": 100e6,             // 100 MPa
            "crack_length": 0.001,       // 1 mm
            "geometry_factor": 1.0       // Simple geometry
        }),
    );

    let request = create_request("materials", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let k_factor = result["value"].as_f64().unwrap();
    assert!(k_factor > 0.0);
    assert_eq!(result["unit"].as_str().unwrap(), "MPa·√m");

    let classification = result["classification"].as_str().unwrap();
    assert!(classification.contains("Safe"));
}

#[test]
fn test_materials_fracture_mechanics_critical() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("materials"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "stress": 500e6,             // Very high stress
            "crack_length": 0.05,        // 50 mm large crack
            "geometry_factor": 1.2       // Complex geometry
        }),
    );

    let request = create_request("materials", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let k_factor = result["value"].as_f64().unwrap();
    assert!(k_factor > 30.0); // Should be high

    let classification = result["classification"].as_str().unwrap();
    assert!(classification.contains("Critical") || classification.contains("Monitor"));
}

// ============================================================================
// Fluid Mechanics Tests - Bernoulli's Equation
// ============================================================================

#[test]
fn test_fluid_mechanics_bernoulli() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("fluid_mechanics"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "density": 1000.0,           // Water
            "pressure_1": 200000.0,      // 200 kPa
            "pressure_2": 100000.0,      // 100 kPa
            "velocity": 5.0,             // 5 m/s inlet
            "height_1": 10.0,            // 10 m elevation
            "height_2": 5.0              // 5 m exit elevation
        }),
    );

    let request = create_request("fluid_mechanics", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let v2 = result["value"].as_f64().unwrap();
    assert!(v2 > 0.0);
    assert_eq!(result["unit"].as_str().unwrap(), "m/s");
}

// ============================================================================
// Fluid Mechanics Tests - Poiseuille's Law
// ============================================================================

#[test]
fn test_fluid_mechanics_poiseuille_laminar() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("fluid_mechanics"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "pressure_1": 10.0,          // Small pressure drop for laminar flow
            "radius": 0.002,             // 2 mm radius (small pipe)
            "length": 1.0,               // 1 m pipe
            "viscosity": 0.001,          // Water viscosity (Pa·s)
            "density": 1000.0            // Water density (kg/m³)
        }),
    );

    let request = create_request("fluid_mechanics", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let flow_rate = result["value"].as_f64().unwrap();
    assert!(flow_rate > 0.0);
    assert_eq!(result["unit"].as_str().unwrap(), "L/s");

    let classification = result["classification"].as_str().unwrap();
    // The classification should indicate laminar flow (Re < 2300)
    assert!(classification.contains("Laminar"));
}

#[test]
fn test_fluid_mechanics_poiseuille_turbulent() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("fluid_mechanics"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "pressure_1": 100000.0,      // Large pressure drop
            "radius": 0.05,              // Large radius
            "length": 1.0,
            "viscosity": 0.001,
            "density": 1000.0
        }),
    );

    let request = create_request("fluid_mechanics", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // High Re number should indicate turbulent flow
    let classification = result["classification"].as_str().unwrap();
    assert!(classification.contains("Turbulent") || classification.contains("invalid"));
}

// ============================================================================
// Fluid Mechanics Tests - Drag Force
// ============================================================================

#[test]
fn test_fluid_mechanics_drag_force() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("fluid_mechanics"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "density": 1.2,                  // Air density
            "velocity": 30.0,                // 30 m/s (~108 km/h)
            "drag_coefficient": 0.3,         // Streamlined car
            "cross_sectional_area": 2.0      // 2 m²
        }),
    );

    let request = create_request("fluid_mechanics", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // F_D = 0.5 * 1.2 * 30² * 0.3 * 2.0 = 324 N
    let drag_force = result["value"].as_f64().unwrap();
    assert!((drag_force - 324.0).abs() < 10.0);
    assert_eq!(result["unit"].as_str().unwrap(), "N");
}

// ============================================================================
// Control Theory Tests - PID Control
// ============================================================================

#[test]
fn test_control_theory_basic_pid() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("control_theory"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "setpoint": 100.0,
            "process_variable": 80.0,
            "kp": 2.0,
            "ki": 0.5,
            "kd": 0.1
        }),
    );

    let request = create_request("control_theory", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Error = 100 - 80 = 20
    // P-term = 2.0 * 20 = 40
    let p_term = result["value"].as_f64().unwrap();
    assert!((p_term - 40.0).abs() < 0.1);
}

#[test]
fn test_control_theory_ziegler_nichols() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("control_theory"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "setpoint": 100.0,
            "process_variable": 90.0,
            "time_constant": 10.0       // Enables Ziegler-Nichols tuning
        }),
    );

    let request = create_request("control_theory", params);
    let response = engineering::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Should return suggested Kp
    let kp_suggested = result["value"].as_f64().unwrap();
    assert!(kp_suggested > 0.0);
    assert_eq!(result["unit"].as_str().unwrap(), "gain");

    let classification = result["classification"].as_str().unwrap();
    assert!(classification.contains("PI"));
}

#[test]
fn test_control_theory_missing_setpoint() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("control_theory"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "process_variable": 80.0  // Missing setpoint
        }),
    );

    let request = create_request("control_theory", params);
    let response = engineering::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Setpoint required"));
}

// ============================================================================
// Error Handling Tests
// ============================================================================

#[test]
fn test_invalid_discipline() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("invalid_discipline"),
    );
    params.insert(
        "parameters".to_string(),
        json!({}),
    );

    let request = create_request("engineering", params);
    let response = engineering::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
}

#[test]
fn test_missing_parameters() {
    let params = HashMap::new(); // Empty parameters

    let request = create_request("engineering", params);
    let response = engineering::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Invalid engineering request"));
}

#[test]
fn test_insufficient_acoustics_parameters() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("acoustics"),
    );
    params.insert(
        "parameters".to_string(),
        json!({}), // No acoustic parameters
    );

    let request = create_request("engineering", params);
    let response = engineering::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Insufficient acoustic parameters"));
}

#[test]
fn test_insufficient_materials_parameters() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("materials"),
    );
    params.insert(
        "parameters".to_string(),
        json!({}),
    );

    let request = create_request("engineering", params);
    let response = engineering::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Insufficient materials parameters"));
}

#[test]
fn test_insufficient_fluid_mechanics_parameters() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("fluid_mechanics"),
    );
    params.insert(
        "parameters".to_string(),
        json!({}),
    );

    let request = create_request("engineering", params);
    let response = engineering::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Insufficient fluid mechanics parameters"));
}

#[test]
fn test_response_structure() {
    let mut params = HashMap::new();
    params.insert(
        "discipline".to_string(),
        json!("acoustics"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "pressure_rms": 0.02
        }),
    );

    let request = create_request("engineering", params);
    let response = engineering::handle(&request);

    // Verify response structure
    assert!(response.success);
    assert_eq!(response.module, "engineering");
    assert_eq!(response.operation, "engineering");
    assert!(response.result.is_some());
    assert!(response.error.is_none());

    // Verify result contains required fields
    let result = get_result(&response);
    assert!(result["value"].is_f64());
    assert!(result["unit"].is_string());
    assert!(result["formula_used"].is_string());
    assert!(result["interpretation"].is_string());
}
