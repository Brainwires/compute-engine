//! Unit tests for thermodynamics API handler
//!
//! Tests all thermodynamics operations including:
//! - Heat transfer: Conduction, Convection, Radiation
//! - Thermal resistance networks (series and parallel)
//! - Entropy calculations (Clausius, Boltzmann, Thermal)
//! - Error handling for invalid inputs and missing parameters

use crate::api::handlers::thermodynamics;
use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};
use std::collections::HashMap;

/// Helper function to create a thermodynamics request
fn create_request(parameters: HashMap<String, Value>) -> ComputationRequest {
    ComputationRequest {
        module: "thermodynamics".to_string(),
        operation: "thermodynamics".to_string(),
        parameters,
    }
}

/// Helper function to extract result value from response
fn get_result(response: &ComputationResponse) -> &Value {
    assert!(response.success, "Expected success response, got error: {:?}", response.error);
    response.result.as_ref().expect("Expected result value")
}

// ============================================================================
// Conduction Tests
// ============================================================================

#[test]
fn test_conduction_with_gradient() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("conduction"));
    params.insert("parameters".to_string(), json!({
        "thermal_conductivity": 0.6, // W/(m·K) - typical for wood
        "area": 2.0,                  // m²
        "temperature_gradient": -50.0 // K/m
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // q = -k·A·(dT/dx) = -0.6 * 2.0 * (-50.0) = 60 W
    let heat_rate = result["value"].as_f64().unwrap();
    assert!((heat_rate - 60.0).abs() < 0.01);
    assert_eq!(result["unit"].as_str().unwrap(), "W");
}

#[test]
fn test_conduction_with_temperatures() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("conduction"));
    params.insert("parameters".to_string(), json!({
        "thermal_conductivity": 0.8,  // W/(m·K)
        "area": 1.5,                   // m²
        "temp_hot": 400.0,             // K
        "temp_cold": 300.0,            // K
        "thickness": 0.1               // m
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // q = k·A·ΔT/L = 0.8 * 1.5 * (400-300) / 0.1 = 1200 W
    let heat_rate = result["value"].as_f64().unwrap();
    assert!((heat_rate - 1200.0).abs() < 0.01);

    // Check thermal resistance
    let additional = result["additional_info"].as_object().unwrap();
    let thermal_r = additional["thermal_resistance"].as_f64().unwrap();
    // R = L/(k·A) = 0.1/(0.8*1.5) ≈ 0.0833 K/W
    assert!((thermal_r - 0.0833).abs() < 0.001);
}

#[test]
fn test_conduction_missing_conductivity() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("conduction"));
    params.insert("parameters".to_string(), json!({
        "area": 1.0
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Thermal conductivity required"));
}

#[test]
fn test_conduction_zero_thickness() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("conduction"));
    params.insert("parameters".to_string(), json!({
        "thermal_conductivity": 0.6,
        "area": 1.0,
        "temp_hot": 300.0,
        "temp_cold": 200.0,
        "thickness": 0.0  // Invalid
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Thickness must be positive"));
}

// ============================================================================
// Convection Tests
// ============================================================================

#[test]
fn test_convection_natural() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("convection"));
    params.insert("parameters".to_string(), json!({
        "heat_transfer_coefficient": 5.0, // W/(m²·K) - natural convection
        "area": 2.0,                       // m²
        "surface_temp": 323.15,            // 50°C
        "fluid_temp": 293.15               // 20°C
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // q = h·A·(Ts - T∞) = 5.0 * 2.0 * (323.15 - 293.15) = 300 W
    let heat_rate = result["value"].as_f64().unwrap();
    assert!((heat_rate - 300.0).abs() < 0.01);

    let interpretation = result["interpretation"].as_str().unwrap();
    assert!(interpretation.contains("Natural convection"));
}

#[test]
fn test_convection_forced() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("convection"));
    params.insert("parameters".to_string(), json!({
        "heat_transfer_coefficient": 50.0, // W/(m²·K) - forced convection
        "area": 1.0,
        "surface_temp": 350.0,
        "fluid_temp": 300.0
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // q = 50.0 * 1.0 * (350 - 300) = 2500 W
    let heat_rate = result["value"].as_f64().unwrap();
    assert!((heat_rate - 2500.0).abs() < 0.01);

    let interpretation = result["interpretation"].as_str().unwrap();
    assert!(interpretation.contains("Forced convection"));
}

#[test]
fn test_convection_phase_change() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("convection"));
    params.insert("parameters".to_string(), json!({
        "heat_transfer_coefficient": 10000.0, // High h for phase change
        "area": 0.5,
        "surface_temp": 373.15,  // 100°C
        "fluid_temp": 373.15     // Same temp (boiling)
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let interpretation = result["interpretation"].as_str().unwrap();
    assert!(interpretation.contains("Phase change") || interpretation.contains("high velocity"));
}

#[test]
fn test_convection_missing_coefficient() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("convection"));
    params.insert("parameters".to_string(), json!({
        "area": 1.0,
        "surface_temp": 300.0,
        "fluid_temp": 280.0
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Heat transfer coefficient required"));
}

// ============================================================================
// Radiation Tests
// ============================================================================

#[test]
fn test_radiation_black_body() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("radiation"));
    params.insert("parameters".to_string(), json!({
        "emissivity": 1.0,       // Black body
        "area": 1.0,             // m²
        "surface_temp_1": 500.0, // K
        "surface_temp_2": 300.0  // K
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // q = σ·ε·A·(T₁⁴ - T₂⁴)
    // σ = 5.67e-8, ε = 1.0, A = 1.0
    // q = 5.67e-8 * 1.0 * 1.0 * (500⁴ - 300⁴) = 3119 W
    let heat_rate = result["value"].as_f64().unwrap();
    assert!(heat_rate > 3000.0 && heat_rate < 3200.0); // Relaxed range for Stefan-Boltzmann

    let interpretation = result["interpretation"].as_str().unwrap();
    assert!(interpretation.contains("Black body"));
}

#[test]
fn test_radiation_gray_body() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("radiation"));
    params.insert("parameters".to_string(), json!({
        "emissivity": 0.7,       // Gray body
        "area": 2.0,
        "surface_temp_1": 400.0,
        "surface_temp_2": 0.0    // Radiation to space
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let heat_rate = result["value"].as_f64().unwrap();
    assert!(heat_rate > 0.0);

    let interpretation = result["interpretation"].as_str().unwrap();
    assert!(interpretation.contains("Gray body"));
}

#[test]
fn test_radiation_reflective_surface() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("radiation"));
    params.insert("parameters".to_string(), json!({
        "emissivity": 0.1,       // Low emissivity (reflective)
        "area": 1.0,
        "surface_temp_1": 350.0
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let interpretation = result["interpretation"].as_str().unwrap();
    assert!(interpretation.contains("Reflective"));
}

#[test]
fn test_radiation_invalid_emissivity() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("radiation"));
    params.insert("parameters".to_string(), json!({
        "emissivity": 1.5,  // Invalid (> 1)
        "area": 1.0,
        "surface_temp_1": 300.0
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Emissivity must be between 0 and 1"));
}

#[test]
fn test_radiation_negative_temperature() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("radiation"));
    params.insert("parameters".to_string(), json!({
        "emissivity": 0.8,
        "area": 1.0,
        "surface_temp_1": -10.0  // Invalid (must be in Kelvin)
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Temperature must be in Kelvin"));
}

// ============================================================================
// Thermal Resistance Tests
// ============================================================================

#[test]
fn test_thermal_resistance_series() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("thermal_resistance"));
    params.insert("parameters".to_string(), json!({
        "resistances": [0.1, 0.2, 0.3],  // K/W
        "configuration": "series"
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let additional = result["additional_info"].as_object().unwrap();
    let total_r = additional["thermal_resistance"].as_f64().unwrap();

    // Series: R_total = 0.1 + 0.2 + 0.3 = 0.6 K/W
    assert!((total_r - 0.6).abs() < 0.001);
    assert_eq!(additional["configuration"].as_str().unwrap(), "series");
}

#[test]
fn test_thermal_resistance_parallel() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("thermal_resistance"));
    params.insert("parameters".to_string(), json!({
        "resistances": [0.2, 0.3, 0.6],  // K/W
        "configuration": "parallel"
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let additional = result["additional_info"].as_object().unwrap();
    let total_r = additional["thermal_resistance"].as_f64().unwrap();

    // Parallel: 1/R_total = 1/0.2 + 1/0.3 + 1/0.6 = 5 + 3.333 + 1.667 = 10
    // R_total = 0.1 K/W
    assert!((total_r - 0.1).abs() < 0.001);
}

#[test]
fn test_thermal_resistance_with_temperatures() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("thermal_resistance"));
    params.insert("parameters".to_string(), json!({
        "resistances": [0.5],
        "configuration": "series",
        "temp_hot": 350.0,   // K
        "temp_cold": 300.0   // K
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // q = ΔT/R = (350 - 300) / 0.5 = 100 W
    let heat_rate = result["value"].as_f64().unwrap();
    assert!((heat_rate - 100.0).abs() < 0.01);
}

#[test]
fn test_thermal_resistance_empty_array() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("thermal_resistance"));
    params.insert("parameters".to_string(), json!({
        "resistances": [],
        "configuration": "series"
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("At least one resistance value required"));
}

#[test]
fn test_thermal_resistance_invalid_configuration() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("thermal_resistance"));
    params.insert("parameters".to_string(), json!({
        "resistances": [0.1, 0.2],
        "configuration": "invalid"
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Configuration must be 'series' or 'parallel'"));
}

// ============================================================================
// Entropy Tests
// ============================================================================

#[test]
fn test_entropy_clausius_positive() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("entropy"));
    params.insert("parameters".to_string(), json!({
        "heat_transfer": 1000.0,  // J (heat added)
        "temperature": 300.0      // K
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // ΔS = Q/T = 1000/300 = 3.333 J/K
    let delta_s = result["value"].as_f64().unwrap();
    assert!((delta_s - 3.333).abs() < 0.01);
    assert_eq!(result["unit"].as_str().unwrap(), "J/K");

    let interpretation = result["interpretation"].as_str().unwrap();
    assert!(interpretation.contains("increases"));
    assert!(interpretation.contains("heat added"));
}

#[test]
fn test_entropy_clausius_negative() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("entropy"));
    params.insert("parameters".to_string(), json!({
        "heat_transfer": -500.0,  // J (heat removed)
        "temperature": 250.0      // K
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // ΔS = -500/250 = -2.0 J/K
    let delta_s = result["value"].as_f64().unwrap();
    assert!((delta_s + 2.0).abs() < 0.01);

    let interpretation = result["interpretation"].as_str().unwrap();
    assert!(interpretation.contains("decreases") || interpretation.contains("removed"));
}

#[test]
fn test_entropy_boltzmann() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("entropy"));
    params.insert("parameters".to_string(), json!({
        "num_microstates": 1e20  // Ω
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // S = k·ln(Ω) = 1.380649e-23 * ln(1e20)
    let entropy = result["value"].as_f64().unwrap();
    let k = 1.380649e-23;
    let expected = k * (1e20_f64).ln();
    assert!((entropy - expected).abs() < 1e-25);

    let formula = result["formula_used"].as_str().unwrap();
    assert!(formula.contains("Boltzmann"));
}

#[test]
fn test_entropy_thermal_heating() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("entropy"));
    params.insert("parameters".to_string(), json!({
        "mass": 2.0,              // kg
        "specific_heat": 4186.0,  // J/(kg·K) - water
        "initial_temp": 293.15,   // 20°C
        "final_temp": 373.15      // 100°C
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // ΔS = m·c·ln(T₂/T₁) = 2.0 * 4186.0 * ln(373.15/293.15)
    let delta_s = result["value"].as_f64().unwrap();
    let expected = 2.0 * 4186.0 * (373.15_f64 / 293.15_f64).ln();
    assert!((delta_s - expected).abs() < 0.1);

    let interpretation = result["interpretation"].as_str().unwrap();
    assert!(interpretation.contains("Heating"));
    assert!(interpretation.contains("increases"));
}

#[test]
fn test_entropy_thermal_cooling() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("entropy"));
    params.insert("parameters".to_string(), json!({
        "mass": 1.0,
        "specific_heat": 900.0,   // J/(kg·K) - aluminum
        "initial_temp": 400.0,    // K
        "final_temp": 300.0       // K (cooling)
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let delta_s = result["value"].as_f64().unwrap();
    // Cooling should give negative entropy change
    assert!(delta_s < 0.0);

    let interpretation = result["interpretation"].as_str().unwrap();
    assert!(interpretation.contains("Cooling"));
}

#[test]
fn test_entropy_negative_temperature() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("entropy"));
    params.insert("parameters".to_string(), json!({
        "heat_transfer": 500.0,
        "temperature": -10.0  // Invalid
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Temperature must be positive"));
}

#[test]
fn test_entropy_insufficient_parameters() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("entropy"));
    params.insert("parameters".to_string(), json!({
        "mass": 1.0  // Not enough for any calculation
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Insufficient parameters"));
}

// ============================================================================
// General Error Handling Tests
// ============================================================================

#[test]
fn test_invalid_json_structure() {
    let mut params = HashMap::new();
    params.insert("invalid_field".to_string(), json!("not_an_operation"));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Invalid thermodynamics request"));
}

#[test]
fn test_response_structure_success() {
    let mut params = HashMap::new();
    params.insert("operation".to_string(), json!("entropy"));
    params.insert("parameters".to_string(), json!({
        "heat_transfer": 100.0,
        "temperature": 300.0
    }));

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    // Verify response structure
    assert!(response.success);
    assert_eq!(response.module, "thermodynamics");
    assert_eq!(response.operation, "thermodynamics");
    assert!(response.result.is_some());
    assert!(response.error.is_none());

    // Verify result contains expected fields
    let result = response.result.as_ref().unwrap();
    assert!(result.get("value").is_some());
    assert!(result.get("unit").is_some());
    assert!(result.get("formula_used").is_some());
    assert!(result.get("interpretation").is_some());
}

#[test]
fn test_response_structure_error() {
    let params = HashMap::new(); // Empty parameters

    let request = create_request(params);
    let response = thermodynamics::handle(&request);

    // Verify error response structure
    assert!(!response.success);
    assert_eq!(response.module, "thermodynamics");
    assert_eq!(response.operation, "thermodynamics");
    assert!(response.result.is_none());
    assert!(response.error.is_some());
}
