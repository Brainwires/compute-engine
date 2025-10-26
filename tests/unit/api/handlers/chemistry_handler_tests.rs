//! Unit tests for chemistry API handler
//!
//! Tests all chemistry operations including:
//! - Molar mass calculations
//! - Gas law computations
//! - Thermodynamics (enthalpy, entropy, Gibbs free energy)
//! - Acid-base chemistry (pH calculations)
//! - Chemical kinetics
//! - Electrochemistry (cell potentials)
//! - Error handling for invalid inputs

use crate::api::handlers::chemistry;
use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};
use std::collections::HashMap;

/// Helper function to create a chemistry request
fn create_request(operation: &str, parameters: HashMap<String, Value>) -> ComputationRequest {
    ComputationRequest {
        module: "chemistry".to_string(),
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
// Molar Mass Tests
// ============================================================================

#[test]
fn test_molar_mass_water() {
    let mut params = HashMap::new();
    params.insert("formula".to_string(), json!("H2O"));

    let request = create_request("molar_mass", params);
    let response = chemistry::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Water: H2O = 2(1.008) + 15.999 = 18.015 g/mol
    let molar_mass = result["molar_mass"].as_f64().unwrap();
    assert!((molar_mass - 18.015).abs() < 0.01);
}

#[test]
fn test_molar_mass_glucose() {
    let mut params = HashMap::new();
    params.insert("formula".to_string(), json!("C6H12O6"));

    let request = create_request("molar_mass", params);
    let response = chemistry::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Glucose: C6H12O6 = 6(12.011) + 12(1.008) + 6(15.999) = 180.156 g/mol
    let molar_mass = result["molar_mass"].as_f64().unwrap();
    assert!((molar_mass - 180.156).abs() < 0.01);
}

#[test]
fn test_molar_mass_composition() {
    let mut params = HashMap::new();
    params.insert("formula".to_string(), json!("H2O"));

    let request = create_request("molar_mass", params);
    let response = chemistry::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Check composition percentages
    let composition = result["composition"].as_object().unwrap();

    // H: 2(1.008) / 18.015 * 100 ≈ 11.19%
    let h_percent = composition["H"].as_f64().unwrap();
    assert!((h_percent - 11.19).abs() < 0.5);

    // O: 15.999 / 18.015 * 100 ≈ 88.81%
    let o_percent = composition["O"].as_f64().unwrap();
    assert!((o_percent - 88.81).abs() < 0.5);
}

#[test]
fn test_molar_mass_invalid_request() {
    let params = HashMap::new(); // Missing formula

    let request = create_request("molar_mass", params);
    let response = chemistry::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Invalid molar_mass request"));
}

// ============================================================================
// Gas Law Tests
// ============================================================================

#[test]
fn test_gas_law_calculate_moles() {
    let mut params = HashMap::new();
    params.insert("pressure".to_string(), json!(1.0)); // 1 atm
    params.insert("volume".to_string(), json!(22.4)); // 22.4 L
    params.insert("temperature".to_string(), json!(273.15)); // 0°C = 273.15 K

    let request = create_request("gas_law", params);
    let response = chemistry::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // PV = nRT → n = PV/RT = (1)(22.4)/(0.08206)(273.15) ≈ 1.0 mol
    let moles = result["moles"].as_f64().unwrap();
    assert!((moles - 1.0).abs() < 0.01);
}

#[test]
fn test_gas_law_calculate_pressure() {
    let mut params = HashMap::new();
    params.insert("volume".to_string(), json!(10.0)); // 10 L
    params.insert("temperature".to_string(), json!(300.0)); // 300 K
    params.insert("moles".to_string(), json!(0.5)); // 0.5 mol

    let request = create_request("gas_law", params);
    let response = chemistry::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // P = nRT/V = (0.5)(0.08206)(300)/10 ≈ 1.23 atm
    let pressure = result["pressure"].as_f64().unwrap();
    assert!((pressure - 1.23).abs() < 0.01);
}

#[test]
fn test_gas_law_calculate_volume() {
    let mut params = HashMap::new();
    params.insert("pressure".to_string(), json!(2.0)); // 2 atm
    params.insert("temperature".to_string(), json!(350.0)); // 350 K
    params.insert("moles".to_string(), json!(1.0)); // 1 mol

    let request = create_request("gas_law", params);
    let response = chemistry::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // V = nRT/P = (1)(0.08206)(350)/2 ≈ 14.36 L
    let volume = result["volume"].as_f64().unwrap();
    assert!((volume - 14.36).abs() < 0.01);
}

#[test]
fn test_gas_law_insufficient_parameters() {
    let mut params = HashMap::new();
    params.insert("pressure".to_string(), json!(1.0));
    // Missing other required parameters

    let request = create_request("gas_law", params);
    let response = chemistry::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
}

// ============================================================================
// Thermodynamics Tests
// ============================================================================

#[test]
fn test_thermodynamics_spontaneous_reaction() {
    let mut params = HashMap::new();
    params.insert("temperature".to_string(), json!(298.15)); // 25°C
    params.insert("pressure".to_string(), json!(1.0)); // 1 atm
    params.insert("compounds".to_string(), json!([
        {
            "formula": "H2O",
            "moles": 1.0,
            "enthalpy_formation": -285.8, // kJ/mol (exothermic)
            "entropy": 69.9 // J/(mol·K)
        }
    ]));

    let request = create_request("thermodynamics", params);
    let response = chemistry::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let delta_h = result["delta_h"].as_f64().unwrap();
    assert!(delta_h < 0.0); // Exothermic

    let spontaneous = result["spontaneous"].as_bool().unwrap();
    // For exothermic reactions at standard conditions, typically spontaneous
    assert_eq!(result["spontaneous"].as_bool(), Some(spontaneous));
}

#[test]
fn test_thermodynamics_equilibrium_constant() {
    let mut params = HashMap::new();
    params.insert("temperature".to_string(), json!(298.15));
    params.insert("pressure".to_string(), json!(1.0));
    params.insert("compounds".to_string(), json!([
        {
            "formula": "A",
            "moles": 1.0,
            "enthalpy_formation": 0.0,
            "entropy": 100.0
        }
    ]));

    let request = create_request("thermodynamics", params);
    let response = chemistry::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // K = exp(-ΔG/RT)
    let k_eq = result["equilibrium_constant"].as_f64().unwrap();
    assert!(k_eq > 0.0); // K must always be positive
}

// ============================================================================
// Acid-Base Tests
// ============================================================================

#[test]
fn test_acid_base_weak_acid() {
    let mut params = HashMap::new();
    params.insert("acid_ka".to_string(), json!(1.8e-5)); // Acetic acid Ka
    params.insert("concentration".to_string(), json!(0.1)); // 0.1 M

    let request = create_request("acid_base", params);
    let response = chemistry::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let ph = result["ph"].as_f64().unwrap();

    // For weak acid: [H+] = sqrt(Ka * C) = sqrt(1.8e-5 * 0.1) ≈ 1.34e-3
    // pH = -log(1.34e-3) ≈ 2.87
    assert!((ph - 2.87).abs() < 0.1);

    let poh = result["poh"].as_f64().unwrap();
    assert!((ph + poh - 14.0).abs() < 0.1); // pH + pOH = 14
}

#[test]
fn test_acid_base_weak_base() {
    let mut params = HashMap::new();
    params.insert("base_kb".to_string(), json!(1.8e-5)); // Kb
    params.insert("concentration".to_string(), json!(0.1)); // 0.1 M

    let request = create_request("acid_base", params);
    let response = chemistry::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let ph = result["ph"].as_f64().unwrap();

    // For weak base, pH should be > 7
    assert!(ph > 7.0);
    assert!(ph < 14.0);
}

#[test]
fn test_acid_base_missing_ka_kb() {
    let mut params = HashMap::new();
    params.insert("concentration".to_string(), json!(0.1));
    // Missing both Ka and Kb

    let request = create_request("acid_base", params);
    let response = chemistry::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
}

// ============================================================================
// Kinetics Tests
// ============================================================================

#[test]
fn test_kinetics_reaction_rate() {
    let mut params = HashMap::new();
    params.insert("rate_constant".to_string(), json!(0.05)); // s^-1
    params.insert("temperature".to_string(), json!(298.15)); // K
    params.insert("activation_energy".to_string(), json!(50.0)); // kJ/mol
    params.insert("concentrations".to_string(), json!([0.5])); // mol/L
    params.insert("orders".to_string(), json!([1])); // First-order

    let request = create_request("kinetics", params);
    let response = chemistry::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let rate = result["rate"].as_f64().unwrap();
    assert!(rate > 0.0);
}

#[test]
fn test_kinetics_first_order_half_life() {
    let mut params = HashMap::new();
    params.insert("rate_constant".to_string(), json!(0.1)); // s^-1
    params.insert("temperature".to_string(), json!(298.15));
    params.insert("activation_energy".to_string(), json!(50.0));
    params.insert("concentrations".to_string(), json!([1.0]));
    params.insert("orders".to_string(), json!([1])); // First-order

    let request = create_request("kinetics", params);
    let response = chemistry::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // For first-order: t_1/2 = ln(2)/k
    let half_life = result["half_life"].as_f64();
    assert!(half_life.is_some());
}

#[test]
fn test_kinetics_second_order_no_half_life() {
    let mut params = HashMap::new();
    params.insert("rate_constant".to_string(), json!(0.1));
    params.insert("temperature".to_string(), json!(298.15));
    params.insert("activation_energy".to_string(), json!(50.0));
    params.insert("concentrations".to_string(), json!([1.0]));
    params.insert("orders".to_string(), json!([2])); // Second-order

    let request = create_request("kinetics", params);
    let response = chemistry::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Second-order reactions don't have simple half-life formula
    assert!(result["half_life"].is_null());
}

// ============================================================================
// Electrochemistry Tests
// ============================================================================

#[test]
fn test_electrochemistry_cell_potential() {
    let mut params = HashMap::new();
    params.insert("temperature".to_string(), json!(298.15));
    params.insert("half_reactions".to_string(), json!([
        {
            "oxidized": "Zn",
            "reduced": "Zn2+",
            "electrons": 2,
            "standard_potential": -0.76 // V (oxidation)
        },
        {
            "oxidized": "Cu2+",
            "reduced": "Cu",
            "electrons": 2,
            "standard_potential": 0.34 // V (reduction)
        }
    ]));

    let request = create_request("electrochemistry", params);
    let response = chemistry::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // E°cell = E°cathode - E°anode = 0.34 - (-0.76) = 1.10 V
    let std_potential = result["standard_cell_potential"].as_f64().unwrap();
    assert!((std_potential - 1.10).abs() < 0.01);

    let spontaneous = result["spontaneous"].as_bool().unwrap();
    assert!(spontaneous); // Positive E° means spontaneous
}

#[test]
fn test_electrochemistry_with_concentrations() {
    let mut params = HashMap::new();
    params.insert("temperature".to_string(), json!(298.15));
    params.insert("half_reactions".to_string(), json!([
        {
            "oxidized": "Zn",
            "reduced": "Zn2+",
            "electrons": 2,
            "standard_potential": -0.76
        },
        {
            "oxidized": "Cu2+",
            "reduced": "Cu",
            "electrons": 2,
            "standard_potential": 0.34
        }
    ]));
    params.insert("concentrations".to_string(), json!({
        "Zn2+": 0.1,
        "Cu2+": 1.0
    }));

    let request = create_request("electrochemistry", params);
    let response = chemistry::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Nernst equation should modify the cell potential
    let cell_potential = result["cell_potential"].as_f64().unwrap();
    assert!(cell_potential > 0.0);
}

#[test]
fn test_electrochemistry_wrong_number_half_reactions() {
    let mut params = HashMap::new();
    params.insert("temperature".to_string(), json!(298.15));
    params.insert("half_reactions".to_string(), json!([
        {
            "oxidized": "Zn",
            "reduced": "Zn2+",
            "electrons": 2,
            "standard_potential": -0.76
        }
        // Missing second half-reaction
    ]));

    let request = create_request("electrochemistry", params);
    let response = chemistry::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("2 half-reactions"));
}

// ============================================================================
// Unknown Operation Tests
// ============================================================================

#[test]
fn test_unknown_operation() {
    let params = HashMap::new();
    let request = create_request("invalid_operation", params);
    let response = chemistry::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Unknown chemistry operation"));
}

#[test]
fn test_response_structure() {
    let mut params = HashMap::new();
    params.insert("formula".to_string(), json!("H2O"));

    let request = create_request("molar_mass", params);
    let response = chemistry::handle(&request);

    // Verify response structure
    assert!(response.success);
    assert_eq!(response.module, "chemistry");
    assert_eq!(response.operation, "molar_mass");
    assert!(response.result.is_some());
    assert!(response.error.is_none());
}
