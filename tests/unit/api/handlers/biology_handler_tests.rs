//! Unit tests for biology API handler
//!
//! Tests all biology operations including:
//! - Enzyme kinetics (Michaelis-Menten, Lineweaver-Burk)
//! - Pharmacokinetics (one-compartment model)
//! - Population genetics (Hardy-Weinberg equilibrium)
//! - Membrane potentials (Goldman-Hodgkin-Katz equation)
//! - Allometric scaling (metabolic rate, surface area, lifespan)
//! - Error handling for invalid inputs

use crate::api::handlers::biology;
use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};
use std::collections::HashMap;

/// Helper function to create a biology request
fn create_request(operation: &str, parameters: HashMap<String, Value>) -> ComputationRequest {
    ComputationRequest {
        module: "biology".to_string(),
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
// Michaelis-Menten Enzyme Kinetics Tests
// ============================================================================

#[test]
fn test_michaelis_menten_basic() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("michaelis_menten"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "vmax": 100.0,
            "km": 10.0,
            "substrate_concentration": 10.0
        }),
    );

    let request = create_request("michaelis_menten", params);
    let response = biology::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // At [S] = Km, v = Vmax/2 = 50.0
    let velocity = result["value"].as_f64().unwrap();
    assert!((velocity - 50.0).abs() < 0.1);
    assert_eq!(result["unit"].as_str().unwrap(), "μmol/(min·mg enzyme)");
}

#[test]
fn test_michaelis_menten_saturation() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("michaelis_menten"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "vmax": 100.0,
            "km": 10.0,
            "substrate_concentration": 150.0 // 15x Km (> 10x)
        }),
    );

    let request = create_request("michaelis_menten", params);
    let response = biology::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // At [S] >> Km, v approaches Vmax
    let velocity = result["value"].as_f64().unwrap();
    assert!(velocity > 90.0); // Should be close to 100

    let interpretation = result["interpretation"].as_str().unwrap();
    assert!(interpretation.contains("saturated"));
}

#[test]
fn test_michaelis_menten_low_substrate() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("michaelis_menten"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "vmax": 100.0,
            "km": 10.0,
            "substrate_concentration": 1.0 // << Km
        }),
    );

    let request = create_request("michaelis_menten", params);
    let response = biology::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let velocity = result["value"].as_f64().unwrap();
    assert!(velocity < 20.0); // Should be much less than Vmax

    let interpretation = result["interpretation"].as_str().unwrap();
    assert!(interpretation.contains("not saturated") || interpretation.contains("Below Km"));
}

#[test]
fn test_michaelis_menten_missing_vmax() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("michaelis_menten"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "km": 10.0,
            "substrate_concentration": 10.0
        }),
    );

    let request = create_request("michaelis_menten", params);
    let response = biology::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Vmax"));
}

#[test]
fn test_michaelis_menten_negative_km() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("michaelis_menten"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "vmax": 100.0,
            "km": -10.0,
            "substrate_concentration": 10.0
        }),
    );

    let request = create_request("michaelis_menten", params);
    let response = biology::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("positive"));
}

// ============================================================================
// Lineweaver-Burk Plot Tests
// ============================================================================

#[test]
fn test_lineweaver_burk_basic() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("lineweaver_burk"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "vmax": 100.0,
            "km": 10.0,
            "substrate_concentration": 5.0
        }),
    );

    let request = create_request("lineweaver_burk", params);
    let response = biology::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Check that reciprocal values are calculated
    let additional = result["additional_data"].as_object().unwrap();
    let x = additional["x_intercept_1_over_s"].as_f64().unwrap();
    let y = additional["y_intercept_1_over_v"].as_f64().unwrap();

    assert!(x > 0.0);
    assert!(y > 0.0);
}

#[test]
fn test_lineweaver_burk_slope() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("lineweaver_burk"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "vmax": 50.0,
            "km": 5.0,
            "substrate_concentration": 10.0
        }),
    );

    let request = create_request("lineweaver_burk", params);
    let response = biology::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let additional = result["additional_data"].as_object().unwrap();

    // Slope should be Km/Vmax = 5.0/50.0 = 0.1
    let slope = additional["slope_km_over_vmax"].as_f64().unwrap();
    assert!((slope - 0.1).abs() < 0.001);

    // Y-intercept should be 1/Vmax = 1/50 = 0.02
    let y_intercept = additional["y_intercept_1_over_vmax"].as_f64().unwrap();
    assert!((y_intercept - 0.02).abs() < 0.001);
}

#[test]
fn test_lineweaver_burk_zero_substrate() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("lineweaver_burk"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "vmax": 100.0,
            "km": 10.0,
            "substrate_concentration": 0.0
        }),
    );

    let request = create_request("lineweaver_burk", params);
    let response = biology::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("positive"));
}

// ============================================================================
// Pharmacokinetics Tests
// ============================================================================

#[test]
fn test_pharmacokinetics_basic() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("pharmacokinetics"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "dose": 500.0, // mg
            "volume_distribution": 50.0, // L
            "bioavailability": 1.0,
            "elimination_rate": 0.1, // 1/h
            "time": 5.0 // hours
        }),
    );

    let request = create_request("pharmacokinetics", params);
    let response = biology::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let concentration = result["value"].as_f64().unwrap();
    assert!(concentration > 0.0);
    assert_eq!(result["unit"].as_str().unwrap(), "mg/L");

    // Check half-life is calculated
    let additional = result["additional_data"].as_object().unwrap();
    let half_life = additional["half_life_hours"].as_f64().unwrap();

    // t½ = 0.693/k = 0.693/0.1 ≈ 6.93 hours
    assert!((half_life - 6.93).abs() < 0.1);
}

#[test]
fn test_pharmacokinetics_partial_bioavailability() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("pharmacokinetics"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "dose": 500.0,
            "volume_distribution": 50.0,
            "bioavailability": 0.7, // 70% bioavailable
            "elimination_rate": 0.1,
            "time": 0.0
        }),
    );

    let request = create_request("pharmacokinetics", params);
    let response = biology::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // At t=0: C = F·Dose/V = 0.7·500/50 = 7.0 mg/L
    let concentration = result["value"].as_f64().unwrap();
    assert!((concentration - 7.0).abs() < 0.1);
}

#[test]
fn test_pharmacokinetics_clearance() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("pharmacokinetics"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "dose": 1000.0,
            "volume_distribution": 100.0,
            "elimination_rate": 0.2,
            "time": 10.0
        }),
    );

    let request = create_request("pharmacokinetics", params);
    let response = biology::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let additional = result["additional_data"].as_object().unwrap();

    // Clearance = k·V = 0.2·100 = 20 L/h
    let clearance = additional["clearance_l_per_h"].as_f64().unwrap();
    assert!((clearance - 20.0).abs() < 0.1);
}

// ============================================================================
// Hardy-Weinberg Equilibrium Tests
// ============================================================================

#[test]
fn test_hardy_weinberg_equal_alleles() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("hardy_weinberg"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "allele_frequency_p": 0.5
        }),
    );

    let request = create_request("hardy_weinberg", params);
    let response = biology::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // p = q = 0.5
    // AA = 0.25, Aa = 0.5, aa = 0.25
    let aa = result["value"].as_f64().unwrap();
    assert!((aa - 0.25).abs() < 0.01);

    let additional = result["additional_data"].as_object().unwrap();
    let aa_het = additional["freq_Aa_heterozygous"].as_f64().unwrap();
    let aa_rec = additional["freq_aa_homozygous_recessive"].as_f64().unwrap();

    assert!((aa_het - 0.5).abs() < 0.01);
    assert!((aa_rec - 0.25).abs() < 0.01);
}

#[test]
fn test_hardy_weinberg_rare_allele() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("hardy_weinberg"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "allele_frequency_p": 0.1
        }),
    );

    let request = create_request("hardy_weinberg", params);
    let response = biology::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // p = 0.1, q = 0.9
    // AA = 0.01, Aa = 0.18, aa = 0.81
    let aa = result["value"].as_f64().unwrap();
    assert!((aa - 0.01).abs() < 0.01);

    let additional = result["additional_data"].as_object().unwrap();
    let aa_rec = additional["freq_aa_homozygous_recessive"].as_f64().unwrap();
    assert!((aa_rec - 0.81).abs() < 0.01);
}

#[test]
fn test_hardy_weinberg_invalid_frequency() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("hardy_weinberg"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "allele_frequency_p": 1.5 // Invalid: > 1.0
        }),
    );

    let request = create_request("hardy_weinberg", params);
    let response = biology::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("between 0 and 1"));
}

// ============================================================================
// Goldman-Hodgkin-Katz Equation Tests
// ============================================================================

#[test]
fn test_goldman_equation_typical_neuron() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("goldman_equation"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "ion_concentrations_inside": [140.0, 5.0, 10.0],  // K+, Na+, Cl-
            "ion_concentrations_outside": [5.0, 145.0, 110.0],
            "permeabilities": [1.0, 0.04, 0.45],
            "temperature": 310.0 // 37°C
        }),
    );

    let request = create_request("goldman_equation", params);
    let response = biology::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let potential = result["value"].as_f64().unwrap();

    // Typical resting potential is around -70 mV
    assert!(potential < 0.0);
    assert!(potential > -100.0);
    assert_eq!(result["unit"].as_str().unwrap(), "mV");
}

#[test]
fn test_goldman_equation_mismatched_arrays() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("goldman_equation"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "ion_concentrations_inside": [140.0, 5.0],
            "ion_concentrations_outside": [5.0, 145.0, 110.0], // Different length
            "permeabilities": [1.0, 0.04, 0.45],
        }),
    );

    let request = create_request("goldman_equation", params);
    let response = biology::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("same length"));
}

// ============================================================================
// Allometric Scaling Tests
// ============================================================================

#[test]
fn test_allometric_metabolic_rate() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("allometric_scaling"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "body_mass": 70.0, // kg (human)
            "scaling_type": "metabolic"
        }),
    );

    let request = create_request("allometric_scaling", params);
    let response = biology::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    // Kleiber's law: BMR = 70·M^0.75
    let bmr = result["value"].as_f64().unwrap();
    assert!(bmr > 1000.0); // Should be reasonable for human
    assert_eq!(result["unit"].as_str().unwrap(), "kcal/day");

    let formula = result["formula_used"].as_str().unwrap();
    assert!(formula.contains("0.75"));
}

#[test]
fn test_allometric_surface_area() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("allometric_scaling"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "body_mass": 70.0,
            "scaling_type": "surface_area"
        }),
    );

    let request = create_request("allometric_scaling", params);
    let response = biology::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let bsa = result["value"].as_f64().unwrap();
    assert!(bsa > 0.0);
    assert_eq!(result["unit"].as_str().unwrap(), "m²");

    let formula = result["formula_used"].as_str().unwrap();
    assert!(formula.contains("0.67"));
}

#[test]
fn test_allometric_lifespan() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("allometric_scaling"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "body_mass": 5000.0, // Elephant
            "scaling_type": "lifespan"
        }),
    );

    let request = create_request("allometric_scaling", params);
    let response = biology::handle(&request);

    assert!(response.success);
    let result = get_result(&response);

    let lifespan = result["value"].as_f64().unwrap();
    assert!(lifespan > 0.0);
    assert_eq!(result["unit"].as_str().unwrap(), "years");

    let formula = result["formula_used"].as_str().unwrap();
    assert!(formula.contains("0.25"));
}

#[test]
fn test_allometric_unknown_type() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("allometric_scaling"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "body_mass": 70.0,
            "scaling_type": "unknown_type"
        }),
    );

    let request = create_request("allometric_scaling", params);
    let response = biology::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Unknown scaling type"));
}

// ============================================================================
// General Error Handling Tests
// ============================================================================

#[test]
fn test_invalid_json_structure() {
    let mut params = HashMap::new();
    params.insert("invalid_field".to_string(), json!("value"));

    let request = create_request("michaelis_menten", params);
    let response = biology::handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
    assert!(response.error.as_ref().unwrap().contains("Invalid biology request"));
}

#[test]
fn test_response_structure() {
    let mut params = HashMap::new();
    params.insert(
        "operation".to_string(),
        json!("hardy_weinberg"),
    );
    params.insert(
        "parameters".to_string(),
        json!({
            "allele_frequency_p": 0.5
        }),
    );

    let request = create_request("hardy_weinberg", params);
    let response = biology::handle(&request);

    // Verify response structure
    assert!(response.success);
    assert_eq!(response.module, "biology");
    assert_eq!(response.operation, "hardy_weinberg");
    assert!(response.result.is_some());
    assert!(response.error.is_none());

    // Verify biology result fields
    let result = get_result(&response);
    assert!(result["value"].is_number());
    assert!(result["unit"].is_string());
    assert!(result["formula_used"].is_string());
    assert!(result["interpretation"].is_string());
}
