//! Unit tests for main API dispatcher
//!
//! Tests the request routing logic in `src/api/handlers/mod.rs` including:
//! - Module routing to correct handlers
//! - Module alias support
//! - Unknown module error handling
//! - All domain categories (Mathematics, Physics, Tools, Specialized, Scientific)

use crate::api::handlers::route_request;
use crate::api::types::ComputationRequest;
use serde_json::json;
use std::collections::HashMap;

fn create_request(module: &str, operation: &str, parameters: serde_json::Value) -> ComputationRequest {
    let params_map: HashMap<String, serde_json::Value> = serde_json::from_value(parameters)
        .expect("Failed to convert parameters to HashMap");

    ComputationRequest {
        module: module.to_string(),
        operation: operation.to_string(),
        parameters: params_map,
    }
}

// ============================================================================
// Mathematics Domain Tests
// ============================================================================

#[test]
fn test_route_tensor_calculus() {
    let params = json!({
        "metric": [[1.0, 0.0], [0.0, 1.0]]
    });
    let request = create_request("tensor_calculus", "christoffel_symbols", params);
    let response = route_request(&request);

    assert_eq!(response.module, "tensor_calculus");
    assert_eq!(response.operation, "christoffel_symbols");
}

#[test]
fn test_route_tensor_alias() {
    let params = json!({
        "metric": [[1.0, 0.0], [0.0, 1.0]]
    });
    let request = create_request("tensor", "christoffel_symbols", params);
    let response = route_request(&request);

    // Handlers preserve the original request module name (the alias)
    assert_eq!(response.module, "tensor");
    assert_eq!(response.operation, "christoffel_symbols");
}

#[test]
fn test_route_advanced_calculus() {
    let params = json!({
        "s": {"real": 2.0}
    });
    let request = create_request("advanced_calculus", "riemann_zeta", params);
    let response = route_request(&request);

    assert_eq!(response.module, "advanced_calculus");
    assert_eq!(response.operation, "riemann_zeta");
}

#[test]
fn test_route_calculus_alias() {
    let params = json!({
        "s": {"real": 2.0}
    });
    let request = create_request("calculus", "riemann_zeta", params);
    let response = route_request(&request);

    assert_eq!(response.module, "calculus");
}

#[test]
fn test_route_linear_algebra() {
    let params = json!({
        "matrix": [[1.0, 2.0], [3.0, 4.0]]
    });
    let request = create_request("linear_algebra", "compute_svd", params);
    let response = route_request(&request);

    assert_eq!(response.module, "linear_algebra");
    assert_eq!(response.operation, "compute_svd");
}

#[test]
fn test_route_linalg_alias() {
    let params = json!({
        "matrix": [[1.0, 2.0], [3.0, 4.0]]
    });
    let request = create_request("linalg", "compute_svd", params);
    let response = route_request(&request);

    assert_eq!(response.module, "linalg");
}

#[test]
fn test_route_special_functions() {
    let params = json!({
        "order": 0,
        "x": 1.0
    });
    let request = create_request("special_functions", "bessel_j", params);
    let response = route_request(&request);

    assert_eq!(response.module, "special_functions");
    assert_eq!(response.operation, "bessel_j");
}

#[test]
fn test_route_special_alias() {
    let params = json!({
        "order": 0,
        "x": 1.0
    });
    let request = create_request("special", "bessel_j", params);
    let response = route_request(&request);

    assert_eq!(response.module, "special");
}

// ============================================================================
// Physics Domain Tests
// ============================================================================

#[test]
fn test_route_fluid_dynamics() {
    let params = json!({
        "density": 1000.0,
        "velocity": 2.0,
        "diameter": 0.1,
        "viscosity": 0.001
    });
    let request = create_request("fluid_dynamics", "reynolds_number", params);
    let response = route_request(&request);

    assert_eq!(response.module, "fluid_dynamics");
    assert_eq!(response.operation, "reynolds_number");
}

#[test]
fn test_route_fluid_alias() {
    let params = json!({
        "density": 1000.0,
        "velocity": 2.0,
        "diameter": 0.1,
        "viscosity": 0.001
    });
    let request = create_request("fluid", "reynolds_number", params);
    let response = route_request(&request);

    assert_eq!(response.module, "fluid");
}

#[test]
fn test_route_quantum_physics() {
    let params = json!({
        "wavelength": 500e-9
    });
    let request = create_request("quantum_physics", "photon_energy", params);
    let response = route_request(&request);

    assert_eq!(response.module, "quantum_physics");
    assert_eq!(response.operation, "photon_energy");
}

#[test]
fn test_route_quantum_alias() {
    let params = json!({
        "wavelength": 500e-9
    });
    let request = create_request("quantum", "photon_energy", params);
    let response = route_request(&request);

    assert_eq!(response.module, "quantum");
}

#[test]
fn test_route_electromagnetism() {
    let params = json!({
        "charge1": 1e-6,
        "charge2": 2e-6,
        "distance": 0.1
    });
    let request = create_request("electromagnetism", "coulomb_force", params);
    let response = route_request(&request);

    assert_eq!(response.module, "electromagnetism");
    assert_eq!(response.operation, "coulomb_force");
}

#[test]
fn test_route_em_alias() {
    let params = json!({
        "charge1": 1e-6,
        "charge2": 2e-6,
        "distance": 0.1
    });
    let request = create_request("em", "coulomb_force", params);
    let response = route_request(&request);

    assert_eq!(response.module, "em");
}

// ============================================================================
// Tools Domain Tests
// ============================================================================

#[test]
fn test_route_signal_processing() {
    let params = json!({
        "signal": [1.0, 2.0, 3.0, 4.0]
    });
    let request = create_request("signal_processing", "fft", params);
    let response = route_request(&request);

    assert_eq!(response.module, "signal_processing");
    assert_eq!(response.operation, "fft");
}

#[test]
fn test_route_signal_alias() {
    let params = json!({
        "signal": [1.0, 2.0, 3.0, 4.0]
    });
    let request = create_request("signal", "fft", params);
    let response = route_request(&request);

    assert_eq!(response.module, "signal");
}

#[test]
fn test_route_dimensional_analysis() {
    let params = json!({
        "quantity": "force",
        "from_unit": "N",
        "to_unit": "lbf"
    });
    let request = create_request("dimensional_analysis", "convert_units", params);
    let response = route_request(&request);

    assert_eq!(response.module, "dimensional_analysis");
    assert_eq!(response.operation, "convert_units");
}

#[test]
fn test_route_dimensional_alias() {
    let params = json!({
        "quantity": "force",
        "from_unit": "N",
        "to_unit": "lbf"
    });
    let request = create_request("dimensional", "convert_units", params);
    let response = route_request(&request);

    assert_eq!(response.module, "dimensional");
}

#[test]
fn test_route_equation_validation() {
    let params = json!({
        "equation": "F = m*a"
    });
    let request = create_request("equation_validation", "validate_equation", params);
    let response = route_request(&request);

    assert_eq!(response.module, "equation_validation");
    assert_eq!(response.operation, "validate_equation");
}

#[test]
fn test_route_validator_alias() {
    let params = json!({
        "equation": "F = m*a"
    });
    let request = create_request("validator", "validate_equation", params);
    let response = route_request(&request);

    assert_eq!(response.module, "validator");
}

#[test]
fn test_route_computational_geometry() {
    let params = json!({
        "points": [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
    });
    let request = create_request("computational_geometry", "convex_hull", params);
    let response = route_request(&request);

    assert_eq!(response.module, "computational_geometry");
    assert_eq!(response.operation, "convex_hull");
}

#[test]
fn test_route_geometry_alias() {
    let params = json!({
        "points": [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
    });
    let request = create_request("geometry", "convex_hull", params);
    let response = route_request(&request);

    assert_eq!(response.module, "geometry");
}

#[test]
fn test_route_numerical_methods() {
    let params = json!({
        "function": "x^2 - 4",
        "initial_guess": 1.0
    });
    let request = create_request("numerical_methods", "newton_raphson", params);
    let response = route_request(&request);

    assert_eq!(response.module, "numerical_methods");
    assert_eq!(response.operation, "newton_raphson");
}

#[test]
fn test_route_numerical_alias() {
    let params = json!({
        "function": "x^2 - 4",
        "initial_guess": 1.0
    });
    let request = create_request("numerical", "newton_raphson", params);
    let response = route_request(&request);

    assert_eq!(response.module, "numerical");
}

#[test]
fn test_route_function_approximator() {
    let params = json!({
        "x_values": [1.0, 2.0, 3.0],
        "y_values": [2.0, 4.0, 6.0]
    });
    let request = create_request("function_approximator", "polynomial_fit", params);
    let response = route_request(&request);

    assert_eq!(response.module, "function_approximator");
    assert_eq!(response.operation, "polynomial_fit");
}

#[test]
fn test_route_approximate_alias() {
    let params = json!({
        "x_values": [1.0, 2.0, 3.0],
        "y_values": [2.0, 4.0, 6.0]
    });
    let request = create_request("approximate", "polynomial_fit", params);
    let response = route_request(&request);

    assert_eq!(response.module, "approximate");
}

// ============================================================================
// Specialized Domain Tests
// ============================================================================

#[test]
fn test_route_stochastic_processes() {
    let params = json!({
        "n_steps": 100,
        "dt": 0.01
    });
    let request = create_request("stochastic_processes", "brownian_motion", params);
    let response = route_request(&request);

    assert_eq!(response.module, "stochastic_processes");
    assert_eq!(response.operation, "brownian_motion");
}

#[test]
fn test_route_stochastic_alias() {
    let params = json!({
        "n_steps": 100,
        "dt": 0.01
    });
    let request = create_request("stochastic", "brownian_motion", params);
    let response = route_request(&request);

    assert_eq!(response.module, "stochastic");
}

#[test]
fn test_route_cryptographic_mathematics() {
    let params = json!({
        "n": 100
    });
    let request = create_request("cryptographic_mathematics", "generate_prime", params);
    let response = route_request(&request);

    assert_eq!(response.module, "cryptographic_mathematics");
    assert_eq!(response.operation, "generate_prime");
}

#[test]
fn test_route_crypto_alias() {
    let params = json!({
        "n": 100
    });
    let request = create_request("crypto", "generate_prime", params);
    let response = route_request(&request);

    assert_eq!(response.module, "crypto");
}

#[test]
fn test_route_symbolic_regression() {
    let params = json!({
        "x_values": [1.0, 2.0, 3.0],
        "y_values": [2.0, 4.0, 6.0]
    });
    let request = create_request("symbolic_regression", "genetic_programming", params);
    let response = route_request(&request);

    assert_eq!(response.module, "symbolic_regression");
    assert_eq!(response.operation, "genetic_programming");
}

#[test]
fn test_route_regression_alias() {
    let params = json!({
        "x_values": [1.0, 2.0, 3.0],
        "y_values": [2.0, 4.0, 6.0]
    });
    let request = create_request("regression", "genetic_programming", params);
    let response = route_request(&request);

    assert_eq!(response.module, "regression");
}

#[test]
fn test_route_statistics() {
    let params = json!({
        "data": [1.0, 2.0, 3.0, 4.0, 5.0]
    });
    let request = create_request("statistics", "mean", params);
    let response = route_request(&request);

    assert_eq!(response.module, "statistics");
    assert_eq!(response.operation, "mean");
}

#[test]
fn test_route_stats_alias() {
    let params = json!({
        "data": [1.0, 2.0, 3.0, 4.0, 5.0]
    });
    let request = create_request("stats", "mean", params);
    let response = route_request(&request);

    assert_eq!(response.module, "stats");
}

#[test]
fn test_route_optimization() {
    let params = json!({
        "function": "x^2",
        "initial_guess": 1.0
    });
    let request = create_request("optimization", "minimize", params);
    let response = route_request(&request);

    assert_eq!(response.module, "optimization");
    assert_eq!(response.operation, "minimize");
}

#[test]
fn test_route_optimize_alias() {
    let params = json!({
        "function": "x^2",
        "initial_guess": 1.0
    });
    let request = create_request("optimize", "minimize", params);
    let response = route_request(&request);

    assert_eq!(response.module, "optimize");
}

#[test]
fn test_route_graph_theory() {
    let params = json!({
        "graph": {"vertices": 4, "edges": [[0, 1], [1, 2], [2, 3]]}
    });
    let request = create_request("graph_theory", "shortest_path", params);
    let response = route_request(&request);

    assert_eq!(response.module, "graph_theory");
    assert_eq!(response.operation, "shortest_path");
}

#[test]
fn test_route_graph_alias() {
    let params = json!({
        "graph": {"vertices": 4, "edges": [[0, 1], [1, 2], [2, 3]]}
    });
    let request = create_request("graph", "shortest_path", params);
    let response = route_request(&request);

    assert_eq!(response.module, "graph");
}

#[test]
fn test_route_information_theory() {
    let params = json!({
        "data": [1, 1, 2, 2, 3, 3, 3]
    });
    let request = create_request("information_theory", "entropy", params);
    let response = route_request(&request);

    assert_eq!(response.module, "information_theory");
    assert_eq!(response.operation, "entropy");
}

#[test]
fn test_route_info_alias() {
    let params = json!({
        "data": [1, 1, 2, 2, 3, 3, 3]
    });
    let request = create_request("info", "entropy", params);
    let response = route_request(&request);

    assert_eq!(response.module, "info");
}

#[test]
fn test_route_chemistry() {
    let params = json!({
        "compound": "H2O"
    });
    let request = create_request("chemistry", "molecular_weight", params);
    let response = route_request(&request);

    assert_eq!(response.module, "chemistry");
    assert_eq!(response.operation, "molecular_weight");
}

#[test]
fn test_route_chem_alias() {
    let params = json!({
        "compound": "H2O"
    });
    let request = create_request("chem", "molecular_weight", params);
    let response = route_request(&request);

    assert_eq!(response.module, "chem");
}

// ============================================================================
// Scientific Formulas Domain Tests (2025 expansion)
// ============================================================================

#[test]
fn test_route_biology() {
    let params = json!({
        "population": 1000,
        "growth_rate": 0.05
    });
    let request = create_request("biology", "population_growth", params);
    let response = route_request(&request);

    assert_eq!(response.module, "biology");
    assert_eq!(response.operation, "population_growth");
}

#[test]
fn test_route_bio_alias() {
    let params = json!({
        "population": 1000,
        "growth_rate": 0.05
    });
    let request = create_request("bio", "population_growth", params);
    let response = route_request(&request);

    assert_eq!(response.module, "bio");
}

#[test]
fn test_route_thermodynamics() {
    let params = json!({
        "pressure": 101325.0,
        "volume": 0.001,
        "temperature": 273.15
    });
    let request = create_request("thermodynamics", "ideal_gas", params);
    let response = route_request(&request);

    assert_eq!(response.module, "thermodynamics");
    assert_eq!(response.operation, "ideal_gas");
}

#[test]
fn test_route_thermo_alias() {
    let params = json!({
        "pressure": 101325.0,
        "volume": 0.001,
        "temperature": 273.15
    });
    let request = create_request("thermo", "ideal_gas", params);
    let response = route_request(&request);

    assert_eq!(response.module, "thermo");
}

#[test]
fn test_route_heat_alias() {
    let params = json!({
        "pressure": 101325.0,
        "volume": 0.001,
        "temperature": 273.15
    });
    let request = create_request("heat", "ideal_gas", params);
    let response = route_request(&request);

    assert_eq!(response.module, "heat");
}

#[test]
fn test_route_optics() {
    let params = json!({
        "wavelength": 500e-9,
        "slit_width": 1e-6
    });
    let request = create_request("optics", "diffraction", params);
    let response = route_request(&request);

    assert_eq!(response.module, "optics");
    assert_eq!(response.operation, "diffraction");
}

#[test]
fn test_route_optical_alias() {
    let params = json!({
        "wavelength": 500e-9,
        "slit_width": 1e-6
    });
    let request = create_request("optical", "diffraction", params);
    let response = route_request(&request);

    assert_eq!(response.module, "optical");
}

#[test]
fn test_route_geophysics() {
    let params = json!({
        "latitude": 45.0
    });
    let request = create_request("geophysics", "gravity_variation", params);
    let response = route_request(&request);

    assert_eq!(response.module, "geophysics");
    assert_eq!(response.operation, "gravity_variation");
}

#[test]
fn test_route_geo_alias() {
    let params = json!({
        "latitude": 45.0
    });
    let request = create_request("geo", "gravity_variation", params);
    let response = route_request(&request);

    assert_eq!(response.module, "geo");
}

#[test]
fn test_route_engineering() {
    let params = json!({
        "stress": 100e6,
        "youngs_modulus": 200e9
    });
    let request = create_request("engineering", "strain", params);
    let response = route_request(&request);

    assert_eq!(response.module, "engineering");
    assert_eq!(response.operation, "strain");
}

#[test]
fn test_route_eng_alias() {
    let params = json!({
        "stress": 100e6,
        "youngs_modulus": 200e9
    });
    let request = create_request("eng", "strain", params);
    let response = route_request(&request);

    assert_eq!(response.module, "eng");
}

#[test]
fn test_route_datetime() {
    let params = json!({
        "year": 2025,
        "month": 10,
        "day": 26
    });
    let request = create_request("datetime", "day_of_week", params);
    let response = route_request(&request);

    assert_eq!(response.module, "datetime");
    assert_eq!(response.operation, "day_of_week");
}

#[test]
fn test_route_date_alias() {
    let params = json!({
        "year": 2025,
        "month": 10,
        "day": 26
    });
    let request = create_request("date", "day_of_week", params);
    let response = route_request(&request);

    assert_eq!(response.module, "date");
}

#[test]
fn test_route_time_alias() {
    let params = json!({
        "year": 2025,
        "month": 10,
        "day": 26
    });
    let request = create_request("time", "day_of_week", params);
    let response = route_request(&request);

    assert_eq!(response.module, "time");
}

// ============================================================================
// Error Handling Tests
// ============================================================================

#[test]
fn test_route_unknown_module() {
    let params = json!({});
    let request = create_request("unknown_module", "some_operation", params);
    let response = route_request(&request);

    assert!(!response.success, "Unknown module should fail");
    assert_eq!(response.module, "unknown_module");
    assert_eq!(response.operation, "some_operation");
    assert!(response.error.is_some(), "Should have error message");

    let error_msg = response.error.unwrap();
    assert!(error_msg.contains("Unknown module"), "Error should mention unknown module");
    assert!(error_msg.contains("unknown_module"), "Error should include module name");
}

#[test]
fn test_route_empty_module_name() {
    let params = json!({});
    let request = create_request("", "some_operation", params);
    let response = route_request(&request);

    assert!(!response.success, "Empty module should fail");
    assert!(response.error.is_some(), "Should have error message");
}

#[test]
fn test_route_case_sensitive_module() {
    // Module names are case-sensitive, so "TENSOR" should not match "tensor"
    let params = json!({
        "metric": [[1.0, 0.0], [0.0, 1.0]]
    });
    let request = create_request("TENSOR", "christoffel_symbols", params);
    let response = route_request(&request);

    assert!(!response.success, "Case mismatch should fail");
    assert!(response.error.is_some(), "Should have error message");
}
