//! Integration tests for the computational engine API

use computational_engine::api::*;
use serde_json::json;

#[test]
fn test_api_linear_algebra() {
    let request = json!({
        "module": "linear_algebra",
        "operation": "compute_svd",
        "parameters": {
            "matrix": [[12.0, -51.0, 4.0], [6.0, 167.0, -68.0], [-4.0, 24.0, -41.0]]
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();

    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_chemistry() {
    let request = json!({
        "module": "chemistry",
        "operation": "molar_mass",
        "parameters": {
            "formula": "H2O"
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();

    assert!(parsed["success"].as_bool().unwrap_or(false));
    if let Some(result) = parsed.get("result") {
        let molar_mass = result["molar_mass"].as_f64().unwrap();
        assert!((molar_mass - 18.015).abs() < 0.1);
    }
}

#[test]
fn test_api_numerical_methods() {
    let request = json!({
        "module": "numerical_methods",
        "operation": "find_root",
        "parameters": {
            "method": "bisection",
            "initial_guess": 0.0,
            "interval": [1.0, 2.0],
            "tolerance": 0.000001,
            "max_iterations": 100,
            "function_type": "polynomial",
            "coefficients": [-2.0, 0.0, 1.0]
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();

    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_electromagnetism() {
    let request = json!({
        "module": "electromagnetism",
        "operation": "em_wave",
        "parameters": {
            "frequency": 1000000000.0,
            "medium": "vacuum"
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();

    assert!(parsed["success"].as_bool().unwrap_or(false));
    if let Some(result) = parsed.get("result") {
        let wavelength = result["wavelength"].as_f64().unwrap();
        assert!((wavelength - 0.3).abs() < 0.1);
    }
}

#[test]
fn test_api_special_functions() {
    let request = json!({
        "module": "special_functions",
        "operation": "gamma",
        "parameters": {
            "x": 5.0,
            "function": "gamma"
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();

    assert!(parsed["success"].as_bool().unwrap_or(false));
    if let Some(result) = parsed.get("result") {
        let value = result["value"].as_f64().unwrap();
        // Γ(5) = 4! = 24
        assert!((value - 24.0).abs() < 0.1);
    }
}

#[test]
fn test_api_unknown_module() {
    let request = json!({
        "module": "nonexistent_module",
        "operation": "some_operation",
        "parameters": {}
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();

    assert!(!parsed["success"].as_bool().unwrap_or(true));
    assert!(parsed.get("error").is_some());
}

#[test]
fn test_list_operations() {
    let ops = list_all_operations();

    // Should have multiple modules
    assert!(ops.len() > 10);

    // Check specific modules exist
    assert!(ops.contains_key("linear_algebra"));
    assert!(ops.contains_key("chemistry"));
    assert!(ops.contains_key("electromagnetism"));
    assert!(ops.contains_key("numerical_methods"));

    // Check operation count
    let total: usize = ops.values().map(|v| v.len()).sum();
    assert!(total >= 190);
}

// Graph Theory tests
#[test]
fn test_api_graph_shortest_path() {
    let request = json!({
        "module": "graph_theory",
        "operation": "shortest_path",
        "parameters": {
            "graph": {
                "nodes": ["A", "B", "C", "D"],
                "edges": [
                    {"from": "A", "to": "B", "weight": 4.0},
                    {"from": "A", "to": "C", "weight": 2.0},
                    {"from": "C", "to": "D", "weight": 5.0},
                    {"from": "B", "to": "D", "weight": 1.0}
                ]
            },
            "source": "A",
            "target": "D"
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_graph_mst() {
    let request = json!({
        "module": "graph_theory",
        "operation": "minimum_spanning_tree",
        "parameters": {
            "graph": {
                "nodes": ["A", "B", "C"],
                "edges": [
                    {"from": "A", "to": "B", "weight": 1.0},
                    {"from": "B", "to": "C", "weight": 2.0},
                    {"from": "A", "to": "C", "weight": 3.0}
                ]
            }
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// Computational Geometry tests
#[test]
fn test_api_convex_hull() {
    let request = json!({
        "module": "computational_geometry",
        "operation": "convex_hull",
        "parameters": {
            "points": [
                {"x": 0.0, "y": 0.0},
                {"x": 1.0, "y": 1.0},
                {"x": 2.0, "y": 0.0},
                {"x": 1.0, "y": -1.0}
            ]
        }
    });

    let response = process_json_request(&request.to_string());
    eprintln!("Convex hull response: {}", response);
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_polygon_area() {
    let request = json!({
        "module": "computational_geometry",
        "operation": "polygon_area",
        "parameters": {
            "vertices": [
                {"x": 0.0, "y": 0.0},
                {"x": 4.0, "y": 0.0},
                {"x": 4.0, "y": 3.0},
                {"x": 0.0, "y": 3.0}
            ]
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
    if let Some(result) = parsed.get("result") {
        let area = result["area"].as_f64().unwrap();
        assert!((area - 12.0).abs() < 0.1);
    }
}

// Information Theory tests
#[test]
fn test_api_shannon_entropy() {
    let request = json!({
        "module": "information_theory",
        "operation": "shannon_entropy",
        "parameters": {
            "data": [0.5, 0.25, 0.25]
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_huffman_coding() {
    let request = json!({
        "module": "information_theory",
        "operation": "huffman_coding",
        "parameters": {
            "symbols": ["A", "B", "C"],
            "frequencies": [5.0, 2.0, 1.0]
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// Additional Chemistry tests
#[test]
fn test_api_chemistry_thermodynamics() {
    let request = json!({
        "module": "chemistry",
        "operation": "thermodynamics",
        "parameters": {
            "temperature": 298.15,
            "pressure": 1.0,
            "compounds": [
                {"formula": "H2O", "moles": 1.0, "enthalpy_formation": -285.8, "entropy": 69.9}
            ]
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_chemistry_gas_law() {
    let request = json!({
        "module": "chemistry",
        "operation": "gas_law",
        "parameters": {
            "pressure": 1.0,
            "volume": 22.4,
            "temperature": 273.15
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// Numerical Methods tests
#[test]
fn test_api_numerical_integration() {
    let request = json!({
        "module": "numerical_methods",
        "operation": "integrate",
        "parameters": {
            "method": "simpson",
            "lower_bound": 0.0,
            "upper_bound": 1.0,
            "num_points": 100,
            "function_type": "polynomial",
            "coefficients": [0.0, 0.0, 1.0]
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_numerical_interpolation() {
    let request = json!({
        "module": "numerical_methods",
        "operation": "interpolate",
        "parameters": {
            "method": "lagrange",
            "x_values": [0.0, 1.0, 2.0],
            "y_values": [0.0, 1.0, 4.0],
            "interpolate_at": [0.5, 1.5]
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// Special Functions tests
#[test]
fn test_api_bessel_function() {
    let request = json!({
        "module": "special_functions",
        "operation": "bessel",
        "parameters": {
            "function_type": "J",
            "order": 0.0,
            "x": 2.0
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_error_function() {
    let request = json!({
        "module": "special_functions",
        "operation": "error_function",
        "parameters": {
            "x": 1.0,
            "function": "erf"
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// Electromagnetism tests
#[test]
fn test_api_em_maxwell() {
    let request = json!({
        "module": "electromagnetism",
        "operation": "maxwell_equations",
        "parameters": {
            "equation": "gauss_electric",
            "charge_density": 1e-6,
            "current_density": [1.0, 0.0, 0.0],
            "electric_field": [100.0, 0.0, 0.0],
            "magnetic_field": [0.0, 0.0, 1e-3],
            "position": [0.0, 0.0, 0.0]
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_em_poynting() {
    let request = json!({
        "module": "electromagnetism",
        "operation": "poynting_vector",
        "parameters": {
            "electric_field": [100.0, 0.0, 0.0],
            "magnetic_field": [0.0, 1e-3, 0.0]
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// Linear Algebra additional tests
#[test]
fn test_api_linear_algebra_eigendecomposition() {
    let request = json!({
        "module": "linear_algebra",
        "operation": "compute_eigendecomposition",
        "parameters": {
            "matrix": [[2.0, 1.0], [1.0, 2.0]]
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_linear_algebra_matrix_rank() {
    let request = json!({
        "module": "linear_algebra",
        "operation": "compute_matrix_rank",
        "parameters": {
            "matrix": [[1.0, 2.0], [2.0, 4.0]]
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// === Additional Comprehensive Tests ===

// More Linear Algebra tests
#[test]
fn test_api_linear_algebra_pca() {
    let request = json!({
        "module": "linear_algebra",
        "operation": "compute_pca",
        "parameters": {
            "matrix": [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]],
            "n_components": 2
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_linear_algebra_pseudoinverse() {
    let request = json!({
        "module": "linear_algebra",
        "operation": "compute_pseudoinverse",
        "parameters": {
            "matrix": [[1.0, 2.0], [3.0, 4.0]]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_linear_algebra_cholesky() {
    let request = json!({
        "module": "linear_algebra",
        "operation": "cholesky",
        "parameters": {
            "matrix": [[4.0, 2.0], [2.0, 3.0]]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// More Graph Theory tests
#[test]
fn test_api_graph_connected_components() {
    let request = json!({
        "module": "graph_theory",
        "operation": "connected_components",
        "parameters": {
            "graph": {
                "nodes": ["A", "B", "C", "D"],
                "edges": [
                    {"from": "A", "to": "B", "weight": null},
                    {"from": "C", "to": "D", "weight": null}
                ]
            }
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_graph_topological_sort() {
    let request = json!({
        "module": "graph_theory",
        "operation": "topological_sort",
        "parameters": {
            "graph": {
                "nodes": ["A", "B", "C"],
                "edges": [
                    {"from": "A", "to": "B", "weight": null},
                    {"from": "B", "to": "C", "weight": null}
                ]
            }
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// More Information Theory tests
#[test]
fn test_api_info_mutual_information() {
    let request = json!({
        "module": "information_theory",
        "operation": "mutual_information",
        "parameters": {
            "x": [1.0, 2.0, 3.0, 4.0],
            "y": [2.0, 4.0, 6.0, 8.0]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_info_conditional_entropy() {
    let request = json!({
        "module": "information_theory",
        "operation": "conditional_entropy",
        "parameters": {
            "x": [1.0, 2.0, 1.0, 2.0],
            "y": [1.0, 1.0, 2.0, 2.0]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// More Chemistry tests
#[test]
fn test_api_chemistry_kinetics() {
    let request = json!({
        "module": "chemistry",
        "operation": "kinetics",
        "parameters": {
            "rate_constant": 1.0,
            "temperature": 298.15,
            "activation_energy": 50.0,
            "concentrations": [1.0],
            "orders": [1]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_chemistry_acid_base() {
    let request = json!({
        "module": "chemistry",
        "operation": "acid_base",
        "parameters": {
            "concentration": 0.1,
            "acid_ka": 1e-5
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_chemistry_electrochemistry() {
    let request = json!({
        "module": "chemistry",
        "operation": "electrochemistry",
        "parameters": {
            "half_reactions": [
                {
                    "oxidized": "Cu2+",
                    "reduced": "Cu",
                    "electrons": 2,
                    "standard_potential": 0.34
                },
                {
                    "oxidized": "Zn2+",
                    "reduced": "Zn",
                    "electrons": 2,
                    "standard_potential": -0.76
                }
            ],
            "temperature": 298.15
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// More Numerical Methods tests
#[test]
fn test_api_numerical_root_finding() {
    let request = json!({
        "module": "numerical_methods",
        "operation": "find_root",
        "parameters": {
            "method": "bisection",
            "initial_guess": 1.0,
            "interval": [0.0, 2.0],
            "tolerance": 0.0001,
            "max_iterations": 100,
            "function_type": "polynomial",
            "coefficients": [-2.0, 0.0, 1.0]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_numerical_ode_solver() {
    let request = json!({
        "module": "numerical_methods",
        "operation": "solve_ode",
        "parameters": {
            "method": "rk4",
            "initial_value": 1.0,
            "t_start": 0.0,
            "t_end": 1.0,
            "step_size": 0.1,
            "function_type": "exponential",
            "function_coefficients": [-1.0]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// More Special Functions tests
#[test]
fn test_api_gamma_function() {
    let request = json!({
        "module": "special_functions",
        "operation": "gamma",
        "parameters": {
            "x": 5.0,
            "function": "gamma"
        }
    });
    let response = process_json_request(&request.to_string());
    eprintln!("Gamma function response: {}", response);
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_elliptic_integral() {
    let request = json!({
        "module": "special_functions",
        "operation": "elliptic",
        "parameters": {
            "k": 0.5,
            "integral_type": "K"
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// More Electromagnetism tests
#[test]
fn test_api_em_waveguide() {
    let request = json!({
        "module": "electromagnetism",
        "operation": "waveguide",
        "parameters": {
            "guide_type": "rectangular",
            "width": 0.02,
            "height": 0.01,
            "mode": "TE10",
            "frequency": 10e9
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_em_transmission_line() {
    let request = json!({
        "module": "electromagnetism",
        "operation": "transmission_line",
        "parameters": {
            "line_type": "coax",
            "z0": 50.0,
            "load_impedance": 75.0,
            "frequency": 1e9,
            "length": 0.1
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_em_antenna() {
    let request = json!({
        "module": "electromagnetism",
        "operation": "antenna",
        "parameters": {
            "antenna_type": "dipole",
            "length": 0.15,
            "frequency": 1e9,
            "distance": 10.0,
            "theta": 1.57,
            "power": 100.0
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// More Computational Geometry tests
#[test]
fn test_api_point_in_polygon() {
    let request = json!({
        "module": "computational_geometry",
        "operation": "point_in_polygon",
        "parameters": {
            "point": {"x": 2.0, "y": 2.0},
            "polygon": [
                {"x": 0.0, "y": 0.0},
                {"x": 4.0, "y": 0.0},
                {"x": 4.0, "y": 4.0},
                {"x": 0.0, "y": 4.0}
            ]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_delaunay_triangulation() {
    let request = json!({
        "module": "computational_geometry",
        "operation": "delaunay_triangulation",
        "parameters": {
            "points": [
                {"x": 0.0, "y": 0.0},
                {"x": 1.0, "y": 0.0},
                {"x": 0.5, "y": 1.0}
            ]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// Debug test to check errors
#[test]
fn debug_failing_tests() {
    // Test gamma function
    let request = json!({
        "module": "special_functions",
        "operation": "gamma",
        "parameters": {"x": 5.0}
    });
    let response = process_json_request(&request.to_string());
    eprintln!("Gamma: {}", response);

    // Test kinetics
    let request2 = json!({
        "module": "chemistry",
        "operation": "kinetics",
        "parameters": {
            "rate_constant": 1.0,
            "temperature": 298.15,
            "activation_energy": 50.0,
            "concentrations": [1.0],
            "orders": [1]
        }
    });
    let response2 = process_json_request(&request2.to_string());
    eprintln!("Kinetics: {}", response2);
}

#[test]
fn debug_remaining_failures() {
    let tests = vec![
        json!({"module": "chemistry", "operation": "acid_base", "parameters": {"ph": 7.0, "calculate": "poh"}}),
        json!({"module": "chemistry", "operation": "electrochemistry", "parameters": {"standard_potential": 0.34, "temperature": 298.15, "reaction_quotient": 1.0, "n_electrons": 2}}),
        json!({"module": "special_functions", "operation": "elliptic_integral", "parameters": {"k": 0.5, "integral_type": "first"}}),
        json!({"module": "electromagnetism", "operation": "waveguide", "parameters": {"width": 0.02, "height": 0.01, "mode": "TE10", "frequency": 10e9}}),
        json!({"module": "electromagnetism", "operation": "transmission_line", "parameters": {"line_type": "coax",
            "z0": 50.0, "load_impedance": 75.0, "frequency": 1e9, "length": 0.1}}),
        json!({"module": "electromagnetism", "operation": "antenna", "parameters": {"antenna_type": "dipole", "length": 0.15, "frequency": 1e9}}),
        json!({"module": "numerical_methods", "operation": "find_root", "parameters": {"method": "bisection", "initial_guess": 1.0,
            "interval": [0.0, 2.0], "tolerance": 0.0001, "max_iterations": 100, "function_type": "polynomial", "coefficients": [-2.0, 0.0, 1.0]}}),
        json!({"module": "numerical_methods", "operation": "solve_ode", "parameters": {"method": "rk4", "initial_value": 1.0, "t_start": 0.0,
            "t_end": 1.0, "step_size": 0.1, "function_type": "exponential",
            "function_coefficients": [-1.0]}}),
    ];

    for (i, test) in tests.iter().enumerate() {
        let response = process_json_request(&test.to_string());
        let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
        if let Some(error) = parsed.get("error") {
            eprintln!("Test {}: {}", i, error.as_str().unwrap_or("unknown"));
        }
    }
}

// ==================== COMPREHENSIVE COVERAGE TESTS ====================
// Adding tests for all previously untested modules to reach 80%+ coverage

// === STATISTICS MODULE TESTS (14 operations) ===
#[test]
fn test_api_statistics_descriptive() {
    let request = json!({
        "module": "statistics",
        "operation": "descriptive",
        "parameters": {
            "data": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_statistics_regression() {
    let request = json!({
        "module": "statistics",
        "operation": "regression",
        "parameters": {
            "x": [1.0, 2.0, 3.0, 4.0, 5.0],
            "y": [2.0, 4.0, 5.0, 4.0, 5.0],
            "degree": 1
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_statistics_hypothesis_test() {
    let request = json!({
        "module": "statistics",
        "operation": "hypothesis_test",
        "parameters": {
            "test_type": "t_test",
            "sample1": [1.0, 2.0, 3.0, 4.0, 5.0],
            "sample2": [2.0, 3.0, 4.0, 5.0, 6.0]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_statistics_distribution() {
    let request = json!({
        "module": "statistics",
        "operation": "distribution",
        "parameters": {
            "distribution_type": "normal",
            "mean": 0.0,
            "std_dev": 1.0,
            "x": 1.0
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// === OPTIMIZATION MODULE TESTS (13 operations) ===
#[test]
fn test_api_optimization_minimize() {
    let request = json!({
        "module": "optimization",
        "operation": "minimize",
        "parameters": {
            "method": "gradient_descent",
            "initial_guess": [0.0, 0.0],
            "learning_rate": 0.01,
            "max_iterations": 1000
        }
    });
    let response = process_json_request(&request.to_string());
    eprintln!("Optimization minimize response: {}", response);
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_optimization_linear_program() {
    let request = json!({
        "module": "optimization",
        "operation": "linear_program",
        "parameters": {
            "objective": [1.0, 2.0],
            "constraints": [[1.0, 1.0]],
            "bounds": [10.0]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_optimization_genetic_algorithm() {
    let request = json!({
        "module": "optimization",
        "operation": "genetic_algorithm",
        "parameters": {
            "population_size": 50,
            "generations": 100,
            "mutation_rate": 0.01
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// === SIGNAL PROCESSING MODULE TESTS (15 operations) ===
#[test]
fn test_api_signal_fft() {
    let request = json!({
        "module": "signal_processing",
        "operation": "fft",
        "parameters": {
            "signal": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
            "sample_rate": 1000.0,
            "window_type": "hanning"
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_signal_filter() {
    let request = json!({
        "module": "signal_processing",
        "operation": "filter",
        "parameters": {
            "signal": [1.0, 2.0, 3.0, 4.0, 5.0],
            "filter_type": "lowpass",
            "cutoff_frequency": 0.5,
            "sample_rate": 1000.0,
            "order": 4
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_signal_convolution() {
    let request = json!({
        "module": "signal_processing",
        "operation": "convolution",
        "parameters": {
            "signal": [1.0, 2.0, 3.0],
            "kernel": [0.5, 0.5]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_signal_spectrogram() {
    let request = json!({
        "module": "signal_processing",
        "operation": "spectrogram",
        "parameters": {
            "signal": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
            "sample_rate": 1000.0,
            "window_size": 4,
            "hop_size": 2,
            "window_type": "hanning"
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// === QUANTUM PHYSICS MODULE TESTS (12 operations) ===
#[test]
fn test_api_quantum_wavefunction() {
    let request = json!({
        "module": "quantum",
        "operation": "wavefunction",
        "parameters": {
            "system": "harmonic_oscillator",
            "n": 0,
            "x": [0.0, 0.1, 0.2]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_quantum_operators() {
    let request = json!({
        "module": "quantum",
        "operation": "operator",
        "parameters": {
            "operator_type": "momentum",
            "dimension": 2
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_quantum_entanglement() {
    let request = json!({
        "module": "quantum",
        "operation": "entanglement",
        "parameters": {
            "state": [[1.0, 0.0], [0.0, 1.0]]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// === STOCHASTIC PROCESSES MODULE TESTS (12 operations) ===
#[test]
fn test_api_stochastic_brownian_motion() {
    let request = json!({
        "module": "stochastic_processes",
        "operation": "brownian_motion",
        "parameters": {
            "time_steps": 100,
            "dt": 0.01,
            "initial_value": 0.0,
            "drift": 0.0,
            "volatility": 1.0
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_stochastic_poisson_process() {
    let request = json!({
        "module": "stochastic_processes",
        "operation": "poisson_process",
        "parameters": {
            "rate": 2.0,
            "time_horizon": 10.0
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_stochastic_markov_chain() {
    let request = json!({
        "module": "stochastic_processes",
        "operation": "markov_chain",
        "parameters": {
            "states": ["A", "B"],
            "transition_matrix": [[0.7, 0.3], [0.4, 0.6]],
            "initial_state": "A",
            "steps": 10
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// === CRYPTOGRAPHY MODULE TESTS (10 operations) ===
#[test]
fn test_api_crypto_rsa() {
    let request = json!({
        "module": "cryptographic_mathematics",
        "operation": "rsa_keygen",
        "parameters": {
            "bit_length": 512
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_crypto_prime_generation() {
    let request = json!({
        "module": "cryptographic_mathematics",
        "operation": "generate_prime",
        "parameters": {
            "bit_length": 256
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_crypto_modular_arithmetic() {
    let request = json!({
        "module": "cryptographic_mathematics",
        "operation": "modular_exponentiation",
        "parameters": {
            "base": "3",
            "exponent": "5",
            "modulus": "7"
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// === FLUID DYNAMICS MODULE TESTS (10 operations) ===
#[test]
fn test_api_fluid_navier_stokes() {
    let request = json!({
        "module": "fluid_dynamics",
        "operation": "navier_stokes",
        "parameters": {
            "velocity": [1.0, 0.0, 0.0],
            "pressure": 101325.0,
            "density": 1.225,
            "viscosity": 1.81e-5
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_fluid_reynolds_number() {
    let request = json!({
        "module": "fluid_dynamics",
        "operation": "reynolds_number",
        "parameters": {
            "velocity": 10.0,
            "length": 1.0,
            "density": 1.225,
            "viscosity": 1.81e-5
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_fluid_bernoulli() {
    let request = json!({
        "module": "fluid_dynamics",
        "operation": "bernoulli",
        "parameters": {
            "v1": 10.0,
            "p1": 101325.0,
            "h1": 0.0,
            "v2": 5.0,
            "h2": 10.0
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// === CALCULUS MODULE TESTS (12 operations) ===
#[test]
fn test_api_calculus_derivative() {
    let request = json!({
        "module": "calculus",
        "operation": "derivative",
        "parameters": {
            "function_type": "polynomial",
            "coefficients": [1.0, 2.0, 3.0],
            "x": 2.0
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_calculus_integral() {
    let request = json!({
        "module": "calculus",
        "operation": "integral",
        "parameters": {
            "function_type": "polynomial",
            "coefficients": [1.0, 0.0, 1.0],
            "lower": 0.0,
            "upper": 1.0
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_calculus_limit() {
    let request = json!({
        "module": "calculus",
        "operation": "limit",
        "parameters": {
            "function_type": "rational",
            "numerator": [1.0, 0.0],
            "denominator": [1.0],
            "point": 0.0
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// === TENSOR CALCULUS MODULE TESTS (15 operations) ===
#[test]
fn test_api_tensor_christoffel() {
    let request = json!({
        "module": "tensor_calculus",
        "operation": "christoffel_symbols",
        "parameters": {
            "metric": [[1.0, 0.0], [0.0, 1.0]],
            "coordinates": ["x", "y"]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_tensor_riemann() {
    let request = json!({
        "module": "tensor_calculus",
        "operation": "riemann_tensor",
        "parameters": {
            "metric": [[1.0, 0.0], [0.0, 1.0]],
            "coordinates": ["x", "y"]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

#[test]
fn test_api_tensor_ricci() {
    let request = json!({
        "module": "tensor_calculus",
        "operation": "ricci_tensor",
        "parameters": {
            "metric": [[1.0, 0.0], [0.0, 1.0]],
            "coordinates": ["x", "y"]
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// === SYMBOLIC REGRESSION MODULE TESTS (8 operations) ===
#[test]
fn test_api_symbolic_fit() {
    let request = json!({
        "module": "symbolic_regression",
        "operation": "fit",
        "parameters": {
            "x": [1.0, 2.0, 3.0, 4.0, 5.0],
            "y": [2.0, 4.0, 6.0, 8.0, 10.0],
            "max_complexity": 5
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// === DIMENSIONAL ANALYSIS MODULE TESTS (5 operations) ===
#[test]
fn test_api_dimensional_check() {
    let request = json!({
        "module": "dimensional_analysis",
        "operation": "check_dimensions",
        "parameters": {
            "equation": "F = m * a",
            "units": {"F": "N", "m": "kg", "a": "m/s^2"}
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// === EQUATION VALIDATION MODULE TESTS (6 operations) ===
#[test]
fn test_api_equation_validate() {
    let request = json!({
        "module": "equation_validation",
        "operation": "validate",
        "parameters": {
            "equation": "E = m * c^2",
            "domain": "physics"
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// === FUNCTION APPROXIMATOR TESTS ===
#[test]
fn test_api_function_approximator() {
    let request = json!({
        "module": "function_approximator",
        "operation": "approximate",
        "parameters": {
            "inputs": [0.0, 1.0, 2.0, 3.0],
            "outputs": [0.0, 2.0, 4.0, 6.0],
            "population_size": 200,
            "generations": 50,
            "complexity": 2.0
        }
    });
    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();
    assert!(parsed["success"].as_bool().unwrap_or(false));
}

// Debug helper test
#[test]
fn debug_check_errors() {
    let tests = vec![
        (
            "statistics",
            "descriptive",
            json!({"data": [1.0, 2.0, 3.0]}),
        ),
        (
            "signal_processing",
            "fft",
            json!({"signal": [1.0, 2.0, 3.0, 4.0]}),
        ),
        (
            "optimization",
            "minimize",
            json!({"initial_guess": [0.0, 0.0]}),
        ),
    ];

    for (module, op, params) in tests {
        let request = json!({"module": module, "operation": op, "parameters": params});
        let response = process_json_request(&request.to_string());
        eprintln!("{}/{}: {}", module, op, response);
    }
}
