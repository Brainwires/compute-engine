//! Unit tests for api/operations.rs module
//!
//! Tests cover:
//! - Basic functionality of list_all_operations()
//! - Presence of all expected modules
//! - Correct operation counts per module
//! - Specific operation existence
//! - Module naming consistency
//! - Data structure integrity

use crate::api::operations::list_all_operations;

#[test]
fn test_list_all_operations_returns_hashmap() {
    let ops = list_all_operations();
    assert!(!ops.is_empty(), "Operations map should not be empty");
}

#[test]
fn test_all_expected_modules_present() {
    let ops = list_all_operations();

    let expected_modules = vec![
        "tensor_calculus",
        "advanced_calculus",
        "fluid_dynamics",
        "signal_processing",
        "stochastic_processes",
        "cryptographic_mathematics",
        "symbolic_regression",
        "dimensional_analysis",
        "equation_validation",
        "computational_geometry",
        "linear_algebra",
        "statistics",
        "optimization",
        "graph_theory",
        "information_theory",
        "electromagnetism",
        "special_functions",
        "numerical_methods",
        "chemistry",
    ];

    for module in expected_modules {
        assert!(
            ops.contains_key(module),
            "Module '{}' should be present in operations map",
            module
        );
    }
}

#[test]
fn test_module_count() {
    let ops = list_all_operations();
    assert_eq!(
        ops.len(),
        19,
        "Should have exactly 19 modules"
    );
}

#[test]
fn test_tensor_calculus_operations() {
    let ops = list_all_operations();
    let tensor_ops = ops.get("tensor_calculus").expect("tensor_calculus module should exist");

    assert_eq!(tensor_ops.len(), 14, "tensor_calculus should have 14 operations");
    assert!(tensor_ops.contains(&"solve_vacuum_einstein".to_string()));
    assert!(tensor_ops.contains(&"compute_riemann_tensor".to_string()));
    assert!(tensor_ops.contains(&"geodesic_equation".to_string()));
}

#[test]
fn test_advanced_calculus_operations() {
    let ops = list_all_operations();
    let calc_ops = ops.get("advanced_calculus").expect("advanced_calculus module should exist");

    assert_eq!(calc_ops.len(), 27, "advanced_calculus should have 27 operations");
    assert!(calc_ops.contains(&"fourier_transform".to_string()));
    assert!(calc_ops.contains(&"laplace_transform".to_string()));
    assert!(calc_ops.contains(&"greens_theorem".to_string()));
}

#[test]
fn test_fluid_dynamics_operations() {
    let ops = list_all_operations();
    let fluid_ops = ops.get("fluid_dynamics").expect("fluid_dynamics module should exist");

    assert_eq!(fluid_ops.len(), 22, "fluid_dynamics should have 22 operations");
    assert!(fluid_ops.contains(&"navier_stokes".to_string()));
    assert!(fluid_ops.contains(&"cavity_flow".to_string()));
    assert!(fluid_ops.contains(&"cfl_condition".to_string()));
}

#[test]
fn test_signal_processing_operations() {
    let ops = list_all_operations();
    let signal_ops = ops.get("signal_processing").expect("signal_processing module should exist");

    assert_eq!(signal_ops.len(), 13, "signal_processing should have 13 operations");
    assert!(signal_ops.contains(&"fft".to_string()));
    assert!(signal_ops.contains(&"wavelet_transform".to_string()));
    assert!(signal_ops.contains(&"autocorrelation".to_string()));
}

#[test]
fn test_stochastic_processes_operations() {
    let ops = list_all_operations();
    let stoch_ops = ops.get("stochastic_processes").expect("stochastic_processes module should exist");

    assert_eq!(stoch_ops.len(), 14, "stochastic_processes should have 14 operations");
    assert!(stoch_ops.contains(&"brownian_motion".to_string()));
    assert!(stoch_ops.contains(&"geometric_brownian_motion".to_string()));
    assert!(stoch_ops.contains(&"heston_model".to_string()));
}

#[test]
fn test_cryptographic_mathematics_operations() {
    let ops = list_all_operations();
    let crypto_ops = ops.get("cryptographic_mathematics").expect("cryptographic_mathematics module should exist");

    assert_eq!(crypto_ops.len(), 10, "cryptographic_mathematics should have 10 operations");
    assert!(crypto_ops.contains(&"generate_prime".to_string()));
    assert!(crypto_ops.contains(&"mod_exp".to_string()));
    assert!(crypto_ops.contains(&"elliptic_curve_point_add".to_string()));
}

#[test]
fn test_symbolic_regression_operations() {
    let ops = list_all_operations();
    let symb_ops = ops.get("symbolic_regression").expect("symbolic_regression module should exist");

    assert_eq!(symb_ops.len(), 8, "symbolic_regression should have 8 operations");
    assert!(symb_ops.contains(&"fit_polynomial".to_string()));
    assert!(symb_ops.contains(&"symbolic_differentiation".to_string()));
    assert!(symb_ops.contains(&"expression_simplification".to_string()));
}

#[test]
fn test_dimensional_analysis_operations() {
    let ops = list_all_operations();
    let dim_ops = ops.get("dimensional_analysis").expect("dimensional_analysis module should exist");

    assert_eq!(dim_ops.len(), 9, "dimensional_analysis should have 9 operations");
    assert!(dim_ops.contains(&"buckingham_pi".to_string()));
    assert!(dim_ops.contains(&"convert_units".to_string()));
    assert!(dim_ops.contains(&"validate_equation_dimensions".to_string()));
}

#[test]
fn test_equation_validation_operations() {
    let ops = list_all_operations();
    let eq_ops = ops.get("equation_validation").expect("equation_validation module should exist");

    assert_eq!(eq_ops.len(), 8, "equation_validation should have 8 operations");
    assert!(eq_ops.contains(&"validate_equation".to_string()));
    assert!(eq_ops.contains(&"check_conservation_laws".to_string()));
    assert!(eq_ops.contains(&"check_symmetries".to_string()));
}

#[test]
fn test_computational_geometry_operations() {
    let ops = list_all_operations();
    let geom_ops = ops.get("computational_geometry").expect("computational_geometry module should exist");

    assert_eq!(geom_ops.len(), 5, "computational_geometry should have 5 operations");
    assert!(geom_ops.contains(&"convex_hull".to_string()));
    assert!(geom_ops.contains(&"delaunay_triangulation".to_string()));
    assert!(geom_ops.contains(&"voronoi_diagram".to_string()));
}

#[test]
fn test_linear_algebra_operations() {
    let ops = list_all_operations();
    let linalg_ops = ops.get("linear_algebra").expect("linear_algebra module should exist");

    assert_eq!(linalg_ops.len(), 10, "linear_algebra should have 10 operations");
    assert!(linalg_ops.contains(&"compute_svd".to_string()));
    assert!(linalg_ops.contains(&"compute_eigendecomposition".to_string()));
    assert!(linalg_ops.contains(&"matrix_exp".to_string()));
}

#[test]
fn test_statistics_operations() {
    let ops = list_all_operations();
    let stats_ops = ops.get("statistics").expect("statistics module should exist");

    assert_eq!(stats_ops.len(), 10, "statistics should have 10 operations");
    assert!(stats_ops.contains(&"hypothesis_test".to_string()));
    assert!(stats_ops.contains(&"mcmc_sampling".to_string()));
    assert!(stats_ops.contains(&"mutual_information".to_string()));
}

#[test]
fn test_optimization_operations() {
    let ops = list_all_operations();
    let opt_ops = ops.get("optimization").expect("optimization module should exist");

    assert_eq!(opt_ops.len(), 3, "optimization should have 3 operations");
    assert!(opt_ops.contains(&"curve_fitting".to_string()));
    assert!(opt_ops.contains(&"gradient_descent".to_string()));
    assert!(opt_ops.contains(&"nelder_mead".to_string()));
}

#[test]
fn test_graph_theory_operations() {
    let ops = list_all_operations();
    let graph_ops = ops.get("graph_theory").expect("graph_theory module should exist");

    assert_eq!(graph_ops.len(), 5, "graph_theory should have 5 operations");
    assert!(graph_ops.contains(&"shortest_path".to_string()));
    assert!(graph_ops.contains(&"minimum_spanning_tree".to_string()));
    assert!(graph_ops.contains(&"topological_sort".to_string()));
}

#[test]
fn test_information_theory_operations() {
    let ops = list_all_operations();
    let info_ops = ops.get("information_theory").expect("information_theory module should exist");

    assert_eq!(info_ops.len(), 7, "information_theory should have 7 operations");
    assert!(info_ops.contains(&"shannon_entropy".to_string()));
    assert!(info_ops.contains(&"huffman_coding".to_string()));
    assert!(info_ops.contains(&"channel_capacity".to_string()));
}

#[test]
fn test_electromagnetism_operations() {
    let ops = list_all_operations();
    let em_ops = ops.get("electromagnetism").expect("electromagnetism module should exist");

    assert_eq!(em_ops.len(), 8, "electromagnetism should have 8 operations");
    assert!(em_ops.contains(&"maxwell_equations".to_string()));
    assert!(em_ops.contains(&"poynting_vector".to_string()));
    assert!(em_ops.contains(&"skin_effect".to_string()));
}

#[test]
fn test_special_functions_operations() {
    let ops = list_all_operations();
    let special_ops = ops.get("special_functions").expect("special_functions module should exist");

    assert_eq!(special_ops.len(), 6, "special_functions should have 6 operations");
    assert!(special_ops.contains(&"bessel".to_string()));
    assert!(special_ops.contains(&"gamma".to_string()));
    assert!(special_ops.contains(&"elliptic".to_string()));
}

#[test]
fn test_numerical_methods_operations() {
    let ops = list_all_operations();
    let num_ops = ops.get("numerical_methods").expect("numerical_methods module should exist");

    assert_eq!(num_ops.len(), 8, "numerical_methods should have 8 operations");
    assert!(num_ops.contains(&"solve_ode".to_string()));
    assert!(num_ops.contains(&"find_root".to_string()));
    assert!(num_ops.contains(&"solve_pde".to_string()));
}

#[test]
fn test_chemistry_operations() {
    let ops = list_all_operations();
    let chem_ops = ops.get("chemistry").expect("chemistry module should exist");

    assert_eq!(chem_ops.len(), 7, "chemistry should have 7 operations");
    assert!(chem_ops.contains(&"balance_equation".to_string()));
    assert!(chem_ops.contains(&"thermodynamics".to_string()));
    assert!(chem_ops.contains(&"molar_mass".to_string()));
}

#[test]
fn test_all_operation_values_are_strings() {
    let ops = list_all_operations();

    for (module, operations) in ops.iter() {
        for op in operations {
            assert!(
                !op.is_empty(),
                "Operation name should not be empty in module '{}'",
                module
            );
            assert!(
                op.chars().all(|c| c.is_alphanumeric() || c == '_'),
                "Operation '{}' in module '{}' should only contain alphanumeric characters and underscores",
                op,
                module
            );
        }
    }
}

#[test]
fn test_no_duplicate_operations_within_modules() {
    let ops = list_all_operations();

    for (module, operations) in ops.iter() {
        let mut seen = std::collections::HashSet::new();
        for op in operations {
            assert!(
                seen.insert(op),
                "Duplicate operation '{}' found in module '{}'",
                op,
                module
            );
        }
    }
}

#[test]
fn test_operations_immutability() {
    // Test that calling list_all_operations() multiple times returns consistent data
    let ops1 = list_all_operations();
    let ops2 = list_all_operations();

    assert_eq!(ops1.len(), ops2.len(), "Consecutive calls should return same number of modules");

    for (key, value1) in ops1.iter() {
        let value2 = ops2.get(key).expect(&format!("Module '{}' should exist in second call", key));
        assert_eq!(value1.len(), value2.len(), "Module '{}' should have same operation count", key);
        assert_eq!(value1, value2, "Module '{}' should have identical operations", key);
    }
}
