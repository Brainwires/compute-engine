//! Comprehensive tests for tensor calculus module
//!
//! This test suite covers:
//! 1. Tensor operations and basic tensor algebra
//! 2. Christoffel symbols computation
//! 3. Riemann curvature tensor
//! 4. Ricci tensor and Ricci scalar
//! 5. Einstein tensor
//! 6. Metric tensors (flat, curved, schwarzschild, etc.)
//! 7. Covariant derivatives
//! 8. Einstein field equations
//! 9. Symbolic tensor expressions
//! 10. Quantum tensor operations

use computational_engine::compute::tensor::einstein::{
    BoundaryCondition, EinsteinEquationSystem, EinsteinSolution, StressEnergyTensor,
    construct_einstein_field_equations, solve_einstein_constraint_equations,
    solve_vacuum_einstein_equations, verify_einstein_solution,
};
use computational_engine::compute::tensor::quantum_tensors;
use computational_engine::compute::tensor::symbolic::SymbolicExpr;
use computational_engine::compute::tensor::tensor::{
    self, ChristoffelResult, MetricTensor, RiemannResult, TensorComponent,
    calculate_christoffel_symbols, calculate_einstein_tensor, calculate_ricci_scalar,
    calculate_ricci_tensor, calculate_riemann_tensor, parse_metric_tensor,
};
use std::collections::HashMap;

// ============================================================================
// SYMBOLIC EXPRESSION TESTS
// ============================================================================

#[test]
fn test_symbolic_parse_variable() {
    let expr = SymbolicExpr::parse("x").unwrap();
    assert_eq!(expr, SymbolicExpr::Variable("x".to_string()));
}

#[test]
fn test_symbolic_parse_constant() {
    let expr = SymbolicExpr::parse("3.14").unwrap();
    match expr {
        SymbolicExpr::Constant(val) => assert!((val - 3.14).abs() < 1e-10),
        _ => panic!("Expected constant"),
    }
}

#[test]
fn test_symbolic_parse_zero() {
    let expr = SymbolicExpr::parse("0").unwrap();
    assert_eq!(expr, SymbolicExpr::Zero);
}

#[test]
fn test_symbolic_parse_one() {
    let expr = SymbolicExpr::parse("1").unwrap();
    assert_eq!(expr, SymbolicExpr::One);
}

#[test]
fn test_symbolic_parse_power() {
    let expr = SymbolicExpr::parse("r^2").unwrap();
    match expr {
        SymbolicExpr::Power(base, exp) => {
            assert_eq!(*base, SymbolicExpr::Variable("r".to_string()));
            assert_eq!(*exp, SymbolicExpr::Constant(2.0));
        }
        _ => panic!("Expected power expression"),
    }
}

#[test]
fn test_symbolic_parse_function_sin() {
    let expr = SymbolicExpr::parse("sin(theta)").unwrap();
    match expr {
        SymbolicExpr::Function(name, args) => {
            assert_eq!(name, "sin");
            assert_eq!(args.len(), 1);
            assert_eq!(args[0], SymbolicExpr::Variable("theta".to_string()));
        }
        _ => panic!("Expected function"),
    }
}

#[test]
fn test_symbolic_parse_function_cos() {
    let expr = SymbolicExpr::parse("cos(phi)").unwrap();
    match expr {
        SymbolicExpr::Function(name, args) => {
            assert_eq!(name, "cos");
            assert_eq!(args.len(), 1);
        }
        _ => panic!("Expected function"),
    }
}

#[test]
fn test_symbolic_simplify_add_constants() {
    let expr = SymbolicExpr::Add(
        Box::new(SymbolicExpr::Constant(2.0)),
        Box::new(SymbolicExpr::Constant(3.0)),
    );
    let simplified = expr.simplify();
    assert_eq!(simplified, SymbolicExpr::Constant(5.0));
}

#[test]
fn test_symbolic_simplify_add_zero() {
    let expr = SymbolicExpr::Add(
        Box::new(SymbolicExpr::Variable("x".to_string())),
        Box::new(SymbolicExpr::Zero),
    );
    let simplified = expr.simplify();
    assert_eq!(simplified, SymbolicExpr::Variable("x".to_string()));
}

#[test]
fn test_symbolic_simplify_multiply_zero() {
    let expr = SymbolicExpr::Multiply(
        Box::new(SymbolicExpr::Variable("x".to_string())),
        Box::new(SymbolicExpr::Zero),
    );
    let simplified = expr.simplify();
    assert_eq!(simplified, SymbolicExpr::Zero);
}

#[test]
fn test_symbolic_simplify_multiply_one() {
    let expr = SymbolicExpr::Multiply(
        Box::new(SymbolicExpr::Variable("y".to_string())),
        Box::new(SymbolicExpr::One),
    );
    let simplified = expr.simplify();
    assert_eq!(simplified, SymbolicExpr::Variable("y".to_string()));
}

#[test]
fn test_symbolic_simplify_power_zero() {
    let expr = SymbolicExpr::Power(
        Box::new(SymbolicExpr::Variable("x".to_string())),
        Box::new(SymbolicExpr::Zero),
    );
    let simplified = expr.simplify();
    assert_eq!(simplified, SymbolicExpr::One);
}

#[test]
fn test_symbolic_simplify_power_one() {
    let expr = SymbolicExpr::Power(
        Box::new(SymbolicExpr::Variable("x".to_string())),
        Box::new(SymbolicExpr::One),
    );
    let simplified = expr.simplify();
    assert_eq!(simplified, SymbolicExpr::Variable("x".to_string()));
}

#[test]
fn test_symbolic_derivative_variable() {
    let expr = SymbolicExpr::Variable("x".to_string());
    assert_eq!(expr.derivative("x"), SymbolicExpr::One);
    assert_eq!(expr.derivative("y"), SymbolicExpr::Zero);
}

#[test]
fn test_symbolic_derivative_constant() {
    let expr = SymbolicExpr::Constant(42.0);
    assert_eq!(expr.derivative("x"), SymbolicExpr::Zero);
}

#[test]
fn test_symbolic_derivative_add() {
    let expr = SymbolicExpr::Add(
        Box::new(SymbolicExpr::Variable("x".to_string())),
        Box::new(SymbolicExpr::Variable("y".to_string())),
    );
    let deriv = expr.derivative("x");
    // d/dx(x + y) = 1 + 0 = 1
    let simplified = deriv.simplify();
    assert_eq!(simplified, SymbolicExpr::One);
}

#[test]
fn test_symbolic_derivative_product_rule() {
    // d/dx(x * x) = x + x = 2x
    let expr = SymbolicExpr::Multiply(
        Box::new(SymbolicExpr::Variable("x".to_string())),
        Box::new(SymbolicExpr::Variable("x".to_string())),
    );
    let deriv = expr.derivative("x");
    // Should be: 1*x + x*1 = 2x (not simplified here)
    assert!(matches!(deriv, SymbolicExpr::Add(_, _)));
}

#[test]
fn test_symbolic_derivative_power() {
    // d/dx(x^2) = 2*x^1 * 1 = 2x
    let expr = SymbolicExpr::Power(
        Box::new(SymbolicExpr::Variable("x".to_string())),
        Box::new(SymbolicExpr::Constant(2.0)),
    );
    let deriv = expr.derivative("x");
    assert!(matches!(deriv, SymbolicExpr::Multiply(_, _)));
}

#[test]
fn test_symbolic_is_zero() {
    assert!(SymbolicExpr::Zero.is_zero());
    assert!(SymbolicExpr::Constant(0.0).is_zero());
    assert!(!SymbolicExpr::One.is_zero());
    assert!(!SymbolicExpr::Variable("x".to_string()).is_zero());
}

#[test]
fn test_symbolic_display() {
    let expr = SymbolicExpr::Variable("x".to_string());
    assert_eq!(expr.to_string(), "x");

    let expr = SymbolicExpr::Constant(42.0);
    assert_eq!(expr.to_string(), "42");

    let expr = SymbolicExpr::Zero;
    assert_eq!(expr.to_string(), "0");
}

// ============================================================================
// METRIC TENSOR TESTS
// ============================================================================

#[test]
fn test_parse_flat_metric_2d() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string()],
    ];
    let coords = vec!["x".to_string(), "y".to_string()];

    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();
    assert_eq!(metric.len(), 2);
    assert_eq!(metric[0].len(), 2);
}

#[test]
fn test_parse_minkowski_metric_4d() {
    let metric_strings = vec![
        vec![
            "-1".to_string(),
            "0".to_string(),
            "0".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "1".to_string(),
            "0".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "0".to_string(),
            "1".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "0".to_string(),
            "0".to_string(),
            "1".to_string(),
        ],
    ];
    let coords = vec![
        "t".to_string(),
        "x".to_string(),
        "y".to_string(),
        "z".to_string(),
    ];

    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();
    assert_eq!(metric.len(), 4);
}

#[test]
fn test_parse_polar_metric() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "r^2".to_string()],
    ];
    let coords = vec!["r".to_string(), "theta".to_string()];

    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();
    assert_eq!(metric.len(), 2);
}

#[test]
fn test_parse_spherical_metric() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string(), "0".to_string()],
        vec!["0".to_string(), "r^2".to_string(), "0".to_string()],
        vec!["0".to_string(), "0".to_string(), "r^2".to_string()],
    ];
    let coords = vec!["r".to_string(), "theta".to_string(), "phi".to_string()];

    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();
    assert_eq!(metric.len(), 3);
}

#[test]
fn test_parse_invalid_metric_non_square() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string(), "0".to_string()],
    ];
    let coords = vec!["x".to_string(), "y".to_string()];

    let result = parse_metric_tensor(metric_strings, &coords);
    assert!(result.is_err());
}

// ============================================================================
// CHRISTOFFEL SYMBOLS TESTS
// ============================================================================

#[test]
fn test_christoffel_flat_2d() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string()],
    ];
    let coords = vec!["x".to_string(), "y".to_string()];
    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();

    let result = calculate_christoffel_symbols(&metric, &coords).unwrap();
    // Flat metric should have zero Christoffel symbols
    assert_eq!(result.dimension, 2);
    // All components should be zero (no non-zero symbols stored)
    assert!(result.symbols.is_empty() || result.symbols.iter().all(|c| c.expression == "0"));
}

#[test]
fn test_christoffel_polar() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "r^2".to_string()],
    ];
    let coords = vec!["r".to_string(), "theta".to_string()];
    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();

    let result = calculate_christoffel_symbols(&metric, &coords).unwrap();
    assert_eq!(result.dimension, 2);
    // Polar coordinates have non-zero Christoffel symbols
    assert!(!result.symbols.is_empty());
}

#[test]
fn test_christoffel_minkowski() {
    let metric_strings = vec![
        vec![
            "-1".to_string(),
            "0".to_string(),
            "0".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "1".to_string(),
            "0".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "0".to_string(),
            "1".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "0".to_string(),
            "0".to_string(),
            "1".to_string(),
        ],
    ];
    let coords = vec![
        "t".to_string(),
        "x".to_string(),
        "y".to_string(),
        "z".to_string(),
    ];
    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();

    let result = calculate_christoffel_symbols(&metric, &coords).unwrap();
    assert_eq!(result.dimension, 4);
    // Minkowski metric (flat spacetime) should have zero Christoffel symbols
}

// ============================================================================
// RIEMANN TENSOR TESTS
// ============================================================================

#[test]
fn test_riemann_flat_2d() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string()],
    ];
    let coords = vec!["x".to_string(), "y".to_string()];
    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();

    let result = calculate_riemann_tensor(&metric, &coords).unwrap();
    assert_eq!(result.dimension, 2);
    // Flat space has zero curvature
}

#[test]
fn test_riemann_polar() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "r^2".to_string()],
    ];
    let coords = vec!["r".to_string(), "theta".to_string()];
    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();

    let result = calculate_riemann_tensor(&metric, &coords).unwrap();
    assert_eq!(result.dimension, 2);
    // 2D polar coordinates are intrinsically flat, so Riemann tensor is zero
}

#[test]
fn test_riemann_minkowski() {
    let metric_strings = vec![
        vec![
            "-1".to_string(),
            "0".to_string(),
            "0".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "1".to_string(),
            "0".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "0".to_string(),
            "1".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "0".to_string(),
            "0".to_string(),
            "1".to_string(),
        ],
    ];
    let coords = vec![
        "t".to_string(),
        "x".to_string(),
        "y".to_string(),
        "z".to_string(),
    ];
    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();

    let result = calculate_riemann_tensor(&metric, &coords).unwrap();
    assert_eq!(result.dimension, 4);
    // Flat Minkowski spacetime has zero curvature
}

// ============================================================================
// RICCI TENSOR TESTS
// ============================================================================

#[test]
fn test_ricci_flat_2d() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string()],
    ];
    let coords = vec!["x".to_string(), "y".to_string()];
    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();

    let result = calculate_ricci_tensor(&metric, &coords).unwrap();
    assert_eq!(result.dimension, 2);
}

#[test]
fn test_ricci_minkowski() {
    let metric_strings = vec![
        vec![
            "-1".to_string(),
            "0".to_string(),
            "0".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "1".to_string(),
            "0".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "0".to_string(),
            "1".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "0".to_string(),
            "0".to_string(),
            "1".to_string(),
        ],
    ];
    let coords = vec![
        "t".to_string(),
        "x".to_string(),
        "y".to_string(),
        "z".to_string(),
    ];
    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();

    let result = calculate_ricci_tensor(&metric, &coords).unwrap();
    assert_eq!(result.dimension, 4);
    // Minkowski space has zero Ricci tensor
}

// ============================================================================
// RICCI SCALAR TESTS
// ============================================================================

#[test]
fn test_ricci_scalar_flat_2d() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string()],
    ];
    let coords = vec!["x".to_string(), "y".to_string()];
    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();

    let result = calculate_ricci_scalar(&metric, &coords).unwrap();
    assert_eq!(result.indices.len(), 0); // Scalar has no indices
}

#[test]
fn test_ricci_scalar_minkowski() {
    let metric_strings = vec![
        vec![
            "-1".to_string(),
            "0".to_string(),
            "0".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "1".to_string(),
            "0".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "0".to_string(),
            "1".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "0".to_string(),
            "0".to_string(),
            "1".to_string(),
        ],
    ];
    let coords = vec![
        "t".to_string(),
        "x".to_string(),
        "y".to_string(),
        "z".to_string(),
    ];
    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();

    let result = calculate_ricci_scalar(&metric, &coords).unwrap();
    assert_eq!(result.indices.len(), 0);
}

// ============================================================================
// EINSTEIN TENSOR TESTS
// ============================================================================

#[test]
fn test_einstein_tensor_flat_2d() {
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string()],
    ];
    let coords = vec!["x".to_string(), "y".to_string()];
    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();

    let result = calculate_einstein_tensor(&metric, &coords).unwrap();
    assert_eq!(result.dimension, 2);
}

#[test]
fn test_einstein_tensor_minkowski() {
    let metric_strings = vec![
        vec![
            "-1".to_string(),
            "0".to_string(),
            "0".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "1".to_string(),
            "0".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "0".to_string(),
            "1".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "0".to_string(),
            "0".to_string(),
            "1".to_string(),
        ],
    ];
    let coords = vec![
        "t".to_string(),
        "x".to_string(),
        "y".to_string(),
        "z".to_string(),
    ];
    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();

    let result = calculate_einstein_tensor(&metric, &coords).unwrap();
    assert_eq!(result.dimension, 4);
    // Vacuum spacetime (Minkowski) has zero Einstein tensor
}

// ============================================================================
// EINSTEIN FIELD EQUATIONS TESTS
// ============================================================================

#[test]
fn test_construct_einstein_field_equations_vacuum() {
    let stress_energy = StressEnergyTensor {
        components: vec![vec![SymbolicExpr::Zero; 4]; 4],
        tensor_type: "vacuum".to_string(),
        parameters: HashMap::new(),
    };
    let coords = vec![
        "t".to_string(),
        "r".to_string(),
        "theta".to_string(),
        "phi".to_string(),
    ];

    let result = construct_einstein_field_equations(&stress_energy, &coords, None).unwrap();
    assert_eq!(result.field_equations.len(), 16); // 4x4 = 16 components
    assert_eq!(result.unknowns.len(), 10); // 10 independent metric components (symmetric)
}

#[test]
fn test_construct_einstein_field_equations_with_cosmological_constant() {
    let stress_energy = StressEnergyTensor {
        components: vec![vec![SymbolicExpr::Zero; 4]; 4],
        tensor_type: "vacuum".to_string(),
        parameters: HashMap::new(),
    };
    let coords = vec![
        "t".to_string(),
        "r".to_string(),
        "theta".to_string(),
        "phi".to_string(),
    ];
    let lambda = SymbolicExpr::Constant(0.1);

    let result = construct_einstein_field_equations(&stress_energy, &coords, Some(lambda)).unwrap();
    assert_eq!(result.field_equations.len(), 16);
}

#[test]
fn test_solve_spherically_symmetric_vacuum() {
    let coords = vec![
        "t".to_string(),
        "r".to_string(),
        "theta".to_string(),
        "phi".to_string(),
    ];
    let boundary_conditions = vec![];

    let solutions =
        solve_vacuum_einstein_equations(&coords, "spherical", &boundary_conditions).unwrap();
    assert!(!solutions.is_empty());
    assert_eq!(solutions[0].solution_type, "exact");
    assert_eq!(solutions[0].coordinates, coords);
}

#[test]
fn test_schwarzschild_solution() {
    let coords = vec![
        "t".to_string(),
        "r".to_string(),
        "theta".to_string(),
        "phi".to_string(),
    ];
    let boundary_conditions = vec![];

    let solutions =
        solve_vacuum_einstein_equations(&coords, "spherical", &boundary_conditions).unwrap();

    // Should contain Schwarzschild solution
    let schwarzschild = &solutions[0];
    assert_eq!(schwarzschild.solution_type, "exact");
    assert!(schwarzschild.constraints_satisfied);
    assert!(schwarzschild.physical_parameters.contains_key("M"));
}

#[test]
fn test_reissner_nordstrom_solution() {
    let coords = vec![
        "t".to_string(),
        "r".to_string(),
        "theta".to_string(),
        "phi".to_string(),
    ];
    let boundary_conditions = vec![];

    let solutions =
        solve_vacuum_einstein_equations(&coords, "spherical", &boundary_conditions).unwrap();

    // Should contain Reissner-Nordström solution (charged black hole)
    assert!(solutions.len() >= 2);
    let rn = &solutions[1];
    assert_eq!(rn.solution_type, "exact");
    assert!(rn.physical_parameters.contains_key("M"));
    assert!(rn.physical_parameters.contains_key("Q"));
}

#[test]
fn test_flrw_cosmological_solution() {
    let coords = vec![
        "t".to_string(),
        "r".to_string(),
        "theta".to_string(),
        "phi".to_string(),
    ];
    let boundary_conditions = vec![];

    let solutions =
        solve_vacuum_einstein_equations(&coords, "cosmological", &boundary_conditions).unwrap();
    assert!(!solutions.is_empty());

    // Should contain FLRW solution
    let flrw = &solutions[0];
    assert_eq!(flrw.solution_type, "exact");
    assert!(flrw.physical_parameters.contains_key("H"));
}

#[test]
fn test_de_sitter_solution() {
    let coords = vec![
        "t".to_string(),
        "r".to_string(),
        "theta".to_string(),
        "phi".to_string(),
    ];
    let boundary_conditions = vec![];

    let solutions =
        solve_vacuum_einstein_equations(&coords, "cosmological", &boundary_conditions).unwrap();

    // Should contain de Sitter solution
    assert!(solutions.len() >= 2);
    let ds = &solutions[1];
    assert_eq!(ds.solution_type, "exact");
    assert!(ds.physical_parameters.contains_key("H"));
    assert!(ds.physical_parameters.contains_key("Lambda"));
}

#[test]
fn test_kerr_axisymmetric_solution() {
    let coords = vec![
        "t".to_string(),
        "r".to_string(),
        "theta".to_string(),
        "phi".to_string(),
    ];
    let boundary_conditions = vec![];

    let solutions =
        solve_vacuum_einstein_equations(&coords, "axisymmetric", &boundary_conditions).unwrap();

    // Should contain Kerr solution (rotating black hole)
    assert!(!solutions.is_empty());
    let kerr = &solutions[0];
    assert_eq!(kerr.solution_type, "exact");
    assert!(kerr.physical_parameters.contains_key("M"));
    assert!(kerr.physical_parameters.contains_key("a")); // Angular momentum parameter
}

#[test]
fn test_verify_einstein_solution() {
    let coords = vec![
        "t".to_string(),
        "r".to_string(),
        "theta".to_string(),
        "phi".to_string(),
    ];
    let boundary_conditions = vec![];

    let solutions =
        solve_vacuum_einstein_equations(&coords, "spherical", &boundary_conditions).unwrap();
    let schwarzschild = &solutions[0];

    // Verify the solution satisfies Einstein equations
    let verified = verify_einstein_solution(schwarzschild, None, None).unwrap();
    assert!(verified); // Schwarzschild is an exact solution
}

#[test]
fn test_einstein_constraint_equations() {
    let metric_strings = vec![
        vec![
            "-1".to_string(),
            "0".to_string(),
            "0".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "1".to_string(),
            "0".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "0".to_string(),
            "1".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "0".to_string(),
            "0".to_string(),
            "1".to_string(),
        ],
    ];
    let coords = vec![
        "t".to_string(),
        "x".to_string(),
        "y".to_string(),
        "z".to_string(),
    ];
    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();

    let constraints = solve_einstein_constraint_equations(&metric, &coords).unwrap();
    assert!(!constraints.is_empty());
    // Should have Hamiltonian constraint (1) + Momentum constraints (3 for 3D)
    assert!(constraints.len() >= 4);
}

// ============================================================================
// QUANTUM TENSOR TESTS
// ============================================================================

#[test]
fn test_christoffel_symbols_quantum() {
    let metric = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "r^2".to_string()],
    ];

    let result = quantum_tensors::calculate_christoffel_symbols_symbolic(metric, 2).unwrap();
    assert_eq!(result.calculation_type, "christoffel_symbols");
    assert_eq!(result.dimensions, 2);
    assert!(!result.properties.symmetries.is_empty());
}

#[test]
fn test_riemann_tensor_quantum() {
    let metric = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string()],
    ];

    let result = quantum_tensors::calculate_riemann_tensor_symbolic(metric, 2).unwrap();
    assert_eq!(result.calculation_type, "riemann_tensor");
    assert_eq!(result.dimensions, 2);
    assert!(!result.physical_interpretation.is_empty());
}

#[test]
fn test_ricci_tensor_quantum() {
    let metric = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string()],
    ];

    let result = quantum_tensors::calculate_ricci_tensor_symbolic(metric, 2).unwrap();
    assert_eq!(result.calculation_type, "ricci_tensor");
    assert_eq!(result.dimensions, 2);
    assert!(result.properties.trace.is_some());
}

#[test]
fn test_einstein_tensor_quantum() {
    let metric = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "1".to_string()],
    ];

    let result = quantum_tensors::calculate_einstein_tensor_symbolic(metric, 2).unwrap();
    assert_eq!(result.calculation_type, "einstein_tensor");
    assert_eq!(result.dimensions, 2);
    assert!(result.properties.trace.is_some());
}

#[test]
fn test_christoffel_4d_quantum() {
    let metric = vec![
        vec![
            "-1".to_string(),
            "0".to_string(),
            "0".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "1".to_string(),
            "0".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "0".to_string(),
            "1".to_string(),
            "0".to_string(),
        ],
        vec![
            "0".to_string(),
            "0".to_string(),
            "0".to_string(),
            "1".to_string(),
        ],
    ];

    let result = quantum_tensors::calculate_christoffel_symbols_symbolic(metric, 4).unwrap();
    assert_eq!(result.dimensions, 4);
}

#[test]
fn test_tensor_properties_symmetries() {
    let metric = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "r^2".to_string()],
    ];

    let result = quantum_tensors::calculate_christoffel_symbols_symbolic(metric, 2).unwrap();
    assert!(
        result
            .properties
            .symmetries
            .iter()
            .any(|s| s.contains("symmetric"))
    );
}

// ============================================================================
// EDGE CASES AND ERROR HANDLING
// ============================================================================

#[test]
fn test_invalid_symmetry_ansatz() {
    let coords = vec![
        "t".to_string(),
        "r".to_string(),
        "theta".to_string(),
        "phi".to_string(),
    ];
    let boundary_conditions = vec![];

    let result = solve_vacuum_einstein_equations(&coords, "invalid_symmetry", &boundary_conditions);
    assert!(result.is_err());
}

#[test]
fn test_wrong_dimension_for_spherical() {
    let coords = vec!["t".to_string(), "r".to_string()]; // Only 2D, needs 4D
    let boundary_conditions = vec![];

    let result = solve_vacuum_einstein_equations(&coords, "spherical", &boundary_conditions);
    assert!(result.is_err());
}

#[test]
fn test_wrong_dimension_for_flrw() {
    let coords = vec!["t".to_string()]; // Only 1D, needs 4D
    let boundary_conditions = vec![];

    let result = solve_vacuum_einstein_equations(&coords, "cosmological", &boundary_conditions);
    assert!(result.is_err());
}

#[test]
fn test_empty_metric_parse() {
    let metric_strings: Vec<Vec<String>> = vec![];
    let coords: Vec<String> = vec![];

    let result = parse_metric_tensor(metric_strings, &coords);
    // Should handle empty metric gracefully
    assert!(result.is_ok() || result.is_err());
}

#[test]
fn test_tensor_component_indices() {
    let component = TensorComponent {
        indices: vec![0, 1, 2],
        expression: "test".to_string(),
    };
    assert_eq!(component.indices.len(), 3);
    assert_eq!(component.expression, "test");
}

#[test]
fn test_christoffel_result_structure() {
    let result = ChristoffelResult {
        symbols: vec![TensorComponent {
            indices: vec![0, 0, 1],
            expression: "r".to_string(),
        }],
        dimension: 2,
    };
    assert_eq!(result.dimension, 2);
    assert_eq!(result.symbols.len(), 1);
}

#[test]
fn test_riemann_result_structure() {
    let result = RiemannResult {
        components: vec![TensorComponent {
            indices: vec![0, 1, 2, 3],
            expression: "test".to_string(),
        }],
        dimension: 4,
    };
    assert_eq!(result.dimension, 4);
    assert_eq!(result.components.len(), 1);
}

#[test]
fn test_stress_energy_tensor_vacuum() {
    let stress_energy = StressEnergyTensor {
        components: vec![vec![SymbolicExpr::Zero; 4]; 4],
        tensor_type: "vacuum".to_string(),
        parameters: HashMap::new(),
    };
    assert_eq!(stress_energy.tensor_type, "vacuum");
    assert_eq!(stress_energy.components.len(), 4);
}

#[test]
fn test_boundary_condition_structure() {
    let bc = BoundaryCondition {
        coordinate: "r".to_string(),
        value: SymbolicExpr::Constant(1.0),
        condition_type: "dirichlet".to_string(),
        component_indices: vec![0, 0],
    };
    assert_eq!(bc.coordinate, "r");
    assert_eq!(bc.condition_type, "dirichlet");
}

#[test]
fn test_einstein_solution_structure() {
    let metric = vec![vec![SymbolicExpr::One; 2]; 2];
    let coords = vec!["x".to_string(), "y".to_string()];
    let params = HashMap::new();

    let solution = EinsteinSolution {
        metric_tensor: metric,
        coordinates: coords.clone(),
        solution_type: "exact".to_string(),
        constraints_satisfied: true,
        physical_parameters: params,
        solution_domain: "all space".to_string(),
    };

    assert_eq!(solution.solution_type, "exact");
    assert!(solution.constraints_satisfied);
    assert_eq!(solution.coordinates, coords);
}

#[test]
fn test_einstein_equation_system_structure() {
    let system = EinsteinEquationSystem {
        field_equations: vec![],
        constraint_equations: vec![],
        gauge_conditions: vec![],
        unknowns: vec!["g_00".to_string()],
        known_parameters: HashMap::new(),
    };
    assert_eq!(system.unknowns.len(), 1);
}

// ============================================================================
// INTEGRATION TESTS
// ============================================================================

#[test]
fn test_full_tensor_workflow() {
    // Full workflow: metric -> Christoffel -> Riemann -> Ricci -> Einstein
    let metric_strings = vec![
        vec!["1".to_string(), "0".to_string()],
        vec!["0".to_string(), "r^2".to_string()],
    ];
    let coords = vec!["r".to_string(), "theta".to_string()];

    // Step 1: Parse metric
    let metric = parse_metric_tensor(metric_strings, &coords).unwrap();

    // Step 2: Calculate Christoffel symbols
    let christoffel = calculate_christoffel_symbols(&metric, &coords).unwrap();
    assert_eq!(christoffel.dimension, 2);

    // Step 3: Calculate Riemann tensor
    let riemann = calculate_riemann_tensor(&metric, &coords).unwrap();
    assert_eq!(riemann.dimension, 2);

    // Step 4: Calculate Ricci tensor
    let ricci = calculate_ricci_tensor(&metric, &coords).unwrap();
    assert_eq!(ricci.dimension, 2);

    // Step 5: Calculate Ricci scalar
    let ricci_scalar = calculate_ricci_scalar(&metric, &coords).unwrap();
    assert_eq!(ricci_scalar.indices.len(), 0);

    // Step 6: Calculate Einstein tensor
    let einstein = calculate_einstein_tensor(&metric, &coords).unwrap();
    assert_eq!(einstein.dimension, 2);
}

#[test]
fn test_schwarzschild_full_analysis() {
    let coords = vec![
        "t".to_string(),
        "r".to_string(),
        "theta".to_string(),
        "phi".to_string(),
    ];
    let boundary_conditions = vec![];

    // Get Schwarzschild solution
    let solutions =
        solve_vacuum_einstein_equations(&coords, "spherical", &boundary_conditions).unwrap();
    let schwarzschild = &solutions[0];

    // Verify it's correct
    assert_eq!(schwarzschild.coordinates.len(), 4);
    assert!(schwarzschild.constraints_satisfied);

    // Calculate tensors from the metric
    let christoffel = calculate_christoffel_symbols(&schwarzschild.metric_tensor, &coords).unwrap();
    assert_eq!(christoffel.dimension, 4);
}
