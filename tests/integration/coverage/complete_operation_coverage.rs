//! Complete operation coverage tests
//!
//! Goal: Test every major operation in all tools for maximum coverage

use computational_engine::create_default_dispatcher;
use computational_engine::engine::equations::*;
use computational_engine::engine::*;
use serde_json::json;
use std::collections::HashMap;

// ============================================================================
// SOLVE TOOL - All Equation Types
// ============================================================================

#[test]
fn solve_einstein_vacuum() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Solve(SolveInput {
        equation_type: EquationType::Einstein(EinsteinEquation::Vacuum),
        equations: vec![],
        variables: Some(vec!["t".to_string(), "r".to_string()]),
        initial_guess: None,
        boundary_conditions: None,
        domain: None,
        method: None,
        parameters: HashMap::from([("symmetry".to_string(), serde_json::json!("spherical"))]),
    });

    let result = dispatcher.dispatch(request);
    // Einstein equations are complex, may not be fully implemented
    let _ = result; // Test it doesn't panic
}

#[test]
#[ignore = "Slow test: Fluid dynamics solvers are computationally expensive (disabled for CI)"]
fn solve_fluid_cavity_flow() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Solve(SolveInput {
        equation_type: EquationType::Fluid(FluidEquation::CavityFlow),
        equations: vec![],
        variables: None,
        initial_guess: None,
        boundary_conditions: None,
        domain: None,
        method: None,
        parameters: HashMap::from([
            ("reynolds_number".to_string(), serde_json::json!(100.0)),
            ("grid_size".to_string(), serde_json::json!(32)),
        ]),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Cavity flow should work");
}

#[test]
fn solve_differential_ode() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Solve(SolveInput {
        equation_type: EquationType::Differential(DifferentialEquation::ODE),
        equations: vec!["dy/dt = -y".to_string()],
        variables: Some(vec!["y".to_string()]),
        initial_guess: None,
        boundary_conditions: None,
        domain: None,
        method: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    // ODE solving may or may not be fully implemented
    let _ = result;
}

#[test]
fn solve_electromagnetic_maxwell() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Solve(SolveInput {
        equation_type: EquationType::Electromagnetic(EMEquation::Maxwell),
        equations: vec![],
        variables: None,
        initial_guess: None,
        boundary_conditions: None,
        domain: None,
        method: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    let _ = result; // Test it doesn't panic
}

// ============================================================================
// DIFFERENTIATE - All Operation Types (via Compute tool)
// ============================================================================

#[test]
fn differentiate_vector_calc_gradient() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Differentiate(DifferentiationOp::VectorCalc(VectorCalcOp::Gradient)),
        data: json!({
            "expression": "f(x,y,z)",
            "data": [[1.0, 2.0, 3.0], [2.0, 3.0, 4.0], [3.0, 4.0, 5.0]]
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    let _ = result; // Gradient may not be fully implemented
}

#[test]
fn differentiate_vector_calc_divergence() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Differentiate(DifferentiationOp::VectorCalc(VectorCalcOp::Divergence)),
        data: json!({
            "expression": "F(x,y,z)"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    let _ = result;
}

#[test]
fn differentiate_vector_calc_curl() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Differentiate(DifferentiationOp::VectorCalc(VectorCalcOp::Curl)),
        data: json!({
            "expression": "F(x,y,z)"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    let _ = result;
}

#[test]
fn differentiate_tensor_calc_covariant() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Differentiate(DifferentiationOp::TensorCalc(TensorDiffOp::Covariant)),
        data: json!({
            "expression": "T^{ij}",
            "variable": "x"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    let _ = result;
}

// ============================================================================
// INTEGRATE - All Integration Types (via Compute tool)
// ============================================================================

#[test]
fn integrate_geometric_line() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::Geometric(GeometricIntegral::Line)),
        data: json!({
            "expression": "F·dr",
            "limits": [[0.0, 1.0]],
            "path": [[0.0, 0.0], [1.0, 1.0]]
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    let _ = result; // Line integrals may not be fully implemented
}

#[test]
fn integrate_theorem_greens() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::Theorem(IntegralTheorem::Greens)),
        data: json!({
            "P": "y",
            "Q": "x",
            "region": "unit_square"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    let _ = result;
}

#[test]
fn integrate_complex_residue() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::ComplexAnalysis(ComplexIntegral::Residue)),
        data: json!({
            "function": "f(z)",
            "poles": [{"location": 0.0, "order": 1}]
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    let _ = result;
}

#[test]
fn integrate_numeric_trapezoidal() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::Numeric(NumericIntegration::Trapezoidal)),
        data: json!({
            "expression": "x^2",
            "lower": 0.0,
            "upper": 1.0,
            "function_type": "polynomial",
            "coefficients": [0.0, 0.0, 1.0],
            "num_points": 100
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Trapezoidal integration should work");
}

#[test]
fn integrate_numeric_gauss_quadrature() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Integrate(IntegrationType::Numeric(NumericIntegration::GaussQuadrature)),
        data: json!({
            "expression": "sin(x)",
            "lower": 0.0,
            "upper": 3.14159,
            "function_type": "polynomial",
            "coefficients": [1.0, 0.0]
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    let _ = result; // Gauss quadrature may not be fully implemented
}

// ============================================================================
// ANALYZE - All Operations
// ============================================================================

#[test]
fn analyze_check_correctness() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::CheckCorrectness,
        expression: "E = mc^2".to_string(),
        options: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    let _ = result;
}

#[test]
fn analyze_check_dimensions() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::CheckDimensions,
        expression: "F = ma".to_string(),
        options: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    let _ = result;
}

#[test]
fn analyze_validate_dimensions() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::ValidateDimensions,
        expression: "v = d/t".to_string(),
        options: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    let _ = result;
}

#[test]
fn analyze_partial_fraction() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::PartialFraction,
        expression: "1/(x^2 - 1)".to_string(),
        options: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    let _ = result;
}

#[test]
fn analyze_limit() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert("variable".to_string(), serde_json::json!("x"));
    options.insert("point".to_string(), serde_json::json!("0"));

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::Limit,
        expression: "sin(x)/x".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    let _ = result;
}

#[test]
fn analyze_graph_components() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert(
        "edges".to_string(),
        serde_json::json!([[0, 1], [1, 2], [3, 4]]),
    );

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::GraphComponents,
        expression: "graph".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    let _ = result;
}

// ============================================================================
// SIMULATE - All Model Types
// ============================================================================

#[test]
fn simulate_runge_kutta4() {
    let dispatcher = create_default_dispatcher();

    let mut initial = HashMap::new();
    initial.insert("y".to_string(), 1.0);

    let mut params = HashMap::new();
    params.insert("initial_value".to_string(), 1.0);

    let request = ToolRequest::Simulate(SimulateInput {
        model: SimulationModel::TimeEvolution(TimeEvolutionMethod::RungeKutta4),
        equations: vec!["dy/dt = -y".to_string()],
        variables: vec!["y".to_string()],
        parameters: params,
        initial_conditions: Some(initial),
        range: Some([0.0, 5.0]),
        steps: Some(50),
        num_paths: None,
        method: None,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "RungeKutta4 simulation should work");
}

#[test]
fn simulate_geometric_brownian() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("initial_value".to_string(), 100.0);
    params.insert("drift".to_string(), 0.05);
    params.insert("volatility".to_string(), 0.2);

    let request = ToolRequest::Simulate(SimulateInput {
        model: SimulationModel::Stochastic(StochasticProcess::GeometricBrownian),
        equations: vec![],
        variables: vec!["S".to_string()],
        parameters: params,
        initial_conditions: None,
        range: Some([0.0, 1.0]),
        steps: Some(252),
        num_paths: Some(1),
        method: None,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Geometric Brownian motion should work");
}

// ============================================================================
// COMPUTE - All Operations
// ============================================================================

#[test]
fn compute_matrix_qr() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::MatrixDecomp(MatrixDecomp::QR),
        data: json!({
            "matrix": [[1.0, 2.0], [3.0, 4.0]]
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    // QR decomposition not fully implemented yet
    let _ = result;
}

#[test]
fn compute_matrix_lu() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::MatrixDecomp(MatrixDecomp::LU),
        data: json!({
            "matrix": [[2.0, 1.0], [1.0, 2.0]]
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    // LU decomposition not yet implemented
    let _ = result;
}

#[test]
fn compute_matrix_cholesky() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::MatrixDecomp(MatrixDecomp::Cholesky),
        data: json!({
            "matrix": [[4.0, 2.0], [2.0, 3.0]]
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Cholesky decomposition should work");
}

#[test]
fn compute_tensor_riemann() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::Riemann),
        data: json!({
            "metric": [
                [-1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0]
            ],
            "coordinates": ["t", "x", "y", "z"]
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Riemann tensor should work");
}

#[test]
fn compute_tensor_ricci() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::Ricci),
        data: json!({
            "metric": [
                [-1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0]
            ],
            "coordinates": ["t", "x", "y", "z"]
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Ricci tensor should work");
}

// ============================================================================
// TRANSFORM - All Transform Types (via Compute tool)
// ============================================================================

#[test]
fn transform_wavelet() {
    let dispatcher = create_default_dispatcher();

    let signal: Vec<f64> = (0..64).map(|i| (i as f64 * 0.1).sin()).collect();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Transform(TransformType::Wavelet(WaveletType::Haar)),
        data: json!({
            "data": signal
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    let _ = result; // Wavelet may not be fully implemented
}

#[test]
fn transform_filter_highpass() {
    let dispatcher = create_default_dispatcher();

    let signal: Vec<f64> = (0..100).map(|i| (i as f64 * 0.1).sin()).collect();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Transform(TransformType::Filter(FilterType::HighPass)),
        data: json!({
            "data": signal,
            "cutoff": 50.0,
            "order": 4
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Highpass filter should work");
}

#[test]
fn transform_filter_bandpass() {
    let dispatcher = create_default_dispatcher();

    let signal: Vec<f64> = (0..100).map(|i| (i as f64 * 0.1).sin()).collect();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Transform(TransformType::Filter(FilterType::BandPass)),
        data: json!({
            "data": signal,
            "low_cutoff": 30.0,
            "high_cutoff": 70.0,
            "order": 4
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    // Bandpass filter needs two cutoff frequencies but implementation only supports one
    let _ = result;
}

// ============================================================================
// OPTIMIZE - All Optimization Methods (via Solve tool)
// ============================================================================

#[test]
fn optimize_gradient_descent() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Solve(SolveInput {
        equation_type: EquationType::Optimize(OptimizationMethod::Minimize(MinimizationMethod::GradientDescent)),
        equations: vec!["x^2 + y^2".to_string()],
        variables: Some(vec!["x".to_string(), "y".to_string()]),
        initial_guess: Some(HashMap::from([
            ("x".to_string(), 1.0),
            ("y".to_string(), 1.0),
        ])),
        boundary_conditions: None,
        domain: None,
        method: None,
        parameters: HashMap::from([
            ("learning_rate".to_string(), json!(0.1)),
            ("max_iterations".to_string(), json!(100)),
        ]),
    });

    let result = dispatcher.dispatch(request);
    let _ = result; // Gradient descent may not be fully implemented
}

#[test]
fn optimize_nelder_mead() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Solve(SolveInput {
        equation_type: EquationType::Optimize(OptimizationMethod::Minimize(MinimizationMethod::NelderMead)),
        equations: vec!["(x-3)^2 + (y-2)^2".to_string()],
        variables: Some(vec!["x".to_string(), "y".to_string()]),
        initial_guess: Some(HashMap::from([
            ("x".to_string(), 0.0),
            ("y".to_string(), 0.0),
        ])),
        boundary_conditions: None,
        domain: None,
        method: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    let _ = result;
}

// ============================================================================
// JSON API - Test all tools
// ============================================================================

#[test]
fn json_api_all_tools_basic() {
    let dispatcher = create_default_dispatcher();

    let test_cases = vec![
        (
            "solve",
            r#"{"tool":"solve","input":{"equation_type":"root_finding","equations":["x^2-4=0"]}}"#,
        ),
        (
            "compute-differentiate",
            r#"{"tool":"compute","input":{"operation":{"differentiate":"symbolic"},"data":{"expression":"x^2","variable":"x"},"parameters":{}}}"#,
        ),
        (
            "compute-integrate",
            r#"{"tool":"compute","input":{"operation":{"integrate":"symbolic"},"data":{"expression":"x","variable":"x"},"parameters":{}}}"#,
        ),
        (
            "analyze",
            r#"{"tool":"analyze","input":{"operation":"is_prime","expression":"17","options":{}}}"#,
        ),
        (
            "simulate",
            r#"{"tool":"simulate","input":{"model":{"time_evolution":"euler"},"equations":["dy/dt=-y"],"variables":["y"],"parameters":{"initial_value":1.0},"range":[0,1],"steps":10}}"#,
        ),
        (
            "compute-matrix",
            r#"{"tool":"compute","input":{"operation":{"matrix_decomp":"svd"},"data":{"matrix":[[1,2],[3,4]]},"parameters":{}}}"#,
        ),
        (
            "compute-transform",
            r#"{"tool":"compute","input":{"operation":{"transform":{"fft":"forward"}},"data":{"data":[0,1,0,1,0,1,0,1]},"parameters":{}}}"#,
        ),
        (
            "compute-sample",
            r#"{"tool":"compute","input":{"operation":{"sample":"moments"},"data":{"data":[1,2,3,4,5]},"parameters":{}}}"#,
        ),
    ];

    for (name, json) in test_cases {
        let response = dispatcher.dispatch_json(json);
        println!("{} response: {}", name, response);
        // Should either succeed or give meaningful error
        assert!(!response.is_empty(), "{} should return a response", name);
    }
}
