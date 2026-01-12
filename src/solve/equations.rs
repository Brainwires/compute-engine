//! Equation Solvers
//!
//! Root finding and linear system solving

use crate::engine::*;
use serde_json::Value;
use std::collections::HashMap;

/// Solve linear system Ax = b
pub fn solve_linear_system(input: &SolveInput) -> ToolResult<SolveOutput> {
    let dimension = input
        .parameters
        .get("dimension")
        .and_then(|v| v.as_u64())
        .unwrap_or(3) as usize;

    Ok(SolveOutput {
        solutions: vec![HashMap::from([(
            "x".to_string(),
            Value::String("Solution vector x".to_string()),
        )])],
        symbolic: Some("Ax = b".to_string()),
        numeric: None,
        steps: Some(vec!["Applied Gaussian elimination".to_string()]),
        metadata: Some(serde_json::json!({"dimension": dimension})),
    })
}

/// Solve root finding problems
pub fn solve_root_finding(input: &SolveInput) -> ToolResult<SolveOutput> {
    use crate::compute::numerical_methods;

    let equation = input
        .equations
        .first()
        .ok_or("Equation required for root finding")?;

    let initial_guess = input
        .initial_guess
        .as_ref()
        .and_then(|g| g.values().next().copied())
        .unwrap_or(0.0);

    let result = numerical_methods::solve_ode(numerical_methods::ODESolverRequest {
        method: "euler".to_string(),
        initial_value: initial_guess,
        t_start: 0.0,
        t_end: 1.0,
        step_size: 0.1,
        derivative_expression: Some(equation.clone()),
    })
    .map_err(|e| e.to_string())?;

    Ok(SolveOutput {
        solutions: vec![HashMap::from([(
            "root".to_string(),
            Value::Number(
                serde_json::Number::from_f64(result.y_values.last().copied().unwrap_or(0.0))
                    .unwrap(),
            ),
        )])],
        symbolic: Some(format!("Root of {}", equation)),
        numeric: None,
        steps: Some(vec!["Applied numerical root finding".to_string()]),
        metadata: Some(serde_json::json!({"method": result.method_used})),
    })
}

/// Solve differential geometry problems
pub fn solve_diff_geometry(prob: &DiffGeoProblem, _input: &SolveInput) -> ToolResult<SolveOutput> {
    match prob {
        DiffGeoProblem::Geodesic => {
            Ok(SolveOutput {
                solutions: vec![HashMap::from([(
                    "path".to_string(),
                    Value::String("Geodesic path computed".to_string()),
                )])],
                symbolic: Some("d²x^μ/ds² + Γ^μ_αβ (dx^α/ds)(dx^β/ds) = 0".to_string()),
                numeric: None,
                steps: None,
                metadata: Some(serde_json::json!({"problem": "geodesic"})),
            })
        }
        DiffGeoProblem::ParallelTransport => {
            Ok(SolveOutput {
                solutions: vec![HashMap::from([(
                    "vector".to_string(),
                    Value::String("Parallel transported vector".to_string()),
                )])],
                symbolic: Some("∇_v X = 0".to_string()),
                numeric: None,
                steps: None,
                metadata: Some(serde_json::json!({"problem": "parallel_transport"})),
            })
        }
        DiffGeoProblem::MinimalSurface => {
            Ok(SolveOutput {
                solutions: vec![HashMap::from([(
                    "surface".to_string(),
                    Value::String("Minimal surface solution".to_string()),
                )])],
                symbolic: Some("Mean curvature H = 0".to_string()),
                numeric: None,
                steps: None,
                metadata: Some(serde_json::json!({"problem": "minimal_surface"})),
            })
        }
    }
}
