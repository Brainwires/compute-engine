//! Differential Equation Solvers
//!
//! Solves differential equations:
//! - ODE (Ordinary Differential Equations)
//! - Initial value problems
//! - PDE (Partial Differential Equations)
//! - Boundary value problems

use crate::engine::*;
use std::collections::HashMap;

/// Solve differential equations
pub fn solve_differential(eq: &DifferentialEquation, input: &SolveInput) -> ToolResult<SolveOutput> {
    use crate::compute::numerical_methods;

    match eq {
        DifferentialEquation::ODE | DifferentialEquation::InitialValue => {
            let initial_value = input
                .initial_guess
                .as_ref()
                .and_then(|ig| ig.values().next().copied())
                .unwrap_or(1.0);

            let range = (0.0, 10.0);

            let result = numerical_methods::solve_ode(numerical_methods::ODESolverRequest {
                method: "rk4".to_string(),
                initial_value,
                t_start: range.0,
                t_end: range.1,
                step_size: (range.1 - range.0) / 100.0,
                derivative_expression: input.equations.first().cloned(),
            })
            .map_err(|e| e.to_string())?;

            Ok(SolveOutput {
                solutions: vec![],
                symbolic: None,
                numeric: Some(vec![HashMap::from([
                    (
                        "t".to_string(),
                        result.t_values.last().copied().unwrap_or(0.0),
                    ),
                    (
                        "y".to_string(),
                        result.y_values.last().copied().unwrap_or(0.0),
                    ),
                ])]),
                steps: Some(vec![format!("Solved ODE using {}", result.method_used)]),
                metadata: Some(serde_json::json!({"steps_taken": result.steps_taken})),
            })
        }
        DifferentialEquation::PDE | DifferentialEquation::BoundaryValue => {
            Err(format!("{:?} equations not yet fully implemented", eq))
        }
    }
}
