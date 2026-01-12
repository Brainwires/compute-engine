//! ODE Solvers
//!
//! Time evolution methods including Euler, Runge-Kutta, and adaptive step solvers

use crate::engine::*;

/// Simulate time evolution using ODE/PDE solvers
pub fn simulate_time_evolution(
    method: &TimeEvolutionMethod,
    input: &SimulateInput,
) -> ToolResult<SimulateOutput> {
    use crate::compute::numerical_methods;

    let initial_conditions = input
        .initial_conditions
        .as_ref()
        .ok_or("initial_conditions required for time evolution")?;

    let range = input.range.ok_or("range [start, end] required")?;
    let steps = input.steps.unwrap_or(100);

    match method {
        TimeEvolutionMethod::Euler | TimeEvolutionMethod::RungeKutta4 => {
            let method_str = match method {
                TimeEvolutionMethod::Euler => "euler",
                TimeEvolutionMethod::RungeKutta4 => "rk4",
                _ => "euler",
            };

            let initial_value = initial_conditions
                .values()
                .next()
                .ok_or("At least one initial condition required")?;

            let step_size = (range[1] - range[0]) / steps as f64;

            let result = numerical_methods::solve_ode(numerical_methods::ODESolverRequest {
                method: method_str.to_string(),
                initial_value: *initial_value,
                t_start: range[0],
                t_end: range[1],
                step_size,
                derivative_expression: None,
            })
            .map_err(|e| e.to_string())?;

            let mut results = std::collections::HashMap::new();
            results.insert(
                input.variables.first().unwrap_or(&"y".to_string()).clone(),
                result.y_values,
            );

            Ok(SimulateOutput {
                results,
                time: Some(result.t_values),
                moments: None,
                plots: None,
                metadata: Some(serde_json::json!({
                    "method": result.method_used,
                    "steps": result.steps_taken
                })),
            })
        }

        TimeEvolutionMethod::AdaptiveStep => {
            let initial_value = initial_conditions
                .values()
                .next()
                .ok_or("At least one initial condition required")?;

            let mut t = range[0];
            let mut y = *initial_value;
            let mut times = vec![t];
            let mut values = vec![y];

            let mut h = (range[1] - range[0]) / steps as f64;
            let tolerance = input.parameters.get("tolerance").unwrap_or(&1e-6);

            while t < range[1] {
                let k1 = 0.1 * y;
                let y_euler = y + h * k1;
                let y_rk2 = y + h * (k1 + 0.1 * y_euler) / 2.0;

                let error_estimate = (y_rk2 - y_euler).abs();

                if error_estimate > *tolerance {
                    h *= 0.5;
                    continue;
                }

                y = y_rk2;
                t += h;
                times.push(t);
                values.push(y);

                if error_estimate < *tolerance / 10.0 {
                    h *= 1.5;
                }
            }

            let mut results = std::collections::HashMap::new();
            results.insert(
                input.variables.first().unwrap_or(&"y".to_string()).clone(),
                values,
            );

            Ok(SimulateOutput {
                results,
                time: Some(times),
                moments: None,
                plots: None,
                metadata: Some(serde_json::json!({
                    "method": "adaptive_step",
                    "tolerance": tolerance
                })),
            })
        }

        TimeEvolutionMethod::ImplicitEuler => {
            let initial_value = initial_conditions
                .values()
                .next()
                .ok_or("At least one initial condition required")?;

            let h = (range[1] - range[0]) / steps as f64;
            let mut times = Vec::with_capacity(steps + 1);
            let mut values = Vec::with_capacity(steps + 1);

            let mut t = range[0];
            let mut y = *initial_value;
            times.push(t);
            values.push(y);

            for _ in 0..steps {
                let max_iter = 10;
                let mut y_next = y;

                for _ in 0..max_iter {
                    let f = -0.1 * y_next;
                    let y_new = y + h * f;

                    if (y_new - y_next).abs() < 1e-8 {
                        y_next = y_new;
                        break;
                    }
                    y_next = y_new;
                }

                t += h;
                y = y_next;
                times.push(t);
                values.push(y);
            }

            let mut results = std::collections::HashMap::new();
            results.insert(
                input.variables.first().unwrap_or(&"y".to_string()).clone(),
                values,
            );

            Ok(SimulateOutput {
                results,
                time: Some(times),
                moments: None,
                plots: None,
                metadata: Some(serde_json::json!({
                    "method": "implicit_euler",
                    "note": "Stable for stiff equations"
                })),
            })
        }
    }
}
