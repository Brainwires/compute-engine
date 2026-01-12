//! Einstein Field Equation Solvers
//!
//! Solves various forms of Einstein's field equations:
//! - Vacuum solutions (Schwarzschild, Kerr)
//! - Solutions with matter sources
//! - Cosmological solutions (FLRW)

use crate::engine::*;
use serde_json::Value;
use std::collections::HashMap;

/// Solve Einstein field equations
pub fn solve_einstein(eq: &EinsteinEquation, input: &SolveInput) -> ToolResult<SolveOutput> {
    use crate::compute::tensor as tensor_calculus;

    // ALWAYS use 4D coordinates for Einstein equations (spacetime)
    let coords = vec![
        "t".to_string(),
        "r".to_string(),
        "theta".to_string(),
        "phi".to_string(),
    ];

    match eq {
        EinsteinEquation::Vacuum => {
            let symmetry = input
                .parameters
                .get("symmetry")
                .and_then(|v| v.as_str())
                .unwrap_or("spherical");

            let solutions =
                tensor_calculus::solve_vacuum_einstein_equations(&coords, symmetry, &[])
                    .map_err(|e| e.to_string())?;

            Ok(SolveOutput {
                solutions: solutions
                    .iter()
                    .map(|sol| {
                        let mut map = HashMap::new();

                        let metric_strings: Vec<Vec<String>> = sol
                            .metric_tensor
                            .iter()
                            .map(|row| row.iter().map(|expr| expr.to_string()).collect())
                            .collect();

                        map.insert("metric".to_string(), serde_json::json!(metric_strings));
                        map.insert(
                            "solution_type".to_string(),
                            Value::String(sol.solution_type.clone()),
                        );

                        let params: HashMap<String, String> = sol
                            .physical_parameters
                            .iter()
                            .map(|(k, v)| (k.clone(), v.to_string()))
                            .collect();
                        map.insert(
                            "physical_parameters".to_string(),
                            serde_json::json!(params),
                        );

                        map.insert("coordinates".to_string(), serde_json::json!(sol.coordinates));

                        map
                    })
                    .collect(),
                symbolic: Some(format!(
                    "Einstein vacuum equations solved with {} symmetry",
                    symmetry
                )),
                numeric: None,
                steps: Some(vec![
                    "Applied symmetry ansatz".to_string(),
                    "Solved field equations".to_string(),
                    "Verified solution".to_string(),
                ]),
                metadata: Some(serde_json::json!({
                    "equation_type": "einstein_vacuum",
                    "symmetry": symmetry,
                    "coordinates": coords
                })),
            })
        }

        EinsteinEquation::WithSource => {
            let _stress_energy = input
                .parameters
                .get("stress_energy")
                .ok_or("stress_energy parameter required for WithSource")?;

            Err("Einstein with source: full implementation pending".to_string())
        }

        EinsteinEquation::Schwarzschild => {
            let solutions =
                tensor_calculus::solve_vacuum_einstein_equations(&coords, "spherical", &[])
                    .map_err(|e| e.to_string())?;

            Ok(SolveOutput {
                solutions: solutions
                    .iter()
                    .map(|sol| {
                        let mut map = HashMap::new();

                        let metric_strings: Vec<Vec<String>> = sol
                            .metric_tensor
                            .iter()
                            .map(|row| row.iter().map(|expr| expr.to_string()).collect())
                            .collect();

                        map.insert("metric".to_string(), serde_json::json!(metric_strings));
                        map.insert(
                            "solution_type".to_string(),
                            Value::String(sol.solution_type.clone()),
                        );

                        let params: HashMap<String, String> = sol
                            .physical_parameters
                            .iter()
                            .map(|(k, v)| (k.clone(), v.to_string()))
                            .collect();
                        map.insert(
                            "physical_parameters".to_string(),
                            serde_json::json!(params),
                        );

                        map.insert("coordinates".to_string(), serde_json::json!(sol.coordinates));

                        map
                    })
                    .collect(),
                symbolic: Some("Schwarzschild solution".to_string()),
                numeric: None,
                steps: None,
                metadata: Some(serde_json::json!({"solution": "schwarzschild"})),
            })
        }

        EinsteinEquation::KerrNewman => {
            let mass = input
                .parameters
                .get("mass")
                .and_then(|v| v.as_f64())
                .ok_or("mass parameter required for Kerr-Newman")?;
            let angular_momentum = input
                .parameters
                .get("angular_momentum")
                .and_then(|v| v.as_f64())
                .unwrap_or(0.0);
            let charge = input
                .parameters
                .get("charge")
                .and_then(|v| v.as_f64())
                .unwrap_or(0.0);

            let m = mass;
            let a = angular_momentum / mass;
            let q = charge;

            let metric_components = serde_json::json!({
                "g_tt": format!("-(r² + a²cos²θ - 2Mr + Q²)/(r² + a²cos²θ)"),
                "g_rr": format!("(r² + a²cos²θ)/(r² - 2Mr + a² + Q²)"),
                "g_θθ": "r² + a²cos²θ",
                "g_φφ": format!("sin²θ[(r² + a²)² - a²sin²θ(r² - 2Mr + a² + Q²)]/(r² + a²cos²θ)"),
                "g_tφ": format!("-2Mra sin²θ/(r² + a²cos²θ)"),
                "description": "Kerr-Newman metric in Boyer-Lindquist coordinates",
                "M": m,
                "a": a,
                "Q": q,
                "horizons": {
                    "r_plus": m + (m.powi(2) - a.powi(2) - q.powi(2)).sqrt(),
                    "r_minus": m - (m.powi(2) - a.powi(2) - q.powi(2)).sqrt()
                }
            });

            Ok(SolveOutput {
                solutions: vec![{
                    let mut map = HashMap::new();
                    map.insert("metric".to_string(), metric_components);
                    map
                }],
                symbolic: Some(format!("Kerr-Newman solution: M={}, a={}, Q={}", m, a, q)),
                numeric: None,
                steps: Some(vec![
                    "Applied Boyer-Lindquist coordinates".to_string(),
                    "Computed metric tensor components".to_string(),
                    "Calculated event horizons".to_string(),
                ]),
                metadata: Some(serde_json::json!({
                    "solution": "kerr_newman",
                    "mass": m,
                    "angular_momentum": a,
                    "charge": q,
                    "type": "rotating_charged_black_hole"
                })),
            })
        }

        EinsteinEquation::FriedmannRobertsonWalker => {
            let solutions =
                tensor_calculus::solve_vacuum_einstein_equations(&coords, "cosmological", &[])
                    .map_err(|e| e.to_string())?;

            Ok(SolveOutput {
                solutions: solutions
                    .iter()
                    .map(|sol| {
                        let mut map = HashMap::new();

                        let metric_strings: Vec<Vec<String>> = sol
                            .metric_tensor
                            .iter()
                            .map(|row| row.iter().map(|expr| expr.to_string()).collect())
                            .collect();

                        map.insert("metric".to_string(), serde_json::json!(metric_strings));
                        map.insert(
                            "solution_type".to_string(),
                            Value::String(sol.solution_type.clone()),
                        );

                        let params: HashMap<String, String> = sol
                            .physical_parameters
                            .iter()
                            .map(|(k, v)| (k.clone(), v.to_string()))
                            .collect();
                        map.insert(
                            "physical_parameters".to_string(),
                            serde_json::json!(params),
                        );

                        map.insert("coordinates".to_string(), serde_json::json!(sol.coordinates));

                        map
                    })
                    .collect(),
                symbolic: Some("Friedmann-Robertson-Walker universe".to_string()),
                numeric: None,
                steps: None,
                metadata: Some(serde_json::json!({"solution": "flrw"})),
            })
        }
    }
}
