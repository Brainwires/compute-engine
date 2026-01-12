//! Optimization Solvers
//!
//! Curve fitting, minimization, and optimization

use crate::engine::*;
use serde_json::Value;
use std::collections::HashMap;

/// Solve linear system Ax = b using Gaussian elimination with partial pivoting
fn solve_linear_system_gauss(a: &[f64], b: &[f64], n: usize) -> Result<Vec<f64>, String> {
    let mut aug = vec![0.0; n * (n + 1)];
    for i in 0..n {
        for j in 0..n {
            aug[i * (n + 1) + j] = a[i * n + j];
        }
        aug[i * (n + 1) + n] = b[i];
    }

    // Forward elimination with partial pivoting
    for col in 0..n {
        let mut max_row = col;
        let mut max_val = aug[col * (n + 1) + col].abs();
        for row in (col + 1)..n {
            let val = aug[row * (n + 1) + col].abs();
            if val > max_val {
                max_val = val;
                max_row = row;
            }
        }

        if max_val < 1e-14 {
            return Err("Singular or nearly singular matrix".to_string());
        }

        if max_row != col {
            for j in 0..=n {
                let tmp = aug[col * (n + 1) + j];
                aug[col * (n + 1) + j] = aug[max_row * (n + 1) + j];
                aug[max_row * (n + 1) + j] = tmp;
            }
        }

        for row in (col + 1)..n {
            let factor = aug[row * (n + 1) + col] / aug[col * (n + 1) + col];
            for j in col..=n {
                aug[row * (n + 1) + j] -= factor * aug[col * (n + 1) + j];
            }
        }
    }

    // Back substitution
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = aug[i * (n + 1) + n];
        for j in (i + 1)..n {
            sum -= aug[i * (n + 1) + j] * x[j];
        }
        x[i] = sum / aug[i * (n + 1) + i];
    }

    Ok(x)
}

/// Solve optimization problems routed through the Solve tool
pub fn solve_optimization(opt_method: &OptimizationMethod, input: &SolveInput) -> ToolResult<SolveOutput> {
    use crate::engine::equations::OptimizationMethod;

    // Get data from parameters.x_data/y_data (preferred) or fall back to initial_guess
    let data: Option<(Vec<f64>, Vec<f64>)> = {
        let x_from_params = input
            .parameters
            .get("x_data")
            .and_then(|v| v.as_array())
            .map(|arr| arr.iter().filter_map(|v| v.as_f64()).collect::<Vec<f64>>());
        let y_from_params = input
            .parameters
            .get("y_data")
            .and_then(|v| v.as_array())
            .map(|arr| arr.iter().filter_map(|v| v.as_f64()).collect::<Vec<f64>>());

        if let (Some(x), Some(y)) = (x_from_params, y_from_params) {
            Some((x, y))
        } else {
            input.initial_guess.as_ref().map(|g| {
                let x_data: Vec<f64> = g.keys().filter_map(|k| k.parse().ok()).collect();
                let y_data: Vec<f64> = g.values().copied().collect();
                (x_data, y_data)
            })
        }
    };

    match opt_method {
        OptimizationMethod::Fit(fit_method) => {
            use crate::engine::equations::FitMethod;
            match fit_method {
                FitMethod::Polynomial => {
                    if let Some((x, y)) = data {
                        let degree = input
                            .parameters
                            .get("degree")
                            .and_then(|v| v.as_u64())
                            .unwrap_or(2) as usize;

                        let n = x.len();
                        let m = degree + 1;
                        if n < m {
                            return Err("Not enough data points for polynomial degree".to_string());
                        }

                        // Build X matrix (Vandermonde)
                        let mut xtx = vec![0.0; m * m];
                        let mut xty = vec![0.0; m];

                        for i in 0..n {
                            let mut xi_pow = vec![1.0; m];
                            for j in 1..m {
                                xi_pow[j] = xi_pow[j - 1] * x[i];
                            }
                            for j in 0..m {
                                xty[j] += xi_pow[j] * y[i];
                                for k in 0..m {
                                    xtx[j * m + k] += xi_pow[j] * xi_pow[k];
                                }
                            }
                        }

                        let coeffs = solve_linear_system_gauss(&xtx, &xty, m)?;

                        Ok(SolveOutput {
                            solutions: vec![{
                                let mut map = HashMap::new();
                                for (i, &c) in coeffs.iter().enumerate() {
                                    map.insert(
                                        format!("a{}", i),
                                        Value::Number(
                                            serde_json::Number::from_f64(c)
                                                .unwrap_or(serde_json::Number::from(0)),
                                        ),
                                    );
                                }
                                map
                            }],
                            symbolic: Some(format!("Polynomial degree {}", degree)),
                            numeric: None,
                            steps: Some(vec!["Polynomial fit completed".to_string()]),
                            metadata: None,
                        })
                    } else {
                        Err("No data provided for polynomial fit".to_string())
                    }
                }
                _ => Ok(SolveOutput {
                    solutions: vec![],
                    symbolic: Some(format!(
                        "Fit method {:?} not yet fully implemented",
                        fit_method
                    )),
                    numeric: None,
                    steps: None,
                    metadata: None,
                }),
            }
        }
        OptimizationMethod::Minimize(min_type) => Ok(SolveOutput {
            solutions: vec![],
            symbolic: Some(format!(
                "Minimization {:?} not yet fully implemented",
                min_type
            )),
            numeric: None,
            steps: None,
            metadata: None,
        }),
        _ => Ok(SolveOutput {
            solutions: vec![],
            symbolic: Some(format!(
                "Optimization method {:?} not yet fully implemented",
                opt_method
            )),
            numeric: None,
            steps: None,
            metadata: None,
        }),
    }
}
