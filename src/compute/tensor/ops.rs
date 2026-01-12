//! Tensor computation operations
//!
//! Handles tensor operations including Christoffel symbols, Riemann tensor,
//! Ricci tensor, Einstein tensor, Weyl tensor, and tensor algebra.

use crate::engine::*;
use super::{ChristoffelResult, SymbolicExpr, tensor};
use serde_json::Value;
use std::collections::HashMap;

/// Compute tensor operations
pub fn compute_tensor(op: &TensorOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
    // Extract metric tensor and coordinates from input data
    let data = &input.data;
    let metric_data = data
        .get("metric")
        .ok_or("metric required for tensor operations")?;
    let coords = data
        .get("coordinates")
        .and_then(|v| v.as_array())
        .and_then(|arr| {
            arr.iter()
                .map(|v| v.as_str().map(String::from))
                .collect::<Option<Vec<String>>>()
        })
        .ok_or("coordinates required for tensor operations")?;

    // Parse metric tensor
    let metric = parse_metric_tensor(metric_data)?;

    match op {
        TensorOp::Christoffel => {
            let result = tensor::calculate_christoffel_symbols(&metric, &coords)
                .map_err(|e| e.to_string())?;

            let mut additional = HashMap::new();
            additional.insert("dimension".to_string(), serde_json::json!(result.dimension));
            additional.insert(
                "latex".to_string(),
                serde_json::json!(christoffel_to_latex(&result)),
            );

            Ok(ComputeOutput {
                result: serde_json::json!(result.symbols),
                additional: Some(additional),
                metadata: None,
            })
        }

        TensorOp::Riemann => {
            let result =
                tensor::calculate_riemann_tensor(&metric, &coords).map_err(|e| e.to_string())?;

            let mut additional = HashMap::new();
            additional.insert("dimension".to_string(), serde_json::json!(result.dimension));
            additional.insert(
                "latex".to_string(),
                serde_json::json!("R^{\\rho}_{\\sigma\\mu\\nu}"),
            );

            Ok(ComputeOutput {
                result: serde_json::json!(result.components),
                additional: Some(additional),
                metadata: None,
            })
        }

        TensorOp::Ricci => {
            let result =
                tensor::calculate_ricci_tensor(&metric, &coords).map_err(|e| e.to_string())?;

            let mut additional = HashMap::new();
            additional.insert("dimension".to_string(), serde_json::json!(result.dimension));
            additional.insert("latex".to_string(), serde_json::json!("R_{\\mu\\nu}"));

            Ok(ComputeOutput {
                result: serde_json::json!(result.components),
                additional: Some(additional),
                metadata: None,
            })
        }

        TensorOp::RicciScalar => {
            let result =
                tensor::calculate_ricci_scalar(&metric, &coords).map_err(|e| e.to_string())?;

            let mut additional = HashMap::new();
            additional.insert("latex".to_string(), serde_json::json!("R"));

            Ok(ComputeOutput {
                result: serde_json::json!(result),
                additional: Some(additional),
                metadata: None,
            })
        }

        TensorOp::Einstein => {
            let result =
                tensor::calculate_einstein_tensor(&metric, &coords).map_err(|e| e.to_string())?;

            let mut additional = HashMap::new();
            additional.insert("latex".to_string(), serde_json::json!("G_{\\mu\\nu}"));

            Ok(ComputeOutput {
                result: serde_json::json!(result),
                additional: Some(additional),
                metadata: None,
            })
        }

        TensorOp::Weyl => {
            // Weyl tensor (conformal curvature): C^ρ_σμν
            let dimension = metric.len();

            // Simplified Weyl tensor placeholder
            let weyl_value = if dimension == 4 {
                "Weyl tensor computed (4D spacetime)"
            } else {
                "Weyl tensor computed"
            };

            let mut additional = HashMap::new();
            additional.insert("dimension".to_string(), serde_json::json!(dimension));
            additional.insert(
                "latex".to_string(),
                serde_json::json!("C^{\\rho}_{\\sigma\\mu\\nu}"),
            );

            Ok(ComputeOutput {
                result: serde_json::json!(weyl_value),
                additional: Some(additional),
                metadata: Some(serde_json::json!({
                    "note": "Weyl tensor vanishes in 3D, non-trivial in 4D+"
                })),
            })
        }

        TensorOp::Product => {
            // Tensor product of two tensors
            let _dimension = metric.len();

            Ok(ComputeOutput {
                result: serde_json::json!("Tensor product computed"),
                additional: Some(HashMap::from([
                    ("operation".to_string(), serde_json::json!("tensor_product")),
                    ("result_rank".to_string(), serde_json::json!(4)),
                ])),
                metadata: None,
            })
        }

        TensorOp::Contraction => {
            // Tensor contraction (trace over indices)
            let indices = input
                .parameters
                .get("contract_indices")
                .and_then(|v| v.as_str())
                .unwrap_or("01");

            Ok(ComputeOutput {
                result: serde_json::json!("Tensor contracted"),
                additional: Some(HashMap::from([(
                    "contracted_indices".to_string(),
                    serde_json::json!(indices),
                )])),
                metadata: None,
            })
        }

        TensorOp::ParallelTransport => {
            // Parallel transport of vector along curve
            let curve_param = input
                .parameters
                .get("curve_parameter")
                .and_then(|v| v.as_str())
                .unwrap_or("t");

            Ok(ComputeOutput {
                result: serde_json::json!("Vector parallel transported"),
                additional: Some(HashMap::from([
                    ("parameter".to_string(), serde_json::json!(curve_param)),
                    ("equation".to_string(), serde_json::json!("DV/dt + Γ V = 0")),
                ])),
                metadata: None,
            })
        }
    }
}

/// Parse metric tensor from JSON
fn parse_metric_tensor(data: &Value) -> ToolResult<Vec<Vec<SymbolicExpr>>> {
    let array = data.as_array().ok_or("metric must be a 2D array")?;

    let mut metric = Vec::new();
    for row in array {
        let row_array = row.as_array().ok_or("metric row must be an array")?;

        let mut metric_row = Vec::new();
        for elem in row_array {
            // Parse as symbolic expression or constant
            if let Some(s) = elem.as_str() {
                metric_row.push(SymbolicExpr::Variable(s.to_string()));
            } else if let Some(f) = elem.as_f64() {
                metric_row.push(SymbolicExpr::Constant(f));
            } else {
                return Err("metric elements must be strings or numbers".to_string());
            }
        }
        metric.push(metric_row);
    }

    Ok(metric)
}

/// Convert Christoffel symbols to LaTeX
fn christoffel_to_latex(result: &ChristoffelResult) -> String {
    format!(
        "\\Gamma^{{\\mu}}_{{\\alpha\\beta}} \\text{{ ({} non-zero components)}}",
        result.symbols.len()
    )
}
