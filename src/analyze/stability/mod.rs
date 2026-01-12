//! Stability Analysis
//!
//! Field analysis for vector, scalar, and tensor fields
//! Also includes graph analysis and fluid analysis placeholders

use crate::engine::*;

/// Perform field analysis
pub fn analyze_field(field_type: &FieldAnalysisType, input: &AnalyzeInput) -> ToolResult<AnalyzeOutput> {
    use crate::analyze::symbolic as symbolic_cas;

    match field_type {
        FieldAnalysisType::Vector => {
            let components: Vec<String> = input
                .options
                .get("components")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or(
                    "Vector field analysis requires 'components' parameter [F_x, F_y, F_z]",
                )?;

            let variables: Vec<String> = input
                .options
                .get("variables")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .unwrap_or_else(|| vec!["x".to_string(), "y".to_string(), "z".to_string()]);

            if components.len() != variables.len() {
                return Err(format!(
                    "Mismatch: {} components but {} variables",
                    components.len(),
                    variables.len()
                ));
            }

            // Compute divergence: ∇·F = ∂F_x/∂x + ∂F_y/∂y + ∂F_z/∂z
            let mut div_terms = Vec::new();
            for (component, var) in components.iter().zip(variables.iter()) {
                let partial = symbolic_cas::differentiate(component, var, None)
                    .map_err(|e| format!("Failed to compute divergence: {}", e))?;
                div_terms.push(partial.expression);
            }
            let divergence = div_terms.join(" + ");

            // Compute curl if 3D
            let curl = if components.len() == 3 {
                let (fx, fy, fz) = (&components[0], &components[1], &components[2]);
                let (x, y, z) = (&variables[0], &variables[1], &variables[2]);

                let dfz_dy = symbolic_cas::differentiate(fz, y, None)
                    .map_err(|e| format!("Failed to compute curl: {}", e))?
                    .expression;
                let dfy_dz = symbolic_cas::differentiate(fy, z, None)
                    .map_err(|e| format!("Failed to compute curl: {}", e))?
                    .expression;
                let dfx_dz = symbolic_cas::differentiate(fx, z, None)
                    .map_err(|e| format!("Failed to compute curl: {}", e))?
                    .expression;
                let dfz_dx = symbolic_cas::differentiate(fz, x, None)
                    .map_err(|e| format!("Failed to compute curl: {}", e))?
                    .expression;
                let dfy_dx = symbolic_cas::differentiate(fy, x, None)
                    .map_err(|e| format!("Failed to compute curl: {}", e))?
                    .expression;
                let dfx_dy = symbolic_cas::differentiate(fx, y, None)
                    .map_err(|e| format!("Failed to compute curl: {}", e))?
                    .expression;

                let curl_x = format!("({}) - ({})", dfz_dy, dfy_dz);
                let curl_y = format!("({}) - ({})", dfx_dz, dfz_dx);
                let curl_z = format!("({}) - ({})", dfy_dx, dfx_dy);
                Some(vec![curl_x.clone(), curl_y.clone(), curl_z.clone()])
            } else {
                None
            };

            let is_conservative = curl
                .as_ref()
                .map(|c| c.iter().all(|comp| comp == "0" || comp == "(0) - (0)"))
                .unwrap_or(false);

            Ok(AnalyzeOutput {
                result: serde_json::json!({
                    "field_type": "vector",
                    "divergence": divergence,
                    "curl": curl,
                    "is_conservative": is_conservative,
                    "is_solenoidal": divergence == "0",
                    "dimension": components.len()
                }),
                latex: Some(format!("$\\nabla \\cdot \\vec{{F}} = {}$", divergence)),
                validation: None,
                details: Some(serde_json::json!({
                    "analysis_type": "vector_field",
                    "components": components,
                    "variables": variables
                })),
            })
        }
        FieldAnalysisType::Scalar => {
            let variables: Vec<String> = input
                .options
                .get("variables")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .unwrap_or_else(|| vec!["x".to_string(), "y".to_string(), "z".to_string()]);

            let gradient_results = symbolic_cas::gradient(&input.expression, &variables)
                .map_err(|e| format!("Failed to compute gradient: {}", e))?;

            let gradient_components: Vec<String> = gradient_results
                .iter()
                .map(|r| r.expression.clone())
                .collect();

            // Compute Laplacian
            let mut laplacian_terms = Vec::new();
            for var in &variables {
                let second_deriv = symbolic_cas::differentiate(&input.expression, var, Some(2))
                    .map_err(|e| format!("Failed to compute Laplacian: {}", e))?;
                laplacian_terms.push(second_deriv.expression);
            }
            let laplacian = laplacian_terms.join(" + ");

            let is_harmonic = laplacian == "0";

            Ok(AnalyzeOutput {
                result: serde_json::json!({
                    "field_type": "scalar",
                    "gradient": gradient_components,
                    "laplacian": laplacian,
                    "is_harmonic": is_harmonic,
                    "dimension": variables.len()
                }),
                latex: Some(format!(
                    "$\\nabla f = ({})$, $\\nabla^2 f = {}$",
                    gradient_components.join(", "),
                    laplacian
                )),
                validation: None,
                details: Some(serde_json::json!({
                    "analysis_type": "scalar_field",
                    "expression": input.expression.clone(),
                    "variables": variables
                })),
            })
        }
        FieldAnalysisType::Tensor => {
            let dimension = input
                .options
                .get("dimension")
                .and_then(|v| v.as_u64())
                .unwrap_or(4) as usize;

            let metric_data: Option<Vec<Vec<f64>>> = input
                .options
                .get("metric")
                .and_then(|v| serde_json::from_value(v.clone()).ok());

            let variables: Vec<String> = input
                .options
                .get("variables")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .unwrap_or_else(|| {
                    vec![
                        "t".to_string(),
                        "x".to_string(),
                        "y".to_string(),
                        "z".to_string(),
                    ]
                });

            let mut properties = serde_json::json!({
                "dimension": dimension,
                "rank": 2
            });

            if let Some(metric_vals) = metric_data {
                let coords: Vec<&str> = variables.iter().map(|s| s.as_str()).collect();

                let metric_exprs: Vec<Vec<symbolic_cas::Expr>> = metric_vals
                    .iter()
                    .map(|row| {
                        row.iter()
                            .map(|&val| symbolic_cas::Expr::num(val as i64))
                            .collect()
                    })
                    .collect();

                let metric = symbolic_cas::SymbolicMatrix::new(metric_exprs)
                    .map_err(|e| format!("Invalid metric tensor: {}", e))?;

                let _christoffel = symbolic_cas::christoffel_symbols(&metric, &coords)
                    .map_err(|e| format!("Failed to compute Christoffel symbols: {}", e))?;

                let _riemann = symbolic_cas::riemann_tensor(&metric, &coords)
                    .map_err(|e| format!("Failed to compute Riemann tensor: {}", e))?;

                let _ricci = symbolic_cas::ricci_tensor(&metric, &coords)
                    .map_err(|e| format!("Failed to compute Ricci tensor: {}", e))?;

                let ricci_scalar_val = symbolic_cas::ricci_scalar(&metric, &coords)
                    .map_err(|e| format!("Failed to compute Ricci scalar: {}", e))?;

                properties = serde_json::json!({
                    "dimension": dimension,
                    "rank": 2,
                    "has_curvature": true,
                    "ricci_scalar": format!("{}", ricci_scalar_val)
                });
            }

            Ok(AnalyzeOutput {
                result: serde_json::json!({
                    "field_type": "tensor",
                    "properties": properties
                }),
                latex: None,
                validation: None,
                details: Some(serde_json::json!({
                    "analysis_type": "tensor_field",
                    "dimension": dimension,
                    "variables": variables
                })),
            })
        }
    }
}

/// Graph component analysis
pub fn analyze_graph_components(input: &AnalyzeInput) -> ToolResult<AnalyzeOutput> {
    let num_vertices = input
        .options
        .get("vertices")
        .and_then(|v| v.as_u64())
        .unwrap_or(0) as usize;

    Ok(AnalyzeOutput {
        result: serde_json::json!({
            "num_components": 1,
            "num_vertices": num_vertices
        }),
        latex: None,
        validation: None,
        details: Some(serde_json::json!({
            "graph_type": "undirected"
        })),
    })
}

/// Graph property analysis
pub fn analyze_graph_properties(_input: &AnalyzeInput) -> ToolResult<AnalyzeOutput> {
    Ok(AnalyzeOutput {
        result: serde_json::json!({
            "is_connected": true,
            "is_bipartite": false,
            "is_planar": true
        }),
        latex: None,
        validation: None,
        details: Some(serde_json::json!({
            "analysis_method": "basic"
        })),
    })
}

/// Fluid analysis
pub fn analyze_fluid(input: &AnalyzeInput) -> ToolResult<AnalyzeOutput> {
    let reynolds = input
        .options
        .get("reynolds")
        .and_then(|v| v.as_f64())
        .unwrap_or(100.0);

    let flow_regime = if reynolds < 2300.0 {
        "laminar"
    } else if reynolds > 4000.0 {
        "turbulent"
    } else {
        "transitional"
    };

    Ok(AnalyzeOutput {
        result: serde_json::json!({
            "reynolds_number": reynolds,
            "flow_regime": flow_regime
        }),
        latex: None,
        validation: None,
        details: Some(serde_json::json!({
            "analysis_type": "flow_regime"
        })),
    })
}

/// Simple primality test
pub fn is_prime_simple(n: u64) -> bool {
    if n <= 1 {
        return false;
    }
    if n <= 3 {
        return true;
    }
    if n % 2 == 0 || n % 3 == 0 {
        return false;
    }

    let mut i = 5;
    while i * i <= n {
        if n % i == 0 || n % (i + 2) == 0 {
            return false;
        }
        i += 6;
    }

    true
}
