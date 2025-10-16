//! Unified Differentiator implementation
//!
//! Routes differentiation requests to vector calculus, tensor calculus, and numeric differentiation modules

use crate::engine::*;

pub struct UnifiedDifferentiator;

impl UnifiedDifferentiator {
    pub fn new() -> Self {
        Self
    }

    /// Perform vector calculus operations (gradient, divergence, curl, laplacian)
    fn differentiate_vector_calc(
        &self,
        op: &VectorCalcOp,
        input: &DifferentiateInput,
    ) -> ToolResult<DifferentiateOutput> {
        use crate::mathematics::symbolic_cas;

        match op {
            VectorCalcOp::Gradient => {
                // Compute symbolic gradient: ∇f = (∂f/∂x, ∂f/∂y, ∂f/∂z, ...)
                let gradient_results = symbolic_cas::gradient(&input.expression, &input.variables)
                    .map_err(|e| format!("Gradient computation failed: {}", e))?;

                let mut derivatives = std::collections::HashMap::new();
                let mut latex_map = std::collections::HashMap::new();

                for (i, result) in gradient_results.iter().enumerate() {
                    if let Some(var) = input.variables.get(i) {
                        derivatives.insert(
                            format!("∂f/∂{}", var),
                            serde_json::json!(result.expression.clone()),
                        );
                        if let Some(latex) = &result.latex {
                            latex_map.insert(format!("∂f/∂{}", var), latex.clone());
                        }
                    }
                }

                // Add full gradient vector representation
                let gradient_components: Vec<String> = gradient_results
                    .iter()
                    .map(|r| r.expression.clone())
                    .collect();
                derivatives.insert(
                    "gradient".to_string(),
                    serde_json::json!(gradient_components),
                );

                Ok(DifferentiateOutput {
                    derivatives,
                    latex: Some(latex_map),
                    metadata: Some(serde_json::json!({
                        "operation": "gradient",
                        "method": "symbolic",
                        "dimensions": input.variables.len()
                    })),
                })
            }

            VectorCalcOp::Divergence => {
                // Divergence of vector field: ∇·F = ∂F_x/∂x + ∂F_y/∂y + ∂F_z/∂z
                // Input expression should be a vector field, parse components
                let components: Vec<String> = input.parameters.get("components")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .or_else(|| {
                        // Try to parse comma-separated components from expression
                        Some(input.expression.split(',').map(|s| s.trim().to_string()).collect())
                    })
                    .ok_or("Divergence requires vector field components (use 'components' parameter or comma-separated expression)")?;

                if components.len() != input.variables.len() {
                    return Err(format!(
                        "Mismatch: {} components but {} variables",
                        components.len(),
                        input.variables.len()
                    ));
                }

                // Compute ∂F_i/∂x_i for each component
                let mut div_terms = Vec::new();
                for (component, var) in components.iter().zip(input.variables.iter()) {
                    let partial = symbolic_cas::differentiate(component, var, None)
                        .map_err(|e| format!("Failed to differentiate component: {}", e))?;
                    div_terms.push(partial.expression);
                }

                // Sum all terms
                let divergence = div_terms.join(" + ");

                let mut derivatives = std::collections::HashMap::new();
                derivatives.insert("divergence".to_string(), serde_json::json!(divergence));

                Ok(DifferentiateOutput {
                    derivatives,
                    latex: Some(std::collections::HashMap::from([(
                        "divergence".to_string(),
                        format!("\\nabla \\cdot \\vec{{F}} = {}", divergence),
                    )])),
                    metadata: Some(serde_json::json!({
                        "operation": "divergence",
                        "method": "symbolic",
                        "components": components
                    })),
                })
            }

            VectorCalcOp::Curl => {
                // Curl of vector field: ∇×F (only defined in 3D)
                let components: Vec<String> = input
                    .parameters
                    .get("components")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("Curl requires vector field components [F_x, F_y, F_z]")?;

                if components.len() != 3 || input.variables.len() != 3 {
                    return Err("Curl is only defined for 3D vector fields".to_string());
                }

                let (fx, fy, fz) = (&components[0], &components[1], &components[2]);
                let (x, y, z) = (
                    &input.variables[0],
                    &input.variables[1],
                    &input.variables[2],
                );

                // Curl = (∂F_z/∂y - ∂F_y/∂z, ∂F_x/∂z - ∂F_z/∂x, ∂F_y/∂x - ∂F_x/∂y)
                let dfz_dy = symbolic_cas::differentiate(fz, y, None)
                    .map_err(|e| format!("Failed to compute ∂F_z/∂y: {}", e))?
                    .expression;
                let dfy_dz = symbolic_cas::differentiate(fy, z, None)
                    .map_err(|e| format!("Failed to compute ∂F_y/∂z: {}", e))?
                    .expression;
                let dfx_dz = symbolic_cas::differentiate(fx, z, None)
                    .map_err(|e| format!("Failed to compute ∂F_x/∂z: {}", e))?
                    .expression;
                let dfz_dx = symbolic_cas::differentiate(fz, x, None)
                    .map_err(|e| format!("Failed to compute ∂F_z/∂x: {}", e))?
                    .expression;
                let dfy_dx = symbolic_cas::differentiate(fy, x, None)
                    .map_err(|e| format!("Failed to compute ∂F_y/∂x: {}", e))?
                    .expression;
                let dfx_dy = symbolic_cas::differentiate(fx, y, None)
                    .map_err(|e| format!("Failed to compute ∂F_x/∂y: {}", e))?
                    .expression;

                let curl_x = format!("({}) - ({})", dfz_dy, dfy_dz);
                let curl_y = format!("({}) - ({})", dfx_dz, dfz_dx);
                let curl_z = format!("({}) - ({})", dfy_dx, dfx_dy);

                let mut derivatives = std::collections::HashMap::new();
                derivatives.insert(
                    "curl".to_string(),
                    serde_json::json!(vec![curl_x.clone(), curl_y.clone(), curl_z.clone()]),
                );
                derivatives.insert("curl_x".to_string(), serde_json::json!(curl_x));
                derivatives.insert("curl_y".to_string(), serde_json::json!(curl_y));
                derivatives.insert("curl_z".to_string(), serde_json::json!(curl_z));

                Ok(DifferentiateOutput {
                    derivatives,
                    latex: Some(std::collections::HashMap::from([(
                        "curl".to_string(),
                        format!(
                            "\\nabla \\times \\vec{{F}} = ({}, {}, {})",
                            curl_x, curl_y, curl_z
                        ),
                    )])),
                    metadata: Some(serde_json::json!({
                        "operation": "curl",
                        "method": "symbolic",
                        "components": components
                    })),
                })
            }

            VectorCalcOp::Laplacian => {
                // Laplacian: ∇²f = ∂²f/∂x² + ∂²f/∂y² + ∂²f/∂z²
                let mut laplacian_terms = Vec::new();

                for var in &input.variables {
                    // Compute second derivative ∂²f/∂var²
                    let second_deriv = symbolic_cas::differentiate(&input.expression, var, Some(2))
                        .map_err(|e| format!("Failed to compute second derivative: {}", e))?;
                    laplacian_terms.push(second_deriv.expression);
                }

                let laplacian = laplacian_terms.join(" + ");

                let mut derivatives = std::collections::HashMap::new();
                derivatives.insert("laplacian".to_string(), serde_json::json!(laplacian));

                Ok(DifferentiateOutput {
                    derivatives,
                    latex: Some(std::collections::HashMap::from([(
                        "laplacian".to_string(),
                        format!("\\nabla^2 f = {}", laplacian),
                    )])),
                    metadata: Some(serde_json::json!({
                        "operation": "laplacian",
                        "method": "symbolic",
                        "dimensions": input.variables.len()
                    })),
                })
            }

            VectorCalcOp::Directional => {
                // Directional derivative: D_v f = ∇f · v = Σ(∂f/∂x_i * v_i)
                let direction: Vec<f64> = input
                    .parameters
                    .get("direction")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("Directional derivative requires 'direction' vector parameter")?;

                if direction.len() != input.variables.len() {
                    return Err(format!(
                        "Direction vector length {} must match variables {}",
                        direction.len(),
                        input.variables.len()
                    ));
                }

                // Compute gradient
                let gradient_results = symbolic_cas::gradient(&input.expression, &input.variables)
                    .map_err(|e| format!("Gradient computation failed: {}", e))?;

                // Compute dot product: ∇f · v
                let mut directional_terms = Vec::new();
                for (i, result) in gradient_results.iter().enumerate() {
                    let v_i = direction[i];
                    if v_i.abs() > 1e-10 {
                        let term = if v_i == 1.0 {
                            result.expression.clone()
                        } else {
                            format!("{} * ({})", v_i, result.expression)
                        };
                        directional_terms.push(term);
                    }
                }

                let directional_deriv = directional_terms.join(" + ");

                let mut derivatives = std::collections::HashMap::new();
                derivatives.insert(
                    "directional".to_string(),
                    serde_json::json!(directional_deriv),
                );

                Ok(DifferentiateOutput {
                    derivatives,
                    latex: Some(std::collections::HashMap::from([(
                        "directional".to_string(),
                        format!("D_{{\\vec{{v}}}} f = {}", directional_deriv),
                    )])),
                    metadata: Some(serde_json::json!({
                        "operation": "directional_derivative",
                        "method": "symbolic",
                        "direction": direction
                    })),
                })
            }
        }
    }

    /// Perform tensor calculus differentiation
    fn differentiate_tensor_calc(
        &self,
        op: &TensorDiffOp,
        input: &DifferentiateInput,
    ) -> ToolResult<DifferentiateOutput> {
        use crate::mathematics::symbolic_cas;

        match op {
            TensorDiffOp::Covariant => {
                // Covariant derivative: ∇_μ T^ν = ∂_μ T^ν + Γ^ν_μλ T^λ
                // Requires: metric tensor, coordinates, tensor field

                // Parse metric tensor from parameters
                let metric_data: Vec<Vec<f64>> = input.parameters.get("metric")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("Covariant derivative requires 'metric' tensor parameter [[g00, g01, ...], [g10, g11, ...], ...]")?;

                let coords: Vec<&str> = input.variables.iter().map(|s| s.as_str()).collect();

                // Build symbolic metric matrix
                let metric_exprs: Vec<Vec<symbolic_cas::Expr>> = metric_data
                    .iter()
                    .map(|row| {
                        row.iter()
                            .map(|&val| symbolic_cas::Expr::num(val as i64))
                            .collect()
                    })
                    .collect();

                let metric = symbolic_cas::SymbolicMatrix::new(metric_exprs)
                    .map_err(|e| format!("Invalid metric tensor: {}", e))?;

                // Compute Christoffel symbols
                let christoffel = symbolic_cas::christoffel_symbols(&metric, &coords)
                    .map_err(|e| format!("Failed to compute Christoffel symbols: {}", e))?;

                // Get tensor components
                let tensor_components: Vec<String> = input
                    .parameters
                    .get("tensor_components")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .unwrap_or_else(|| vec![input.expression.clone()]);

                // For a vector field (rank-1 contravariant tensor), compute ∇_μ V^ν
                let dimension = coords.len();
                let mut covariant_derivs = Vec::new();

                for mu in 0..dimension {
                    for nu in 0..dimension {
                        if nu < tensor_components.len() {
                            // Partial derivative term: ∂_μ V^ν
                            let partial = symbolic_cas::differentiate(
                                &tensor_components[nu],
                                coords[mu],
                                None,
                            )
                            .map_err(|e| format!("Failed to compute partial derivative: {}", e))?;

                            // Connection term: Γ^ν_μλ V^λ (sum over λ)
                            let mut connection_terms = Vec::new();
                            for lambda in 0..dimension {
                                if lambda < tensor_components.len() {
                                    let gamma =
                                        christoffel.get(&[nu, mu, lambda]).map_err(|e| {
                                            format!("Failed to get Christoffel symbol: {}", e)
                                        })?;

                                    let gamma_str = format!("{}", gamma);
                                    let component = &tensor_components[lambda];

                                    // Γ^ν_μλ * V^λ
                                    if gamma_str != "0" {
                                        connection_terms
                                            .push(format!("({}) * ({})", gamma_str, component));
                                    }
                                }
                            }

                            let connection_sum = if connection_terms.is_empty() {
                                "0".to_string()
                            } else {
                                connection_terms.join(" + ")
                            };

                            let covariant_deriv = if connection_sum == "0" {
                                partial.expression
                            } else {
                                format!("({}) + ({})", partial.expression, connection_sum)
                            };

                            covariant_derivs.push((mu, nu, covariant_deriv));
                        }
                    }
                }

                let mut derivatives = std::collections::HashMap::new();
                for (mu, nu, deriv) in &covariant_derivs {
                    derivatives.insert(
                        format!("∇_{}_V^{}", coords[*mu], nu),
                        serde_json::json!(deriv),
                    );
                }

                // Add full covariant derivative tensor representation
                let full_derivs: Vec<String> =
                    covariant_derivs.iter().map(|(_, _, d)| d.clone()).collect();
                derivatives.insert("covariant".to_string(), serde_json::json!(full_derivs));

                Ok(DifferentiateOutput {
                    derivatives,
                    latex: Some(std::collections::HashMap::from([(
                        "covariant".to_string(),
                        format!(
                            "\\nabla_\\mu V^\\nu = \\partial_\\mu V^\\nu + \\Gamma^\\nu_{{\\mu\\lambda}} V^\\lambda"
                        ),
                    )])),
                    metadata: Some(serde_json::json!({
                        "operation": "covariant_derivative",
                        "method": "symbolic_tensor_calculus",
                        "dimension": dimension,
                        "christoffel_computed": true
                    })),
                })
            }

            TensorDiffOp::Lie => {
                // Lie derivative: L_X Y along vector field X
                // For scalar: L_X f = X · ∇f = X^μ ∂_μ f
                // For vector: L_X Y^μ = X^ν ∂_ν Y^μ - Y^ν ∂_ν X^μ (Lie bracket)

                let vector_field: Vec<String> = input
                    .parameters
                    .get("vector_field")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("Lie derivative requires 'vector_field' parameter [X^0, X^1, ...]")?;

                if vector_field.len() != input.variables.len() {
                    return Err(format!(
                        "Vector field dimension {} must match variables {}",
                        vector_field.len(),
                        input.variables.len()
                    ));
                }

                // Check if we're operating on a scalar or vector
                let target_components: Option<Vec<String>> = input
                    .parameters
                    .get("target_components")
                    .and_then(|v| serde_json::from_value(v.clone()).ok());

                let mut derivatives = std::collections::HashMap::new();

                match target_components {
                    None => {
                        // Scalar function: L_X f = X · ∇f = Σ X^μ ∂_μ f
                        let mut lie_terms = Vec::new();

                        for (i, var) in input.variables.iter().enumerate() {
                            let partial = symbolic_cas::differentiate(&input.expression, var, None)
                                .map_err(|e| format!("Failed to compute partial: {}", e))?;

                            let term = format!("({}) * ({})", vector_field[i], partial.expression);
                            lie_terms.push(term);
                        }

                        let lie_deriv = lie_terms.join(" + ");
                        derivatives.insert("lie".to_string(), serde_json::json!(lie_deriv));

                        Ok(DifferentiateOutput {
                            derivatives,
                            latex: Some(std::collections::HashMap::from([(
                                "lie".to_string(),
                                format!(
                                    "\\mathcal{{L}}_X f = X^\\mu \\partial_\\mu f = {}",
                                    lie_deriv
                                ),
                            )])),
                            metadata: Some(serde_json::json!({
                                "operation": "lie_derivative",
                                "type": "scalar",
                                "method": "directional_derivative"
                            })),
                        })
                    }
                    Some(target) => {
                        // Vector field: L_X Y^μ = X^ν ∂_ν Y^μ - Y^ν ∂_ν X^μ
                        let dimension = input.variables.len();
                        let mut lie_components = Vec::new();

                        for mu in 0..dimension {
                            if mu >= target.len() {
                                break;
                            }

                            // First term: X^ν ∂_ν Y^μ (sum over ν)
                            let mut first_terms = Vec::new();
                            for (nu, var) in input.variables.iter().enumerate() {
                                let partial_y = symbolic_cas::differentiate(&target[mu], var, None)
                                    .map_err(|e| format!("Failed to compute ∂_ν Y^μ: {}", e))?;

                                first_terms.push(format!(
                                    "({}) * ({})",
                                    vector_field[nu], partial_y.expression
                                ));
                            }

                            // Second term: -Y^ν ∂_ν X^μ (sum over ν)
                            let mut second_terms = Vec::new();
                            for (nu, var) in input.variables.iter().enumerate() {
                                let partial_x =
                                    symbolic_cas::differentiate(&vector_field[mu], var, None)
                                        .map_err(|e| format!("Failed to compute ∂_ν X^μ: {}", e))?;

                                second_terms
                                    .push(format!("({}) * ({})", target[nu], partial_x.expression));
                            }

                            let first_sum = first_terms.join(" + ");
                            let second_sum = second_terms.join(" + ");
                            let lie_component = format!("({}) - ({})", first_sum, second_sum);

                            lie_components.push(lie_component.clone());
                            derivatives
                                .insert(format!("L_X Y^{}", mu), serde_json::json!(lie_component));
                        }

                        derivatives.insert("lie".to_string(), serde_json::json!(lie_components));

                        Ok(DifferentiateOutput {
                            derivatives,
                            latex: Some(std::collections::HashMap::from([
                                ("lie".to_string(), "\\mathcal{L}_X Y^\\mu = X^\\nu \\partial_\\nu Y^\\mu - Y^\\nu \\partial_\\nu X^\\mu".to_string())
                            ])),
                            metadata: Some(serde_json::json!({
                                "operation": "lie_derivative",
                                "type": "vector_field",
                                "method": "lie_bracket",
                                "dimension": dimension
                            })),
                        })
                    }
                }
            }

            TensorDiffOp::ExteriorDerivative => {
                // Exterior derivative: d(ω) for differential forms
                // For 0-form (scalar): df = Σ(∂f/∂x_i) dx_i (this is the gradient as a 1-form)
                // For 1-form: d(Σ f_i dx_i) = Σ_i Σ_j (∂f_i/∂x_j) dx_j ∧ dx_i
                // Property: d² = 0

                let form_degree = input
                    .parameters
                    .get("form_degree")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(0) as usize;

                let mut derivatives = std::collections::HashMap::new();

                match form_degree {
                    0 => {
                        // 0-form (scalar function): df = gradient as 1-form
                        let gradient_results =
                            symbolic_cas::gradient(&input.expression, &input.variables)
                                .map_err(|e| format!("Failed to compute gradient: {}", e))?;

                        let mut form_terms = Vec::new();
                        for (i, result) in gradient_results.iter().enumerate() {
                            if let Some(var) = input.variables.get(i) {
                                form_terms.push(format!("({}) d{}", result.expression, var));
                            }
                        }

                        let exterior_deriv = form_terms.join(" + ");
                        derivatives
                            .insert("exterior".to_string(), serde_json::json!(exterior_deriv));

                        Ok(DifferentiateOutput {
                            derivatives,
                            latex: Some(std::collections::HashMap::from([(
                                "exterior".to_string(),
                                format!("df = \\sum_i \\frac{{\\partial f}}{{\\partial x_i}} dx_i"),
                            )])),
                            metadata: Some(serde_json::json!({
                                "operation": "exterior_derivative",
                                "form_degree": 0,
                                "result_degree": 1,
                                "representation": "1-form",
                                "property": "d² = 0"
                            })),
                        })
                    }
                    1 => {
                        // 1-form: ω = Σ f_i dx_i
                        // dω = Σ_i<j (∂f_j/∂x_i - ∂f_i/∂x_j) dx_i ∧ dx_j

                        let form_components: Vec<String> = input.parameters.get("form_components")
                            .and_then(|v| serde_json::from_value(v.clone()).ok())
                            .ok_or("1-form exterior derivative requires 'form_components' [f_1, f_2, ...]")?;

                        let n = input.variables.len();
                        let mut wedge_terms = Vec::new();

                        for i in 0..n {
                            for j in (i + 1)..n {
                                if i < form_components.len() && j < form_components.len() {
                                    // ∂f_j/∂x_i
                                    let df_j_dx_i = symbolic_cas::differentiate(
                                        &form_components[j],
                                        &input.variables[i],
                                        None,
                                    )
                                    .map_err(|e| format!("Failed to compute ∂f_j/∂x_i: {}", e))?;

                                    // ∂f_i/∂x_j
                                    let df_i_dx_j = symbolic_cas::differentiate(
                                        &form_components[i],
                                        &input.variables[j],
                                        None,
                                    )
                                    .map_err(|e| format!("Failed to compute ∂f_i/∂x_j: {}", e))?;

                                    let wedge_coeff = format!(
                                        "({}) - ({})",
                                        df_j_dx_i.expression, df_i_dx_j.expression
                                    );
                                    let wedge_term = format!(
                                        "({}) d{} ∧ d{}",
                                        wedge_coeff, input.variables[i], input.variables[j]
                                    );
                                    wedge_terms.push(wedge_term);
                                }
                            }
                        }

                        let exterior_deriv = wedge_terms.join(" + ");
                        derivatives
                            .insert("exterior".to_string(), serde_json::json!(exterior_deriv));

                        Ok(DifferentiateOutput {
                            derivatives,
                            latex: Some(std::collections::HashMap::from([(
                                "exterior".to_string(),
                                format!(
                                    "d\\omega = \\sum_{{i<j}} \\left(\\frac{{\\partial f_j}}{{\\partial x_i}} - \\frac{{\\partial f_i}}{{\\partial x_j}}\\right) dx_i \\wedge dx_j"
                                ),
                            )])),
                            metadata: Some(serde_json::json!({
                                "operation": "exterior_derivative",
                                "form_degree": 1,
                                "result_degree": 2,
                                "representation": "2-form",
                                "property": "d² = 0",
                                "antisymmetric": true
                            })),
                        })
                    }
                    _ => {
                        // Higher degree forms - return symbolic representation
                        let exterior_deriv = format!("d({})", input.expression);
                        derivatives
                            .insert("exterior".to_string(), serde_json::json!(exterior_deriv));

                        Ok(DifferentiateOutput {
                            derivatives,
                            latex: Some(std::collections::HashMap::from([(
                                "exterior".to_string(),
                                format!(
                                    "d(\\omega_{}) \\text{{ where }} \\omega \\text{{ is a {}-form}}",
                                    form_degree, form_degree
                                ),
                            )])),
                            metadata: Some(serde_json::json!({
                                "operation": "exterior_derivative",
                                "form_degree": form_degree,
                                "result_degree": form_degree + 1,
                                "property": "d² = 0",
                                "note": "Symbolic representation - full computation requires explicit wedge products"
                            })),
                        })
                    }
                }
            }
        }
    }

    /// Perform numeric differentiation using finite differences
    fn differentiate_numeric(&self, input: &DifferentiateInput) -> ToolResult<DifferentiateOutput> {
        use crate::tools::numerical_methods;

        // Extract x and y values from parameters or expression
        let x_values: Vec<f64> = input
            .parameters
            .get("x_values")
            .and_then(|v| serde_json::from_value(v.clone()).ok())
            .ok_or("x_values required for numeric differentiation")?;

        let y_values: Vec<f64> = input
            .parameters
            .get("y_values")
            .and_then(|v| serde_json::from_value(v.clone()).ok())
            .ok_or("y_values required for numeric differentiation")?;

        let order = input
            .order
            .as_ref()
            .and_then(|orders| orders.first())
            .copied()
            .unwrap_or(1);

        let method = input
            .parameters
            .get("method")
            .and_then(|v| v.as_str())
            .unwrap_or("central");

        let result = numerical_methods::differentiate(numerical_methods::DifferentiationRequest {
            method: method.to_string(),
            x_values: x_values.clone(),
            y_values: y_values.clone(),
            order,
        })
        .map_err(|e| e.to_string())?;

        let mut derivatives = std::collections::HashMap::new();
        let var_name = input.variables.first().unwrap_or(&"x".to_string()).clone();

        derivatives.insert(var_name.clone(), serde_json::json!(result.derivatives));

        Ok(DifferentiateOutput {
            derivatives,
            latex: Some(std::collections::HashMap::from([(
                var_name,
                format!("d{}/dx", input.expression),
            )])),
            metadata: Some(serde_json::json!({
                "method": result.method_used,
                "order": order,
                "points": x_values.len()
            })),
        })
    }

    /// Perform symbolic differentiation
    fn differentiate_symbolic(
        &self,
        input: &DifferentiateInput,
    ) -> ToolResult<DifferentiateOutput> {
        use crate::mathematics::symbolic_cas;

        let expr = input.expression.trim();
        let mut derivatives = std::collections::HashMap::new();
        let mut latex_map = std::collections::HashMap::new();

        // Differentiate with respect to each variable
        for var in &input.variables {
            let order = input
                .parameters
                .get("order")
                .and_then(|v| v.as_u64())
                .map(|o| o as usize);

            // Use symbolic CAS for differentiation
            let result = symbolic_cas::differentiate(expr, var, order)
                .map_err(|e| format!("Symbolic differentiation failed: {}", e))?;

            derivatives.insert(var.clone(), serde_json::json!(result.expression));
            if let Some(latex) = result.latex {
                latex_map.insert(var.clone(), latex);
            }
        }

        Ok(DifferentiateOutput {
            derivatives,
            latex: if latex_map.is_empty() {
                None
            } else {
                Some(latex_map)
            },
            metadata: Some(serde_json::json!({
                "operation": "symbolic_differentiation",
                "method": "symbolic_cas"
            })),
        })
    }

    /// Perform variational calculus
    fn differentiate_variational(
        &self,
        input: &DifferentiateInput,
    ) -> ToolResult<DifferentiateOutput> {
        use crate::mathematics::symbolic_cas;

        // Variational calculus - Euler-Lagrange equations: d/dx(∂L/∂y') - ∂L/∂y = 0
        // Lagrangian L depends on (x, y, y') where y' = dy/dx

        let lagrangian = input.expression.clone();

        // Get the independent variable (typically 'x' or 't')
        let independent_var = input
            .parameters
            .get("independent_var")
            .and_then(|v| v.as_str())
            .unwrap_or("x");

        // Get the dependent variable (typically 'y')
        let dependent_var = input
            .parameters
            .get("dependent_var")
            .and_then(|v| v.as_str())
            .unwrap_or("y");

        // Derivative variable (typically 'y_prime' or 'dy')
        let derivative_var = input
            .parameters
            .get("derivative_var")
            .and_then(|v| v.as_str())
            .unwrap_or("y_prime");

        let mut derivatives = std::collections::HashMap::new();

        // Compute ∂L/∂y
        let d_l_d_y = symbolic_cas::differentiate(&lagrangian, dependent_var, None)
            .map_err(|e| format!("Failed to compute ∂L/∂y: {}", e))?;

        derivatives.insert(
            "∂L/∂y".to_string(),
            serde_json::json!(d_l_d_y.expression.clone()),
        );

        // Compute ∂L/∂y'
        let d_l_d_yprime = symbolic_cas::differentiate(&lagrangian, derivative_var, None)
            .map_err(|e| format!("Failed to compute ∂L/∂y': {}", e))?;

        derivatives.insert(
            "∂L/∂y'".to_string(),
            serde_json::json!(d_l_d_yprime.expression.clone()),
        );

        // For the full Euler-Lagrange equation, we need d/dx(∂L/∂y')
        // This requires treating ∂L/∂y' as a function of x (via chain rule)
        // d/dx(∂L/∂y') = ∂²L/∂y∂x + ∂²L/∂y²·y' + ∂²L/∂y'∂y'·y''

        // For now, compute the partial derivative with respect to x
        let d_dx_partial =
            symbolic_cas::differentiate(&d_l_d_yprime.expression, independent_var, None);

        let d_dx_term = match d_dx_partial {
            Ok(result) => result.expression,
            Err(_) => format!("d/dx({})", d_l_d_yprime.expression),
        };

        derivatives.insert(
            "d/dx(∂L/∂y')".to_string(),
            serde_json::json!(d_dx_term.clone()),
        );

        // Form the Euler-Lagrange equation: d/dx(∂L/∂y') - ∂L/∂y = 0
        let euler_lagrange = format!("({}) - ({}) = 0", d_dx_term, d_l_d_y.expression);
        derivatives.insert(
            "euler_lagrange".to_string(),
            serde_json::json!(euler_lagrange.clone()),
        );

        // Compute second derivatives for more complete analysis
        // ∂²L/∂y²
        let d2_l_d_y2 = symbolic_cas::differentiate(&d_l_d_y.expression, dependent_var, None);
        if let Ok(result) = d2_l_d_y2 {
            derivatives.insert("∂²L/∂y²".to_string(), serde_json::json!(result.expression));
        }

        // ∂²L/∂y'²
        let d2_l_d_yprime2 =
            symbolic_cas::differentiate(&d_l_d_yprime.expression, derivative_var, None);
        if let Ok(result) = d2_l_d_yprime2 {
            derivatives.insert("∂²L/∂y'²".to_string(), serde_json::json!(result.expression));
        }

        // ∂²L/∂y∂y' (mixed partial)
        let d2_l_d_y_d_yprime =
            symbolic_cas::differentiate(&d_l_d_y.expression, derivative_var, None);
        if let Ok(result) = d2_l_d_y_d_yprime {
            derivatives.insert(
                "∂²L/∂y∂y'".to_string(),
                serde_json::json!(result.expression),
            );
        }

        Ok(DifferentiateOutput {
            derivatives,
            latex: Some(std::collections::HashMap::from([
                (
                    "euler_lagrange".to_string(),
                    format!(
                        "\\frac{{d}}{{d{}}}\\left(\\frac{{\\partial L}}{{\\partial y'}}\\right) - \\frac{{\\partial L}}{{\\partial {}}} = 0",
                        independent_var, dependent_var
                    ),
                ),
                (
                    "∂L/∂y".to_string(),
                    format!(
                        "\\frac{{\\partial L}}{{\\partial {}}} = {}",
                        dependent_var, d_l_d_y.expression
                    ),
                ),
                (
                    "∂L/∂y'".to_string(),
                    format!(
                        "\\frac{{\\partial L}}{{\\partial y'}} = {}",
                        d_l_d_yprime.expression
                    ),
                ),
            ])),
            metadata: Some(serde_json::json!({
                "operation": "variational_calculus",
                "lagrangian": lagrangian,
                "independent_variable": independent_var,
                "dependent_variable": dependent_var,
                "method": "euler_lagrange_equation",
                "note": "Full d/dx(∂L/∂y') requires chain rule expansion with second derivatives"
            })),
        })
    }

    /// Perform differential forms operations
    fn differentiate_forms(&self, input: &DifferentiateInput) -> ToolResult<DifferentiateOutput> {
        // Differential forms - simplified interface that routes to exterior derivative
        // For full exterior derivative computation, use TensorCalc(ExteriorDerivative)

        // This provides a simplified interface for common differential forms operations
        // Auto-detect form degree from expression structure

        let form_degree = input
            .parameters
            .get("form_degree")
            .and_then(|v| v.as_u64())
            .unwrap_or(0) as usize;

        // If form_degree is 0 and no explicit form_components, treat as scalar function
        if form_degree == 0 && !input.parameters.contains_key("form_components") {
            // Route to exterior derivative implementation for 0-forms (scalar functions)
            let modified_input = DifferentiateInput {
                operation: DifferentiationOp::TensorCalc(TensorDiffOp::ExteriorDerivative),
                expression: input.expression.clone(),
                variables: input.variables.clone(),
                order: input.order.clone(),
                evaluate_at: input.evaluate_at.clone(),
                parameters: input.parameters.clone(),
            };
            return self
                .differentiate_tensor_calc(&TensorDiffOp::ExteriorDerivative, &modified_input);
        }

        // For 1-forms and higher, provide guidance
        if form_degree >= 1 {
            let mut derivatives = std::collections::HashMap::new();

            // Check if form_components are provided
            if input.parameters.contains_key("form_components") {
                // Route to full exterior derivative implementation
                let modified_input = DifferentiateInput {
                    operation: DifferentiationOp::TensorCalc(TensorDiffOp::ExteriorDerivative),
                    expression: input.expression.clone(),
                    variables: input.variables.clone(),
                    order: input.order.clone(),
                    evaluate_at: input.evaluate_at.clone(),
                    parameters: input.parameters.clone(),
                };
                return self
                    .differentiate_tensor_calc(&TensorDiffOp::ExteriorDerivative, &modified_input);
            }

            // Otherwise provide symbolic representation
            let exterior_deriv = format!("d({})", input.expression);
            derivatives.insert(
                "exterior_derivative".to_string(),
                serde_json::json!(exterior_deriv.clone()),
            );

            return Ok(DifferentiateOutput {
                derivatives,
                latex: Some(std::collections::HashMap::from([(
                    "exterior_derivative".to_string(),
                    format!(
                        "d(\\omega) \\text{{ where }} \\omega \\text{{ is a {}-form}}",
                        form_degree
                    ),
                )])),
                metadata: Some(serde_json::json!({
                    "operation": "differential_forms",
                    "form_degree": form_degree,
                    "result_degree": form_degree + 1,
                    "property": "d² = 0",
                    "note": "Provide 'form_components' parameter for explicit computation. Use TensorCalc(ExteriorDerivative) for full control."
                })),
            });
        }

        // Default case - provide symbolic representation
        use crate::mathematics::symbolic_cas;

        let mut derivatives = std::collections::HashMap::new();

        // Compute gradient as 1-form representation
        let gradient_results = symbolic_cas::gradient(&input.expression, &input.variables)
            .map_err(|e| format!("Failed to compute gradient: {}", e))?;

        let mut form_terms = Vec::new();
        for (i, result) in gradient_results.iter().enumerate() {
            if let Some(var) = input.variables.get(i) {
                form_terms.push(format!("({}) d{}", result.expression, var));
            }
        }

        let differential = form_terms.join(" + ");
        derivatives.insert(
            "differential".to_string(),
            serde_json::json!(differential.clone()),
        );
        derivatives.insert(
            "exterior_derivative".to_string(),
            serde_json::json!(differential.clone()),
        );

        Ok(DifferentiateOutput {
            derivatives,
            latex: Some(std::collections::HashMap::from([(
                "differential".to_string(),
                format!("df = \\sum_i \\frac{{\\partial f}}{{\\partial x_i}} dx_i"),
            )])),
            metadata: Some(serde_json::json!({
                "operation": "differential_forms",
                "form_degree": 0,
                "result_degree": 1,
                "representation": "1-form (differential)",
                "property": "d² = 0",
                "note": "Computed as gradient expressed as 1-form"
            })),
        })
    }
}

impl Differentiate for UnifiedDifferentiator {
    fn differentiate(&self, input: &DifferentiateInput) -> ToolResult<DifferentiateOutput> {
        match &input.operation {
            DifferentiationOp::VectorCalc(op) => self.differentiate_vector_calc(op, input),

            DifferentiationOp::TensorCalc(op) => self.differentiate_tensor_calc(op, input),

            DifferentiationOp::Numeric => self.differentiate_numeric(input),

            DifferentiationOp::Symbolic => self.differentiate_symbolic(input),

            DifferentiationOp::Variational => self.differentiate_variational(input),

            DifferentiationOp::DifferentialForms => self.differentiate_forms(input),
        }
    }
}
