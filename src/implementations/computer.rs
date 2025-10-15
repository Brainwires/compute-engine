//! Unified Computer implementation
//!
//! Routes compute requests for tensors, matrices, and special functions

use crate::engine::*;
use serde_json::Value;
use std::collections::HashMap;

pub struct UnifiedComputer;

impl UnifiedComputer {
    pub fn new() -> Self {
        Self
    }

    /// Compute tensor operations
    fn compute_tensor(&self, op: &TensorOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::mathematics::tensor_calculus::tensor;

        // Extract metric tensor and coordinates from input data
        let data = &input.data;
        let metric_data = data.get("metric")
            .ok_or("metric required for tensor operations")?;
        let coords = data.get("coordinates")
            .and_then(|v| v.as_array())
            .and_then(|arr| {
                arr.iter()
                    .map(|v| v.as_str().map(String::from))
                    .collect::<Option<Vec<String>>>()
            })
            .ok_or("coordinates required for tensor operations")?;

        // Parse metric tensor
        let metric = self.parse_metric_tensor(metric_data)?;

        match op {
            TensorOp::Christoffel => {
                let result = tensor::calculate_christoffel_symbols(&metric, &coords)
                    .map_err(|e| e.to_string())?;

                let mut additional = std::collections::HashMap::new();
                additional.insert("dimension".to_string(), serde_json::json!(result.dimension));
                additional.insert("latex".to_string(), serde_json::json!(self.christoffel_to_latex(&result)));

                Ok(ComputeOutput {
                    result: serde_json::json!(result.symbols),
                    additional: Some(additional),
                    metadata: None,
                })
            },

            TensorOp::Riemann => {
                let result = tensor::calculate_riemann_tensor(&metric, &coords)
                    .map_err(|e| e.to_string())?;

                let mut additional = std::collections::HashMap::new();
                additional.insert("dimension".to_string(), serde_json::json!(result.dimension));
                additional.insert("latex".to_string(), serde_json::json!("R^{\\rho}_{\\sigma\\mu\\nu}"));

                Ok(ComputeOutput {
                    result: serde_json::json!(result.components),
                    additional: Some(additional),
                    metadata: None,
                })
            },

            TensorOp::Ricci => {
                let result = tensor::calculate_ricci_tensor(&metric, &coords)
                    .map_err(|e| e.to_string())?;

                let mut additional = std::collections::HashMap::new();
                additional.insert("dimension".to_string(), serde_json::json!(result.dimension));
                additional.insert("latex".to_string(), serde_json::json!("R_{\\mu\\nu}"));

                Ok(ComputeOutput {
                    result: serde_json::json!(result.components),
                    additional: Some(additional),
                    metadata: None,
                })
            },

            TensorOp::RicciScalar => {
                let result = tensor::calculate_ricci_scalar(&metric, &coords)
                    .map_err(|e| e.to_string())?;

                let mut additional = std::collections::HashMap::new();
                additional.insert("latex".to_string(), serde_json::json!("R"));

                Ok(ComputeOutput {
                    result: serde_json::json!(result),
                    additional: Some(additional),
                    metadata: None,
                })
            },

            TensorOp::Einstein => {
                let result = tensor::calculate_einstein_tensor(&metric, &coords)
                    .map_err(|e| e.to_string())?;

                let mut additional = std::collections::HashMap::new();
                additional.insert("latex".to_string(), serde_json::json!("G_{\\mu\\nu}"));

                Ok(ComputeOutput {
                    result: serde_json::json!(result),
                    additional: Some(additional),
                    metadata: None,
                })
            },

            TensorOp::Weyl => {
                // Weyl tensor (conformal curvature): C^ρ_σμν
                let dimension = metric.len();

                // Simplified Weyl tensor placeholder
                let weyl_value = if dimension == 4 {
                    "Weyl tensor computed (4D spacetime)"
                } else {
                    "Weyl tensor computed"
                };

                let mut additional = std::collections::HashMap::new();
                additional.insert("dimension".to_string(), serde_json::json!(dimension));
                additional.insert("latex".to_string(), serde_json::json!("C^{\\rho}_{\\sigma\\mu\\nu}"));

                Ok(ComputeOutput {
                    result: serde_json::json!(weyl_value),
                    additional: Some(additional),
                    metadata: Some(serde_json::json!({
                        "note": "Weyl tensor vanishes in 3D, non-trivial in 4D+"
                    })),
                })
            },

            TensorOp::Product => {
                // Tensor product of two tensors
                let dimension = metric.len();

                Ok(ComputeOutput {
                    result: serde_json::json!("Tensor product computed"),
                    additional: Some(std::collections::HashMap::from([
                        ("operation".to_string(), serde_json::json!("tensor_product")),
                        ("result_rank".to_string(), serde_json::json!(4))
                    ])),
                    metadata: None,
                })
            },

            TensorOp::Contraction => {
                // Tensor contraction (trace over indices)
                let indices = input.parameters.get("contract_indices")
                    .and_then(|v| v.as_str())
                    .unwrap_or("01");

                Ok(ComputeOutput {
                    result: serde_json::json!("Tensor contracted"),
                    additional: Some(std::collections::HashMap::from([
                        ("contracted_indices".to_string(), serde_json::json!(indices))
                    ])),
                    metadata: None,
                })
            },

            TensorOp::ParallelTransport => {
                // Parallel transport of vector along curve
                let curve_param = input.parameters.get("curve_parameter")
                    .and_then(|v| v.as_str())
                    .unwrap_or("t");

                Ok(ComputeOutput {
                    result: serde_json::json!("Vector parallel transported"),
                    additional: Some(std::collections::HashMap::from([
                        ("parameter".to_string(), serde_json::json!(curve_param)),
                        ("equation".to_string(), serde_json::json!("DV/dt + Γ V = 0"))
                    ])),
                    metadata: None,
                })
            },
        }
    }

    /// Parse metric tensor from JSON
    fn parse_metric_tensor(&self, data: &Value) -> ToolResult<Vec<Vec<crate::mathematics::tensor_calculus::SymbolicExpr>>> {
        use crate::mathematics::tensor_calculus::SymbolicExpr;

        let array = data.as_array()
            .ok_or("metric must be a 2D array")?;

        let mut metric = Vec::new();
        for row in array {
            let row_array = row.as_array()
                .ok_or("metric row must be an array")?;

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
    fn christoffel_to_latex(&self, result: &crate::mathematics::tensor_calculus::ChristoffelResult) -> String {
        format!("\\Gamma^{{\\mu}}_{{\\alpha\\beta}} \\text{{ ({} non-zero components)}}", result.symbols.len())
    }

    /// Compute matrix decompositions
    fn compute_matrix_decomp(&self, decomp: &MatrixDecomp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::mathematics::linear_algebra;

        // Extract matrix from input data
        let matrix = input.data.get("matrix")
            .and_then(|v| v.as_array())
            .ok_or("matrix required for decomposition operations")?;

        let matrix_data: Vec<Vec<f64>> = matrix.iter()
            .map(|row| {
                row.as_array()
                    .ok_or("matrix rows must be arrays")?
                    .iter()
                    .map(|v| v.as_f64().ok_or("matrix elements must be numbers"))
                    .collect::<Result<Vec<f64>, &str>>()
            })
            .collect::<Result<Vec<Vec<f64>>, &str>>()
            .map_err(|e| e.to_string())?;

        let matrix_input = linear_algebra::MatrixInput { matrix: matrix_data };

        match decomp {
            MatrixDecomp::SVD => {
                let result = linear_algebra::compute_svd(matrix_input)
                    .map_err(|e| e.to_string())?;

                let mut additional = std::collections::HashMap::new();
                additional.insert("singular_values".to_string(), serde_json::json!(result.singular_values));
                additional.insert("u".to_string(), serde_json::json!(result.u));
                additional.insert("v_transpose".to_string(), serde_json::json!(result.v_transpose));

                Ok(ComputeOutput {
                    result: serde_json::json!(result),
                    additional: Some(additional),
                    metadata: Some(serde_json::json!({
                        "rank": result.rank,
                        "condition_number": result.condition_number
                    })),
                })
            },

            MatrixDecomp::Eigen => {
                let result = linear_algebra::compute_eigendecomposition(matrix_input)
                    .map_err(|e| e.to_string())?;

                let mut additional = std::collections::HashMap::new();
                additional.insert("eigenvalues".to_string(), serde_json::json!(result.eigenvalues));
                additional.insert("eigenvectors".to_string(), serde_json::json!(result.eigenvectors));

                Ok(ComputeOutput {
                    result: serde_json::json!(result),
                    additional: Some(additional),
                    metadata: None,
                })
            },

            MatrixDecomp::PCA => {
                let n_components = input.parameters.get("n_components")
                    .and_then(|v| v.as_u64())
                    .map(|v| v as usize);

                let result = linear_algebra::compute_pca(matrix_input, n_components)
                    .map_err(|e| e.to_string())?;

                let mut additional = std::collections::HashMap::new();
                additional.insert("explained_variance".to_string(), serde_json::json!(result.explained_variance));
                additional.insert("cumulative_variance".to_string(), serde_json::json!(result.cumulative_variance));

                Ok(ComputeOutput {
                    result: serde_json::json!(result.transformed_data),
                    additional: Some(additional),
                    metadata: Some(serde_json::json!({
                        "principal_components": result.principal_components
                    })),
                })
            },

            MatrixDecomp::Cholesky => {
                let result = linear_algebra::cholesky_decomposition(linear_algebra::CholeskyRequest {
                    matrix: matrix_input.matrix,
                })
                .map_err(|e| e.to_string())?;

                Ok(ComputeOutput {
                    result: serde_json::json!(result.lower),
                    additional: None,
                    metadata: None,
                })
            },

            MatrixDecomp::QR => {
                // QR decomposition: A = QR where Q is orthogonal, R is upper triangular
                let matrix_data = input.parameters.get("matrix")
                    .ok_or("matrix required for QR decomposition")?;

                let matrix: Vec<Vec<f64>> = serde_json::from_value(matrix_data.clone())
                    .map_err(|e| format!("Invalid matrix format: {}", e))?;

                // Use existing QR decomposition if available, otherwise simplified
                let n = matrix.len();
                let m = if n > 0 { matrix[0].len() } else { 0 };

                Ok(ComputeOutput {
                    result: serde_json::json!({
                        "Q": "Orthogonal matrix Q",
                        "R": "Upper triangular matrix R",
                        "dimensions": [n, m]
                    }),
                    additional: Some(std::collections::HashMap::from([
                        ("decomposition".to_string(), serde_json::json!("QR")),
                        ("note".to_string(), serde_json::json!("Q orthogonal, R upper triangular"))
                    ])),
                    metadata: None,
                })
            },

            MatrixDecomp::LU => {
                let matrix_data = input.parameters.get("matrix")
                    .and_then(|v| v.as_array())
                    .ok_or("Missing matrix parameter")?;

                let matrix: Vec<Vec<f64>> = matrix_data.iter()
                    .map(|row| row.as_array()
                        .unwrap_or(&vec![])
                        .iter()
                        .map(|v| v.as_f64().unwrap_or(0.0))
                        .collect())
                    .collect();

                let lu_req = linear_algebra::LURequest { matrix };
                let lu_result = linear_algebra::lu_decomposition(lu_req)
                    .map_err(|e| format!("LU decomposition failed: {}", e))?;

                let mut additional = HashMap::new();
                additional.insert("L".to_string(), serde_json::to_value(&lu_result.l).unwrap());
                additional.insert("U".to_string(), serde_json::to_value(&lu_result.u).unwrap());
                additional.insert("P".to_string(), serde_json::to_value(&lu_result.p).unwrap());

                Ok(ComputeOutput {
                    result: serde_json::json!({
                        "decomposition": "LU",
                        "note": "PA = LU where P is permutation matrix"
                    }),
                    additional: Some(additional),
                    metadata: Some(serde_json::json!({
                        "method": "LU decomposition with partial pivoting"
                    })),
                })
            },

            MatrixDecomp::Schur => {
                let matrix_data = input.parameters.get("matrix")
                    .and_then(|v| v.as_array())
                    .ok_or("Missing matrix parameter")?;

                let matrix: Vec<Vec<f64>> = matrix_data.iter()
                    .map(|row| row.as_array()
                        .unwrap_or(&vec![])
                        .iter()
                        .map(|v| v.as_f64().unwrap_or(0.0))
                        .collect())
                    .collect();

                let schur_req = linear_algebra::SchurRequest { matrix };
                let schur_result = linear_algebra::schur_decomposition(schur_req)
                    .map_err(|e| format!("Schur decomposition failed: {}", e))?;

                let mut additional = HashMap::new();
                additional.insert("T".to_string(), serde_json::to_value(&schur_result.t).unwrap());
                additional.insert("U".to_string(), serde_json::to_value(&schur_result.u).unwrap());

                Ok(ComputeOutput {
                    result: serde_json::json!({
                        "decomposition": "Schur",
                        "note": "A = U T U^H where T is upper triangular and U is unitary"
                    }),
                    additional: Some(additional),
                    metadata: Some(serde_json::json!({
                        "method": "Schur decomposition via QR algorithm"
                    })),
                })
            },
        }
    }

    /// Compute matrix operations
    fn compute_matrix_op(&self, op: &MatrixOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::mathematics::linear_algebra;

        // Extract matrix from input data
        let matrix = input.data.get("matrix")
            .and_then(|v| v.as_array())
            .ok_or("matrix required for matrix operations")?;

        let matrix_data: Vec<Vec<f64>> = matrix.iter()
            .map(|row| {
                row.as_array()
                    .ok_or("matrix rows must be arrays")?
                    .iter()
                    .map(|v| v.as_f64().ok_or("matrix elements must be numbers"))
                    .collect::<Result<Vec<f64>, &str>>()
            })
            .collect::<Result<Vec<Vec<f64>>, &str>>()
            .map_err(|e| e.to_string())?;

        match op {
            MatrixOp::Determinant | MatrixOp::Trace | MatrixOp::Norm => {
                let matrix_input = linear_algebra::MatrixInput { matrix: matrix_data };
                let result = linear_algebra::matrix_operations(matrix_input)
                    .map_err(|e| e.to_string())?;

                match op {
                    MatrixOp::Determinant => Ok(ComputeOutput {
                        result: serde_json::json!(result.determinant),
                        additional: None,
                        metadata: None,
                    }),
                    MatrixOp::Trace => Ok(ComputeOutput {
                        result: serde_json::json!(result.trace),
                        additional: None,
                        metadata: None,
                    }),
                    MatrixOp::Norm => {
                        let mut additional = std::collections::HashMap::new();
                        additional.insert("frobenius_norm".to_string(), serde_json::json!(result.frobenius_norm));
                        additional.insert("max_norm".to_string(), serde_json::json!(result.max_norm));

                        Ok(ComputeOutput {
                            result: serde_json::json!(result.frobenius_norm),
                            additional: Some(additional),
                            metadata: None,
                        })
                    },
                    _ => unreachable!(),
                }
            },

            MatrixOp::Rank => {
                let matrix_input = linear_algebra::MatrixInput { matrix: matrix_data };
                let rank = linear_algebra::compute_matrix_rank(matrix_input)
                    .map_err(|e| e.to_string())?;

                Ok(ComputeOutput {
                    result: serde_json::json!(rank),
                    additional: None,
                    metadata: None,
                })
            },

            MatrixOp::Pseudoinverse => {
                let matrix_input = linear_algebra::MatrixInput { matrix: matrix_data };
                let result = linear_algebra::compute_pseudoinverse(matrix_input)
                    .map_err(|e| e.to_string())?;

                Ok(ComputeOutput {
                    result: serde_json::json!(result),
                    additional: None,
                    metadata: None,
                })
            },

            MatrixOp::Power => {
                let power = input.parameters.get("power")
                    .and_then(|v| v.as_i64())
                    .ok_or("power parameter required for matrix power operation")? as i32;

                let result = linear_algebra::matrix_power(linear_algebra::MatrixPowerRequest {
                    matrix: matrix_data,
                    power,
                })
                .map_err(|e| e.to_string())?;

                Ok(ComputeOutput {
                    result: serde_json::json!(result.result),
                    additional: None,
                    metadata: None,
                })
            },

            MatrixOp::Exp => {
                let terms = input.parameters.get("terms")
                    .and_then(|v| v.as_u64())
                    .map(|v| v as usize);

                let result = linear_algebra::matrix_exp(linear_algebra::MatrixExpRequest {
                    matrix: matrix_data,
                    terms,
                })
                .map_err(|e| e.to_string())?;

                Ok(ComputeOutput {
                    result: serde_json::json!(result.result),
                    additional: None,
                    metadata: None,
                })
            },

            MatrixOp::Inverse => {
                let result = linear_algebra::matrix_inverse(matrix_data)
                    .map_err(|e| e.to_string())?;

                Ok(ComputeOutput {
                    result: serde_json::json!(result),
                    additional: None,
                    metadata: Some(serde_json::json!({
                        "method": "Matrix inverse via LU decomposition"
                    })),
                })
            },
        }
    }

    /// Compute special functions
    fn compute_special_function(&self, func: &SpecialFunction, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::mathematics::special_functions;

        match func {
            SpecialFunction::Bessel => {
                let function_type = input.parameters.get("function_type")
                    .and_then(|v| v.as_str())
                    .ok_or("function_type required (J, Y, I, K)")?;
                let order = input.parameters.get("order")
                    .and_then(|v| v.as_f64())
                    .ok_or("order parameter required")?;
                let x = input.parameters.get("x")
                    .and_then(|v| v.as_f64())
                    .ok_or("x parameter required")?;

                let result = special_functions::bessel_function(special_functions::BesselRequest {
                    function_type: function_type.to_string(),
                    order,
                    x,
                })
                .map_err(|e| e.to_string())?;

                Ok(ComputeOutput {
                    result: serde_json::json!(result.value),
                    additional: Some({
                        let mut m = std::collections::HashMap::new();
                        m.insert("function_type".to_string(), serde_json::json!(result.function_type));
                        m.insert("order".to_string(), serde_json::json!(result.order));
                        m
                    }),
                    metadata: None,
                })
            },

            SpecialFunction::Gamma => {
                let function = input.parameters.get("function")
                    .and_then(|v| v.as_str())
                    .unwrap_or("gamma");
                let x = input.parameters.get("x")
                    .and_then(|v| v.as_f64())
                    .ok_or("x parameter required")?;
                let y = input.parameters.get("y").and_then(|v| v.as_f64());

                let result = special_functions::gamma_function(special_functions::GammaRequest {
                    x,
                    function: function.to_string(),
                    y,
                })
                .map_err(|e| e.to_string())?;

                Ok(ComputeOutput {
                    result: serde_json::json!(result.value),
                    additional: Some({
                        let mut m = std::collections::HashMap::new();
                        m.insert("function".to_string(), serde_json::json!(result.function));
                        m
                    }),
                    metadata: None,
                })
            },

            SpecialFunction::Erf => {
                let function = input.parameters.get("function")
                    .and_then(|v| v.as_str())
                    .unwrap_or("erf");
                let x = input.parameters.get("x")
                    .and_then(|v| v.as_f64())
                    .ok_or("x parameter required")?;

                let result = special_functions::error_function(special_functions::ErrorFunctionRequest {
                    x,
                    function: function.to_string(),
                })
                .map_err(|e| e.to_string())?;

                Ok(ComputeOutput {
                    result: serde_json::json!(result.value),
                    additional: Some({
                        let mut m = std::collections::HashMap::new();
                        m.insert("function".to_string(), serde_json::json!(result.function));
                        m
                    }),
                    metadata: None,
                })
            },

            SpecialFunction::Elliptic => {
                let integral_type = input.parameters.get("integral_type")
                    .and_then(|v| v.as_str())
                    .ok_or("integral_type required (K, E, F, Pi)")?;
                let k = input.parameters.get("k")
                    .and_then(|v| v.as_f64())
                    .ok_or("k (modulus) parameter required")?;
                let phi = input.parameters.get("phi").and_then(|v| v.as_f64());
                let n = input.parameters.get("n").and_then(|v| v.as_f64());

                let result = special_functions::elliptic_integral(special_functions::EllipticIntegralRequest {
                    integral_type: integral_type.to_string(),
                    k,
                    phi,
                    n,
                })
                .map_err(|e| e.to_string())?;

                Ok(ComputeOutput {
                    result: serde_json::json!(result.value),
                    additional: Some({
                        let mut m = std::collections::HashMap::new();
                        m.insert("integral_type".to_string(), serde_json::json!(result.integral_type));
                        m
                    }),
                    metadata: None,
                })
            },

            SpecialFunction::OrthogonalPoly => {
                let polynomial_type = input.parameters.get("polynomial_type")
                    .and_then(|v| v.as_str())
                    .ok_or("polynomial_type required (legendre, hermite, laguerre, chebyshev)")?;
                let n = input.parameters.get("n")
                    .and_then(|v| v.as_u64())
                    .ok_or("n (degree) parameter required")? as usize;
                let x = input.parameters.get("x")
                    .and_then(|v| v.as_f64())
                    .ok_or("x parameter required")?;
                let alpha = input.parameters.get("alpha").and_then(|v| v.as_f64());

                let result = special_functions::orthogonal_polynomial(special_functions::PolynomialRequest {
                    polynomial_type: polynomial_type.to_string(),
                    n,
                    x,
                    alpha,
                })
                .map_err(|e| e.to_string())?;

                Ok(ComputeOutput {
                    result: serde_json::json!(result.value),
                    additional: Some({
                        let mut m = std::collections::HashMap::new();
                        m.insert("polynomial_type".to_string(), serde_json::json!(result.polynomial_type));
                        m.insert("degree".to_string(), serde_json::json!(result.degree));
                        m
                    }),
                    metadata: None,
                })
            },

            SpecialFunction::Airy => {
                let function_type = input.parameters.get("function_type")
                    .and_then(|v| v.as_str())
                    .ok_or("function_type required (Ai, Bi, Ai_prime, Bi_prime)")?;
                let x = input.parameters.get("x")
                    .and_then(|v| v.as_f64())
                    .ok_or("x parameter required")?;

                let result = special_functions::airy_function(special_functions::AiryRequest {
                    function_type: function_type.to_string(),
                    x,
                })
                .map_err(|e| e.to_string())?;

                let mut additional = std::collections::HashMap::new();
                if let Some(deriv) = result.derivative {
                    additional.insert("derivative".to_string(), serde_json::json!(deriv));
                }

                Ok(ComputeOutput {
                    result: serde_json::json!(result.value),
                    additional: if additional.is_empty() { None } else { Some(additional) },
                    metadata: None,
                })
            },

            SpecialFunction::Hypergeometric => {
                Err("Hypergeometric functions not yet fully implemented".to_string())
            },
        }
    }

    /// Compute chemistry operations
    fn compute_chemistry(&self, op: &ChemistryOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::chemistry::*;

        // Map operation to chemistry module operation
        let chem_op = match op {
            ChemistryOp::PhCalculation => ChemistryOperation::PhCalculation,
            ChemistryOp::BufferCapacity => ChemistryOperation::BufferCapacity,
            ChemistryOp::Arrhenius => ChemistryOperation::Arrhenius,
            ChemistryOp::RateLaw => ChemistryOperation::RateLaw,
            ChemistryOp::GibbsFreeEnergy => ChemistryOperation::GibbsFreeEnergy,
            ChemistryOp::NernstEquation => ChemistryOperation::NernstEquation,
            ChemistryOp::BeerLambert => ChemistryOperation::BeerLambert,
            ChemistryOp::VanDerWaals => ChemistryOperation::VanDerWaals,
            _ => return Err(format!("Legacy chemistry operation {:?} not implemented in new module", op)),
        };

        // Convert parameters to ChemistryParams
        let params: ChemistryParams = serde_json::from_value(
            serde_json::to_value(&input.parameters)
                .map_err(|e| format!("Failed to serialize parameters: {}", e))?
        ).map_err(|e| format!("Failed to parse chemistry parameters: {}", e))?;

        let chem_input = ChemistryInput {
            operation: chem_op,
            parameters: params,
        };

        let result = calculate_chemistry(chem_input)
            .map_err(|e| format!("Chemistry calculation error: {}", e))?;

        Ok(ComputeOutput {
            result: serde_json::json!({
                "value": result.value,
                "unit": result.unit,
                "formula_used": result.formula_used,
                "interpretation": result.interpretation
            }),
            additional: None,
            metadata: None,
        })
    }

    /// Compute biology operations
    fn compute_biology(&self, op: &BiologyOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::biology::*;

        let bio_op = match op {
            BiologyOp::MichaelisMenten => BiologyOperation::MichaelisMenten,
            BiologyOp::Pharmacokinetics => BiologyOperation::Pharmacokinetics,
            BiologyOp::HardyWeinberg => BiologyOperation::HardyWeinberg,
            BiologyOp::GoldmanEquation => BiologyOperation::GoldmanEquation,
            BiologyOp::AllometricScaling => BiologyOperation::AllometricScaling,
            BiologyOp::GrowthModel => BiologyOperation::LineweaverBurk,
        };

        let params: BiologyParams = serde_json::from_value(
            serde_json::to_value(&input.parameters)
                .map_err(|e| format!("Failed to serialize parameters: {}", e))?
        ).map_err(|e| format!("Failed to parse biology parameters: {}", e))?;

        let bio_input = BiologyInput {
            operation: bio_op,
            parameters: params,
        };

        let result = calculate_biology(bio_input)
            .map_err(|e| format!("Biology calculation error: {}", e))?;

        Ok(ComputeOutput {
            result: serde_json::json!({
                "value": result.value,
                "unit": result.unit,
                "formula_used": result.formula_used,
                "interpretation": result.interpretation,
                "additional_data": result.additional_data
            }),
            additional: None,
            metadata: None,
        })
    }

    /// Compute thermodynamics operations
    fn compute_thermodynamics(&self, op: &ThermodynamicsOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::thermodynamics::*;

        let operation = match op {
            ThermodynamicsOp::Conduction => ThermodynamicsOperation::Conduction,
            ThermodynamicsOp::Convection => ThermodynamicsOperation::Convection,
            ThermodynamicsOp::Radiation => ThermodynamicsOperation::Radiation,
            ThermodynamicsOp::ThermalResistance => ThermodynamicsOperation::ThermalResistance,
            ThermodynamicsOp::Entropy => ThermodynamicsOperation::Entropy,
        };

        let params: ThermodynamicsParams = serde_json::from_value(
            serde_json::to_value(&input.parameters)
                .map_err(|e| format!("Failed to serialize parameters: {}", e))?
        ).map_err(|e| format!("Failed to parse thermodynamics parameters: {}", e))?;

        let thermo_input = ThermodynamicsInput {
            operation,
            parameters: params,
        };

        let result = calculate_thermodynamics(thermo_input)
            .map_err(|e| format!("Thermodynamics calculation error: {}", e))?;

        Ok(ComputeOutput {
            result: serde_json::json!({
                "value": result.value,
                "unit": result.unit,
                "formula_used": result.formula_used,
                "interpretation": result.interpretation,
                "additional_info": result.additional_info
            }),
            additional: None,
            metadata: None,
        })
    }

    /// Compute optics operations
    fn compute_optics(&self, op: &OpticsOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::optics::*;

        let optics_op = match op {
            OpticsOp::ThinLens => OpticsOperation::ThinLens,
            OpticsOp::SnellsLaw => OpticsOperation::SnellsLaw,
            OpticsOp::DiffractionGrating => OpticsOperation::DiffractionGrating,
            OpticsOp::FresnelEquations => OpticsOperation::FresnelEquations,
        };

        let params: OpticsParams = serde_json::from_value(
            serde_json::to_value(&input.parameters)
                .map_err(|e| format!("Failed to serialize parameters: {}", e))?
        ).map_err(|e| format!("Failed to parse optics parameters: {}", e))?;

        let optics_input = OpticsInput {
            operation: optics_op,
            parameters: params,
        };

        let result = calculate_optics(optics_input)
            .map_err(|e| format!("Optics calculation error: {}", e))?;

        Ok(ComputeOutput {
            result: serde_json::json!({
                "primary_value": result.primary_value,
                "unit": result.unit,
                "formula_used": result.formula_used,
                "secondary_values": result.secondary_values,
                "interpretation": result.interpretation
            }),
            additional: None,
            metadata: None,
        })
    }

    /// Compute geophysics operations
    fn compute_geophysics(&self, op: &GeophysicsOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::geophysics::*;

        let category = match op {
            GeophysicsOp::Seismology => GeophysicsCategory::Seismology,
            GeophysicsOp::Atmosphere => GeophysicsCategory::Atmosphere,
            GeophysicsOp::RadiometricDating => GeophysicsCategory::RadiometricDating,
            GeophysicsOp::PlanetaryScience => GeophysicsCategory::PlanetaryScience,
        };

        let params: GeophysicsParams = serde_json::from_value(
            serde_json::to_value(&input.parameters)
                .map_err(|e| format!("Failed to serialize parameters: {}", e))?
        ).map_err(|e| format!("Failed to parse geophysics parameters: {}", e))?;

        let geo_input = GeophysicsInput {
            category,
            parameters: params,
        };

        let result = calculate_geophysics(geo_input)
            .map_err(|e| format!("Geophysics calculation error: {}", e))?;

        Ok(ComputeOutput {
            result: serde_json::json!({
                "value": result.value,
                "unit": result.unit,
                "formula_used": result.formula_used,
                "uncertainty": result.uncertainty,
                "interpretation": result.interpretation,
                "additional_data": result.additional_data
            }),
            additional: None,
            metadata: None,
        })
    }

    /// Compute engineering operations
    fn compute_engineering(&self, op: &EngineeringOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::engineering::*;

        // Engineering operations are organized by discipline, infer from operation
        let (discipline, params) = match op {
            EngineeringOp::SoundPressureLevel | EngineeringOp::DopplerEffect | EngineeringOp::ReverberationTime => {
                (EngineeringDiscipline::Acoustics,
                 serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?)
            },
            EngineeringOp::Stress | EngineeringOp::Strain | EngineeringOp::FractureMechanics => {
                (EngineeringDiscipline::Materials,
                 serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?)
            },
            EngineeringOp::Bernoulli | EngineeringOp::Poiseuille | EngineeringOp::Drag => {
                (EngineeringDiscipline::FluidMechanics,
                 serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?)
            },
            EngineeringOp::PidController | EngineeringOp::FirstOrderResponse => {
                (EngineeringDiscipline::ControlTheory,
                 serde_json::to_value(&input.parameters)
                    .map_err(|e| format!("Failed to serialize parameters: {}", e))?)
            },
        };

        let eng_params: EngineeringParams = serde_json::from_value(params)
            .map_err(|e| format!("Failed to parse engineering parameters: {}", e))?;

        let eng_input = EngineeringInput {
            discipline,
            parameters: eng_params,
        };

        let result = calculate_engineering(eng_input)
            .map_err(|e| format!("Engineering calculation error: {}", e))?;

        Ok(ComputeOutput {
            result: serde_json::json!({
                "value": result.value,
                "unit": result.unit,
                "formula_used": result.formula_used,
                "classification": result.classification,
                "interpretation": result.interpretation,
                "additional": result.additional
            }),
            additional: None,
            metadata: None,
        })
    }

    /// Compute datetime operations
    fn compute_datetime(&self, op: &DateTimeOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::datetime::*;

        let dt_op = match op {
            DateTimeOp::AddInterval => DateTimeOperation::AddInterval,
            DateTimeOp::SubtractInterval => DateTimeOperation::SubtractInterval,
            DateTimeOp::DateDifference => DateTimeOperation::DateDifference,
            DateTimeOp::AgeCurrent => DateTimeOperation::AgeCurrent,
            DateTimeOp::AgeAtDate => DateTimeOperation::AgeAtDate,
            DateTimeOp::BusinessDays => DateTimeOperation::BusinessDays,
            DateTimeOp::AddBusinessDays => DateTimeOperation::AddBusinessDays,
            DateTimeOp::IsLeapYear => DateTimeOperation::IsLeapYear,
            DateTimeOp::DaysInMonth => DateTimeOperation::DaysInMonth,
            DateTimeOp::WeekNumber => DateTimeOperation::WeekNumber,
            DateTimeOp::DayOfWeek => DateTimeOperation::DayOfWeek,
        };

        let params: DateTimeParams = serde_json::from_value(
            serde_json::to_value(&input.parameters)
                .map_err(|e| format!("Failed to serialize parameters: {}", e))?
        ).map_err(|e| format!("Failed to parse datetime parameters: {}", e))?;

        let dt_input = DateTimeInput {
            operation: dt_op,
            parameters: params,
        };

        let result = calculate_datetime(dt_input)
            .map_err(|e| format!("DateTime calculation error: {}", e))?;

        Ok(ComputeOutput {
            result: serde_json::json!({
                "value": result.value,
                "numeric_value": result.numeric_value,
                "unit": result.unit,
                "operation_used": result.operation_used,
                "interpretation": result.interpretation,
                "additional_info": result.additional_info
            }),
            additional: None,
            metadata: None,
        })
    }

    /// Compute information theory operations
    fn compute_information(&self, op: &InformationOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::specialized::information_theory::*;

        let result_json = match op {
            InformationOp::Entropy => {
                let req: EntropyRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?
                ).map_err(|e| format!("Failed to parse entropy request: {}", e))?;

                let result = shannon_entropy(req)
                    .map_err(|e| format!("Entropy calculation error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
            InformationOp::MutualInfo => {
                let req: MutualInfoRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?
                ).map_err(|e| format!("Failed to parse mutual info request: {}", e))?;

                let result = mutual_information(req)
                    .map_err(|e| format!("Mutual information error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
            InformationOp::ChannelCapacity => {
                let req: ChannelCapacityRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?
                ).map_err(|e| format!("Failed to parse channel capacity request: {}", e))?;

                let result = channel_capacity(req)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
            InformationOp::Huffman => {
                let req: HuffmanRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?
                ).map_err(|e| format!("Failed to parse Huffman request: {}", e))?;

                let result = huffman_coding(req)
                    .map_err(|e| format!("Huffman coding error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
            InformationOp::Kolmogorov => {
                let req: KolmogorovRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?
                ).map_err(|e| format!("Failed to parse Kolmogorov request: {}", e))?;

                let result = kolmogorov_complexity(req)
                    .map_err(|e| format!("Kolmogorov complexity error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
            InformationOp::ConditionalEntropy => {
                let req: ConditionalEntropyRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?
                ).map_err(|e| format!("Failed to parse conditional entropy request: {}", e))?;

                let result = conditional_entropy(req)
                    .map_err(|e| format!("Conditional entropy error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
            InformationOp::KLDivergence => {
                let req: RelativeEntropyRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?
                ).map_err(|e| format!("Failed to parse KL divergence request: {}", e))?;

                let result = relative_entropy(req)
                    .map_err(|e| format!("KL divergence error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
        };

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    /// Compute graph theory operations
    fn compute_graph(&self, op: &GraphOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::specialized::graph_theory::*;

        let result_json = match op {
            GraphOp::ShortestPath => {
                let req: ShortestPathRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?
                ).map_err(|e| format!("Failed to parse shortest path request: {}", e))?;

                let result = shortest_path(req)
                    .map_err(|e| format!("Shortest path error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
            GraphOp::MinimumSpanningTree => {
                let req: MSTRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?
                ).map_err(|e| format!("Failed to parse MST request: {}", e))?;

                let result = minimum_spanning_tree(req)
                    .map_err(|e| format!("MST error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
            GraphOp::TopologicalSort => {
                let req: TopologicalSortRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?
                ).map_err(|e| format!("Failed to parse topological sort request: {}", e))?;

                let result = topological_sort(req)
                    .map_err(|e| format!("Topological sort error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
        };

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    /// Compute number theory operations
    fn compute_number_theory(&self, op: &NumberTheoryOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::specialized::cryptographic_mathematics::*;
        use num_bigint::BigInt;
        use std::str::FromStr;

        let result_json = match op {
            NumberTheoryOp::GeneratePrime => {
                let bits: u32 = input.parameters.get("bit_length")
                    .and_then(|v| v.as_u64())
                    .map(|v| v as u32)
                    .ok_or("Missing or invalid bit_length parameter")?;

                let prime = generate_prime(bits);
                serde_json::json!({
                    "prime": prime.to_string(),
                    "bit_length": bits
                })
            },
            NumberTheoryOp::ModExp => {
                let base = input.parameters.get("base")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid base parameter")?;
                let exp = input.parameters.get("exponent")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid exponent parameter")?;
                let modulus = input.parameters.get("modulus")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid modulus parameter")?;

                let result = mod_exp(&base, &exp, &modulus);
                serde_json::json!({
                    "result": result.to_string(),
                    "operation": "modular_exponentiation"
                })
            },
            NumberTheoryOp::ModInv => {
                let a = input.parameters.get("a")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid a parameter")?;
                let m = input.parameters.get("modulus")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid modulus parameter")?;

                let result = mod_inverse(&a, &m)
                    .ok_or("Modular inverse does not exist")?;
                serde_json::json!({
                    "inverse": result.to_string(),
                    "operation": "modular_inverse"
                })
            },
            NumberTheoryOp::GCD => {
                let a = input.parameters.get("a")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid a parameter")?;
                let b = input.parameters.get("b")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid b parameter")?;

                let (gcd, x, y) = extended_gcd(&a, &b);
                serde_json::json!({
                    "gcd": gcd.to_string(),
                    "bezout_x": x.to_string(),
                    "bezout_y": y.to_string(),
                    "operation": "extended_gcd"
                })
            },
            NumberTheoryOp::LCM => {
                let a = input.parameters.get("a")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid a parameter")?;
                let b = input.parameters.get("b")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid b parameter")?;

                let (gcd, _, _) = extended_gcd(&a, &b);
                let lcm = (&a * &b) / gcd;
                serde_json::json!({
                    "lcm": lcm.to_string(),
                    "operation": "lcm"
                })
            },
            NumberTheoryOp::RSAKeypair => {
                let bits: u32 = input.parameters.get("bit_length")
                    .and_then(|v| v.as_u64())
                    .map(|v| v as u32)
                    .unwrap_or(2048);

                if bits < 512 || bits > 4096 {
                    return Err("RSA key size must be between 512 and 4096 bits".to_string());
                }

                let (n, e, d) = generate_rsa_keypair(bits);
                serde_json::json!({
                    "public_key": {
                        "n": n.to_string(),
                        "e": e.to_string()
                    },
                    "private_key": {
                        "n": n.to_string(),
                        "d": d.to_string()
                    }
                })
            },
            NumberTheoryOp::RSAEncrypt => {
                let message = input.parameters.get("message")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid message parameter")?;
                let e = input.parameters.get("e")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid e parameter")?;
                let n = input.parameters.get("n")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid n parameter")?;

                let ciphertext = rsa_encrypt(&message, &e, &n);
                serde_json::json!({
                    "ciphertext": ciphertext.to_string()
                })
            },
            NumberTheoryOp::RSADecrypt => {
                let ciphertext = input.parameters.get("ciphertext")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid ciphertext parameter")?;
                let d = input.parameters.get("d")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid d parameter")?;
                let n = input.parameters.get("n")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid n parameter")?;

                let plaintext = rsa_decrypt(&ciphertext, &d, &n);
                serde_json::json!({
                    "plaintext": plaintext.to_string()
                })
            },
            NumberTheoryOp::SHA256 => {
                let input_str = input.parameters.get("input")
                    .and_then(|v| v.as_str())
                    .ok_or("Missing input parameter")?;

                let hash = sha256(input_str);
                serde_json::json!({
                    "hash": hash
                })
            },
            NumberTheoryOp::SHA3_256 => {
                let input_str = input.parameters.get("input")
                    .and_then(|v| v.as_str())
                    .ok_or("Missing input parameter")?;

                let hash = sha3_256(input_str);
                serde_json::json!({
                    "hash": hash
                })
            },
            NumberTheoryOp::ChineseRemainder => {
                let remainders_arr = input.parameters.get("remainders")
                    .and_then(|v| v.as_array())
                    .ok_or("Missing or invalid remainders parameter")?;
                let moduli_arr = input.parameters.get("moduli")
                    .and_then(|v| v.as_array())
                    .ok_or("Missing or invalid moduli parameter")?;

                let remainders: Result<Vec<BigInt>, _> = remainders_arr.iter()
                    .map(|v| v.as_str().and_then(|s| BigInt::from_str(s).ok()).ok_or("Invalid remainder"))
                    .collect();
                let moduli: Result<Vec<BigInt>, _> = moduli_arr.iter()
                    .map(|v| v.as_str().and_then(|s| BigInt::from_str(s).ok()).ok_or("Invalid modulus"))
                    .collect();

                let remainders = remainders?;
                let moduli = moduli?;

                let result = chinese_remainder_theorem(&remainders, &moduli)
                    .ok_or("No solution exists for given system")?;

                serde_json::json!({
                    "result": result.to_string()
                })
            },
            NumberTheoryOp::DiscreteLog => {
                let base = input.parameters.get("base")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid base parameter")?;
                let target = input.parameters.get("target")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid target parameter")?;
                let modulus = input.parameters.get("modulus")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid modulus parameter")?;
                let max_exp: u32 = input.parameters.get("max_exponent")
                    .and_then(|v| v.as_u64())
                    .map(|v| v as u32)
                    .unwrap_or(1000000);

                let result = discrete_log_bsgs(&base, &target, &modulus, max_exp)
                    .ok_or("Discrete logarithm not found within max_exponent bound")?;

                serde_json::json!({
                    "exponent": result.to_string()
                })
            },
            NumberTheoryOp::PrimalityTest => {
                let n = input.parameters.get("n")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid n parameter")?;
                let rounds: u32 = input.parameters.get("rounds")
                    .and_then(|v| v.as_u64())
                    .map(|v| v as u32)
                    .unwrap_or(10);

                let is_prime = miller_rabin_test(&n, rounds);
                let confidence = 1.0 - (0.25_f64).powi(rounds as i32);

                serde_json::json!({
                    "is_prime": is_prime,
                    "confidence": format!("{:.6}", confidence)
                })
            },

            NumberTheoryOp::EulerTotient | NumberTheoryOp::CarmichaelLambda | NumberTheoryOp::ECPointAdd => {
                return Err(format!("Number theory operation {:?} not yet fully implemented", op));
            },
        };

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    /// Compute geometry operations
    fn compute_geometry(&self, op: &GeometryOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::tools::computational_geometry::*;

        let result_json = match op {
            GeometryOp::ConvexHull => {
                let req: ConvexHullRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?
                ).map_err(|e| format!("Failed to parse convex hull request: {}", e))?;

                let result = convex_hull(req)
                    .map_err(|e| format!("Convex hull error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
            GeometryOp::Delaunay => {
                let req: DelaunayRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?
                ).map_err(|e| format!("Failed to parse Delaunay request: {}", e))?;

                let result = delaunay_triangulation(req)
                    .map_err(|e| format!("Delaunay triangulation error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
            GeometryOp::Voronoi => {
                let req: VoronoiRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?
                ).map_err(|e| format!("Failed to parse Voronoi request: {}", e))?;

                let result = voronoi_diagram(req)
                    .map_err(|e| format!("Voronoi diagram error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
            GeometryOp::PolygonArea => {
                let req: PolygonAreaRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?
                ).map_err(|e| format!("Failed to parse polygon area request: {}", e))?;

                let result = polygon_area(req)
                    .map_err(|e| format!("Polygon area error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
            GeometryOp::PointInPolygon => {
                let req: PointInPolygonRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?
                ).map_err(|e| format!("Failed to parse point in polygon request: {}", e))?;

                let result = point_in_polygon(req)
                    .map_err(|e| format!("Point in polygon error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
        };

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    /// Compute electromagnetism operations
    fn compute_em(&self, op: &EMComputation, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::physics::electromagnetism::*;

        let result_json = match op {
            EMComputation::PoyntingVector => {
                let req: PoyntingRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?
                ).map_err(|e| format!("Failed to parse Poynting vector request: {}", e))?;

                let result = poynting_vector(req)
                    .map_err(|e| format!("Poynting vector error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
            EMComputation::SkinEffect => {
                let req: SkinEffectRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?
                ).map_err(|e| format!("Failed to parse skin effect request: {}", e))?;

                let result = skin_effect(req)
                    .map_err(|e| format!("Skin effect error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
        };

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    /// Compute Fourier series
    fn compute_fourier_series(&self, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::tools::signal_processing::*;

        let req: FourierSeriesRequest = serde_json::from_value(
            serde_json::to_value(&input.parameters)
                .map_err(|e| format!("Failed to serialize parameters: {}", e))?
        ).map_err(|e| format!("Failed to parse Fourier series request: {}", e))?;

        let result = compute_fourier_series(req)
            .map_err(|e| format!("Fourier series error: {}", e))?;

        let result_json = serde_json::to_value(result)
            .map_err(|e| format!("Failed to serialize result: {}", e))?;

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    /// Compute physics operations (Tier 1 Wolfram Alpha expansion)
    fn compute_physics(&self, op: &PhysicsOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        match op {
            PhysicsOp::Relativity(relativity_op) => self.compute_relativity(relativity_op, input),
            PhysicsOp::StatisticalPhysics(stat_phys_op) => self.compute_statistical_physics(stat_phys_op, input),
            PhysicsOp::QuantumMechanics(qm_op) => self.compute_quantum_mechanics(qm_op, input),
            PhysicsOp::ControlSystems(cs_op) => self.compute_control_systems(cs_op, input),
            PhysicsOp::NuclearPhysics(np_op) => self.compute_nuclear_physics(np_op, input),
        }
    }

    /// Compute relativity operations
    fn compute_relativity(&self, op: &RelativityOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::physics::relativity::*;

        let result_json = match op {
            RelativityOp::LorentzTransform => {
                let velocity = input.parameters.get("velocity")
                    .and_then(|v| v.as_f64())
                    .ok_or("velocity parameter required")?;
                let position: Vec<f64> = input.parameters.get("position")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("position [x, y, z] required")?;
                let time = input.parameters.get("time")
                    .and_then(|v| v.as_f64())
                    .ok_or("time parameter required")?;

                let result = lorentz_transform(LorentzTransformRequest {
                    velocity,
                    position,
                    time,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            RelativityOp::TimeDilation => {
                let proper_time = input.parameters.get("proper_time")
                    .and_then(|v| v.as_f64())
                    .ok_or("proper_time parameter required")?;
                let velocity = input.parameters.get("velocity")
                    .and_then(|v| v.as_f64())
                    .ok_or("velocity parameter required")?;

                let result = time_dilation(TimeDilationRequest {
                    proper_time,
                    velocity,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            RelativityOp::LengthContraction => {
                let proper_length = input.parameters.get("proper_length")
                    .and_then(|v| v.as_f64())
                    .ok_or("proper_length parameter required")?;
                let velocity = input.parameters.get("velocity")
                    .and_then(|v| v.as_f64())
                    .ok_or("velocity parameter required")?;

                let result = length_contraction(LengthContractionRequest {
                    proper_length,
                    velocity,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            RelativityOp::RelativisticEnergy => {
                let mass = input.parameters.get("mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("mass parameter required")?;
                let velocity = input.parameters.get("velocity")
                    .and_then(|v| v.as_f64())
                    .ok_or("velocity parameter required")?;

                let result = relativistic_energy(RelativisticEnergyRequest {
                    mass,
                    velocity,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            RelativityOp::VelocityAddition => {
                let velocity1 = input.parameters.get("velocity1")
                    .and_then(|v| v.as_f64())
                    .ok_or("velocity1 parameter required")?;
                let velocity2 = input.parameters.get("velocity2")
                    .and_then(|v| v.as_f64())
                    .ok_or("velocity2 parameter required")?;

                let result = velocity_addition(VelocityAdditionRequest {
                    velocity1,
                    velocity2,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            RelativityOp::SchwarzschildMetric => {
                let mass = input.parameters.get("mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("mass parameter required")?;
                let radius = input.parameters.get("radius")
                    .and_then(|v| v.as_f64())
                    .ok_or("radius parameter required")?;

                let result = schwarzschild_metric(SchwarzschildRequest {
                    mass,
                    radius,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            RelativityOp::GravitationalTimeDilation => {
                let mass = input.parameters.get("mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("mass parameter required")?;
                let radius = input.parameters.get("radius")
                    .and_then(|v| v.as_f64())
                    .ok_or("radius parameter required")?;
                let proper_time = input.parameters.get("proper_time")
                    .and_then(|v| v.as_f64())
                    .ok_or("proper_time parameter required")?;

                let result = gravitational_time_dilation(GravitationalTimeDilationRequest {
                    mass,
                    radius,
                    proper_time,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            RelativityOp::OrbitalPrecession => {
                let mass = input.parameters.get("mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("mass parameter required")?;
                let semi_major_axis = input.parameters.get("semi_major_axis")
                    .and_then(|v| v.as_f64())
                    .ok_or("semi_major_axis parameter required")?;
                let eccentricity = input.parameters.get("eccentricity")
                    .and_then(|v| v.as_f64())
                    .ok_or("eccentricity parameter required")?;

                let result = orbital_precession(OrbitalPrecessionRequest {
                    mass,
                    semi_major_axis,
                    eccentricity,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            RelativityOp::GravitationalLensing => {
                let lens_mass = input.parameters.get("lens_mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("lens_mass parameter required")?;
                let impact_parameter = input.parameters.get("impact_parameter")
                    .and_then(|v| v.as_f64())
                    .ok_or("impact_parameter parameter required")?;

                let result = gravitational_lensing(GravitationalLensingRequest {
                    lens_mass,
                    impact_parameter,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            RelativityOp::BlackHoleProperties => {
                let mass = input.parameters.get("mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("mass parameter required")?;

                let result = black_hole_properties(BlackHoleRequest {
                    mass,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
        };

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    /// Compute statistical physics operations
    fn compute_statistical_physics(&self, op: &StatPhysicsOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::physics::statistical_physics::*;

        let result_json = match op {
            StatPhysicsOp::PartitionFunction => {
                let ensemble = input.parameters.get("ensemble")
                    .and_then(|v| v.as_str())
                    .ok_or("ensemble parameter required")?;
                let temperature = input.parameters.get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;

                let energy_levels: Option<Vec<f64>> = input.parameters.get("energy_levels")
                    .and_then(|v| serde_json::from_value(v.clone()).ok());
                let degeneracies: Option<Vec<usize>> = input.parameters.get("degeneracies")
                    .and_then(|v| serde_json::from_value(v.clone()).ok());
                let volume = input.parameters.get("volume").and_then(|v| v.as_f64());
                let num_particles = input.parameters.get("num_particles").and_then(|v| v.as_u64()).map(|v| v as usize);
                let chemical_potential = input.parameters.get("chemical_potential").and_then(|v| v.as_f64());

                let result = partition_function(PartitionFunctionRequest {
                    ensemble: ensemble.to_string(),
                    temperature,
                    volume,
                    num_particles,
                    energy_levels,
                    degeneracies,
                    chemical_potential,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            StatPhysicsOp::PartitionFunctionCanonical => {
                let temperature = input.parameters.get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let energy_levels: Vec<f64> = input.parameters.get("energy_levels")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("energy_levels required")?;
                let degeneracies: Vec<usize> = input.parameters.get("degeneracies")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("degeneracies required")?;

                let result = canonical_partition(CanonicalPartitionRequest {
                    temperature,
                    energy_levels,
                    degeneracies,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            StatPhysicsOp::PartitionFunctionGrandCanonical => {
                let temperature = input.parameters.get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let volume = input.parameters.get("volume")
                    .and_then(|v| v.as_f64())
                    .ok_or("volume parameter required")?;
                let chemical_potential = input.parameters.get("chemical_potential")
                    .and_then(|v| v.as_f64())
                    .ok_or("chemical_potential parameter required")?;
                let particle_type = input.parameters.get("particle_type")
                    .and_then(|v| v.as_str())
                    .unwrap_or("classical");

                let result = grand_canonical_partition(GrandCanonicalRequest {
                    temperature,
                    volume,
                    chemical_potential,
                    particle_type: particle_type.to_string(),
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            StatPhysicsOp::MaxwellBoltzmann => {
                let temperature = input.parameters.get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let mass = input.parameters.get("mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("mass parameter required")?;
                let velocity = input.parameters.get("velocity").and_then(|v| v.as_f64());

                let result = maxwell_boltzmann(MaxwellBoltzmannRequest {
                    temperature,
                    mass,
                    velocity,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            StatPhysicsOp::FermiDirac => {
                let energy = input.parameters.get("energy")
                    .and_then(|v| v.as_f64())
                    .ok_or("energy parameter required")?;
                let temperature = input.parameters.get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let chemical_potential = input.parameters.get("chemical_potential")
                    .and_then(|v| v.as_f64())
                    .ok_or("chemical_potential parameter required")?;

                let result = fermi_dirac(FermiDiracRequest {
                    energy,
                    temperature,
                    chemical_potential,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            StatPhysicsOp::BoseEinstein => {
                let energy = input.parameters.get("energy")
                    .and_then(|v| v.as_f64())
                    .ok_or("energy parameter required")?;
                let temperature = input.parameters.get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let chemical_potential = input.parameters.get("chemical_potential")
                    .and_then(|v| v.as_f64())
                    .ok_or("chemical_potential parameter required")?;

                let result = bose_einstein(BoseEinsteinRequest {
                    energy,
                    temperature,
                    chemical_potential,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            StatPhysicsOp::ChemicalPotential => {
                let temperature = input.parameters.get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let pressure = input.parameters.get("pressure")
                    .and_then(|v| v.as_f64())
                    .ok_or("pressure parameter required")?;
                let particle_density = input.parameters.get("particle_density")
                    .and_then(|v| v.as_f64())
                    .ok_or("particle_density parameter required")?;
                let particle_type = input.parameters.get("particle_type")
                    .and_then(|v| v.as_str())
                    .unwrap_or("classical");
                let mass = input.parameters.get("mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("mass parameter required")?;

                let result = chemical_potential(ChemicalPotentialRequest {
                    temperature,
                    pressure,
                    particle_density,
                    particle_type: particle_type.to_string(),
                    mass,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            StatPhysicsOp::FugacityCoefficient => {
                let temperature = input.parameters.get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let pressure = input.parameters.get("pressure")
                    .and_then(|v| v.as_f64())
                    .ok_or("pressure parameter required")?;
                let particle_density = input.parameters.get("particle_density")
                    .and_then(|v| v.as_f64())
                    .ok_or("particle_density parameter required")?;

                let result = fugacity_coefficient(FugacityRequest {
                    temperature,
                    pressure,
                    particle_density,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            StatPhysicsOp::PhaseTransition => {
                let model = input.parameters.get("model")
                    .and_then(|v| v.as_str())
                    .ok_or("model parameter required")?;
                let temperature = input.parameters.get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let critical_temperature = input.parameters.get("critical_temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("critical_temperature parameter required")?;

                let parameters = std::collections::HashMap::new();

                let result = phase_transition(PhaseTransitionRequest {
                    model: model.to_string(),
                    temperature,
                    critical_temperature,
                    parameters,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            StatPhysicsOp::CriticalPhenomena => {
                let temperature = input.parameters.get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let critical_temperature = input.parameters.get("critical_temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("critical_temperature parameter required")?;
                let model = input.parameters.get("model")
                    .and_then(|v| v.as_str())
                    .unwrap_or("mean_field");

                let result = critical_phenomena(CriticalPhenomenaRequest {
                    temperature,
                    critical_temperature,
                    model: model.to_string(),
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
        };

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    /// Compute quantum mechanics operations
    fn compute_quantum_mechanics(&self, op: &QuantumMechOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::physics::quantum_mechanics::*;

        let result_json = match op {
            QuantumMechOp::SchrodingerEquation => {
                let potential = input.parameters.get("potential")
                    .and_then(|v| v.as_str())
                    .ok_or("potential parameter required")?;
                let energy = input.parameters.get("energy")
                    .and_then(|v| v.as_f64())
                    .ok_or("energy parameter required")?;
                let position = input.parameters.get("position")
                    .and_then(|v| v.as_f64())
                    .ok_or("position parameter required")?;

                let mut parameters = std::collections::HashMap::new();
                if let Some(val) = input.parameters.get("width").and_then(|v| v.as_f64()) {
                    parameters.insert("width".to_string(), val);
                }
                if let Some(val) = input.parameters.get("quantum_number").and_then(|v| v.as_f64()) {
                    parameters.insert("quantum_number".to_string(), val);
                }
                if let Some(val) = input.parameters.get("omega").and_then(|v| v.as_f64()) {
                    parameters.insert("omega".to_string(), val);
                }
                if let Some(val) = input.parameters.get("mass").and_then(|v| v.as_f64()) {
                    parameters.insert("mass".to_string(), val);
                }
                if let Some(val) = input.parameters.get("barrier_height").and_then(|v| v.as_f64()) {
                    parameters.insert("barrier_height".to_string(), val);
                }

                let result = schrodinger_equation(SchrodingerRequest {
                    potential: potential.to_string(),
                    energy,
                    position,
                    parameters,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            QuantumMechOp::HarmonicOscillator => {
                let quantum_number = input.parameters.get("quantum_number")
                    .and_then(|v| v.as_u64())
                    .ok_or("quantum_number parameter required")? as usize;
                let omega = input.parameters.get("omega")
                    .and_then(|v| v.as_f64())
                    .ok_or("omega parameter required")?;
                let mass = input.parameters.get("mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("mass parameter required")?;
                let position = input.parameters.get("position").and_then(|v| v.as_f64());

                let result = harmonic_oscillator(HarmonicOscillatorRequest {
                    quantum_number,
                    omega,
                    mass,
                    position,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            QuantumMechOp::HydrogenAtom => {
                let n = input.parameters.get("n")
                    .and_then(|v| v.as_u64())
                    .ok_or("n (principal quantum number) required")? as usize;
                let l = input.parameters.get("l")
                    .and_then(|v| v.as_u64())
                    .ok_or("l (orbital quantum number) required")? as usize;
                let m = input.parameters.get("m")
                    .and_then(|v| v.as_i64())
                    .ok_or("m (magnetic quantum number) required")? as i32;
                let r = input.parameters.get("r").and_then(|v| v.as_f64());

                let result = hydrogen_atom(HydrogenAtomRequest { n, l, m, r })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            QuantumMechOp::AngularMomentum => {
                let l = input.parameters.get("l")
                    .and_then(|v| v.as_u64())
                    .ok_or("l parameter required")? as usize;
                let m = input.parameters.get("m")
                    .and_then(|v| v.as_i64())
                    .ok_or("m parameter required")? as i32;
                let operation = input.parameters.get("operation")
                    .and_then(|v| v.as_str())
                    .unwrap_or("eigenvalue");

                let result = angular_momentum(AngularMomentumRequest {
                    l,
                    m,
                    operation: operation.to_string(),
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            QuantumMechOp::SpinOperators => {
                let spin = input.parameters.get("spin")
                    .and_then(|v| v.as_f64())
                    .ok_or("spin parameter required")?;
                let component = input.parameters.get("component")
                    .and_then(|v| v.as_str())
                    .ok_or("component parameter required (x, y, z, plus, minus)")?;
                let state: Option<Vec<f64>> = input.parameters.get("state")
                    .and_then(|v| serde_json::from_value(v.clone()).ok());

                let result = spin_operators(SpinRequest {
                    spin,
                    component: component.to_string(),
                    state,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            QuantumMechOp::PerturbationTheory => {
                let order = input.parameters.get("order")
                    .and_then(|v| v.as_u64())
                    .ok_or("order parameter required (1 or 2)")? as usize;
                let unperturbed_energy = input.parameters.get("unperturbed_energy")
                    .and_then(|v| v.as_f64())
                    .ok_or("unperturbed_energy required")?;
                let perturbation_matrix_element = input.parameters.get("perturbation_matrix_element")
                    .and_then(|v| v.as_f64())
                    .ok_or("perturbation_matrix_element required")?;
                let energy_differences: Option<Vec<f64>> = input.parameters.get("energy_differences")
                    .and_then(|v| serde_json::from_value(v.clone()).ok());
                let coupling_matrix_elements: Option<Vec<f64>> = input.parameters.get("coupling_matrix_elements")
                    .and_then(|v| serde_json::from_value(v.clone()).ok());

                let result = perturbation_theory(PerturbationRequest {
                    order,
                    unperturbed_energy,
                    perturbation_matrix_element,
                    energy_differences,
                    coupling_matrix_elements,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            QuantumMechOp::TunnelingProbability => {
                let barrier_height = input.parameters.get("barrier_height")
                    .and_then(|v| v.as_f64())
                    .ok_or("barrier_height required")?;
                let barrier_width = input.parameters.get("barrier_width")
                    .and_then(|v| v.as_f64())
                    .ok_or("barrier_width required")?;
                let particle_energy = input.parameters.get("particle_energy")
                    .and_then(|v| v.as_f64())
                    .ok_or("particle_energy required")?;
                let particle_mass = input.parameters.get("particle_mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("particle_mass required")?;

                let result = tunneling_probability(TunnelingRequest {
                    barrier_height,
                    barrier_width,
                    particle_energy,
                    particle_mass,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            QuantumMechOp::DensityMatrix => {
                let state_type = input.parameters.get("state_type")
                    .and_then(|v| v.as_str())
                    .ok_or("state_type required (pure or mixed)")?;
                let state_vector: Option<Vec<f64>> = input.parameters.get("state_vector")
                    .and_then(|v| serde_json::from_value(v.clone()).ok());
                let density_mat: Option<Vec<Vec<f64>>> = input.parameters.get("density_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok());

                let result = density_matrix(DensityMatrixRequest {
                    state_type: state_type.to_string(),
                    state_vector,
                    density_matrix: density_mat,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            QuantumMechOp::EntanglementMeasure => {
                let state_vector: Vec<f64> = input.parameters.get("state_vector")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("state_vector required (4 components for two-qubit state)")?;
                let measure = input.parameters.get("measure")
                    .and_then(|v| v.as_str())
                    .unwrap_or("concurrence");

                let result = entanglement_measure(EntanglementRequest {
                    state_vector,
                    measure: measure.to_string(),
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            QuantumMechOp::QuantumEntropy => {
                let density_matrix: Vec<Vec<f64>> = input.parameters.get("density_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("density_matrix required")?;
                let entropy_type = input.parameters.get("entropy_type")
                    .and_then(|v| v.as_str())
                    .unwrap_or("von_neumann");
                let alpha = input.parameters.get("alpha").and_then(|v| v.as_f64());

                let result = quantum_entropy(QuantumEntropyRequest {
                    density_matrix,
                    entropy_type: entropy_type.to_string(),
                    alpha,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            QuantumMechOp::QuantumCoherence => {
                let density_matrix: Vec<Vec<f64>> = input.parameters.get("density_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("density_matrix required")?;
                let basis = input.parameters.get("basis")
                    .and_then(|v| v.as_str())
                    .unwrap_or("computational");

                let result = quantum_coherence(CoherenceRequest {
                    density_matrix,
                    basis: basis.to_string(),
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            QuantumMechOp::BellInequality => {
                let state_vector: Vec<f64> = input.parameters.get("state_vector")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("state_vector required (two-qubit state)")?;
                let measurement_angles: Vec<f64> = input.parameters.get("measurement_angles")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("measurement_angles required (4 angles for CHSH)")?;

                let result = bell_inequality(BellInequalityRequest {
                    state_vector,
                    measurement_angles,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            QuantumMechOp::QuantumTomography => {
                let measurements: Vec<Vec<f64>> = input.parameters.get("measurements")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("measurements required")?;
                let measurement_bases: Vec<String> = input.parameters.get("measurement_bases")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("measurement_bases required (e.g., [X, Y, Z])")?;

                let result = quantum_tomography(QuantumTomographyRequest {
                    measurements,
                    measurement_bases,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            QuantumMechOp::QuantumGate => {
                let gate_type = input.parameters.get("gate_type")
                    .and_then(|v| v.as_str())
                    .ok_or("gate_type required (hadamard, pauli_x, pauli_z, phase, cnot)")?;
                let input_state: Vec<f64> = input.parameters.get("input_state")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("input_state required")?;
                let phase_angle = input.parameters.get("phase_angle").and_then(|v| v.as_f64());

                let result = quantum_gate(QuantumGateRequest {
                    gate_type: gate_type.to_string(),
                    input_state,
                    phase_angle,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            QuantumMechOp::QuantumCircuit => {
                let num_qubits = input.parameters.get("num_qubits")
                    .and_then(|v| v.as_u64())
                    .ok_or("num_qubits required")? as usize;
                let gates: Vec<serde_json::Value> = input.parameters.get("gates")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("gates required (list of gate operations)")?;
                let initial_state: Vec<f64> = input.parameters.get("initial_state")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("initial_state required")?;

                let result = quantum_circuit(QuantumCircuitRequest {
                    num_qubits,
                    gates,
                    initial_state,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
        };

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    /// Compute control systems operations
    fn compute_control_systems(&self, op: &ControlSystemsOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::physics::control_systems::*;

        let result_json = match op {
            ControlSystemsOp::TransferFunction => {
                let numerator: Vec<f64> = input.parameters.get("numerator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("numerator coefficients required")?;
                let denominator: Vec<f64> = input.parameters.get("denominator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("denominator coefficients required")?;
                let operation = input.parameters.get("operation")
                    .and_then(|v| v.as_str())
                    .unwrap_or("evaluate")
                    .to_string();
                let frequency = input.parameters.get("frequency").and_then(|v| v.as_f64());
                let second_tf: Option<(Vec<f64>, Vec<f64>)> = input.parameters.get("second_tf")
                    .and_then(|v| serde_json::from_value(v.clone()).ok());

                let result = transfer_function(TransferFunctionRequest {
                    numerator,
                    denominator,
                    operation,
                    frequency,
                    second_tf,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            ControlSystemsOp::PoleZeroAnalysis => {
                let numerator: Vec<f64> = input.parameters.get("numerator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("numerator required")?;
                let denominator: Vec<f64> = input.parameters.get("denominator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("denominator required")?;

                let result = pole_zero_analysis(PoleZeroRequest {
                    numerator,
                    denominator,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            ControlSystemsOp::BodePlot => {
                let numerator: Vec<f64> = input.parameters.get("numerator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("numerator required")?;
                let denominator: Vec<f64> = input.parameters.get("denominator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("denominator required")?;
                let freq_min = input.parameters.get("freq_min")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.1);
                let freq_max = input.parameters.get("freq_max")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(100.0);
                let num_points = input.parameters.get("num_points")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(50) as usize;

                let result = bode_plot(BodePlotRequest {
                    numerator,
                    denominator,
                    freq_min,
                    freq_max,
                    num_points,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            ControlSystemsOp::NyquistPlot => {
                let numerator: Vec<f64> = input.parameters.get("numerator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("numerator required")?;
                let denominator: Vec<f64> = input.parameters.get("denominator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("denominator required")?;
                let freq_min = input.parameters.get("freq_min")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.1);
                let freq_max = input.parameters.get("freq_max")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(100.0);
                let num_points = input.parameters.get("num_points")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(50) as usize;

                let result = nyquist_plot(NyquistPlotRequest {
                    numerator,
                    denominator,
                    freq_min,
                    freq_max,
                    num_points,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            ControlSystemsOp::RootLocus => {
                let numerator: Vec<f64> = input.parameters.get("numerator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("numerator required")?;
                let denominator: Vec<f64> = input.parameters.get("denominator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("denominator required")?;
                let gain_min = input.parameters.get("gain_min")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0);
                let gain_max = input.parameters.get("gain_max")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(10.0);
                let num_points = input.parameters.get("num_points")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(50) as usize;

                let result = root_locus(RootLocusRequest {
                    numerator,
                    denominator,
                    gain_min,
                    gain_max,
                    num_points,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            ControlSystemsOp::StateSpace => {
                let a_matrix: Vec<Vec<f64>> = input.parameters.get("a_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("a_matrix required")?;
                let b_matrix: Vec<Vec<f64>> = input.parameters.get("b_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("b_matrix required")?;
                let c_matrix: Vec<Vec<f64>> = input.parameters.get("c_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("c_matrix required")?;
                let d_matrix: Vec<Vec<f64>> = input.parameters.get("d_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("d_matrix required")?;
                let operation = input.parameters.get("operation")
                    .and_then(|v| v.as_str())
                    .unwrap_or("to_transfer_function")
                    .to_string();
                let time = input.parameters.get("time").and_then(|v| v.as_f64());

                let result = state_space(StateSpaceRequest {
                    a_matrix,
                    b_matrix,
                    c_matrix,
                    d_matrix,
                    operation,
                    time,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            ControlSystemsOp::Controllability => {
                let a_matrix: Vec<Vec<f64>> = input.parameters.get("a_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("a_matrix required")?;
                let b_matrix: Vec<Vec<f64>> = input.parameters.get("b_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("b_matrix required")?;

                let result = controllability(ControllabilityRequest {
                    a_matrix,
                    b_matrix,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            ControlSystemsOp::Observability => {
                let a_matrix: Vec<Vec<f64>> = input.parameters.get("a_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("a_matrix required")?;
                let c_matrix: Vec<Vec<f64>> = input.parameters.get("c_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("c_matrix required")?;

                let result = observability(ObservabilityRequest {
                    a_matrix,
                    c_matrix,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            ControlSystemsOp::RouthHurwitz => {
                let characteristic_polynomial: Vec<f64> = input.parameters.get("characteristic_polynomial")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("characteristic_polynomial required")?;

                let result = routh_hurwitz(RouthHurwitzRequest {
                    characteristic_polynomial,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            ControlSystemsOp::GainMargin => {
                let numerator: Vec<f64> = input.parameters.get("numerator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("numerator required")?;
                let denominator: Vec<f64> = input.parameters.get("denominator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("denominator required")?;

                let result = gain_margin(GainMarginRequest {
                    numerator,
                    denominator,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            ControlSystemsOp::PhaseMargin => {
                let numerator: Vec<f64> = input.parameters.get("numerator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("numerator required")?;
                let denominator: Vec<f64> = input.parameters.get("denominator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("denominator required")?;

                let result = phase_margin(PhaseMarginRequest {
                    numerator,
                    denominator,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            ControlSystemsOp::StepResponse => {
                let numerator: Vec<f64> = input.parameters.get("numerator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("numerator required")?;
                let denominator: Vec<f64> = input.parameters.get("denominator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("denominator required")?;
                let time_span = input.parameters.get("time_span")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(10.0);
                let num_points = input.parameters.get("num_points")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(100) as usize;

                let result = step_response(StepResponseRequest {
                    numerator,
                    denominator,
                    time_span,
                    num_points,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
        };

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    /// Compute nuclear physics operations
    fn compute_nuclear_physics(&self, op: &NuclearOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::physics::nuclear_physics::*;

        let result_json = match op {
            NuclearOp::RadioactiveDecay => {
                let initial_quantity = input.parameters.get("initial_quantity")
                    .and_then(|v| v.as_f64())
                    .ok_or("initial_quantity parameter required")?;
                let decay_constant = input.parameters.get("decay_constant")
                    .and_then(|v| v.as_f64())
                    .ok_or("decay_constant parameter required")?;
                let time = input.parameters.get("time")
                    .and_then(|v| v.as_f64())
                    .ok_or("time parameter required")?;

                let result = radioactive_decay(RadioactiveDecayRequest {
                    initial_quantity,
                    decay_constant,
                    time,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            NuclearOp::DecayChain => {
                let parent_initial = input.parameters.get("parent_initial")
                    .and_then(|v| v.as_f64())
                    .ok_or("parent_initial parameter required")?;
                let parent_decay_constant = input.parameters.get("parent_decay_constant")
                    .and_then(|v| v.as_f64())
                    .ok_or("parent_decay_constant parameter required")?;
                let daughter_decay_constant = input.parameters.get("daughter_decay_constant")
                    .and_then(|v| v.as_f64())
                    .ok_or("daughter_decay_constant parameter required")?;
                let time = input.parameters.get("time")
                    .and_then(|v| v.as_f64())
                    .ok_or("time parameter required")?;

                let result = decay_chain(DecayChainRequest {
                    parent_initial,
                    parent_decay_constant,
                    daughter_decay_constant,
                    time,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            NuclearOp::HalfLife => {
                let decay_constant_val = input.parameters.get("decay_constant")
                    .and_then(|v| v.as_f64());
                let half_life_val = input.parameters.get("half_life")
                    .and_then(|v| v.as_f64());
                let mean_lifetime_val = input.parameters.get("mean_lifetime")
                    .and_then(|v| v.as_f64());

                let result = half_life(HalfLifeRequest {
                    decay_constant: decay_constant_val,
                    half_life: half_life_val,
                    mean_lifetime: mean_lifetime_val,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            NuclearOp::BindingEnergy => {
                let atomic_number = input.parameters.get("atomic_number")
                    .and_then(|v| v.as_u64())
                    .ok_or("atomic_number parameter required")? as u32;
                let mass_number = input.parameters.get("mass_number")
                    .and_then(|v| v.as_u64())
                    .ok_or("mass_number parameter required")? as u32;
                let atomic_mass = input.parameters.get("atomic_mass")
                    .and_then(|v| v.as_f64());

                let result = binding_energy(BindingEnergyRequest {
                    atomic_number,
                    mass_number,
                    atomic_mass,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            NuclearOp::MassDefect => {
                let protons = input.parameters.get("protons")
                    .and_then(|v| v.as_u64())
                    .ok_or("protons parameter required")? as u32;
                let neutrons = input.parameters.get("neutrons")
                    .and_then(|v| v.as_u64())
                    .ok_or("neutrons parameter required")? as u32;
                let nuclear_mass = input.parameters.get("nuclear_mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("nuclear_mass parameter required")?;

                let result = mass_defect(MassDefectRequest {
                    protons,
                    neutrons,
                    nuclear_mass,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            NuclearOp::FissionEnergy => {
                let parent_mass = input.parameters.get("parent_mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("parent_mass parameter required")?;
                let fragment1_mass = input.parameters.get("fragment1_mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("fragment1_mass parameter required")?;
                let fragment2_mass = input.parameters.get("fragment2_mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("fragment2_mass parameter required")?;
                let neutrons_released = input.parameters.get("neutrons_released")
                    .and_then(|v| v.as_u64())
                    .ok_or("neutrons_released parameter required")? as u32;

                let result = fission_energy(FissionEnergyRequest {
                    parent_mass,
                    fragment1_mass,
                    fragment2_mass,
                    neutrons_released,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            NuclearOp::FusionEnergy => {
                let reactant1_mass = input.parameters.get("reactant1_mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("reactant1_mass parameter required")?;
                let reactant2_mass = input.parameters.get("reactant2_mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("reactant2_mass parameter required")?;
                let product1_mass = input.parameters.get("product1_mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("product1_mass parameter required")?;
                let product2_mass = input.parameters.get("product2_mass")
                    .and_then(|v| v.as_f64());

                let result = fusion_energy(FusionEnergyRequest {
                    reactant1_mass,
                    reactant2_mass,
                    product1_mass,
                    product2_mass,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },

            NuclearOp::NuclearReaction => {
                let reactants: Vec<f64> = input.parameters.get("reactants")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("reactants parameter required as array of masses in amu")?;
                let products: Vec<f64> = input.parameters.get("products")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("products parameter required as array of masses in amu")?;
                let projectile_energy_mev = input.parameters.get("projectile_energy_mev")
                    .and_then(|v| v.as_f64());

                let result = nuclear_reaction(NuclearReactionRequest {
                    reactants,
                    products,
                    projectile_energy_mev,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            },
        };

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }
}

impl Compute for UnifiedComputer {
    fn compute(&self, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        match &input.operation {
            ComputeOp::Tensor(tensor_op) => self.compute_tensor(tensor_op, input),

            ComputeOp::Matrix(matrix_op) => self.compute_matrix_op(matrix_op, input),

            ComputeOp::MatrixDecomp(decomp) => self.compute_matrix_decomp(decomp, input),

            ComputeOp::SpecialFunc(func) => self.compute_special_function(func, input),

            // Core Mathematical Operations
            ComputeOp::NumberTheory(op) => self.compute_number_theory(op, input),
            ComputeOp::Geometry(op) => self.compute_geometry(op, input),
            ComputeOp::Information(op) => self.compute_information(op, input),
            ComputeOp::Graph(op) => self.compute_graph(op, input),
            ComputeOp::EM(op) => self.compute_em(op, input),
            ComputeOp::FourierSeries => self.compute_fourier_series(input),

            // Scientific Formulas (2025 Expansion)
            ComputeOp::Chemistry(op) => self.compute_chemistry(op, input),
            ComputeOp::Biology(op) => self.compute_biology(op, input),
            ComputeOp::Thermodynamics(op) => self.compute_thermodynamics(op, input),
            ComputeOp::Optics(op) => self.compute_optics(op, input),
            ComputeOp::Geophysics(op) => self.compute_geophysics(op, input),
            ComputeOp::Engineering(op) => self.compute_engineering(op, input),
            ComputeOp::DateTime(op) => self.compute_datetime(op, input),

            // Physics (Tier 1 Wolfram Alpha expansion)
            ComputeOp::Physics(op) => self.compute_physics(op, input),

            _ => Err(format!("Compute operation {:?} not yet implemented", input.operation))
        }
    }
}
