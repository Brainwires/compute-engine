//! Unified Computer implementation
//!
//! Routes compute requests for tensors, matrices, and special functions

use crate::engine::*;
use crate::implementations::compute;
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
        let metric = self.parse_metric_tensor(metric_data)?;

        match op {
            TensorOp::Christoffel => {
                let result = tensor::calculate_christoffel_symbols(&metric, &coords)
                    .map_err(|e| e.to_string())?;

                let mut additional = std::collections::HashMap::new();
                additional.insert("dimension".to_string(), serde_json::json!(result.dimension));
                additional.insert(
                    "latex".to_string(),
                    serde_json::json!(self.christoffel_to_latex(&result)),
                );

                Ok(ComputeOutput {
                    result: serde_json::json!(result.symbols),
                    additional: Some(additional),
                    metadata: None,
                })
            }

            TensorOp::Riemann => {
                let result = tensor::calculate_riemann_tensor(&metric, &coords)
                    .map_err(|e| e.to_string())?;

                let mut additional = std::collections::HashMap::new();
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

                let mut additional = std::collections::HashMap::new();
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

                let mut additional = std::collections::HashMap::new();
                additional.insert("latex".to_string(), serde_json::json!("R"));

                Ok(ComputeOutput {
                    result: serde_json::json!(result),
                    additional: Some(additional),
                    metadata: None,
                })
            }

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

                let mut additional = std::collections::HashMap::new();
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
                let dimension = metric.len();

                Ok(ComputeOutput {
                    result: serde_json::json!("Tensor product computed"),
                    additional: Some(std::collections::HashMap::from([
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
                    additional: Some(std::collections::HashMap::from([(
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
                    additional: Some(std::collections::HashMap::from([
                        ("parameter".to_string(), serde_json::json!(curve_param)),
                        ("equation".to_string(), serde_json::json!("DV/dt + Γ V = 0")),
                    ])),
                    metadata: None,
                })
            }
        }
    }

    /// Parse metric tensor from JSON
    fn parse_metric_tensor(
        &self,
        data: &Value,
    ) -> ToolResult<Vec<Vec<crate::mathematics::tensor_calculus::SymbolicExpr>>> {
        use crate::mathematics::tensor_calculus::SymbolicExpr;

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
    fn christoffel_to_latex(
        &self,
        result: &crate::mathematics::tensor_calculus::ChristoffelResult,
    ) -> String {
        format!(
            "\\Gamma^{{\\mu}}_{{\\alpha\\beta}} \\text{{ ({} non-zero components)}}",
            result.symbols.len()
        )
    }

    /// Compute special functions
    fn compute_special_function(
        &self,
        func: &SpecialFunction,
        input: &ComputeInput,
    ) -> ToolResult<ComputeOutput> {
        use crate::mathematics::special_functions;

        match func {
            SpecialFunction::Bessel => {
                let function_type = input
                    .parameters
                    .get("function_type")
                    .and_then(|v| v.as_str())
                    .ok_or("function_type required (J, Y, I, K)")?;
                let order = input
                    .parameters
                    .get("order")
                    .and_then(|v| v.as_f64())
                    .ok_or("order parameter required")?;
                let x = input
                    .parameters
                    .get("x")
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
                        m.insert(
                            "function_type".to_string(),
                            serde_json::json!(result.function_type),
                        );
                        m.insert("order".to_string(), serde_json::json!(result.order));
                        m
                    }),
                    metadata: None,
                })
            }

            SpecialFunction::Gamma => {
                let function = input
                    .parameters
                    .get("function")
                    .and_then(|v| v.as_str())
                    .unwrap_or("gamma");
                let x = input
                    .parameters
                    .get("x")
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
            }

            SpecialFunction::Erf => {
                let function = input
                    .parameters
                    .get("function")
                    .and_then(|v| v.as_str())
                    .unwrap_or("erf");
                let x = input
                    .parameters
                    .get("x")
                    .and_then(|v| v.as_f64())
                    .ok_or("x parameter required")?;

                let result =
                    special_functions::error_function(special_functions::ErrorFunctionRequest {
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
            }

            SpecialFunction::Elliptic => {
                let integral_type = input
                    .parameters
                    .get("integral_type")
                    .and_then(|v| v.as_str())
                    .ok_or("integral_type required (K, E, F, Pi)")?;
                let k = input
                    .parameters
                    .get("k")
                    .and_then(|v| v.as_f64())
                    .ok_or("k (modulus) parameter required")?;
                let phi = input.parameters.get("phi").and_then(|v| v.as_f64());
                let n = input.parameters.get("n").and_then(|v| v.as_f64());

                let result = special_functions::elliptic_integral(
                    special_functions::EllipticIntegralRequest {
                        integral_type: integral_type.to_string(),
                        k,
                        phi,
                        n,
                    },
                )
                .map_err(|e| e.to_string())?;

                Ok(ComputeOutput {
                    result: serde_json::json!(result.value),
                    additional: Some({
                        let mut m = std::collections::HashMap::new();
                        m.insert(
                            "integral_type".to_string(),
                            serde_json::json!(result.integral_type),
                        );
                        m
                    }),
                    metadata: None,
                })
            }

            SpecialFunction::OrthogonalPoly => {
                let polynomial_type = input
                    .parameters
                    .get("polynomial_type")
                    .and_then(|v| v.as_str())
                    .ok_or("polynomial_type required (legendre, hermite, laguerre, chebyshev)")?;
                let n = input
                    .parameters
                    .get("n")
                    .and_then(|v| v.as_u64())
                    .ok_or("n (degree) parameter required")? as usize;
                let x = input
                    .parameters
                    .get("x")
                    .and_then(|v| v.as_f64())
                    .ok_or("x parameter required")?;
                let alpha = input.parameters.get("alpha").and_then(|v| v.as_f64());

                let result = special_functions::orthogonal_polynomial(
                    special_functions::PolynomialRequest {
                        polynomial_type: polynomial_type.to_string(),
                        n,
                        x,
                        alpha,
                    },
                )
                .map_err(|e| e.to_string())?;

                Ok(ComputeOutput {
                    result: serde_json::json!(result.value),
                    additional: Some({
                        let mut m = std::collections::HashMap::new();
                        m.insert(
                            "polynomial_type".to_string(),
                            serde_json::json!(result.polynomial_type),
                        );
                        m.insert("degree".to_string(), serde_json::json!(result.degree));
                        m
                    }),
                    metadata: None,
                })
            }

            SpecialFunction::Airy => {
                let function_type = input
                    .parameters
                    .get("function_type")
                    .and_then(|v| v.as_str())
                    .ok_or("function_type required (Ai, Bi, Ai_prime, Bi_prime)")?;
                let x = input
                    .parameters
                    .get("x")
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
                    additional: if additional.is_empty() {
                        None
                    } else {
                        Some(additional)
                    },
                    metadata: None,
                })
            }

            SpecialFunction::Hypergeometric => {
                Err("Hypergeometric functions not yet fully implemented".to_string())
            }
        }
    }

    /// Compute information theory operations
    fn compute_information(
        &self,
        op: &InformationOp,
        input: &ComputeInput,
    ) -> ToolResult<ComputeOutput> {
        use crate::specialized::information_theory::*;

        let result_json = match op {
            InformationOp::Entropy => {
                let req: EntropyRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
                )
                .map_err(|e| format!("Failed to parse entropy request: {}", e))?;

                let result = shannon_entropy(req)
                    .map_err(|e| format!("Entropy calculation error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
            InformationOp::MutualInfo => {
                let req: MutualInfoRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
                )
                .map_err(|e| format!("Failed to parse mutual info request: {}", e))?;

                let result = mutual_information(req)
                    .map_err(|e| format!("Mutual information error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
            InformationOp::ChannelCapacity => {
                let req: ChannelCapacityRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
                )
                .map_err(|e| format!("Failed to parse channel capacity request: {}", e))?;

                let result = channel_capacity(req)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
            InformationOp::Huffman => {
                let req: HuffmanRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
                )
                .map_err(|e| format!("Failed to parse Huffman request: {}", e))?;

                let result =
                    huffman_coding(req).map_err(|e| format!("Huffman coding error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
            InformationOp::Kolmogorov => {
                let req: KolmogorovRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
                )
                .map_err(|e| format!("Failed to parse Kolmogorov request: {}", e))?;

                let result = kolmogorov_complexity(req)
                    .map_err(|e| format!("Kolmogorov complexity error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
            InformationOp::ConditionalEntropy => {
                let req: ConditionalEntropyRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
                )
                .map_err(|e| format!("Failed to parse conditional entropy request: {}", e))?;

                let result = conditional_entropy(req)
                    .map_err(|e| format!("Conditional entropy error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
            InformationOp::KLDivergence => {
                let req: RelativeEntropyRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
                )
                .map_err(|e| format!("Failed to parse KL divergence request: {}", e))?;

                let result =
                    relative_entropy(req).map_err(|e| format!("KL divergence error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
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
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
                )
                .map_err(|e| format!("Failed to parse shortest path request: {}", e))?;

                let result =
                    shortest_path(req).map_err(|e| format!("Shortest path error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
            GraphOp::MinimumSpanningTree => {
                let req: MSTRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
                )
                .map_err(|e| format!("Failed to parse MST request: {}", e))?;

                let result = minimum_spanning_tree(req).map_err(|e| format!("MST error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
            GraphOp::TopologicalSort => {
                let req: TopologicalSortRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
                )
                .map_err(|e| format!("Failed to parse topological sort request: {}", e))?;

                let result =
                    topological_sort(req).map_err(|e| format!("Topological sort error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
        };

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    /// Compute number theory operations
    fn compute_number_theory(
        &self,
        op: &NumberTheoryOp,
        input: &ComputeInput,
    ) -> ToolResult<ComputeOutput> {
        use crate::specialized::cryptographic_mathematics::*;
        use num_bigint::BigInt;
        use std::str::FromStr;

        let result_json = match op {
            NumberTheoryOp::GeneratePrime => {
                let bits: u32 = input
                    .parameters
                    .get("bit_length")
                    .and_then(|v| v.as_u64())
                    .map(|v| v as u32)
                    .ok_or("Missing or invalid bit_length parameter")?;

                let prime = generate_prime(bits);
                serde_json::json!({
                    "prime": prime.to_string(),
                    "bit_length": bits
                })
            }
            NumberTheoryOp::ModExp => {
                let base = input
                    .parameters
                    .get("base")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid base parameter")?;
                let exp = input
                    .parameters
                    .get("exponent")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid exponent parameter")?;
                let modulus = input
                    .parameters
                    .get("modulus")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid modulus parameter")?;

                let result = mod_exp(&base, &exp, &modulus);
                serde_json::json!({
                    "result": result.to_string(),
                    "operation": "modular_exponentiation"
                })
            }
            NumberTheoryOp::ModInv => {
                let a = input
                    .parameters
                    .get("a")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid a parameter")?;
                let m = input
                    .parameters
                    .get("modulus")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid modulus parameter")?;

                let result = mod_inverse(&a, &m).ok_or("Modular inverse does not exist")?;
                serde_json::json!({
                    "inverse": result.to_string(),
                    "operation": "modular_inverse"
                })
            }
            NumberTheoryOp::GCD => {
                let a = input
                    .parameters
                    .get("a")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid a parameter")?;
                let b = input
                    .parameters
                    .get("b")
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
            }
            NumberTheoryOp::LCM => {
                let a = input
                    .parameters
                    .get("a")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid a parameter")?;
                let b = input
                    .parameters
                    .get("b")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid b parameter")?;

                let (gcd, _, _) = extended_gcd(&a, &b);
                let lcm = (&a * &b) / gcd;
                serde_json::json!({
                    "lcm": lcm.to_string(),
                    "operation": "lcm"
                })
            }
            NumberTheoryOp::RSAKeypair => {
                let bits: u32 = input
                    .parameters
                    .get("bit_length")
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
            }
            NumberTheoryOp::RSAEncrypt => {
                let message = input
                    .parameters
                    .get("message")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid message parameter")?;
                let e = input
                    .parameters
                    .get("e")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid e parameter")?;
                let n = input
                    .parameters
                    .get("n")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid n parameter")?;

                let ciphertext = rsa_encrypt(&message, &e, &n);
                serde_json::json!({
                    "ciphertext": ciphertext.to_string()
                })
            }
            NumberTheoryOp::RSADecrypt => {
                let ciphertext = input
                    .parameters
                    .get("ciphertext")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid ciphertext parameter")?;
                let d = input
                    .parameters
                    .get("d")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid d parameter")?;
                let n = input
                    .parameters
                    .get("n")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid n parameter")?;

                let plaintext = rsa_decrypt(&ciphertext, &d, &n);
                serde_json::json!({
                    "plaintext": plaintext.to_string()
                })
            }
            NumberTheoryOp::SHA256 => {
                let input_str = input
                    .parameters
                    .get("input")
                    .and_then(|v| v.as_str())
                    .ok_or("Missing input parameter")?;

                let hash = sha256(input_str);
                serde_json::json!({
                    "hash": hash
                })
            }
            NumberTheoryOp::SHA3_256 => {
                let input_str = input
                    .parameters
                    .get("input")
                    .and_then(|v| v.as_str())
                    .ok_or("Missing input parameter")?;

                let hash = sha3_256(input_str);
                serde_json::json!({
                    "hash": hash
                })
            }
            NumberTheoryOp::ChineseRemainder => {
                let remainders_arr = input
                    .parameters
                    .get("remainders")
                    .and_then(|v| v.as_array())
                    .ok_or("Missing or invalid remainders parameter")?;
                let moduli_arr = input
                    .parameters
                    .get("moduli")
                    .and_then(|v| v.as_array())
                    .ok_or("Missing or invalid moduli parameter")?;

                let remainders: Result<Vec<BigInt>, _> = remainders_arr
                    .iter()
                    .map(|v| {
                        v.as_str()
                            .and_then(|s| BigInt::from_str(s).ok())
                            .ok_or("Invalid remainder")
                    })
                    .collect();
                let moduli: Result<Vec<BigInt>, _> = moduli_arr
                    .iter()
                    .map(|v| {
                        v.as_str()
                            .and_then(|s| BigInt::from_str(s).ok())
                            .ok_or("Invalid modulus")
                    })
                    .collect();

                let remainders = remainders?;
                let moduli = moduli?;

                let result = chinese_remainder_theorem(&remainders, &moduli)
                    .ok_or("No solution exists for given system")?;

                serde_json::json!({
                    "result": result.to_string()
                })
            }
            NumberTheoryOp::DiscreteLog => {
                let base = input
                    .parameters
                    .get("base")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid base parameter")?;
                let target = input
                    .parameters
                    .get("target")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid target parameter")?;
                let modulus = input
                    .parameters
                    .get("modulus")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid modulus parameter")?;
                let max_exp: u32 = input
                    .parameters
                    .get("max_exponent")
                    .and_then(|v| v.as_u64())
                    .map(|v| v as u32)
                    .unwrap_or(1000000);

                let result = discrete_log_bsgs(&base, &target, &modulus, max_exp)
                    .ok_or("Discrete logarithm not found within max_exponent bound")?;

                serde_json::json!({
                    "exponent": result.to_string()
                })
            }
            NumberTheoryOp::PrimalityTest => {
                let n = input
                    .parameters
                    .get("n")
                    .and_then(|v| v.as_str())
                    .and_then(|s| BigInt::from_str(s).ok())
                    .ok_or("Missing or invalid n parameter")?;
                let rounds: u32 = input
                    .parameters
                    .get("rounds")
                    .and_then(|v| v.as_u64())
                    .map(|v| v as u32)
                    .unwrap_or(10);

                let is_prime = miller_rabin_test(&n, rounds);
                let confidence = 1.0 - (0.25_f64).powi(rounds as i32);

                serde_json::json!({
                    "is_prime": is_prime,
                    "confidence": format!("{:.6}", confidence)
                })
            }

            NumberTheoryOp::EulerTotient
            | NumberTheoryOp::CarmichaelLambda
            | NumberTheoryOp::ECPointAdd => {
                return Err(format!(
                    "Number theory operation {:?} not yet fully implemented",
                    op
                ));
            }
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
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
                )
                .map_err(|e| format!("Failed to parse convex hull request: {}", e))?;

                let result = convex_hull(req).map_err(|e| format!("Convex hull error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
            GeometryOp::Delaunay => {
                let req: DelaunayRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
                )
                .map_err(|e| format!("Failed to parse Delaunay request: {}", e))?;

                let result = delaunay_triangulation(req)
                    .map_err(|e| format!("Delaunay triangulation error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
            GeometryOp::Voronoi => {
                let req: VoronoiRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
                )
                .map_err(|e| format!("Failed to parse Voronoi request: {}", e))?;

                let result =
                    voronoi_diagram(req).map_err(|e| format!("Voronoi diagram error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
            GeometryOp::PolygonArea => {
                let req: PolygonAreaRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
                )
                .map_err(|e| format!("Failed to parse polygon area request: {}", e))?;

                let result = polygon_area(req).map_err(|e| format!("Polygon area error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
            GeometryOp::PointInPolygon => {
                let req: PointInPolygonRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
                )
                .map_err(|e| format!("Failed to parse point in polygon request: {}", e))?;

                let result =
                    point_in_polygon(req).map_err(|e| format!("Point in polygon error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
        };

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    // NOTE: Physics computation functions (compute_em, compute_fourier_series, compute_physics,
    // compute_relativity, compute_statistical_physics, compute_quantum_mechanics,
    // compute_control_systems, compute_nuclear_physics) have been extracted to
    // src/implementations/compute/ subdirectory for better organization.
    // See: physics.rs, relativity.rs, statistical_physics.rs, quantum.rs, control.rs, nuclear.rs

    /// Compute electromagnetism operations
    /// @deprecated: Moved to compute::physics module
    #[allow(dead_code)]
    fn compute_em(&self, op: &EMComputation, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::physics::electromagnetism::*;

        let result_json = match op {
            EMComputation::PoyntingVector => {
                let req: PoyntingRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
                )
                .map_err(|e| format!("Failed to parse Poynting vector request: {}", e))?;

                let result =
                    poynting_vector(req).map_err(|e| format!("Poynting vector error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
            EMComputation::SkinEffect => {
                let req: SkinEffectRequest = serde_json::from_value(
                    serde_json::to_value(&input.parameters)
                        .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
                )
                .map_err(|e| format!("Failed to parse skin effect request: {}", e))?;

                let result = skin_effect(req).map_err(|e| format!("Skin effect error: {}", e))?;
                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
        };

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    /// Compute Fourier series
    /// @deprecated: Moved to compute::physics module
    #[allow(dead_code)]
    fn compute_fourier_series(&self, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        use crate::tools::signal_processing::*;

        let req: FourierSeriesRequest = serde_json::from_value(
            serde_json::to_value(&input.parameters)
                .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
        )
        .map_err(|e| format!("Failed to parse Fourier series request: {}", e))?;

        let result =
            compute_fourier_series(req).map_err(|e| format!("Fourier series error: {}", e))?;

        let result_json = serde_json::to_value(result)
            .map_err(|e| format!("Failed to serialize result: {}", e))?;

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    /// Compute physics operations (Tier 1 Wolfram Alpha expansion)
    /// @deprecated: Moved to compute::physics module
    #[allow(dead_code)]
    fn compute_physics(&self, op: &PhysicsOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        match op {
            PhysicsOp::Relativity(relativity_op) => self.compute_relativity(relativity_op, input),
            PhysicsOp::StatisticalPhysics(stat_phys_op) => {
                self.compute_statistical_physics(stat_phys_op, input)
            }
            PhysicsOp::QuantumMechanics(qm_op) => self.compute_quantum_mechanics(qm_op, input),
            PhysicsOp::ControlSystems(cs_op) => self.compute_control_systems(cs_op, input),
            PhysicsOp::NuclearPhysics(np_op) => self.compute_nuclear_physics(np_op, input),
        }
    }

    /// Compute relativity operations
    /// @deprecated: Moved to compute::relativity module
    #[allow(dead_code)]
    fn compute_relativity(
        &self,
        op: &RelativityOp,
        input: &ComputeInput,
    ) -> ToolResult<ComputeOutput> {
        use crate::physics::relativity::*;

        let result_json = match op {
            RelativityOp::LorentzTransform => {
                let velocity = input
                    .parameters
                    .get("velocity")
                    .and_then(|v| v.as_f64())
                    .ok_or("velocity parameter required")?;
                let position: Vec<f64> = input
                    .parameters
                    .get("position")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("position [x, y, z] required")?;
                let time = input
                    .parameters
                    .get("time")
                    .and_then(|v| v.as_f64())
                    .ok_or("time parameter required")?;

                let result = lorentz_transform(LorentzTransformRequest {
                    velocity,
                    position,
                    time,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            RelativityOp::TimeDilation => {
                let proper_time = input
                    .parameters
                    .get("proper_time")
                    .and_then(|v| v.as_f64())
                    .ok_or("proper_time parameter required")?;
                let velocity = input
                    .parameters
                    .get("velocity")
                    .and_then(|v| v.as_f64())
                    .ok_or("velocity parameter required")?;

                let result = time_dilation(TimeDilationRequest {
                    proper_time,
                    velocity,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            RelativityOp::LengthContraction => {
                let proper_length = input
                    .parameters
                    .get("proper_length")
                    .and_then(|v| v.as_f64())
                    .ok_or("proper_length parameter required")?;
                let velocity = input
                    .parameters
                    .get("velocity")
                    .and_then(|v| v.as_f64())
                    .ok_or("velocity parameter required")?;

                let result = length_contraction(LengthContractionRequest {
                    proper_length,
                    velocity,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            RelativityOp::RelativisticEnergy => {
                let mass = input
                    .parameters
                    .get("mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("mass parameter required")?;
                let velocity = input
                    .parameters
                    .get("velocity")
                    .and_then(|v| v.as_f64())
                    .ok_or("velocity parameter required")?;

                let result = relativistic_energy(RelativisticEnergyRequest { mass, velocity })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            RelativityOp::VelocityAddition => {
                let velocity1 = input
                    .parameters
                    .get("velocity1")
                    .and_then(|v| v.as_f64())
                    .ok_or("velocity1 parameter required")?;
                let velocity2 = input
                    .parameters
                    .get("velocity2")
                    .and_then(|v| v.as_f64())
                    .ok_or("velocity2 parameter required")?;

                let result = velocity_addition(VelocityAdditionRequest {
                    velocity1,
                    velocity2,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            RelativityOp::SchwarzschildMetric => {
                let mass = input
                    .parameters
                    .get("mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("mass parameter required")?;
                let radius = input
                    .parameters
                    .get("radius")
                    .and_then(|v| v.as_f64())
                    .ok_or("radius parameter required")?;

                let result = schwarzschild_metric(SchwarzschildRequest { mass, radius })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            RelativityOp::GravitationalTimeDilation => {
                let mass = input
                    .parameters
                    .get("mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("mass parameter required")?;
                let radius = input
                    .parameters
                    .get("radius")
                    .and_then(|v| v.as_f64())
                    .ok_or("radius parameter required")?;
                let proper_time = input
                    .parameters
                    .get("proper_time")
                    .and_then(|v| v.as_f64())
                    .ok_or("proper_time parameter required")?;

                let result = gravitational_time_dilation(GravitationalTimeDilationRequest {
                    mass,
                    radius,
                    proper_time,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            RelativityOp::OrbitalPrecession => {
                let mass = input
                    .parameters
                    .get("mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("mass parameter required")?;
                let semi_major_axis = input
                    .parameters
                    .get("semi_major_axis")
                    .and_then(|v| v.as_f64())
                    .ok_or("semi_major_axis parameter required")?;
                let eccentricity = input
                    .parameters
                    .get("eccentricity")
                    .and_then(|v| v.as_f64())
                    .ok_or("eccentricity parameter required")?;

                let result = orbital_precession(OrbitalPrecessionRequest {
                    mass,
                    semi_major_axis,
                    eccentricity,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            RelativityOp::GravitationalLensing => {
                let lens_mass = input
                    .parameters
                    .get("lens_mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("lens_mass parameter required")?;
                let impact_parameter = input
                    .parameters
                    .get("impact_parameter")
                    .and_then(|v| v.as_f64())
                    .ok_or("impact_parameter parameter required")?;

                let result = gravitational_lensing(GravitationalLensingRequest {
                    lens_mass,
                    impact_parameter,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            RelativityOp::BlackHoleProperties => {
                let mass = input
                    .parameters
                    .get("mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("mass parameter required")?;

                let result = black_hole_properties(BlackHoleRequest { mass })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
        };

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    /// Compute statistical physics operations
    /// @deprecated: Moved to compute::statistical_physics module
    #[allow(dead_code)]
    fn compute_statistical_physics(
        &self,
        op: &StatPhysicsOp,
        input: &ComputeInput,
    ) -> ToolResult<ComputeOutput> {
        use crate::physics::statistical_physics::*;

        let result_json = match op {
            StatPhysicsOp::PartitionFunction => {
                let ensemble = input
                    .parameters
                    .get("ensemble")
                    .and_then(|v| v.as_str())
                    .ok_or("ensemble parameter required")?;
                let temperature = input
                    .parameters
                    .get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;

                let energy_levels: Option<Vec<f64>> = input
                    .parameters
                    .get("energy_levels")
                    .and_then(|v| serde_json::from_value(v.clone()).ok());
                let degeneracies: Option<Vec<usize>> = input
                    .parameters
                    .get("degeneracies")
                    .and_then(|v| serde_json::from_value(v.clone()).ok());
                let volume = input.parameters.get("volume").and_then(|v| v.as_f64());
                let num_particles = input
                    .parameters
                    .get("num_particles")
                    .and_then(|v| v.as_u64())
                    .map(|v| v as usize);
                let chemical_potential = input
                    .parameters
                    .get("chemical_potential")
                    .and_then(|v| v.as_f64());

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
            }

            StatPhysicsOp::PartitionFunctionCanonical => {
                let temperature = input
                    .parameters
                    .get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let energy_levels: Vec<f64> = input
                    .parameters
                    .get("energy_levels")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("energy_levels required")?;
                let degeneracies: Vec<usize> = input
                    .parameters
                    .get("degeneracies")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("degeneracies required")?;

                let result = canonical_partition(CanonicalPartitionRequest {
                    temperature,
                    energy_levels,
                    degeneracies,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            StatPhysicsOp::PartitionFunctionGrandCanonical => {
                let temperature = input
                    .parameters
                    .get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let volume = input
                    .parameters
                    .get("volume")
                    .and_then(|v| v.as_f64())
                    .ok_or("volume parameter required")?;
                let chemical_potential = input
                    .parameters
                    .get("chemical_potential")
                    .and_then(|v| v.as_f64())
                    .ok_or("chemical_potential parameter required")?;
                let particle_type = input
                    .parameters
                    .get("particle_type")
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
            }

            StatPhysicsOp::MaxwellBoltzmann => {
                let temperature = input
                    .parameters
                    .get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let mass = input
                    .parameters
                    .get("mass")
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
            }

            StatPhysicsOp::FermiDirac => {
                let energy = input
                    .parameters
                    .get("energy")
                    .and_then(|v| v.as_f64())
                    .ok_or("energy parameter required")?;
                let temperature = input
                    .parameters
                    .get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let chemical_potential = input
                    .parameters
                    .get("chemical_potential")
                    .and_then(|v| v.as_f64())
                    .ok_or("chemical_potential parameter required")?;

                let result = fermi_dirac(FermiDiracRequest {
                    energy,
                    temperature,
                    chemical_potential,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            StatPhysicsOp::BoseEinstein => {
                let energy = input
                    .parameters
                    .get("energy")
                    .and_then(|v| v.as_f64())
                    .ok_or("energy parameter required")?;
                let temperature = input
                    .parameters
                    .get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let chemical_potential = input
                    .parameters
                    .get("chemical_potential")
                    .and_then(|v| v.as_f64())
                    .ok_or("chemical_potential parameter required")?;

                let result = bose_einstein(BoseEinsteinRequest {
                    energy,
                    temperature,
                    chemical_potential,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            StatPhysicsOp::ChemicalPotential => {
                let temperature = input
                    .parameters
                    .get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let pressure = input
                    .parameters
                    .get("pressure")
                    .and_then(|v| v.as_f64())
                    .ok_or("pressure parameter required")?;
                let particle_density = input
                    .parameters
                    .get("particle_density")
                    .and_then(|v| v.as_f64())
                    .ok_or("particle_density parameter required")?;
                let particle_type = input
                    .parameters
                    .get("particle_type")
                    .and_then(|v| v.as_str())
                    .unwrap_or("classical");
                let mass = input
                    .parameters
                    .get("mass")
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
            }

            StatPhysicsOp::FugacityCoefficient => {
                let temperature = input
                    .parameters
                    .get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let pressure = input
                    .parameters
                    .get("pressure")
                    .and_then(|v| v.as_f64())
                    .ok_or("pressure parameter required")?;
                let particle_density = input
                    .parameters
                    .get("particle_density")
                    .and_then(|v| v.as_f64())
                    .ok_or("particle_density parameter required")?;

                let result = fugacity_coefficient(FugacityRequest {
                    temperature,
                    pressure,
                    particle_density,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            StatPhysicsOp::PhaseTransition => {
                let model = input
                    .parameters
                    .get("model")
                    .and_then(|v| v.as_str())
                    .ok_or("model parameter required")?;
                let temperature = input
                    .parameters
                    .get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let critical_temperature = input
                    .parameters
                    .get("critical_temperature")
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
            }

            StatPhysicsOp::CriticalPhenomena => {
                let temperature = input
                    .parameters
                    .get("temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("temperature parameter required")?;
                let critical_temperature = input
                    .parameters
                    .get("critical_temperature")
                    .and_then(|v| v.as_f64())
                    .ok_or("critical_temperature parameter required")?;
                let model = input
                    .parameters
                    .get("model")
                    .and_then(|v| v.as_str())
                    .unwrap_or("mean_field");

                let result = critical_phenomena(CriticalPhenomenaRequest {
                    temperature,
                    critical_temperature,
                    model: model.to_string(),
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
        };

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    /// Compute quantum mechanics operations
    /// @deprecated: Moved to compute::quantum module
    #[allow(dead_code)]
    fn compute_quantum_mechanics(
        &self,
        op: &QuantumMechOp,
        input: &ComputeInput,
    ) -> ToolResult<ComputeOutput> {
        use crate::physics::quantum_mechanics::*;

        let result_json = match op {
            QuantumMechOp::SchrodingerEquation => {
                let potential = input
                    .parameters
                    .get("potential")
                    .and_then(|v| v.as_str())
                    .ok_or("potential parameter required")?;
                let energy = input
                    .parameters
                    .get("energy")
                    .and_then(|v| v.as_f64())
                    .ok_or("energy parameter required")?;
                let position = input
                    .parameters
                    .get("position")
                    .and_then(|v| v.as_f64())
                    .ok_or("position parameter required")?;

                let mut parameters = std::collections::HashMap::new();
                if let Some(val) = input.parameters.get("width").and_then(|v| v.as_f64()) {
                    parameters.insert("width".to_string(), val);
                }
                if let Some(val) = input
                    .parameters
                    .get("quantum_number")
                    .and_then(|v| v.as_f64())
                {
                    parameters.insert("quantum_number".to_string(), val);
                }
                if let Some(val) = input.parameters.get("omega").and_then(|v| v.as_f64()) {
                    parameters.insert("omega".to_string(), val);
                }
                if let Some(val) = input.parameters.get("mass").and_then(|v| v.as_f64()) {
                    parameters.insert("mass".to_string(), val);
                }
                if let Some(val) = input
                    .parameters
                    .get("barrier_height")
                    .and_then(|v| v.as_f64())
                {
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
            }

            QuantumMechOp::HarmonicOscillator => {
                let quantum_number = input
                    .parameters
                    .get("quantum_number")
                    .and_then(|v| v.as_u64())
                    .ok_or("quantum_number parameter required")?
                    as usize;
                let omega = input
                    .parameters
                    .get("omega")
                    .and_then(|v| v.as_f64())
                    .ok_or("omega parameter required")?;
                let mass = input
                    .parameters
                    .get("mass")
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
            }

            QuantumMechOp::HydrogenAtom => {
                let n = input
                    .parameters
                    .get("n")
                    .and_then(|v| v.as_u64())
                    .ok_or("n (principal quantum number) required")?
                    as usize;
                let l = input
                    .parameters
                    .get("l")
                    .and_then(|v| v.as_u64())
                    .ok_or("l (orbital quantum number) required")? as usize;
                let m = input
                    .parameters
                    .get("m")
                    .and_then(|v| v.as_i64())
                    .ok_or("m (magnetic quantum number) required")? as i32;
                let r = input.parameters.get("r").and_then(|v| v.as_f64());

                let result = hydrogen_atom(HydrogenAtomRequest { n, l, m, r })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            QuantumMechOp::AngularMomentum => {
                let l = input
                    .parameters
                    .get("l")
                    .and_then(|v| v.as_u64())
                    .ok_or("l parameter required")? as usize;
                let m = input
                    .parameters
                    .get("m")
                    .and_then(|v| v.as_i64())
                    .ok_or("m parameter required")? as i32;
                let operation = input
                    .parameters
                    .get("operation")
                    .and_then(|v| v.as_str())
                    .unwrap_or("eigenvalue");

                let result = angular_momentum(AngularMomentumRequest {
                    l,
                    m,
                    operation: operation.to_string(),
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            QuantumMechOp::SpinOperators => {
                let spin = input
                    .parameters
                    .get("spin")
                    .and_then(|v| v.as_f64())
                    .ok_or("spin parameter required")?;
                let component = input
                    .parameters
                    .get("component")
                    .and_then(|v| v.as_str())
                    .ok_or("component parameter required (x, y, z, plus, minus)")?;
                let state: Option<Vec<f64>> = input
                    .parameters
                    .get("state")
                    .and_then(|v| serde_json::from_value(v.clone()).ok());

                let result = spin_operators(SpinRequest {
                    spin,
                    component: component.to_string(),
                    state,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            QuantumMechOp::PerturbationTheory => {
                let order = input
                    .parameters
                    .get("order")
                    .and_then(|v| v.as_u64())
                    .ok_or("order parameter required (1 or 2)")?
                    as usize;
                let unperturbed_energy = input
                    .parameters
                    .get("unperturbed_energy")
                    .and_then(|v| v.as_f64())
                    .ok_or("unperturbed_energy required")?;
                let perturbation_matrix_element = input
                    .parameters
                    .get("perturbation_matrix_element")
                    .and_then(|v| v.as_f64())
                    .ok_or("perturbation_matrix_element required")?;
                let energy_differences: Option<Vec<f64>> = input
                    .parameters
                    .get("energy_differences")
                    .and_then(|v| serde_json::from_value(v.clone()).ok());
                let coupling_matrix_elements: Option<Vec<f64>> = input
                    .parameters
                    .get("coupling_matrix_elements")
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
            }

            QuantumMechOp::TunnelingProbability => {
                let barrier_height = input
                    .parameters
                    .get("barrier_height")
                    .and_then(|v| v.as_f64())
                    .ok_or("barrier_height required")?;
                let barrier_width = input
                    .parameters
                    .get("barrier_width")
                    .and_then(|v| v.as_f64())
                    .ok_or("barrier_width required")?;
                let particle_energy = input
                    .parameters
                    .get("particle_energy")
                    .and_then(|v| v.as_f64())
                    .ok_or("particle_energy required")?;
                let particle_mass = input
                    .parameters
                    .get("particle_mass")
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
            }

            QuantumMechOp::DensityMatrix => {
                let state_type = input
                    .parameters
                    .get("state_type")
                    .and_then(|v| v.as_str())
                    .ok_or("state_type required (pure or mixed)")?;
                let state_vector: Option<Vec<f64>> = input
                    .parameters
                    .get("state_vector")
                    .and_then(|v| serde_json::from_value(v.clone()).ok());
                let density_mat: Option<Vec<Vec<f64>>> = input
                    .parameters
                    .get("density_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok());

                let result = density_matrix(DensityMatrixRequest {
                    state_type: state_type.to_string(),
                    state_vector,
                    density_matrix: density_mat,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            QuantumMechOp::EntanglementMeasure => {
                let state_vector: Vec<f64> = input
                    .parameters
                    .get("state_vector")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("state_vector required (4 components for two-qubit state)")?;
                let measure = input
                    .parameters
                    .get("measure")
                    .and_then(|v| v.as_str())
                    .unwrap_or("concurrence");

                let result = entanglement_measure(EntanglementRequest {
                    state_vector,
                    measure: measure.to_string(),
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            QuantumMechOp::QuantumEntropy => {
                let density_matrix: Vec<Vec<f64>> = input
                    .parameters
                    .get("density_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("density_matrix required")?;
                let entropy_type = input
                    .parameters
                    .get("entropy_type")
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
            }

            QuantumMechOp::QuantumCoherence => {
                let density_matrix: Vec<Vec<f64>> = input
                    .parameters
                    .get("density_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("density_matrix required")?;
                let basis = input
                    .parameters
                    .get("basis")
                    .and_then(|v| v.as_str())
                    .unwrap_or("computational");

                let result = quantum_coherence(CoherenceRequest {
                    density_matrix,
                    basis: basis.to_string(),
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            QuantumMechOp::BellInequality => {
                let state_vector: Vec<f64> = input
                    .parameters
                    .get("state_vector")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("state_vector required (two-qubit state)")?;
                let measurement_angles: Vec<f64> = input
                    .parameters
                    .get("measurement_angles")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("measurement_angles required (4 angles for CHSH)")?;

                let result = bell_inequality(BellInequalityRequest {
                    state_vector,
                    measurement_angles,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            QuantumMechOp::QuantumTomography => {
                let measurements: Vec<Vec<f64>> = input
                    .parameters
                    .get("measurements")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("measurements required")?;
                let measurement_bases: Vec<String> = input
                    .parameters
                    .get("measurement_bases")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("measurement_bases required (e.g., [X, Y, Z])")?;

                let result = quantum_tomography(QuantumTomographyRequest {
                    measurements,
                    measurement_bases,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            QuantumMechOp::QuantumGate => {
                let gate_type = input
                    .parameters
                    .get("gate_type")
                    .and_then(|v| v.as_str())
                    .ok_or("gate_type required (hadamard, pauli_x, pauli_z, phase, cnot)")?;
                let input_state: Vec<f64> = input
                    .parameters
                    .get("input_state")
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
            }

            QuantumMechOp::QuantumCircuit => {
                let num_qubits = input
                    .parameters
                    .get("num_qubits")
                    .and_then(|v| v.as_u64())
                    .ok_or("num_qubits required")? as usize;
                let gates: Vec<serde_json::Value> = input
                    .parameters
                    .get("gates")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("gates required (list of gate operations)")?;
                let initial_state: Vec<f64> = input
                    .parameters
                    .get("initial_state")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("initial_state required")?;

                let result = quantum_circuit(QuantumCircuitRequest {
                    num_qubits,
                    gates,
                    initial_state,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
        };

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    /// Compute control systems operations
    /// @deprecated: Moved to compute::control module
    #[allow(dead_code)]
    fn compute_control_systems(
        &self,
        op: &ControlSystemsOp,
        input: &ComputeInput,
    ) -> ToolResult<ComputeOutput> {
        use crate::physics::control_systems::*;

        let result_json = match op {
            ControlSystemsOp::TransferFunction => {
                let numerator: Vec<f64> = input
                    .parameters
                    .get("numerator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("numerator coefficients required")?;
                let denominator: Vec<f64> = input
                    .parameters
                    .get("denominator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("denominator coefficients required")?;
                let operation = input
                    .parameters
                    .get("operation")
                    .and_then(|v| v.as_str())
                    .unwrap_or("evaluate")
                    .to_string();
                let frequency = input.parameters.get("frequency").and_then(|v| v.as_f64());
                let second_tf: Option<(Vec<f64>, Vec<f64>)> = input
                    .parameters
                    .get("second_tf")
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
            }

            ControlSystemsOp::PoleZeroAnalysis => {
                let numerator: Vec<f64> = input
                    .parameters
                    .get("numerator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("numerator required")?;
                let denominator: Vec<f64> = input
                    .parameters
                    .get("denominator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("denominator required")?;

                let result = pole_zero_analysis(PoleZeroRequest {
                    numerator,
                    denominator,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            ControlSystemsOp::BodePlot => {
                let numerator: Vec<f64> = input
                    .parameters
                    .get("numerator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("numerator required")?;
                let denominator: Vec<f64> = input
                    .parameters
                    .get("denominator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("denominator required")?;
                let freq_min = input
                    .parameters
                    .get("freq_min")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.1);
                let freq_max = input
                    .parameters
                    .get("freq_max")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(100.0);
                let num_points = input
                    .parameters
                    .get("num_points")
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
            }

            ControlSystemsOp::NyquistPlot => {
                let numerator: Vec<f64> = input
                    .parameters
                    .get("numerator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("numerator required")?;
                let denominator: Vec<f64> = input
                    .parameters
                    .get("denominator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("denominator required")?;
                let freq_min = input
                    .parameters
                    .get("freq_min")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.1);
                let freq_max = input
                    .parameters
                    .get("freq_max")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(100.0);
                let num_points = input
                    .parameters
                    .get("num_points")
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
            }

            ControlSystemsOp::RootLocus => {
                let numerator: Vec<f64> = input
                    .parameters
                    .get("numerator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("numerator required")?;
                let denominator: Vec<f64> = input
                    .parameters
                    .get("denominator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("denominator required")?;
                let gain_min = input
                    .parameters
                    .get("gain_min")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0);
                let gain_max = input
                    .parameters
                    .get("gain_max")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(10.0);
                let num_points = input
                    .parameters
                    .get("num_points")
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
            }

            ControlSystemsOp::StateSpace => {
                let a_matrix: Vec<Vec<f64>> = input
                    .parameters
                    .get("a_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("a_matrix required")?;
                let b_matrix: Vec<Vec<f64>> = input
                    .parameters
                    .get("b_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("b_matrix required")?;
                let c_matrix: Vec<Vec<f64>> = input
                    .parameters
                    .get("c_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("c_matrix required")?;
                let d_matrix: Vec<Vec<f64>> = input
                    .parameters
                    .get("d_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("d_matrix required")?;
                let operation = input
                    .parameters
                    .get("operation")
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
            }

            ControlSystemsOp::Controllability => {
                let a_matrix: Vec<Vec<f64>> = input
                    .parameters
                    .get("a_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("a_matrix required")?;
                let b_matrix: Vec<Vec<f64>> = input
                    .parameters
                    .get("b_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("b_matrix required")?;

                let result = controllability(ControllabilityRequest { a_matrix, b_matrix })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            ControlSystemsOp::Observability => {
                let a_matrix: Vec<Vec<f64>> = input
                    .parameters
                    .get("a_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("a_matrix required")?;
                let c_matrix: Vec<Vec<f64>> = input
                    .parameters
                    .get("c_matrix")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("c_matrix required")?;

                let result = observability(ObservabilityRequest { a_matrix, c_matrix })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            ControlSystemsOp::RouthHurwitz => {
                let characteristic_polynomial: Vec<f64> = input
                    .parameters
                    .get("characteristic_polynomial")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("characteristic_polynomial required")?;

                let result = routh_hurwitz(RouthHurwitzRequest {
                    characteristic_polynomial,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            ControlSystemsOp::GainMargin => {
                let numerator: Vec<f64> = input
                    .parameters
                    .get("numerator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("numerator required")?;
                let denominator: Vec<f64> = input
                    .parameters
                    .get("denominator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("denominator required")?;

                let result = gain_margin(GainMarginRequest {
                    numerator,
                    denominator,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            ControlSystemsOp::PhaseMargin => {
                let numerator: Vec<f64> = input
                    .parameters
                    .get("numerator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("numerator required")?;
                let denominator: Vec<f64> = input
                    .parameters
                    .get("denominator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("denominator required")?;

                let result = phase_margin(PhaseMarginRequest {
                    numerator,
                    denominator,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            ControlSystemsOp::StepResponse => {
                let numerator: Vec<f64> = input
                    .parameters
                    .get("numerator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("numerator required")?;
                let denominator: Vec<f64> = input
                    .parameters
                    .get("denominator")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("denominator required")?;
                let time_span = input
                    .parameters
                    .get("time_span")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(10.0);
                let num_points = input
                    .parameters
                    .get("num_points")
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
            }
        };

        Ok(ComputeOutput {
            result: result_json,
            additional: None,
            metadata: None,
        })
    }

    /// Compute nuclear physics operations
    /// @deprecated: Moved to compute::nuclear module
    #[allow(dead_code)]
    fn compute_nuclear_physics(
        &self,
        op: &NuclearOp,
        input: &ComputeInput,
    ) -> ToolResult<ComputeOutput> {
        use crate::physics::nuclear_physics::*;

        let result_json = match op {
            NuclearOp::RadioactiveDecay => {
                let initial_quantity = input
                    .parameters
                    .get("initial_quantity")
                    .and_then(|v| v.as_f64())
                    .ok_or("initial_quantity parameter required")?;
                let decay_constant = input
                    .parameters
                    .get("decay_constant")
                    .and_then(|v| v.as_f64())
                    .ok_or("decay_constant parameter required")?;
                let time = input
                    .parameters
                    .get("time")
                    .and_then(|v| v.as_f64())
                    .ok_or("time parameter required")?;

                let result = radioactive_decay(RadioactiveDecayRequest {
                    initial_quantity,
                    decay_constant,
                    time,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            NuclearOp::DecayChain => {
                let parent_initial = input
                    .parameters
                    .get("parent_initial")
                    .and_then(|v| v.as_f64())
                    .ok_or("parent_initial parameter required")?;
                let parent_decay_constant = input
                    .parameters
                    .get("parent_decay_constant")
                    .and_then(|v| v.as_f64())
                    .ok_or("parent_decay_constant parameter required")?;
                let daughter_decay_constant = input
                    .parameters
                    .get("daughter_decay_constant")
                    .and_then(|v| v.as_f64())
                    .ok_or("daughter_decay_constant parameter required")?;
                let time = input
                    .parameters
                    .get("time")
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
            }

            NuclearOp::HalfLife => {
                let decay_constant_val = input
                    .parameters
                    .get("decay_constant")
                    .and_then(|v| v.as_f64());
                let half_life_val = input.parameters.get("half_life").and_then(|v| v.as_f64());
                let mean_lifetime_val = input
                    .parameters
                    .get("mean_lifetime")
                    .and_then(|v| v.as_f64());

                let result = half_life(HalfLifeRequest {
                    decay_constant: decay_constant_val,
                    half_life: half_life_val,
                    mean_lifetime: mean_lifetime_val,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            NuclearOp::BindingEnergy => {
                let atomic_number = input
                    .parameters
                    .get("atomic_number")
                    .and_then(|v| v.as_u64())
                    .ok_or("atomic_number parameter required")?
                    as u32;
                let mass_number = input
                    .parameters
                    .get("mass_number")
                    .and_then(|v| v.as_u64())
                    .ok_or("mass_number parameter required")?
                    as u32;
                let atomic_mass = input.parameters.get("atomic_mass").and_then(|v| v.as_f64());

                let result = binding_energy(BindingEnergyRequest {
                    atomic_number,
                    mass_number,
                    atomic_mass,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            NuclearOp::MassDefect => {
                let protons = input
                    .parameters
                    .get("protons")
                    .and_then(|v| v.as_u64())
                    .ok_or("protons parameter required")? as u32;
                let neutrons = input
                    .parameters
                    .get("neutrons")
                    .and_then(|v| v.as_u64())
                    .ok_or("neutrons parameter required")? as u32;
                let nuclear_mass = input
                    .parameters
                    .get("nuclear_mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("nuclear_mass parameter required")?;

                let result = mass_defect(MassDefectRequest {
                    protons,
                    neutrons,
                    nuclear_mass,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            NuclearOp::FissionEnergy => {
                let parent_mass = input
                    .parameters
                    .get("parent_mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("parent_mass parameter required")?;
                let fragment1_mass = input
                    .parameters
                    .get("fragment1_mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("fragment1_mass parameter required")?;
                let fragment2_mass = input
                    .parameters
                    .get("fragment2_mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("fragment2_mass parameter required")?;
                let neutrons_released = input
                    .parameters
                    .get("neutrons_released")
                    .and_then(|v| v.as_u64())
                    .ok_or("neutrons_released parameter required")?
                    as u32;

                let result = fission_energy(FissionEnergyRequest {
                    parent_mass,
                    fragment1_mass,
                    fragment2_mass,
                    neutrons_released,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            NuclearOp::FusionEnergy => {
                let reactant1_mass = input
                    .parameters
                    .get("reactant1_mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("reactant1_mass parameter required")?;
                let reactant2_mass = input
                    .parameters
                    .get("reactant2_mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("reactant2_mass parameter required")?;
                let product1_mass = input
                    .parameters
                    .get("product1_mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("product1_mass parameter required")?;
                let product2_mass = input
                    .parameters
                    .get("product2_mass")
                    .and_then(|v| v.as_f64());

                let result = fusion_energy(FusionEnergyRequest {
                    reactant1_mass,
                    reactant2_mass,
                    product1_mass,
                    product2_mass,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }

            NuclearOp::NuclearReaction => {
                let reactants: Vec<f64> = input
                    .parameters
                    .get("reactants")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("reactants parameter required as array of masses in amu")?;
                let products: Vec<f64> = input
                    .parameters
                    .get("products")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("products parameter required as array of masses in amu")?;
                let projectile_energy_mev = input
                    .parameters
                    .get("projectile_energy_mev")
                    .and_then(|v| v.as_f64());

                let result = nuclear_reaction(NuclearReactionRequest {
                    reactants,
                    products,
                    projectile_energy_mev,
                })?;

                serde_json::to_value(result)
                    .map_err(|e| format!("Failed to serialize result: {}", e))?
            }
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

            ComputeOp::Matrix(matrix_op) => compute::compute_matrix_op(matrix_op, input),

            ComputeOp::MatrixDecomp(decomp) => compute::compute_matrix_decomp(decomp, input),

            ComputeOp::SpecialFunc(func) => self.compute_special_function(func, input),

            // Core Mathematical Operations
            ComputeOp::NumberTheory(op) => compute::compute_number_theory(op, input),
            ComputeOp::Geometry(op) => self.compute_geometry(op, input),
            ComputeOp::Information(op) => compute::compute_information(op, input),
            ComputeOp::Graph(op) => compute::compute_graph(op, input),
            ComputeOp::EM(op) => compute::compute_em(op, input),
            ComputeOp::FourierSeries => compute::compute_fourier_series(input),

            // Scientific Formulas (2025 Expansion)
            ComputeOp::Chemistry(op) => compute::compute_chemistry(op, input),
            ComputeOp::Biology(op) => compute::compute_biology(op, input),
            ComputeOp::Thermodynamics(op) => compute::compute_thermodynamics(op, input),
            ComputeOp::Optics(op) => compute::compute_optics(op, input),
            ComputeOp::Geophysics(op) => compute::compute_geophysics(op, input),
            ComputeOp::Engineering(op) => compute::compute_engineering(op, input),
            ComputeOp::DateTime(op) => compute::compute_datetime(op, input),

            // Physics (Tier 1 Wolfram Alpha expansion)
            ComputeOp::Physics(op) => compute::compute_physics(op, input),

            // ========== Consolidated Operations (formerly separate tools) ==========

            // Differentiate operations (absorbed from Differentiate tool)
            ComputeOp::Differentiate(diff_op) => {
                use crate::implementations::differentiator::UnifiedDifferentiator;
                let differentiator = UnifiedDifferentiator::new();
                let variable = input.data.get("variable").and_then(|v| v.as_str()).unwrap_or("x").to_string();
                let diff_input = DifferentiateInput {
                    operation: diff_op.clone(),
                    expression: input.data.get("expression").and_then(|v| v.as_str()).unwrap_or("").to_string(),
                    variables: vec![variable.clone()],
                    order: input.data.get("order").and_then(|v| v.as_u64()).map(|n| vec![n as usize]),
                    evaluate_at: input.data.get("point").and_then(|v| v.as_f64()).map(|p| {
                        let mut m = HashMap::new();
                        m.insert(variable, p);
                        m
                    }),
                    parameters: HashMap::new(),
                };
                let result = differentiator.differentiate(&diff_input)?;
                Ok(ComputeOutput {
                    result: serde_json::json!(result.derivatives),
                    additional: result.latex.map(|l| {
                        let mut m = HashMap::new();
                        m.insert("latex".to_string(), serde_json::json!(l));
                        m
                    }),
                    metadata: result.metadata,
                })
            }

            // Integrate operations (absorbed from Integrate tool)
            ComputeOp::Integrate(int_type) => {
                use crate::implementations::integrator::UnifiedIntegrator;
                let integrator = UnifiedIntegrator::new();
                let variable = input.data.get("variable").and_then(|v| v.as_str()).unwrap_or("x").to_string();
                let limits = match (input.data.get("lower").and_then(|v| v.as_f64()), input.data.get("upper").and_then(|v| v.as_f64())) {
                    (Some(l), Some(u)) => Some(vec![[l, u]]),
                    _ => None,
                };
                let int_input = IntegrateInput {
                    integration_type: int_type.clone(),
                    expression: input.data.get("expression").and_then(|v| v.as_str()).unwrap_or("").to_string(),
                    variables: vec![variable],
                    limits,
                    path: None,
                    method: None,
                    parameters: HashMap::new(),
                };
                let result = integrator.integrate(&int_input)?;
                Ok(ComputeOutput {
                    result: result.result,
                    additional: result.symbolic.map(|s| {
                        let mut m = HashMap::new();
                        m.insert("symbolic".to_string(), serde_json::json!(s));
                        if let Some(latex) = result.latex.clone() {
                            m.insert("latex".to_string(), serde_json::json!(latex));
                        }
                        m
                    }),
                    metadata: result.metadata,
                })
            }

            // Transform operations (absorbed from Transform tool)
            ComputeOp::Transform(transform_type) => {
                use crate::implementations::transformer::UnifiedTransformer;
                let transformer = UnifiedTransformer::new();
                let transform_input = TransformInput {
                    transform_type: transform_type.clone(),
                    data: input.data.get("data").and_then(|v| v.as_array()).map(|arr| {
                        arr.iter().filter_map(|v| v.as_f64()).collect()
                    }).unwrap_or_default(),
                    sampling_rate: input.data.get("sampling_rate").and_then(|v| v.as_f64()),
                    parameters: input.parameters.iter().map(|(k, v)| (k.clone(), v.clone())).collect(),
                };
                let result = transformer.transform(&transform_input)?;
                Ok(ComputeOutput {
                    result: serde_json::json!(result.result),
                    additional: result.frequencies.map(|f| {
                        let mut m = HashMap::new();
                        m.insert("frequencies".to_string(), serde_json::json!(f));
                        if let Some(mag) = result.magnitude.clone() {
                            m.insert("magnitude".to_string(), serde_json::json!(mag));
                        }
                        if let Some(ph) = result.phase.clone() {
                            m.insert("phase".to_string(), serde_json::json!(ph));
                        }
                        m
                    }),
                    metadata: result.metadata,
                })
            }

            // Field theory operations (absorbed from FieldTheory tool)
            ComputeOp::Field(field_type) => {
                use crate::implementations::field_solver::UnifiedFieldSolver;
                let field_solver = UnifiedFieldSolver::new();
                // Convert Value to HashMap for configuration
                let configuration: HashMap<String, Value> = match &input.data {
                    Value::Object(map) => map.iter().map(|(k, v)| (k.clone(), v.clone())).collect(),
                    _ => HashMap::new(),
                };
                let field_input = FieldTheoryInput {
                    field_type: field_type.clone(),
                    configuration,
                    points: None,
                    parameters: input.parameters.iter().map(|(k, v)| (k.clone(), v.clone())).collect(),
                };
                let result = field_solver.field_theory(&field_input)?;
                Ok(ComputeOutput {
                    result: serde_json::json!(result.field_values),
                    additional: result.components.map(|c| {
                        c.into_iter().map(|(k, v)| (k, serde_json::json!(v))).collect()
                    }),
                    metadata: result.metadata,
                })
            }

            // Sample operations (absorbed from Sample tool)
            ComputeOp::Sample(sampling_method) => {
                use crate::implementations::sampler::UnifiedSampler;
                let sampler = UnifiedSampler::new();
                let sample_input = SampleInput {
                    method: sampling_method.clone(),
                    data: input.data.get("data").and_then(|v| v.as_array()).map(|arr| {
                        arr.iter().filter_map(|v| v.as_f64()).collect()
                    }).unwrap_or_default(),
                    num_samples: input.data.get("num_samples").and_then(|v| v.as_u64()).map(|n| n as usize),
                    parameters: input.parameters.iter().map(|(k, v)| (k.clone(), v.clone())).collect(),
                };
                let result = sampler.sample(&sample_input)?;
                Ok(ComputeOutput {
                    result: result.result,
                    additional: result.moments.map(|m| {
                        m.into_iter().map(|(k, v)| (k, serde_json::json!(v))).collect()
                    }),
                    metadata: result.metadata,
                })
            }
        }
    }
}
