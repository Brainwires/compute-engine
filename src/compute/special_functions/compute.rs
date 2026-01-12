//! Special functions computation
//!
//! Computes Bessel, Gamma, Error, Elliptic, Orthogonal polynomials, Airy, Hypergeometric functions

use crate::engine::*;
use std::collections::HashMap;

use super::{
    airy_function, AiryRequest,
    bessel_function, BesselRequest,
    elliptic_integral, EllipticIntegralRequest,
    error_function, ErrorFunctionRequest,
    gamma_function, GammaRequest,
    orthogonal_polynomial, PolynomialRequest,
};

/// Compute special functions
pub fn compute_special_function(func: &SpecialFunction, input: &ComputeInput) -> ToolResult<ComputeOutput> {

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

            let result = bessel_function(BesselRequest {
                function_type: function_type.to_string(),
                order,
                x,
            })
            .map_err(|e| e.to_string())?;

            Ok(ComputeOutput {
                result: serde_json::json!(result.value),
                additional: Some({
                    let mut m = HashMap::new();
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

            let result = gamma_function(GammaRequest {
                x,
                function: function.to_string(),
                y,
            })
            .map_err(|e| e.to_string())?;

            Ok(ComputeOutput {
                result: serde_json::json!(result.value),
                additional: Some({
                    let mut m = HashMap::new();
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
                error_function(ErrorFunctionRequest {
                    x,
                    function: function.to_string(),
                })
                .map_err(|e| e.to_string())?;

            Ok(ComputeOutput {
                result: serde_json::json!(result.value),
                additional: Some({
                    let mut m = HashMap::new();
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

            let result = elliptic_integral(
                EllipticIntegralRequest {
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
                    let mut m = HashMap::new();
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

            let result = orthogonal_polynomial(
                PolynomialRequest {
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
                    let mut m = HashMap::new();
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

            let result = airy_function(AiryRequest {
                function_type: function_type.to_string(),
                x,
            })
            .map_err(|e| e.to_string())?;

            let mut additional = HashMap::new();
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
