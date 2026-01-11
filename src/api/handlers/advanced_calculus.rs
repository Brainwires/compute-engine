//! Advanced calculus operations handler

use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::json;

pub fn handle(request: &ComputationRequest) -> ComputationResponse {
    match request.operation.as_str() {
        "derivative" => handle_derivative(request),
        "integral" => handle_integral(request),
        "limit" => handle_limit(request),
        _ => {
            return ComputationResponse::error(
                request.module.clone(),
                request.operation.clone(),
                format!("Unknown operation: {}", request.operation),
            );
        }
    }
}

fn handle_derivative(request: &ComputationRequest) -> ComputationResponse {
    // Simple polynomial derivative: d/dx (a₀ + a₁x + a₂x² + ...) = a₁ + 2a₂x + 3a₃x² + ...
    if let Some(coeffs) = request.parameters.get("coefficients") {
        if let Some(coeffs_arr) = coeffs.as_array() {
            let coefficients: Vec<f64> = coeffs_arr.iter().filter_map(|v| v.as_f64()).collect();

            // Compute derivative coefficients
            let mut deriv_coeffs = Vec::new();
            for (i, &coeff) in coefficients.iter().enumerate().skip(1) {
                deriv_coeffs.push(coeff * i as f64);
            }

            // If x is provided, evaluate at that point
            if let Some(x_val) = request.parameters.get("x").and_then(|v| v.as_f64()) {
                let mut result = 0.0;
                for (i, &coeff) in deriv_coeffs.iter().enumerate() {
                    result += coeff * x_val.powi(i as i32);
                }

                return ComputationResponse::success(
                    request.module.clone(),
                    request.operation.clone(),
                    json!({
                        "derivative": result,
                        "x": x_val,
                        "derivative_coefficients": deriv_coeffs
                    }),
                );
            }

            return ComputationResponse::success(
                request.module.clone(),
                request.operation.clone(),
                json!({ "derivative_coefficients": deriv_coeffs }),
            );
        }
    }

    ComputationResponse::error(
        request.module.clone(),
        request.operation.clone(),
        "Invalid parameters for derivative".to_string(),
    )
}

fn handle_integral(request: &ComputationRequest) -> ComputationResponse {
    // Simple polynomial integral: ∫(a₀ + a₁x + a₂x² + ...)dx = a₀x + a₁x²/2 + a₂x³/3 + ... + C
    if let Some(coeffs) = request.parameters.get("coefficients") {
        if let Some(coeffs_arr) = coeffs.as_array() {
            let coefficients: Vec<f64> = coeffs_arr.iter().filter_map(|v| v.as_f64()).collect();

            // Compute integral coefficients (indefinite integral, C=0)
            let mut integral_coeffs = vec![0.0]; // Constant term (C)
            for (i, &coeff) in coefficients.iter().enumerate() {
                integral_coeffs.push(coeff / (i + 1) as f64);
            }

            // If bounds are provided, compute definite integral
            if let (Some(a), Some(b)) = (
                request.parameters.get("lower").and_then(|v| v.as_f64()),
                request.parameters.get("upper").and_then(|v| v.as_f64()),
            ) {
                let eval_at = |x: f64| -> f64 {
                    let mut sum = 0.0;
                    for (i, &coeff) in integral_coeffs.iter().enumerate() {
                        sum += coeff * x.powi(i as i32);
                    }
                    sum
                };

                let result = eval_at(b) - eval_at(a);

                return ComputationResponse::success(
                    request.module.clone(),
                    request.operation.clone(),
                    json!({
                        "integral": result,
                        "lower": a,
                        "upper": b,
                        "antiderivative_coefficients": integral_coeffs
                    }),
                );
            }

            return ComputationResponse::success(
                request.module.clone(),
                request.operation.clone(),
                json!({ "antiderivative_coefficients": integral_coeffs }),
            );
        }
    }

    ComputationResponse::error(
        request.module.clone(),
        request.operation.clone(),
        "Invalid parameters for integral".to_string(),
    )
}

fn handle_limit(request: &ComputationRequest) -> ComputationResponse {
    // Handle rational functions: P(x)/Q(x)
    if let Some(function_type) = request.parameters.get("function_type").and_then(|v| v.as_str()) {
        if function_type == "rational" {
            if let (Some(num_arr), Some(den_arr), Some(point)) = (
                request.parameters.get("numerator").and_then(|v| v.as_array()),
                request.parameters.get("denominator").and_then(|v| v.as_array()),
                request.parameters.get("point").and_then(|v| v.as_f64()),
            ) {
                let numerator: Vec<f64> = num_arr.iter().filter_map(|v| v.as_f64()).collect();
                let denominator: Vec<f64> = den_arr.iter().filter_map(|v| v.as_f64()).collect();

                let eval_poly = |coeffs: &[f64], x: f64| -> f64 {
                    coeffs.iter().enumerate().fold(0.0, |sum, (i, &coeff)| sum + coeff * x.powi(i as i32))
                };

                let num_val = eval_poly(&numerator, point);
                let den_val = eval_poly(&denominator, point);

                if den_val.abs() < 1e-10 {
                    // Indeterminate form - use L'Hôpital's rule or return infinity
                    if num_val.abs() < 1e-10 {
                        // 0/0 - apply L'Hôpital's rule (derivative of num / derivative of den)
                        let num_deriv: Vec<f64> = numerator.iter().enumerate().skip(1)
                            .map(|(i, &c)| c * i as f64).collect();
                        let den_deriv: Vec<f64> = denominator.iter().enumerate().skip(1)
                            .map(|(i, &c)| c * i as f64).collect();

                        let num_deriv_val = eval_poly(&num_deriv, point);
                        let den_deriv_val = eval_poly(&den_deriv, point);

                        if den_deriv_val.abs() > 1e-10 {
                            let result = num_deriv_val / den_deriv_val;
                            return ComputationResponse::success(
                                request.module.clone(),
                                request.operation.clone(),
                                json!({
                                    "limit": result,
                                    "point": point,
                                    "method": "lhopital"
                                }),
                            );
                        }
                    }
                    return ComputationResponse::success(
                        request.module.clone(),
                        request.operation.clone(),
                        json!({
                            "limit": if num_val > 0.0 { f64::INFINITY } else { f64::NEG_INFINITY },
                            "point": point,
                            "infinite": true
                        }),
                    );
                }

                let result = num_val / den_val;
                return ComputationResponse::success(
                    request.module.clone(),
                    request.operation.clone(),
                    json!({
                        "limit": result,
                        "point": point
                    }),
                );
            }
        }
    }

    // Fallback for polynomial limits
    if let Some(coeffs) = request.parameters.get("coefficients") {
        if let Some(coeffs_arr) = coeffs.as_array() {
            let coefficients: Vec<f64> = coeffs_arr.iter().filter_map(|v| v.as_f64()).collect();

            if let Some(x_val) = request.parameters.get("x").and_then(|v| v.as_f64()) {
                let mut result = 0.0;
                for (i, &coeff) in coefficients.iter().enumerate() {
                    result += coeff * x_val.powi(i as i32);
                }

                return ComputationResponse::success(
                    request.module.clone(),
                    request.operation.clone(),
                    json!({
                        "limit": result,
                        "x": x_val
                    }),
                );
            }
        }
    }

    ComputationResponse::error(
        request.module.clone(),
        request.operation.clone(),
        "Invalid parameters for limit".to_string(),
    )
}

#[cfg(test)]
#[path = "../../../tests/unit/api/handlers/advanced_calculus_handler_tests.rs"]
mod tests;
