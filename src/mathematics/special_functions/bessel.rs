//! Bessel Functions
//!
//! Implements Bessel functions of the first and second kind (J, Y),
//! and modified Bessel functions (I, K).

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// Request/Response types

#[derive(Debug, Deserialize)]
pub struct BesselRequest {
    pub function_type: String, // "J", "Y", "I", "K"
    pub order: f64,
    pub x: f64,
}

#[derive(Debug, Serialize)]
pub struct BesselResult {
    pub value: f64,
    pub function_type: String,
    pub order: f64,
}

// Factorial helper
fn factorial(n: usize) -> f64 {
    (1..=n).map(|i| i as f64).product()
}

/// Main Bessel function dispatcher
pub fn bessel_function(request: BesselRequest) -> Result<BesselResult, String> {
    let value = match request.function_type.as_str() {
        "J" => bessel_j(request.order, request.x)?,
        "Y" => bessel_y(request.order, request.x)?,
        "I" => bessel_i(request.order, request.x)?,
        "K" => bessel_k(request.order, request.x)?,
        _ => {
            return Err(format!(
                "Unknown Bessel function: {}",
                request.function_type
            ));
        }
    };

    Ok(BesselResult {
        value,
        function_type: request.function_type,
        order: request.order,
    })
}

fn bessel_j(nu: f64, x: f64) -> Result<f64, String> {
    // Bessel J using series expansion (simplified)
    if nu.fract() != 0.0 && nu != nu.floor() {
        return Err("Only integer orders supported in this implementation".to_string());
    }

    let n = nu as usize;
    let mut sum = 0.0;

    for k in 0..50_usize {
        let term = (-1.0_f64).powi(k as i32) * (x / 2.0).powi((2 * k + n) as i32)
            / (factorial(k) * factorial(k + n));
        sum += term;
        if term.abs() < 1e-10 {
            break;
        }
    }

    Ok(sum)
}

fn bessel_y(nu: f64, x: f64) -> Result<f64, String> {
    // Bessel Y (Neumann function) using relation to J
    // Y_n(x) = (J_n(x) * cos(nπ) - J_{-n}(x)) / sin(nπ)
    // For integer n, use limit formula: Y_n(x) = lim_{ν→n} Y_ν(x)

    if x <= 0.0 {
        return Err("Bessel Y requires x > 0".to_string());
    }

    if nu.fract() != 0.0 {
        return Err("Only integer orders supported in this implementation".to_string());
    }

    let n = nu as i32;

    // For integer orders, use series expansion:
    // Y_n(x) = (2/π) * [ln(x/2) * J_n(x) + sum_k - sum_j]
    // where sum_k involves psi function terms

    let j_n = bessel_j(nu, x)?;
    let log_term = (x / 2.0).ln();

    // Series approximation for small x
    if x < 2.0 {
        // Use series expansion
        let mut sum = 0.0;

        // Leading term with logarithm
        let result = (2.0 / PI) * (log_term + 0.5772156649) * j_n; // 0.5772... is Euler's gamma

        // Additional series terms for better accuracy
        let n_abs = n.abs() as usize;
        for k in 0..n_abs {
            let term = factorial(n_abs - k - 1) / factorial(k)
                * (2.0 / x).powi(n_abs as i32 - 2 * k as i32);
            sum += term;
        }

        let correction = -(1.0 / PI) * sum;
        Ok(result + correction)
    } else {
        // For larger x, use asymptotic approximation
        // Y_n(x) ≈ sqrt(2/(πx)) * sin(x - nπ/2 - π/4)
        let phase = x - (n as f64) * PI / 2.0 - PI / 4.0;
        Ok((2.0 / (PI * x)).sqrt() * phase.sin())
    }
}

fn bessel_i(nu: f64, x: f64) -> Result<f64, String> {
    // Modified Bessel I
    if nu.fract() != 0.0 {
        return Err("Only integer orders supported".to_string());
    }

    let n = nu as usize;
    let mut sum = 0.0;

    for k in 0..50_usize {
        let term = (x / 2.0).powi((2 * k + n) as i32) / (factorial(k) * factorial(k + n));
        sum += term;
        if term.abs() < 1e-10 {
            break;
        }
    }

    Ok(sum)
}

fn bessel_k(nu: f64, x: f64) -> Result<f64, String> {
    // Modified Bessel K (simplified)
    let i = bessel_i(nu, x)?;
    Ok(PI / (2.0 * x * i)) // Simplified approximation
}

