//! Gamma Function Module
//!
//! Implements the gamma function and related functions:
//! - Gamma function (Γ)
//! - Log-gamma function (ln Γ)
//! - Digamma function (ψ, derivative of log-gamma)
//! - Beta function (B)

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// Request/Response types

#[derive(Debug, Deserialize)]
pub struct GammaRequest {
    pub x: f64,
    pub function: String, // "gamma", "log_gamma", "digamma", "beta"
    pub y: Option<f64>,   // For beta function
}

#[derive(Debug, Serialize)]
pub struct GammaResult {
    pub value: f64,
    pub function: String,
}

/// Main entry point for gamma function operations
pub fn gamma_function(request: GammaRequest) -> Result<GammaResult, String> {
    let value = match request.function.as_str() {
        "gamma" => gamma(request.x)?,
        "log_gamma" | "lgamma" => log_gamma(request.x)?,
        "digamma" => digamma(request.x)?,
        "beta" => {
            let y = request.y.ok_or("Beta function requires two arguments")?;
            beta(request.x, y)?
        }
        _ => return Err(format!("Unknown gamma function: {}", request.function)),
    };

    Ok(GammaResult {
        value,
        function: request.function,
    })
}

/// Compute gamma function using Lanczos approximation
pub fn gamma(z: f64) -> Result<f64, String> {
    // Lanczos approximation
    if z < 0.5 {
        // Reflection formula: Γ(z)Γ(1-z) = π/sin(πz)
        return Ok(PI / ((PI * z).sin() * gamma(1.0 - z)?));
    }

    const G: f64 = 7.0;
    const COEF: [f64; 9] = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];

    let z = z - 1.0;
    let mut x = COEF[0];

    for i in 1..9 {
        x += COEF[i] / (z + i as f64);
    }

    let t = z + G + 0.5;
    Ok((2.0 * PI).sqrt() * t.powf(z + 0.5) * (-t).exp() * x)
}

/// Compute natural logarithm of gamma function
pub fn log_gamma(z: f64) -> Result<f64, String> {
    Ok(gamma(z)?.ln())
}

/// Compute digamma function (psi function, derivative of log-gamma)
pub fn digamma(z: f64) -> Result<f64, String> {
    // Asymptotic series approximation
    let mut result = z.ln();
    result -= 1.0 / (2.0 * z);
    result -= 1.0 / (12.0 * z * z);
    Ok(result)
}

/// Compute beta function: B(a,b) = Γ(a)Γ(b)/Γ(a+b)
pub fn beta(a: f64, b: f64) -> Result<f64, String> {
    Ok(gamma(a)? * gamma(b)? / gamma(a + b)?)
}

