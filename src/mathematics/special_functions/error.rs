//! Error Function Module
//!
//! Implements error functions and related functions:
//! - Error function (erf)
//! - Complementary error function (erfc)
//! - Scaled complementary error function (erfcx)
//! - Imaginary error function (erfi)

use serde::{Deserialize, Serialize};

// Request/Response types

#[derive(Debug, Deserialize)]
pub struct ErrorFunctionRequest {
    pub x: f64,
    pub function: String, // "erf", "erfc", "erfcx", "erfi"
}

#[derive(Debug, Serialize)]
pub struct ErrorFunctionResult {
    pub value: f64,
    pub function: String,
}

/// Main entry point for error function operations
pub fn error_function(request: ErrorFunctionRequest) -> Result<ErrorFunctionResult, String> {
    let value = match request.function.as_str() {
        "erf" => erf(request.x),
        "erfc" => erfc(request.x),
        "erfcx" => erfcx(request.x),
        "erfi" => erfi(request.x),
        _ => return Err(format!("Unknown error function: {}", request.function)),
    };

    Ok(ErrorFunctionResult {
        value,
        function: request.function,
    })
}

/// Compute error function using series expansion
pub fn erf(x: f64) -> f64 {
    // Error function using series expansion
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;
    let p = 0.3275911;

    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();

    let t = 1.0 / (1.0 + p * x);
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();

    sign * y
}

/// Compute complementary error function
pub fn erfc(x: f64) -> f64 {
    1.0 - erf(x)
}

/// Compute scaled complementary error function: erfcx(x) = exp(x²) * erfc(x)
pub fn erfcx(x: f64) -> f64 {
    (x * x).exp() * erfc(x)
}

/// Compute imaginary error function: erfi(x) = -i * erf(ix)
pub fn erfi(x: f64) -> f64 {
    -erf(-x)
}

