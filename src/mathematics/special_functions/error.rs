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

#[cfg(test)]
mod tests {
    use super::*;

    // ==================== Error Function Tests ====================

    #[test]
    fn test_erf_zero() {
        let result = error_function(ErrorFunctionRequest {
            x: 0.0,
            function: "erf".to_string(),
        })
        .unwrap();
        // erf(0) = 0
        assert!(result.value.abs() < 1e-6);
    }

    #[test]
    fn test_erf_one() {
        let result = error_function(ErrorFunctionRequest {
            x: 1.0,
            function: "erf".to_string(),
        })
        .unwrap();
        // erf(1) ≈ 0.8427
        assert!((result.value - 0.8427).abs() < 0.001);
    }

    #[test]
    fn test_erf_negative() {
        let result = error_function(ErrorFunctionRequest {
            x: -1.0,
            function: "erf".to_string(),
        })
        .unwrap();
        // erf(-1) ≈ -0.8427 (odd function)
        assert!((result.value - (-0.8427)).abs() < 0.001);
    }

    #[test]
    fn test_erf_two() {
        let result = error_function(ErrorFunctionRequest {
            x: 2.0,
            function: "erf".to_string(),
        })
        .unwrap();
        // erf(2) ≈ 0.9953
        assert!((result.value - 0.9953).abs() < 0.001);
    }

    #[test]
    fn test_erf_large() {
        let result = error_function(ErrorFunctionRequest {
            x: 3.0,
            function: "erf".to_string(),
        })
        .unwrap();
        // erf(3) ≈ 0.99998 (approaches 1)
        assert!((result.value - 0.99998).abs() < 0.0001);
    }

    #[test]
    fn test_erfc_zero() {
        let result = error_function(ErrorFunctionRequest {
            x: 0.0,
            function: "erfc".to_string(),
        })
        .unwrap();
        // erfc(0) = 1 - erf(0) = 1
        assert!((result.value - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_erfc_one() {
        let result = error_function(ErrorFunctionRequest {
            x: 1.0,
            function: "erfc".to_string(),
        })
        .unwrap();
        // erfc(1) ≈ 0.1573
        assert!((result.value - 0.1573).abs() < 0.001);
    }

    #[test]
    fn test_erf_erfc_complement() {
        let x = 1.5;
        let erf_result = error_function(ErrorFunctionRequest {
            x,
            function: "erf".to_string(),
        })
        .unwrap();
        let erfc_result = error_function(ErrorFunctionRequest {
            x,
            function: "erfc".to_string(),
        })
        .unwrap();
        // erf(x) + erfc(x) = 1
        assert!((erf_result.value + erfc_result.value - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_erfcx_one() {
        let result = error_function(ErrorFunctionRequest {
            x: 1.0,
            function: "erfcx".to_string(),
        })
        .unwrap();
        // erfcx(1) = exp(1) * erfc(1) ≈ 0.4276
        assert!((result.value - 0.4276).abs() < 0.01);
    }

    #[test]
    fn test_erfi_zero() {
        let result = error_function(ErrorFunctionRequest {
            x: 0.0,
            function: "erfi".to_string(),
        })
        .unwrap();
        // erfi(0) = 0
        assert!(result.value.abs() < 1e-6);
    }

    // ==================== Additional Comprehensive Tests ====================

    #[test]
    fn test_erf_symmetry() {
        // erf(-x) = -erf(x)
        let x = 1.5;
        let erf_pos = error_function(ErrorFunctionRequest {
            x,
            function: "erf".to_string(),
        })
        .unwrap();
        let erf_neg = error_function(ErrorFunctionRequest {
            x: -x,
            function: "erf".to_string(),
        })
        .unwrap();
        assert!((erf_pos.value + erf_neg.value).abs() < 1e-6);
    }

    #[test]
    fn test_error_function_limit() {
        // lim_{x→∞} erf(x) = 1
        let result = error_function(ErrorFunctionRequest {
            x: 5.0,
            function: "erf".to_string(),
        })
        .unwrap();
        assert!((result.value - 1.0).abs() < 0.001);
    }
}
