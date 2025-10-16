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

#[cfg(test)]
mod tests {
    use super::*;

    // ==================== Gamma Function Tests ====================

    #[test]
    fn test_gamma_integer() {
        let result = gamma_function(GammaRequest {
            x: 5.0,
            function: "gamma".to_string(),
            y: None,
        })
        .unwrap();
        // Γ(5) = 4! = 24
        assert!((result.value - 24.0).abs() < 0.01);
    }

    #[test]
    fn test_gamma_half() {
        let result = gamma_function(GammaRequest {
            x: 0.5,
            function: "gamma".to_string(),
            y: None,
        })
        .unwrap();
        // Γ(0.5) = √π ≈ 1.77245
        assert!((result.value - 1.77245).abs() < 0.01);
    }

    #[test]
    fn test_gamma_one() {
        let result = gamma_function(GammaRequest {
            x: 1.0,
            function: "gamma".to_string(),
            y: None,
        })
        .unwrap();
        // Γ(1) = 0! = 1
        assert!((result.value - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_gamma_two() {
        let result = gamma_function(GammaRequest {
            x: 2.0,
            function: "gamma".to_string(),
            y: None,
        })
        .unwrap();
        // Γ(2) = 1! = 1
        assert!((result.value - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_gamma_three() {
        let result = gamma_function(GammaRequest {
            x: 3.0,
            function: "gamma".to_string(),
            y: None,
        })
        .unwrap();
        // Γ(3) = 2! = 2
        assert!((result.value - 2.0).abs() < 0.01);
    }

    #[test]
    fn test_gamma_seven_halves() {
        let result = gamma_function(GammaRequest {
            x: 3.5,
            function: "gamma".to_string(),
            y: None,
        })
        .unwrap();
        // Γ(3.5) = 2.5 * 1.5 * 0.5 * Γ(0.5) ≈ 3.32335
        assert!((result.value - 3.32335).abs() < 0.01);
    }

    #[test]
    fn test_log_gamma() {
        let result = gamma_function(GammaRequest {
            x: 5.0,
            function: "log_gamma".to_string(),
            y: None,
        })
        .unwrap();
        // ln(Γ(5)) = ln(24) ≈ 3.17805
        assert!((result.value - 3.17805).abs() < 0.01);
    }

    #[test]
    fn test_digamma() {
        let result = gamma_function(GammaRequest {
            x: 1.0,
            function: "digamma".to_string(),
            y: None,
        })
        .unwrap();
        // ψ(1) ≈ -0.5772 (negative Euler-Mascheroni constant)
        assert!((result.value - (-0.5772)).abs() < 0.1);
    }

    #[test]
    fn test_beta_function() {
        let result = gamma_function(GammaRequest {
            x: 2.0,
            function: "beta".to_string(),
            y: Some(3.0),
        })
        .unwrap();
        // B(2,3) = Γ(2)Γ(3)/Γ(5) = 1*2/24 = 1/12 ≈ 0.0833
        assert!((result.value - 0.0833).abs() < 0.01);
    }

    #[test]
    fn test_beta_symmetric() {
        let result1 = gamma_function(GammaRequest {
            x: 3.0,
            function: "beta".to_string(),
            y: Some(5.0),
        })
        .unwrap();
        let result2 = gamma_function(GammaRequest {
            x: 5.0,
            function: "beta".to_string(),
            y: Some(3.0),
        })
        .unwrap();
        // B(3,5) = B(5,3) (symmetry)
        assert!((result1.value - result2.value).abs() < 1e-6);
    }

    #[test]
    fn test_beta_halves() {
        let result = gamma_function(GammaRequest {
            x: 0.5,
            function: "beta".to_string(),
            y: Some(0.5),
        })
        .unwrap();
        // B(0.5,0.5) = π
        assert!((result.value - PI).abs() < 0.01);
    }

    // ==================== Additional Comprehensive Tests ====================

    #[test]
    fn test_gamma_reflection_formula() {
        // Γ(z)Γ(1-z) = π/sin(πz)
        let z = 0.3;
        let gamma_z = gamma_function(GammaRequest {
            x: z,
            function: "gamma".to_string(),
            y: None,
        })
        .unwrap();
        let gamma_1_minus_z = gamma_function(GammaRequest {
            x: 1.0 - z,
            function: "gamma".to_string(),
            y: None,
        })
        .unwrap();
        let product = gamma_z.value * gamma_1_minus_z.value;
        let expected = PI / (PI * z).sin();
        assert!((product - expected).abs() < 0.1);
    }

    #[test]
    fn test_beta_integral_relation() {
        // B(a,b) = ∫[0 to 1] t^(a-1) * (1-t)^(b-1) dt
        // For a=2, b=2: B(2,2) = 1/6
        let result = gamma_function(GammaRequest {
            x: 2.0,
            function: "beta".to_string(),
            y: Some(2.0),
        })
        .unwrap();
        assert!((result.value - 1.0 / 6.0).abs() < 0.01);
    }

    #[test]
    fn test_gamma_factorial_sequence() {
        // Γ(n+1) = n! for positive integers
        fn factorial(n: usize) -> f64 {
            (1..=n).map(|i| i as f64).product()
        }

        for n in 1..8 {
            let gamma_n_plus_1 = gamma_function(GammaRequest {
                x: (n + 1) as f64,
                function: "gamma".to_string(),
                y: None,
            })
            .unwrap();
            let factorial_n = factorial(n);
            assert!((gamma_n_plus_1.value - factorial_n).abs() < 0.01);
        }
    }
}
