//! Orthogonal Polynomials
//!
//! Implements classical orthogonal polynomials:
//! - Legendre polynomials P_n(x)
//! - Hermite polynomials H_n(x) (physicist's convention)
//! - Laguerre polynomials L_n(x)
//! - Chebyshev polynomials T_n(x) (first kind)

use serde::{Deserialize, Serialize};

/// Request structure for orthogonal polynomial computations
#[derive(Debug, Deserialize)]
pub struct PolynomialRequest {
    pub polynomial_type: String, // "legendre", "hermite", "laguerre", "chebyshev"
    pub n: usize, // degree
    pub x: f64,
    pub alpha: Option<f64>, // For generalized polynomials
}

/// Result structure for orthogonal polynomial computations
#[derive(Debug, Serialize)]
pub struct PolynomialResult {
    pub value: f64,
    pub polynomial_type: String,
    pub degree: usize,
}

/// Compute orthogonal polynomial based on the request type
///
/// # Arguments
/// * `request` - PolynomialRequest specifying the polynomial type, degree, and evaluation point
///
/// # Returns
/// * `Result<PolynomialResult, String>` - The computed value or error message
///
/// # Examples
/// ```
/// use computational_engine::compute::special_functions::{
///     orthogonal_polynomial, PolynomialRequest
/// };
///
/// let request = PolynomialRequest {
///     polynomial_type: "legendre".to_string(),
///     n: 2,
///     x: 0.5,
///     alpha: None,
/// };
/// let result = orthogonal_polynomial(request).unwrap();
/// // P_2(0.5) = (3*0.5^2 - 1)/2 = -0.125
/// assert!((result.value - (-0.125)).abs() < 1e-6);
/// ```
pub fn orthogonal_polynomial(request: PolynomialRequest) -> Result<PolynomialResult, String> {
    let value = match request.polynomial_type.as_str() {
        "legendre" => legendre_p(request.n, request.x)?,
        "hermite" => hermite_h(request.n, request.x)?,
        "laguerre" => laguerre_l(request.n, request.x)?,
        "chebyshev" | "chebyshev_t" => chebyshev_t(request.n, request.x)?,
        _ => return Err(format!("Unknown polynomial: {}", request.polynomial_type)),
    };

    Ok(PolynomialResult {
        value,
        polynomial_type: request.polynomial_type,
        degree: request.n,
    })
}

/// Legendre polynomial P_n(x)
///
/// Computes the Legendre polynomial of degree n at point x using
/// the three-term recurrence relation:
/// P_0(x) = 1
/// P_1(x) = x
/// P_{n+1}(x) = [(2n+1)xP_n(x) - nP_{n-1}(x)] / (n+1)
///
/// # Arguments
/// * `n` - Degree of the polynomial
/// * `x` - Evaluation point
///
/// # Returns
/// * `Result<f64, String>` - The value of P_n(x)
pub fn legendre_p(n: usize, x: f64) -> Result<f64, String> {
    // Legendre polynomial using recurrence
    if n == 0 {
        return Ok(1.0);
    }
    if n == 1 {
        return Ok(x);
    }

    let mut p0 = 1.0;
    let mut p1 = x;

    for k in 2..=n {
        let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
        p0 = p1;
        p1 = p2;
    }

    Ok(p1)
}

/// Hermite polynomial H_n(x) (physicist's convention)
///
/// Computes the Hermite polynomial of degree n at point x using
/// the three-term recurrence relation:
/// H_0(x) = 1
/// H_1(x) = 2x
/// H_{n+1}(x) = 2xH_n(x) - 2nH_{n-1}(x)
///
/// # Arguments
/// * `n` - Degree of the polynomial
/// * `x` - Evaluation point
///
/// # Returns
/// * `Result<f64, String>` - The value of H_n(x)
pub fn hermite_h(n: usize, x: f64) -> Result<f64, String> {
    // Hermite polynomial (physicist's)
    if n == 0 {
        return Ok(1.0);
    }
    if n == 1 {
        return Ok(2.0 * x);
    }

    let mut h0 = 1.0;
    let mut h1 = 2.0 * x;

    for k in 2..=n {
        let h2 = 2.0 * x * h1 - 2.0 * (k - 1) as f64 * h0;
        h0 = h1;
        h1 = h2;
    }

    Ok(h1)
}

/// Laguerre polynomial L_n(x)
///
/// Computes the Laguerre polynomial of degree n at point x using
/// the three-term recurrence relation:
/// L_0(x) = 1
/// L_1(x) = 1 - x
/// L_{n+1}(x) = [(2n+1-x)L_n(x) - nL_{n-1}(x)] / (n+1)
///
/// # Arguments
/// * `n` - Degree of the polynomial
/// * `x` - Evaluation point
///
/// # Returns
/// * `Result<f64, String>` - The value of L_n(x)
pub fn laguerre_l(n: usize, x: f64) -> Result<f64, String> {
    // Laguerre polynomial using recurrence relation:
    // L_{k+1}(x) = [(2k+1-x)L_k(x) - k*L_{k-1}(x)] / (k+1)
    if n == 0 {
        return Ok(1.0);
    }
    if n == 1 {
        return Ok(1.0 - x);
    }

    let mut l0 = 1.0;
    let mut l1 = 1.0 - x;

    for k in 1..n {
        let k_f64 = k as f64;
        let l2 = ((2.0 * k_f64 + 1.0 - x) * l1 - k_f64 * l0) / (k_f64 + 1.0);
        l0 = l1;
        l1 = l2;
    }

    Ok(l1)
}

/// Chebyshev polynomial T_n(x) (first kind)
///
/// Computes the Chebyshev polynomial of the first kind of degree n at point x.
/// Uses the three-term recurrence relation:
/// T_0(x) = 1
/// T_1(x) = x
/// T_{n+1}(x) = 2xT_n(x) - T_{n-1}(x)
///
/// These polynomials satisfy T_n(cos(θ)) = cos(nθ).
///
/// # Arguments
/// * `n` - Degree of the polynomial
/// * `x` - Evaluation point
///
/// # Returns
/// * `Result<f64, String>` - The value of T_n(x)
pub fn chebyshev_t(n: usize, x: f64) -> Result<f64, String> {
    // Chebyshev polynomial of first kind
    if n == 0 {
        return Ok(1.0);
    }
    if n == 1 {
        return Ok(x);
    }

    let mut t0 = 1.0;
    let mut t1 = x;

    for _ in 2..=n {
        let t2 = 2.0 * x * t1 - t0;
        t0 = t1;
        t1 = t2;
    }

    Ok(t1)
}

