//! Elliptic Integrals
//!
//! Implements elliptic integrals commonly used in physics and engineering:
//! - Complete elliptic integral of the first kind (K)
//! - Complete elliptic integral of the second kind (E)
//! - Incomplete elliptic integral of the first kind (F)
//! - Incomplete elliptic integral of the third kind (Pi)

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Request structure for elliptic integral computations
#[derive(Debug, Deserialize)]
pub struct EllipticIntegralRequest {
    pub integral_type: String, // "K", "E", "F", "Pi"
    pub k: f64, // modulus
    pub phi: Option<f64>, // For incomplete integrals
    pub n: Option<f64>, // For Pi integral
}

/// Result structure for elliptic integral computations
#[derive(Debug, Serialize)]
pub struct EllipticIntegralResult {
    pub value: f64,
    pub integral_type: String,
}

/// Compute elliptic integral based on the request type
///
/// # Arguments
/// * `request` - EllipticIntegralRequest specifying the integral type and parameters
///
/// # Returns
/// * `Result<EllipticIntegralResult, String>` - The computed value or error message
///
/// # Examples
/// ```
/// use computational_engine::compute::special_functions::{
///     elliptic_integral, EllipticIntegralRequest
/// };
///
/// let request = EllipticIntegralRequest {
///     integral_type: "K".to_string(),
///     k: 0.5,
///     phi: None,
///     n: None,
/// };
/// let result = elliptic_integral(request).unwrap();
/// assert!((result.value - 1.685750).abs() < 0.01);
/// ```
pub fn elliptic_integral(request: EllipticIntegralRequest) -> Result<EllipticIntegralResult, String> {
    let value = match request.integral_type.as_str() {
        "K" => elliptic_k(request.k)?,
        "E" => elliptic_e(request.k)?,
        "F" => {
            let phi = request.phi.ok_or("Incomplete elliptic integral F requires phi parameter")?;
            elliptic_f(phi, request.k)?
        },
        "Pi" => {
            let phi = request.phi.ok_or("Elliptic integral Pi requires phi parameter")?;
            let n = request.n.ok_or("Elliptic integral Pi requires n parameter")?;
            elliptic_pi(phi, n, request.k)?
        },
        _ => return Err(format!("Elliptic integral {} not implemented", request.integral_type)),
    };

    Ok(EllipticIntegralResult {
        value,
        integral_type: request.integral_type,
    })
}

/// Complete elliptic integral of the first kind K(k)
///
/// Uses the Arithmetic-Geometric Mean (AGM) algorithm for efficient computation.
///
/// # Arguments
/// * `k` - Elliptic modulus (must satisfy |k| < 1)
///
/// # Returns
/// * `Result<f64, String>` - The value of K(k) or error if |k| >= 1
pub fn elliptic_k(k: f64) -> Result<f64, String> {
    // Complete elliptic integral of first kind (AGM algorithm)
    if k.abs() >= 1.0 {
        return Err("Modulus must be |k| < 1".to_string());
    }

    let mut a = 1.0;
    let mut g = (1.0 - k * k).sqrt();

    for _ in 0..10 {
        let a_new = (a + g) / 2.0;
        let g_new = (a * g).sqrt();
        if (a - g).abs() < 1e-10 {
            break;
        }
        a = a_new;
        g = g_new;
    }

    Ok(PI / (2.0 * a))
}

/// Complete elliptic integral of the second kind E(k)
///
/// # Arguments
/// * `k` - Elliptic modulus (must satisfy |k| < 1)
///
/// # Returns
/// * `Result<f64, String>` - The value of E(k) or error if |k| >= 1
pub fn elliptic_e(k: f64) -> Result<f64, String> {
    // Complete elliptic integral of second kind
    if k.abs() >= 1.0 {
        return Err("Modulus must be |k| < 1".to_string());
    }

    // Approximate using polynomial
    let m = k * k;
    Ok((1.0 - m / 4.0 - 3.0 * m * m / 64.0) * elliptic_k(k)?)
}

/// Incomplete elliptic integral of the first kind F(φ, k)
///
/// Computes F(φ,k) = ∫[0 to φ] dθ / sqrt(1 - k²sin²θ)
/// using numerical integration with Simpson's rule.
///
/// # Arguments
/// * `phi` - Upper limit of integration
/// * `k` - Elliptic modulus (must satisfy |k| < 1)
///
/// # Returns
/// * `Result<f64, String>` - The value of F(φ,k) or error if |k| >= 1
pub fn elliptic_f(phi: f64, k: f64) -> Result<f64, String> {
    // Incomplete elliptic integral of the first kind
    // F(φ,k) = ∫[0 to φ] dθ / sqrt(1 - k²sin²θ)

    if k.abs() >= 1.0 {
        return Err("Modulus must be |k| < 1".to_string());
    }

    // Use numerical integration with adaptive Simpson's rule
    let n = 1000; // number of integration steps
    let h = phi / n as f64;
    let k2 = k * k;

    let mut sum = 0.0;

    // Simpson's 1/3 rule
    for i in 0..=n {
        let theta = i as f64 * h;
        let sin_theta = theta.sin();
        let integrand = 1.0 / (1.0 - k2 * sin_theta * sin_theta).sqrt();

        let coefficient = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };

        sum += coefficient * integrand;
    }

    Ok(sum * h / 3.0)
}

/// Incomplete elliptic integral of the third kind Π(φ, n, k)
///
/// Computes Π(φ,n,k) = ∫[0 to φ] dθ / [(1 - n·sin²θ) · sqrt(1 - k²sin²θ)]
/// using numerical integration with Simpson's rule.
///
/// # Arguments
/// * `phi` - Upper limit of integration
/// * `n` - Characteristic parameter (must satisfy n < 1)
/// * `k` - Elliptic modulus (must satisfy |k| < 1)
///
/// # Returns
/// * `Result<f64, String>` - The value of Π(φ,n,k) or error if constraints violated
pub fn elliptic_pi(phi: f64, n: f64, k: f64) -> Result<f64, String> {
    // Incomplete elliptic integral of the third kind
    // Π(φ,n,k) = ∫[0 to φ] dθ / [(1 - n*sin²θ) * sqrt(1 - k²sin²θ)]

    if k.abs() >= 1.0 {
        return Err("Modulus must be |k| < 1".to_string());
    }

    if n >= 1.0 {
        return Err("Parameter n must be n < 1 for convergence".to_string());
    }

    // Use numerical integration
    let num_steps = 1000;
    let h = phi / num_steps as f64;
    let k2 = k * k;

    let mut sum = 0.0;

    // Simpson's 1/3 rule
    for i in 0..=num_steps {
        let theta = i as f64 * h;
        let sin_theta = theta.sin();
        let sin2_theta = sin_theta * sin_theta;

        let denominator_part1 = 1.0 - n * sin2_theta;
        let denominator_part2 = (1.0 - k2 * sin2_theta).sqrt();

        if denominator_part1.abs() < 1e-10 || denominator_part2.abs() < 1e-10 {
            return Err("Singular point encountered in integration".to_string());
        }

        let integrand = 1.0 / (denominator_part1 * denominator_part2);

        let coefficient = if i == 0 || i == num_steps {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };

        sum += coefficient * integrand;
    }

    Ok(sum * h / 3.0)
}

