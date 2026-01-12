//! Special Functions Module
//!
//! Implements special mathematical functions commonly used in physics and engineering:
//! - Bessel functions (J, Y, I, K)
//! - Gamma and Beta functions
//! - Error functions (erf, erfc, erfcx)
//! - Elliptic integrals
//! - Orthogonal polynomials (Legendre, Hermite, Laguerre, Chebyshev)
//! - Airy functions
//! - Hypergeometric functions (planned)

mod airy;
mod bessel;
mod elliptic;
mod error;
mod gamma;
mod polynomials;
mod compute;

pub use airy::{airy_function, AiryRequest, AiryResult};
pub use bessel::{bessel_function, BesselRequest, BesselResult};
pub use elliptic::{
    elliptic_e, elliptic_f, elliptic_integral, elliptic_k, elliptic_pi, EllipticIntegralRequest,
    EllipticIntegralResult,
};
pub use error::{error_function, ErrorFunctionRequest, ErrorFunctionResult};
pub use gamma::{gamma_function, GammaRequest, GammaResult};
pub use polynomials::{
    chebyshev_t, hermite_h, laguerre_l, legendre_p, orthogonal_polynomial, PolynomialRequest,
    PolynomialResult,
};
pub use compute::compute_special_function;

use serde::{Deserialize, Serialize};

// ==================== Hypergeometric Functions ====================
// TODO: Extract to separate module

#[derive(Debug, Deserialize)]
pub struct HypergeometricRequest {
    pub function_type: String, // "1F1", "2F1", "0F1"
    pub a: Vec<f64>,
    pub b: Vec<f64>,
    pub z: f64,
}

#[derive(Debug, Serialize)]
pub struct HypergeometricResult {
    pub value: f64,
    pub function_type: String,
}

// Test module - part of this module, can access private functions
#[cfg(test)]
#[path = "../../../tests/unit/compute/special_functions_tests.rs"]
mod tests;

// Advanced calculus tests (comprehensive unit tests for all special functions)
#[cfg(test)]
#[path = "../../../tests/unit/compute/calculus/advanced_calculus_tests.rs"]
mod advanced_calculus_tests;
