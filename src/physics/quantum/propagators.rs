//! QFT Propagators and Green's Functions

use num_complex::Complex64;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Fine structure constant
pub const ALPHA_EM: f64 = 1.0 / 137.036;

/// Scalar propagator (Klein-Gordon): i/(p² - m² + iε)
pub fn scalar_propagator(p_squared: f64, mass: f64) -> Complex64 {
    let mass_sq = mass * mass;
    let epsilon = 1e-10;
    let denominator = Complex64::new(p_squared - mass_sq, epsilon);
    Complex64::i() / denominator
}

/// Fermion propagator (Dirac): i(γ·p + m)/(p² - m² + iε)
pub fn fermion_propagator_trace(p_squared: f64, mass: f64) -> Complex64 {
    let mass_sq = mass * mass;
    let epsilon = 1e-10;
    let denominator = Complex64::new(p_squared - mass_sq, epsilon);
    Complex64::new(mass, 0.0) / denominator
}

/// Photon propagator (Feynman gauge): -i g_μν / q²
pub fn photon_propagator(q_squared: f64) -> Complex64 {
    let epsilon = 1e-10;
    let denominator = Complex64::new(q_squared, epsilon);
    -Complex64::i() / denominator
}

/// Massive vector boson propagator (W/Z): -i(g_μν - q_μq_ν/m²)/(q² - m² + iε)
pub fn vector_boson_propagator(q_squared: f64, mass: f64) -> Complex64 {
    let mass_sq = mass * mass;
    let epsilon = 1e-10;
    let denominator = Complex64::new(q_squared - mass_sq, epsilon);
    -Complex64::i() / denominator
}

/// Yukawa potential (scalar exchange): V(r) = -g² exp(-mr)/(4π r)
pub fn yukawa_potential(r: f64, mass: f64, coupling: f64) -> f64 {
    if r < 1e-10 {
        return 0.0;
    }
    let g_sq = coupling * coupling;
    -(g_sq / (4.0 * PI * r)) * (-mass * r).exp()
}

/// Coulomb potential (photon exchange, m→0 limit): V(r) = -α/(4π r)
pub fn coulomb_potential(r: f64, alpha: f64) -> f64 {
    if r < 1e-10 {
        return 0.0;
    }
    -alpha / (4.0 * PI * r)
}

