//! Friedmann Equations - Evolution of the Universe

use super::{CosmologyParams, Redshift, C, G, PI};
use serde::{Deserialize, Serialize};

/// Hubble parameter H(z) as function of redshift
/// H²(z) = H₀² [Ω_m(1+z)³ + Ω_r(1+z)⁴ + Ω_λ + Ω_k(1+z)²]
pub fn hubble_parameter(z: Redshift, params: &CosmologyParams) -> f64 {
    let h0 = params.hubble_constant_si();
    let one_plus_z = 1.0 + z.0;

    let term_m = params.omega_m * one_plus_z.powi(3);
    let term_r = params.omega_r * one_plus_z.powi(4);
    let term_lambda = params.omega_lambda;
    let term_k = params.omega_k() * one_plus_z.powi(2);

    h0 * (term_m + term_r + term_lambda + term_k).sqrt()
}

/// Deceleration parameter q = -ä·a/ȧ²
/// q(z) = [Ω_m(1+z)³ + Ω_r(1+z)⁴ - 2Ω_λ] / [2H²(z)/H₀²]
pub fn deceleration_parameter(z: Redshift, params: &CosmologyParams) -> f64 {
    let one_plus_z = 1.0 + z.0;
    let h_ratio_sq = (hubble_parameter(z, params) / params.hubble_constant_si()).powi(2);

    let numerator = params.omega_m * one_plus_z.powi(3)
        + params.omega_r * one_plus_z.powi(4)
        - 2.0 * params.omega_lambda;

    numerator / (2.0 * h_ratio_sq)
}

/// Lookback time: time between emission at redshift z and today
pub fn lookback_time(z: Redshift, params: &CosmologyParams) -> f64 {
    // Numerical integration of dt/dz = -(1+z)/(H(z))
    let n_steps = 100;
    let dz = z.0 / n_steps as f64;
    let mut time = 0.0;

    for i in 0..n_steps {
        let z_i = z.0 - (i as f64 + 0.5) * dz;
        let z_point = Redshift(z_i);
        let h_z = hubble_parameter(z_point, params);

        time += dz / ((1.0 + z_i) * h_z);
    }

    time
}

/// Comoving distance (line-of-sight distance in expanding universe)
pub fn comoving_distance(z: Redshift, params: &CosmologyParams) -> f64 {
    // D_c = c ∫₀^z dz'/H(z')
    let n_steps = 100;
    let dz = z.0 / n_steps as f64;
    let mut distance = 0.0;

    for i in 0..n_steps {
        let z_i = (i as f64 + 0.5) * dz;
        let z_point = Redshift(z_i);
        let h_z = hubble_parameter(z_point, params);

        distance += C * dz / h_z;
    }

    distance
}

/// Luminosity distance (apparent brightness distance)
/// D_L = (1+z) D_c
pub fn luminosity_distance(z: Redshift, params: &CosmologyParams) -> f64 {
    let d_c = comoving_distance(z, params);
    (1.0 + z.0) * d_c
}

/// Angular diameter distance (apparent size distance)
/// D_A = D_c / (1+z)
pub fn angular_diameter_distance(z: Redshift, params: &CosmologyParams) -> f64 {
    let d_c = comoving_distance(z, params);
    d_c / (1.0 + z.0)
}

/// Cosmic event horizon (maximum distance light can travel)
pub fn event_horizon(params: &CosmologyParams) -> f64 {
    // Approximate for flat ΛCDM
    let h0 = params.hubble_constant_si();
    C / h0 * (params.omega_lambda.sqrt())
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UniverseEvolution {
    pub redshifts: Vec<f64>,
    pub scale_factors: Vec<f64>,
    pub hubble_params: Vec<f64>,
    pub ages: Vec<f64>,
}

/// Simulate universe evolution from z_max to z=0
pub fn evolve_universe(z_max: f64, n_steps: usize, params: &CosmologyParams) -> UniverseEvolution {
    let dz = z_max / n_steps as f64;

    let mut redshifts = vec![];
    let mut scale_factors = vec![];
    let mut hubble_params = vec![];
    let mut ages = vec![];

    for i in 0..=n_steps {
        let z = z_max - i as f64 * dz;
        let z_obj = Redshift(z);

        redshifts.push(z);
        scale_factors.push(z_obj.to_scale_factor());
        hubble_params.push(hubble_parameter(z_obj, params));
        ages.push(z_obj.age(params));
    }

    UniverseEvolution {
        redshifts,
        scale_factors,
        hubble_params,
        ages,
    }
}

