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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hubble_parameter_today() {
        let params = CosmologyParams::planck_2018();
        let h_today = hubble_parameter(Redshift(0.0), &params);
        let h0 = params.hubble_constant_si();

        // Should equal H0 today
        assert!((h_today - h0).abs() < 1e-10);
    }

    #[test]
    fn test_hubble_parameter_evolution() {
        let params = CosmologyParams::planck_2018();

        let h_today = hubble_parameter(Redshift(0.0), &params);
        let h_early = hubble_parameter(Redshift(1000.0), &params);

        // H(z) should be larger in the early universe
        assert!(h_early > h_today);
    }

    #[test]
    fn test_deceleration_parameter() {
        let params = CosmologyParams::planck_2018();

        let q_today = deceleration_parameter(Redshift(0.0), &params);

        // With dark energy, universe is accelerating (q < 0)
        assert!(q_today < 0.0);
    }

    #[test]
    fn test_lookback_time() {
        let params = CosmologyParams::planck_2018();

        let t_lookback = lookback_time(Redshift(1.0), &params);

        // Lookback time to z=1 should be several billion years
        assert!(t_lookback > 1e17); // >3 billion years in seconds
        assert!(t_lookback.is_finite());
    }

    #[test]
    fn test_comoving_distance() {
        let params = CosmologyParams::planck_2018();

        let d_c = comoving_distance(Redshift(1.0), &params);

        // Comoving distance to z=1 should be a few Gpc
        assert!(d_c > 0.0);
        assert!(d_c.is_finite());
    }

    #[test]
    fn test_luminosity_distance() {
        let params = CosmologyParams::planck_2018();

        let z = Redshift(1.0);
        let d_c = comoving_distance(z, &params);
        let d_l = luminosity_distance(z, &params);

        // D_L = (1+z) D_c
        assert!((d_l - 2.0 * d_c).abs() < 1e-6);
    }

    #[test]
    fn test_angular_diameter_distance() {
        let params = CosmologyParams::planck_2018();

        let z = Redshift(1.0);
        let d_c = comoving_distance(z, &params);
        let d_a = angular_diameter_distance(z, &params);

        // D_A = D_c / (1+z)
        assert!((d_a * 2.0 - d_c).abs() < 1e-6);
    }

    #[test]
    fn test_distance_relationship() {
        let params = CosmologyParams::planck_2018();
        let z = Redshift(2.0);

        let d_l = luminosity_distance(z, &params);
        let d_a = angular_diameter_distance(z, &params);

        // D_L = (1+z)² D_A
        let expected_ratio = (1.0 + z.0).powi(2);
        let actual_ratio = d_l / d_a;

        assert!((actual_ratio - expected_ratio).abs() < 1e-6);
    }

    #[test]
    fn test_evolve_universe() {
        let params = CosmologyParams::planck_2018();
        let evolution = evolve_universe(10.0, 100, &params);

        assert_eq!(evolution.redshifts.len(), 101);
        assert_eq!(evolution.scale_factors.len(), 101);

        // Scale factor should decrease going back in time
        assert!(evolution.scale_factors[0] < *evolution.scale_factors.last().unwrap());

        // Hubble parameter should increase going back in time
        assert!(evolution.hubble_params[0] > *evolution.hubble_params.last().unwrap());
    }
}
