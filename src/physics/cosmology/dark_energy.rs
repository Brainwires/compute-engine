//! Dark Energy and Equation of State

use super::{CosmologyParams, Redshift};
use serde::{Deserialize, Serialize};

/// Dark energy equation of state: w = P/ρ
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum DarkEnergyModel {
    /// Cosmological constant: w = -1
    LambdaCDM,
    /// Quintessence: -1 < w < -1/3
    Quintessence { w0: f64 },
    /// Phantom energy: w < -1
    Phantom { w0: f64 },
    /// Time-varying: w(z) = w0 + wa z/(1+z)
    CPL { w0: f64, wa: f64 },
}

impl DarkEnergyModel {
    /// Equation of state parameter at given redshift
    pub fn w(&self, z: f64) -> f64 {
        match self {
            DarkEnergyModel::LambdaCDM => -1.0,
            DarkEnergyModel::Quintessence { w0 } => *w0,
            DarkEnergyModel::Phantom { w0 } => *w0,
            DarkEnergyModel::CPL { w0, wa } => {
                // Chevallier-Polarski-Linder parametrization
                w0 + wa * z / (1.0 + z)
            }
        }
    }

    /// Dark energy density evolution: ρ_DE(z) = ρ_DE,0 (1+z)^(3(1+w))
    pub fn density_evolution(&self, z: Redshift) -> f64 {
        match self {
            DarkEnergyModel::LambdaCDM => {
                // Constant density
                1.0
            }
            DarkEnergyModel::Quintessence { w0 } | DarkEnergyModel::Phantom { w0 } => {
                (1.0 + z.0).powf(3.0 * (1.0 + w0))
            }
            DarkEnergyModel::CPL { w0, wa } => {
                // Numerical integration needed for CPL
                // Simplified approximation
                let w_eff = w0 + wa * z.0 / (2.0 * (1.0 + z.0));
                (1.0 + z.0).powf(3.0 * (1.0 + w_eff))
            }
        }
    }
}

/// Modified Hubble parameter with dark energy model
pub fn hubble_parameter_de(
    z: Redshift,
    params: &CosmologyParams,
    de_model: &DarkEnergyModel,
) -> f64 {
    let h0 = params.hubble_constant_si();
    let one_plus_z = 1.0 + z.0;

    let term_m = params.omega_m * one_plus_z.powi(3);
    let term_r = params.omega_r * one_plus_z.powi(4);
    let term_de = params.omega_lambda * de_model.density_evolution(z);
    let term_k = params.omega_k() * one_plus_z.powi(2);

    h0 * (term_m + term_r + term_de + term_k).sqrt()
}

/// Rip time for phantom energy (Big Rip singularity)
pub fn big_rip_time(params: &CosmologyParams, w: f64) -> Option<f64> {
    if w >= -1.0 {
        return None; // No Big Rip for w >= -1
    }

    let h0 = params.hubble_constant_si();
    let factor = 2.0 / (3.0 * (1.0 + w) * h0 * params.omega_lambda.sqrt());

    Some(factor.abs())
}

/// Future scale factor evolution
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FutureEvolution {
    pub times: Vec<f64>,      // Time in seconds from now
    pub scale_factors: Vec<f64>,
    pub hubble_params: Vec<f64>,
}

/// Predict future evolution of universe
pub fn future_evolution(
    duration: f64,
    n_steps: usize,
    params: &CosmologyParams,
    de_model: &DarkEnergyModel,
) -> FutureEvolution {
    let dt = duration / n_steps as f64;
    let h0 = params.hubble_constant_si();

    let mut times = vec![];
    let mut scale_factors = vec![1.0]; // Start at a=1 today
    let mut hubble_params = vec![];

    for i in 0..=n_steps {
        let t = i as f64 * dt;
        times.push(t);

        // Simple forward integration
        if i > 0 {
            let a_prev = scale_factors[i - 1];
            let z_prev = Redshift(1.0 / a_prev - 1.0);
            let h = hubble_parameter_de(z_prev, params, de_model);

            // da/dt = a * H
            let a_new = a_prev + a_prev * h * dt;
            scale_factors.push(a_new);
        }

        let a = scale_factors[i];
        let z = Redshift((1.0 / a - 1.0).max(0.0));
        let h = hubble_parameter_de(z, params, de_model);
        hubble_params.push(h);
    }

    FutureEvolution {
        times,
        scale_factors,
        hubble_params,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lambda_cdm_eos() {
        let model = DarkEnergyModel::LambdaCDM;

        assert_eq!(model.w(0.0), -1.0);
        assert_eq!(model.w(1.0), -1.0);
        assert_eq!(model.w(1000.0), -1.0);
    }

    #[test]
    fn test_quintessence_eos() {
        let model = DarkEnergyModel::Quintessence { w0: -0.8 };

        let w = model.w(0.0);
        assert!(w > -1.0 && w < -1.0/3.0);
    }

    #[test]
    fn test_phantom_eos() {
        let model = DarkEnergyModel::Phantom { w0: -1.2 };

        let w = model.w(0.0);
        assert!(w < -1.0);
    }

    #[test]
    fn test_cpl_evolution() {
        let model = DarkEnergyModel::CPL { w0: -0.9, wa: -0.1 };

        let w_today = model.w(0.0);
        let w_past = model.w(1.0);

        // w should vary with redshift
        assert!((w_today - (-0.9)).abs() < 0.01);
        assert!(w_past != w_today);
    }

    #[test]
    fn test_lambda_density_constant() {
        let model = DarkEnergyModel::LambdaCDM;

        let rho_today = model.density_evolution(Redshift(0.0));
        let rho_past = model.density_evolution(Redshift(1.0));

        // Cosmological constant has constant density
        assert_eq!(rho_today, 1.0);
        assert_eq!(rho_past, 1.0);
    }

    #[test]
    fn test_quintessence_density_evolution() {
        let model = DarkEnergyModel::Quintessence { w0: -0.8 };

        let rho_today = model.density_evolution(Redshift(0.0));
        let rho_past = model.density_evolution(Redshift(1.0));

        // Quintessence density should increase going back in time
        assert!(rho_past > rho_today);
    }

    #[test]
    fn test_big_rip_phantom() {
        let params = CosmologyParams::planck_2018();

        // Phantom energy (w < -1) leads to Big Rip
        let t_rip = big_rip_time(&params, -1.2);
        assert!(t_rip.is_some());
        assert!(t_rip.unwrap() > 0.0);

        // Lambda CDM (w = -1) has no Big Rip
        let t_rip_lambda = big_rip_time(&params, -1.0);
        assert!(t_rip_lambda.is_none());
    }

    #[test]
    fn test_hubble_parameter_de() {
        let params = CosmologyParams::planck_2018();
        let model = DarkEnergyModel::LambdaCDM;

        let h_today = hubble_parameter_de(Redshift(0.0), &params, &model);
        let h0 = params.hubble_constant_si();

        // Should equal H0 today
        assert!((h_today - h0).abs() / h0 < 0.01);
    }

    #[test]
    fn test_future_evolution() {
        let params = CosmologyParams::planck_2018();
        let model = DarkEnergyModel::LambdaCDM;

        // Simulate 1 billion years into future
        let evolution = future_evolution(3.15e16, 100, &params, &model);

        assert_eq!(evolution.times.len(), 101);
        assert_eq!(evolution.scale_factors.len(), 101);

        // Scale factor should increase
        assert!(evolution.scale_factors.last().unwrap() > &1.0);

        // All values should be finite
        for a in &evolution.scale_factors {
            assert!(a.is_finite());
            assert!(*a > 0.0);
        }
    }

    #[test]
    fn test_accelerating_expansion() {
        let params = CosmologyParams::planck_2018();
        let model = DarkEnergyModel::LambdaCDM;

        let evolution = future_evolution(1e17, 100, &params, &model);

        // Check that expansion is accelerating (da/dt increasing)
        let early_rate = (evolution.scale_factors[10] - evolution.scale_factors[0]) / 10.0;
        let late_rate = (evolution.scale_factors[100] - evolution.scale_factors[90]) / 10.0;

        assert!(late_rate > early_rate);
    }
}
