//! Scattering Cross Sections

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Fine structure constant
pub const ALPHA_EM: f64 = 1.0 / 137.036;

/// Convert GeV⁻² to barns
const GEV2_TO_BARN: f64 = 0.3894;

/// Total cross section for e⁺e⁻ → μ⁺μ⁻
pub fn ee_to_mumu_cross_section(s: f64, m_mu: f64) -> f64 {
    let m_sq = m_mu * m_mu;

    if s < 4.0 * m_sq {
        return 0.0;
    }

    let beta = (1.0 - 4.0 * m_sq / s).sqrt();
    let factor = 1.0 + m_sq / s;

    let sigma_gev2 = (4.0 * PI * ALPHA_EM * ALPHA_EM) / (3.0 * s) * factor * beta;
    sigma_gev2 * GEV2_TO_BARN
}

/// Compton scattering cross section: γe → γe (Klein-Nishina formula)
pub fn compton_scattering_cross_section(omega: f64, m_e: f64) -> f64 {
    let x = omega / m_e;

    if x < 1e-3 {
        // Thomson limit
        let r_e = 2.818e-15;
        return (8.0 * PI / 3.0) * r_e * r_e * 1e28;
    }

    // Full Klein-Nishina
    let r_e = 2.818e-15;
    let r_e_sq = r_e * r_e;

    let factor1 = (1.0 + x) / (x * x * x);
    let factor2 = (2.0 * x * (1.0 + x)) / (1.0 + 2.0 * x) - ((1.0 + 2.0 * x).ln());
    let factor3 = ((1.0 + 2.0 * x).ln()) / (2.0 * x);
    let factor4 = (1.0 + 3.0 * x) / ((1.0 + 2.0 * x).powi(2));

    let sigma = 2.0 * PI * r_e_sq * (factor1 * (factor2 + factor3) - factor4);
    sigma * 1e28
}

/// Rutherford scattering differential cross section
pub fn rutherford_scattering_differential(theta: f64, z1: f64, z2: f64, energy: f64) -> f64 {
    let half_theta = theta / 2.0;
    let sin_half = half_theta.sin();

    if sin_half.abs() < 1e-10 {
        return 0.0;
    }

    let numerator = (z1 * z2 * ALPHA_EM).powi(2);
    let denominator = 16.0 * energy * energy * sin_half.powi(4);

    numerator / denominator
}

/// Breit-Wigner resonance cross section
pub fn breit_wigner_cross_section(
    s: f64,
    resonance_mass: f64,
    resonance_width: f64,
    peak_cross_section: f64,
) -> f64 {
    let m_sq = resonance_mass * resonance_mass;
    let gamma = resonance_width;

    let denominator = (s - m_sq).powi(2) + (resonance_mass * gamma).powi(2);

    peak_cross_section * (m_sq * gamma * gamma) / denominator
}

/// Event rate from cross section and luminosity
pub fn event_rate(cross_section_barns: f64, luminosity_cm2s: f64) -> f64 {
    let sigma_cm2 = cross_section_barns * 1e-24;
    sigma_cm2 * luminosity_cm2s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ee_to_mumu() {
        let s = 10000.0; // MeV²
        let m_mu = 105.66; // MeV
        let sigma = ee_to_mumu_cross_section(s, m_mu);
        assert!(sigma >= 0.0);
        assert!(sigma.is_finite());
    }

    #[test]
    fn test_compton_thomson_limit() {
        let omega = 0.001; // MeV
        let m_e = 0.511; // MeV
        let sigma = compton_scattering_cross_section(omega, m_e);
        assert!(sigma > 0.0);
        assert!(sigma.is_finite());
    }

    #[test]
    fn test_rutherford_scattering() {
        let theta = PI / 4.0;
        let dsigma = rutherford_scattering_differential(theta, 1.0, 79.0, 5.0);
        assert!(dsigma > 0.0);
    }

    #[test]
    fn test_breit_wigner() {
        let m_z = 91187.6;
        let gamma_z = 2495.0;
        let peak_sigma = 1.0;

        let s_on = m_z * m_z;
        let sigma_on = breit_wigner_cross_section(s_on, m_z, gamma_z, peak_sigma);

        let s_off = (m_z + 10000.0).powi(2);
        let sigma_off = breit_wigner_cross_section(s_off, m_z, gamma_z, peak_sigma);

        assert!(sigma_on > sigma_off);
        assert!(sigma_on.is_finite());
    }

    #[test]
    fn test_event_rate() {
        let sigma = 1.0;
        let lumi = 1e34;
        let rate = event_rate(sigma, lumi);
        assert!(rate > 1e9 && rate < 1e11);
    }
}
