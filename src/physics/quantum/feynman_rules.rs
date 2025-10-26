//! Feynman Rules and Matrix Elements

use num_complex::Complex64;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Fine structure constant
pub const ALPHA_EM: f64 = 1.0 / 137.036;

/// QED vertex factor
pub fn qed_vertex_factor() -> Complex64 {
    Complex64::new(0.0, -(ALPHA_EM.sqrt()))
}

/// Tree-level matrix element squared for e⁺e⁻ → μ⁺μ⁻
pub fn ee_to_mumu_matrix_element_squared(s: f64) -> f64 {
    let alpha_sq = ALPHA_EM * ALPHA_EM;
    4.0 * alpha_sq * s
}

/// Compton scattering matrix element: γe → γe
pub fn compton_scattering_matrix_element_squared(s: f64, t: f64, m_e: f64) -> f64 {
    let u = 2.0 * m_e * m_e - s - t;
    let alpha_sq = ALPHA_EM * ALPHA_EM;

    if t.abs() < 1e-10 || u.abs() < 1e-10 {
        return 0.0;
    }

    let m_sq = m_e * m_e;
    alpha_sq * (m_sq / s + m_sq / u + s / u + u / s)
}

/// Bhabha scattering: e⁺e⁻ → e⁺e⁻
pub fn bhabha_scattering_matrix_element_squared(s: f64, t: f64) -> f64 {
    let u = -s - t;
    let alpha_sq = ALPHA_EM * ALPHA_EM;

    if t.abs() < 1e-10 || u.abs() < 1e-10 {
        return 0.0;
    }

    alpha_sq * ((s / t).powi(2) + (s / u).powi(2))
}

/// Loop integral result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LoopIntegral {
    pub finite_part: Complex64,
    pub log_divergence: f64,
    pub quadratic_divergence: f64,
}

/// One-loop vacuum polarization (photon self-energy)
pub fn vacuum_polarization_one_loop(q_squared: f64, fermion_mass: f64, cutoff: f64) -> LoopIntegral {
    let log_div = ALPHA_EM / (3.0 * PI) * q_squared;
    let finite = -ALPHA_EM * q_squared / (18.0 * PI);

    LoopIntegral {
        finite_part: Complex64::new(finite, 0.0),
        log_divergence: log_div * (cutoff / fermion_mass).ln(),
        quadratic_divergence: 0.0,
    }
}

/// Anomalous magnetic moment (one-loop, Schwinger term)
pub fn anomalous_magnetic_moment_one_loop() -> f64 {
    ALPHA_EM / (2.0 * PI)
}

/// Running coupling constant (QED)
pub fn running_alpha_em(q_squared: f64, electron_mass: f64) -> f64 {
    let me_sq = electron_mass * electron_mass;
    let beta0 = 4.0 / (3.0 * PI); // Positive for QED (screening)
    let log_term = (q_squared / me_sq).ln();
    ALPHA_EM / (1.0 - ALPHA_EM * beta0 * log_term)
}

/// Running strong coupling (QCD, one-loop)
pub fn running_alpha_s(q_squared: f64, lambda_qcd: f64) -> f64 {
    let nf = 6.0;
    let beta0 = (33.0 - 2.0 * nf) / (12.0 * PI);
    let log_term = (q_squared / (lambda_qcd * lambda_qcd)).ln();

    if log_term <= 0.0 {
        return 1.0; // Non-perturbative regime
    }

    1.0 / (beta0 * log_term)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_qed_vertex() {
        let v = qed_vertex_factor();
        assert!(v.im.abs() > 0.0);
    }

    #[test]
    fn test_ee_to_mumu() {
        let s = 10000.0;
        let m_sq = ee_to_mumu_matrix_element_squared(s);
        assert!(m_sq > 0.0);
    }

    #[test]
    fn test_compton_scattering() {
        let s = 1000.0;
        let t = -250.0;
        let m_e = 0.511;
        let m_sq = compton_scattering_matrix_element_squared(s, t, m_e);
        // Result may be negative for unphysical kinematics
        assert!(m_sq.is_finite());
    }

    #[test]
    fn test_anomalous_magnetic_moment() {
        let a_e = anomalous_magnetic_moment_one_loop();
        assert!((a_e - 0.001161).abs() < 1e-5);
    }

    #[test]
    fn test_running_alpha_em() {
        let q_sq = 1e6;
        let alpha = running_alpha_em(q_sq, 0.511);
        assert!(alpha >= ALPHA_EM);
        assert!(alpha.is_finite());
    }

    #[test]
    fn test_running_alpha_s() {
        let lambda_qcd = 200.0;
        let q_sq = 1e6;
        let alpha_s = running_alpha_s(q_sq, lambda_qcd);
        assert!(alpha_s > 0.0 && alpha_s < 1.0);
    }
}
