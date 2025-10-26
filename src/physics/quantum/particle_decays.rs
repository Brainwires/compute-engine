//! Particle Decay Widths and Lifetimes

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Physical constants
pub const HBAR: f64 = 1.054571817e-34; // J·s
pub const ALPHA_EM: f64 = 1.0 / 137.036;

/// Muon decay width: μ → e ν_μ ν̄_e
pub fn muon_decay_width(m_mu: f64) -> f64 {
    let g_f = 1.1663787e-5; // Fermi constant in GeV⁻²
    let m_gev = m_mu / 1000.0;
    let width_gev = (g_f * g_f * m_gev.powi(5)) / (192.0 * PI.powi(3));
    width_gev * 1000.0
}

/// Particle lifetime from decay width
pub fn lifetime_from_width(width_mev: f64) -> f64 {
    let width_joules = width_mev * 1.602176634e-13;
    HBAR / width_joules
}

/// W boson decay width (total)
pub fn w_boson_decay_width(m_w: f64) -> f64 {
    let g_f = 1.1663787e-5 * 1000.0; // MeV⁻²
    (3.0 * g_f * m_w.powi(3)) / (2.0 * 2.0_f64.sqrt() * PI)
}

/// Z boson decay width (total)
pub fn z_boson_decay_width(m_z: f64) -> f64 {
    let g_f = 1.1663787e-5 * 1000.0; // MeV⁻²
    (2.0_f64.sqrt() * g_f * m_z.powi(3)) / (6.0 * PI)
}

/// Higgs decay to fermions: H → f f̄
pub fn higgs_to_fermion_width(m_h: f64, m_f: f64, n_c: f64) -> f64 {
    let g_f = 1.1663787e-5 * 1000.0; // MeV⁻²

    if m_f > m_h / 2.0 {
        return 0.0;
    }

    let ratio = 4.0 * m_f * m_f / (m_h * m_h);
    if ratio >= 1.0 {
        return 0.0;
    }

    let suppression = (1.0 - ratio).powf(1.5);
    (n_c * g_f * m_f * m_f * m_h) / (4.0 * 2.0_f64.sqrt() * PI) * suppression
}

/// Top quark decay width: t → Wb
pub fn top_quark_decay_width(m_t: f64, m_w: f64) -> f64 {
    let g_f = 1.1663787e-5 * 1000.0; // MeV⁻²

    let r_w = m_w / m_t;
    let r_w_sq = r_w * r_w;

    let factor1 = (1.0 - r_w_sq).powi(2);
    let factor2 = 1.0 + 2.0 * r_w_sq;

    (g_f * m_t.powi(3)) / (8.0 * PI * 2.0_f64.sqrt()) * factor1 * factor2
}

/// Branching ratio
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BranchingRatio {
    pub channel: String,
    pub partial_width: f64,
    pub branching_ratio: f64,
}

/// Calculate branching ratios
pub fn calculate_branching_ratios(channels: Vec<(&str, f64)>) -> Vec<BranchingRatio> {
    let total_width: f64 = channels.iter().map(|(_, w)| w).sum();

    channels
        .into_iter()
        .map(|(name, width)| BranchingRatio {
            channel: name.to_string(),
            partial_width: width,
            branching_ratio: width / total_width,
        })
        .collect()
}

