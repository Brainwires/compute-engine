//! Plasma Waves and Instabilities

use super::{PlasmaParams, E, M_E};
use serde::{Deserialize, Serialize};

/// Wave dispersion relation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WaveDispersion {
    pub omega: f64,      // Angular frequency (rad/s)
    pub k: f64,          // Wave number (m⁻¹)
    pub v_phase: f64,    // Phase velocity (m/s)
    pub v_group: f64,    // Group velocity (m/s)
}

/// Langmuir wave (electron plasma oscillation)
/// ω² = ω_pe² + 3k²v_th,e²
pub fn langmuir_wave(k: f64, params: &PlasmaParams) -> WaveDispersion {
    let omega_pe = super::plasma_frequency(params.n_e);
    let v_th = super::thermal_velocity(params.t_e, M_E);

    let omega_sq = omega_pe * omega_pe + 3.0 * k * k * v_th * v_th;
    let omega = omega_sq.sqrt();

    let v_phase = omega / k;
    let v_group = 3.0 * v_th * v_th * k / omega;

    WaveDispersion {
        omega,
        k,
        v_phase,
        v_group,
    }
}

/// Alfvén wave (MHD wave along field lines)
/// ω = k∥ v_A
pub fn alfven_wave(k_parallel: f64, params: &PlasmaParams) -> WaveDispersion {
    let rho = params.n_i * params.ion_mass;
    let v_a = super::alfven_velocity(params.b_field, rho);

    let omega = k_parallel * v_a;

    WaveDispersion {
        omega,
        k: k_parallel,
        v_phase: v_a,
        v_group: v_a,
    }
}

/// Ion acoustic wave
/// ω = k √(k_B T_e / m_i) for T_i << T_e
pub fn ion_acoustic_wave(k: f64, params: &PlasmaParams) -> WaveDispersion {
    let c_s = (params.t_e * E / params.ion_mass).sqrt(); // Sound speed

    let omega = k * c_s;

    WaveDispersion {
        omega,
        k,
        v_phase: c_s,
        v_group: c_s,
    }
}

/// Upper hybrid frequency: ω_uh = √(ω_pe² + ω_ce²)
pub fn upper_hybrid_frequency(params: &PlasmaParams) -> f64 {
    let omega_pe = super::plasma_frequency(params.n_e);
    let omega_ce = super::electron_cyclotron_frequency(params.b_field);

    (omega_pe * omega_pe + omega_ce * omega_ce).sqrt()
}

/// Lower hybrid frequency: ω_lh ≈ √(ω_ci · ω_ce)
pub fn lower_hybrid_frequency(params: &PlasmaParams) -> f64 {
    let omega_ci = super::ion_cyclotron_frequency(params.b_field, params.z_ion, params.ion_mass);
    let omega_ce = super::electron_cyclotron_frequency(params.b_field);

    (omega_ci * omega_ce).sqrt()
}

/// Two-stream instability growth rate
/// γ = √(3)/2 * ω_pe * (n_b/n_0)^(1/3)
pub fn two_stream_growth_rate(n_beam: f64, n_background: f64) -> f64 {
    let omega_pe = super::plasma_frequency(n_background);
    let ratio = n_beam / n_background;

    (3.0_f64).sqrt() / 2.0 * omega_pe * ratio.powf(1.0 / 3.0)
}

/// Weibel instability (temperature anisotropy)
pub fn weibel_growth_rate(params: &PlasmaParams, t_perp: f64, t_parallel: f64) -> f64 {
    let omega_pe = super::plasma_frequency(params.n_e);
    let anisotropy = (t_perp - t_parallel) / t_parallel;

    if anisotropy > 0.0 {
        omega_pe * anisotropy.sqrt()
    } else {
        0.0
    }
}

/// Electromagnetic wave in plasma (dispersion relation)
/// ω² = ω_pe² + k²c²
pub fn em_wave_in_plasma(k: f64, params: &PlasmaParams) -> WaveDispersion {
    let omega_pe = super::plasma_frequency(params.n_e);
    let c = super::C;

    let omega_sq = omega_pe * omega_pe + k * k * c * c;
    let omega = omega_sq.sqrt();

    let v_phase = omega / k;
    let v_group = k * c * c / omega;

    WaveDispersion {
        omega,
        k,
        v_phase,
        v_group,
    }
}

/// Cutoff frequency (wave cannot propagate below this)
pub fn cutoff_frequency(params: &PlasmaParams) -> f64 {
    super::plasma_frequency(params.n_e)
}

