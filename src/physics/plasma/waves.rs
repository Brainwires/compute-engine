//! Plasma Waves and Instabilities

use super::{PlasmaParams, E, M_E, EPSILON_0};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::plasma::{plasma_frequency, electron_cyclotron_frequency, C};

    #[test]
    fn test_langmuir_wave() {
        let params = PlasmaParams::tokamak();
        let k = 1e3; // m⁻¹

        let wave = langmuir_wave(k, &params);

        // Frequency should be close to plasma frequency
        let omega_pe = plasma_frequency(params.n_e);
        assert!(wave.omega >= omega_pe);
        assert!(wave.v_phase > 0.0);
    }

    #[test]
    fn test_alfven_wave() {
        let params = PlasmaParams::tokamak();
        let k = 1e3;

        let wave = alfven_wave(k, &params);

        // Alfvén wave is non-dispersive
        assert!((wave.v_phase - wave.v_group).abs() < 1e-6);
        assert!(wave.omega > 0.0);
    }

    #[test]
    fn test_ion_acoustic_wave() {
        let params = PlasmaParams::tokamak();
        let k = 1e3;

        let wave = ion_acoustic_wave(k, &params);

        // Ion acoustic wave is also non-dispersive
        assert!((wave.v_phase - wave.v_group).abs() / wave.v_phase < 0.01);
        assert!(wave.omega > 0.0);
    }

    #[test]
    fn test_hybrid_frequencies() {
        let params = PlasmaParams::tokamak();

        let omega_uh = upper_hybrid_frequency(&params);
        let omega_lh = lower_hybrid_frequency(&params);

        let omega_pe = plasma_frequency(params.n_e);
        let omega_ce = electron_cyclotron_frequency(params.b_field);

        // Upper hybrid should be higher than both
        assert!(omega_uh > omega_pe);
        assert!(omega_uh > omega_ce);

        // Lower hybrid should be lower
        assert!(omega_lh < omega_uh);
    }

    #[test]
    fn test_two_stream_instability() {
        let n_b = 1e18;
        let n_0 = 1e20;

        let gamma = two_stream_growth_rate(n_b, n_0);

        // Growth rate should be positive
        assert!(gamma > 0.0);
        assert!(gamma.is_finite());
    }

    #[test]
    fn test_weibel_instability() {
        let params = PlasmaParams::tokamak();
        let t_perp = 20000.0; // Hot perpendicular
        let t_parallel = 10000.0;

        let gamma = weibel_growth_rate(&params, t_perp, t_parallel);

        // Should be unstable (positive growth rate)
        assert!(gamma > 0.0);
    }

    #[test]
    fn test_weibel_stable() {
        let params = PlasmaParams::tokamak();
        let t_perp = 10000.0;
        let t_parallel = 20000.0; // Hot parallel

        let gamma = weibel_growth_rate(&params, t_perp, t_parallel);

        // Should be stable (zero growth rate)
        assert_eq!(gamma, 0.0);
    }

    #[test]
    fn test_em_wave_dispersion() {
        let params = PlasmaParams::tokamak();
        let k = 1e6;

        let wave = em_wave_in_plasma(k, &params);

        // At high k, should approach light speed
        assert!(wave.v_phase.is_finite());
        assert!(wave.v_group.is_finite());
        assert!(wave.v_group < C);
    }

    #[test]
    fn test_cutoff_frequency() {
        let params = PlasmaParams::tokamak();
        let f_cutoff = cutoff_frequency(&params);

        // Should equal plasma frequency
        let omega_pe = plasma_frequency(params.n_e);
        assert_eq!(f_cutoff, omega_pe);
    }
}
