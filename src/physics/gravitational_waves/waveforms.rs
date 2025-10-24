//! Gravitational Waveform Generation
//!
//! Post-Newtonian approximations for binary inspiral waveforms

use super::{BinarySystem, Waveform, C, G, M_SUN, PI};
use serde::{Deserialize, Serialize};

/// Generate inspiral waveform using TaylorF2 approximation
pub fn taylor_f2_waveform(
    binary: &BinarySystem,
    f_min: f64,
    f_max: f64,
    sample_rate: f64,
) -> Waveform {
    let dt = 1.0 / sample_rate;
    let m_chirp_si = binary.chirp_mass() * M_SUN;
    let distance_si = binary.distance * 1e6 * super::PC; // Mpc to meters

    let mut times = vec![];
    let mut h_plus = vec![];
    let mut h_cross = vec![];
    let mut frequencies = vec![];

    let mut f = f_min;
    let mut t = 0.0;

    while f < f_max && f.is_finite() {
        // Post-Newtonian phase evolution
        let psi = phase_at_frequency(binary, f);

        // Amplitude: h₀ = (G M_chirp / c² d)^(5/6) (π f)^(2/3)
        let amp = (G * m_chirp_si / (C * C * distance_si)).powf(5.0 / 6.0)
            * (PI * f).powf(2.0 / 3.0);

        // Polarizations (assuming face-on for simplicity)
        let cos_psi = psi.cos();
        let sin_psi = psi.sin();

        h_plus.push(amp * (1.0 + binary.inclination.cos().powi(2)) * cos_psi);
        h_cross.push(amp * 2.0 * binary.inclination.cos() * sin_psi);

        times.push(t);
        frequencies.push(f);

        // Frequency evolution: df/dt = (96/5) π^(8/3) (G M_chirp / c³)^(5/3) f^(11/3)
        let df_dt = (96.0 / 5.0) * PI.powf(8.0 / 3.0)
            * (G * m_chirp_si / C.powi(3)).powf(5.0 / 3.0)
            * f.powf(11.0 / 3.0);

        f += df_dt * dt;
        t += dt;

        // Safety limit
        if times.len() > 1_000_000 {
            break;
        }
    }

    Waveform { times, h_plus, h_cross, frequency: frequencies }
}

/// Post-Newtonian phase at given frequency (2PN approximation)
fn phase_at_frequency(binary: &BinarySystem, f: f64) -> f64 {
    let m_chirp = binary.chirp_mass() * M_SUN;
    let eta = binary.symmetric_mass_ratio();

    let v = (PI * G * binary.total_mass() * M_SUN * f / C.powi(3)).powf(1.0 / 3.0);
    let v2 = v * v;
    let v3 = v2 * v;
    let v4 = v2 * v2;

    // Newtonian term
    let phi_n = binary.coalescence_phase - 2.0 / (eta * v.powf(5.0));

    // PN corrections
    let phi_1pn = phi_n * (1.0 + (3715.0 / 756.0 + 55.0 / 9.0 * eta) * v2);
    let phi_2pn = phi_1pn + (-16.0 * PI * v3);

    phi_2pn
}

/// Generate ringdown waveform (post-merger)
pub fn ringdown_waveform(
    binary: &BinarySystem,
    duration: f64,
    sample_rate: f64,
) -> Waveform {
    let dt = 1.0 / sample_rate;
    let n_samples = (duration * sample_rate) as usize;

    let m_final = binary.total_mass() * M_SUN;
    let distance_si = binary.distance * 1e6 * super::PC;

    // Ringdown frequency and damping time (approximations)
    let f_ringdown = C.powi(3) / (2.0 * PI * G * m_final) / 6.28; // ~1/(2π r_s)
    let tau = 4.0 * G * m_final / C.powi(3); // Damping time ~M

    let mut times = vec![];
    let mut h_plus = vec![];
    let mut h_cross = vec![];
    let mut frequencies = vec![];

    // Initial amplitude
    let amp0 = G * m_final / (C * C * distance_si);

    for i in 0..n_samples {
        let t = i as f64 * dt;
        let exp_decay = (-t / tau).exp();
        let osc = (2.0 * PI * f_ringdown * t).cos();

        let h = amp0 * exp_decay * osc;

        times.push(t);
        h_plus.push(h);
        h_cross.push(h * 0.5); // Simplified polarization
        frequencies.push(f_ringdown);
    }

    Waveform { times, h_plus, h_cross, frequency: frequencies }
}

/// Simplified IMRPhenomD-like waveform (Inspiral-Merger-Ringdown)
pub fn imr_waveform(binary: &BinarySystem, sample_rate: f64) -> Waveform {
    // Generate inspiral up to ISCO
    let f_isco = binary.isco_frequency();
    let f_min = 10.0; // Start at 10 Hz

    let mut inspiral = taylor_f2_waveform(binary, f_min, f_isco * 0.9, sample_rate);

    // Generate ringdown
    let ringdown = ringdown_waveform(binary, 0.1, sample_rate); // 100ms

    // Concatenate (simplified - real IMR has intermediate merger phase)
    let time_offset = inspiral.times.last().copied().unwrap_or(0.0);
    inspiral.times.extend(ringdown.times.iter().map(|t| t + time_offset));
    inspiral.h_plus.extend(ringdown.h_plus);
    inspiral.h_cross.extend(ringdown.h_cross);
    inspiral.frequency.extend(ringdown.frequency);

    inspiral
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::gravitational_waves::BinaryType;

    #[test]
    fn test_taylor_f2_waveform() {
        let binary = BinarySystem {
            mass1: 30.0,
            mass2: 30.0,
            distance: 410.0,
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryBlackHole,
        };

        let waveform = taylor_f2_waveform(&binary, 20.0, 100.0, 2048.0);

        assert!(waveform.len() > 0);
        assert_eq!(waveform.times.len(), waveform.h_plus.len());
        assert_eq!(waveform.times.len(), waveform.h_cross.len());

        // Check that frequency increases monotonically
        for i in 1..waveform.frequency.len().min(100) {
            assert!(waveform.frequency[i] >= waveform.frequency[i - 1]);
        }
    }

    #[test]
    fn test_ringdown_waveform() {
        let binary = BinarySystem {
            mass1: 30.0,
            mass2: 30.0,
            distance: 410.0,
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryBlackHole,
        };

        let waveform = ringdown_waveform(&binary, 0.05, 2048.0);

        assert!(waveform.len() > 0);

        // Ringdown should decay exponentially
        let early_amp = waveform.h_plus[10].abs();
        let late_amp = waveform.h_plus[waveform.len() - 10].abs();
        assert!(late_amp < early_amp);
    }

    #[test]
    fn test_imr_waveform() {
        let binary = BinarySystem {
            mass1: 36.0,
            mass2: 29.0,
            distance: 410.0,
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryBlackHole,
        };

        let waveform = imr_waveform(&binary, 2048.0);

        assert!(waveform.len() > 100);

        // Peak strain should be finite and non-zero
        let peak = waveform.peak_strain();
        assert!(peak > 0.0 && peak.is_finite());
    }

    #[test]
    fn test_waveform_strain_amplitude() {
        let binary = BinarySystem {
            mass1: 30.0,
            mass2: 30.0,
            distance: 100.0, // Closer -> larger strain
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryBlackHole,
        };

        let waveform = taylor_f2_waveform(&binary, 30.0, 50.0, 1024.0);

        // Strain should be finite
        for &h in &waveform.h_plus {
            assert!(h.is_finite());
        }

        // At least some non-zero values
        let has_signal = waveform.h_plus.iter().any(|&h| h.abs() > 1e-30);
        assert!(has_signal);
    }
}
