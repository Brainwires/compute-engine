//! Filter Design
//!
//! Design calculations for passive filters (RC, RL, RLC)

use std::f64::consts::PI;

/// Design RC low-pass filter
///
/// # Arguments
/// * `cutoff_frequency` - Desired cutoff frequency in Hz
/// * `resistance` - Resistor value in ohms (or None to calculate)
/// * `capacitance` - Capacitor value in farads (or None to calculate)
///
/// # Returns
/// (R, C) tuple for the filter
pub fn design_rc_lowpass(
    cutoff_frequency: f64,
    resistance: Option<f64>,
    capacitance: Option<f64>,
) -> (f64, f64) {
    match (resistance, capacitance) {
        (Some(r), None) => {
            let c = 1.0 / (2.0 * PI * cutoff_frequency * r);
            (r, c)
        }
        (None, Some(c)) => {
            let r = 1.0 / (2.0 * PI * cutoff_frequency * c);
            (r, c)
        }
        _ => {
            // Default: Use 10kΩ
            let r = 10_000.0;
            let c = 1.0 / (2.0 * PI * cutoff_frequency * r);
            (r, c)
        }
    }
}

/// Design RC high-pass filter
pub fn design_rc_highpass(
    cutoff_frequency: f64,
    resistance: Option<f64>,
    capacitance: Option<f64>,
) -> (f64, f64) {
    // Same calculation as low-pass, just different topology
    design_rc_lowpass(cutoff_frequency, resistance, capacitance)
}

/// Design second-order Butterworth low-pass filter
///
/// # Arguments
/// * `cutoff_frequency` - Cutoff frequency in Hz
/// * `load_resistance` - Load resistance in ohms
///
/// # Returns
/// (R1, R2, C1, C2) for Sallen-Key topology
pub fn design_butterworth_lowpass(cutoff_frequency: f64, load_resistance: f64) -> (f64, f64, f64, f64) {
    let omega = 2.0 * PI * cutoff_frequency;

    // For Butterworth: R1 = R2 = R, C1 = C2 = C
    // ω_c = 1/(RC)
    let r = load_resistance;
    let c = 1.0 / (omega * r);

    (r, r, c, c)
}

/// Design band-pass filter using RLC circuit
///
/// # Arguments
/// * `center_frequency` - Center frequency in Hz
/// * `bandwidth` - Bandwidth in Hz
/// * `resistance` - Series resistance in ohms
///
/// # Returns
/// (R, L, C) for series RLC band-pass
pub fn design_rlc_bandpass(
    center_frequency: f64,
    bandwidth: f64,
    resistance: f64,
) -> (f64, f64, f64) {
    let omega_0 = 2.0 * PI * center_frequency;
    let q = center_frequency / bandwidth;

    // Q = ω₀L/R = 1/(ω₀RC)
    let l = (q * resistance) / omega_0;
    let c = 1.0 / (omega_0 * omega_0 * l);

    (resistance, l, c)
}

/// Design notch (band-stop) filter
pub fn design_notch_filter(
    notch_frequency: f64,
    q_factor: f64,
    resistance: f64,
) -> (f64, f64, f64) {
    let omega_0 = 2.0 * PI * notch_frequency;

    let l = (q_factor * resistance) / omega_0;
    let c = 1.0 / (omega_0 * omega_0 * l);

    (resistance, l, c)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rc_lowpass_design() {
        let (r, c) = design_rc_lowpass(1000.0, Some(10_000.0), None);
        assert_eq!(r, 10_000.0);
        // f_c = 1/(2πRC) => C = 1/(2π*1000*10000) ≈ 15.9 nF
        assert!((c - 15.9e-9).abs() < 1e-9);
    }

    #[test]
    fn test_butterworth() {
        let (r1, r2, c1, c2) = design_butterworth_lowpass(1000.0, 10_000.0);
        assert_eq!(r1, r2); // Butterworth has equal resistors
        assert_eq!(c1, c2); // And equal capacitors
    }
}
