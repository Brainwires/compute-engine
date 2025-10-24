//! AC Circuit Analysis
//!
//! AC circuit calculations including phasors, impedance, RLC circuits,
//! resonance, and frequency response.

use num_complex::Complex64;
use std::f64::consts::PI;

/// Calculate RMS (Root Mean Square) value from peak amplitude
///
/// RMS = Peak / √2
///
/// # Arguments
/// * `peak` - Peak amplitude
///
/// # Returns
/// RMS value
pub fn peak_to_rms(peak: f64) -> f64 {
    peak / 2_f64.sqrt()
}

/// Calculate peak amplitude from RMS value
///
/// Peak = RMS * √2
pub fn rms_to_peak(rms: f64) -> f64 {
    rms * 2_f64.sqrt()
}

/// Calculate capacitive reactance: X_C = 1 / (2πfC)
///
/// # Arguments
/// * `frequency` - Frequency in hertz (Hz)
/// * `capacitance` - Capacitance in farads (F)
///
/// # Returns
/// Capacitive reactance in ohms (Ω)
///
/// # Example
/// ```
/// use computational_engine::electrical::ac_analysis::capacitive_reactance;
/// let xc = capacitive_reactance(60.0, 10e-6); // 60Hz, 10μF
/// assert!((xc - 265.26).abs() < 0.01);
/// ```
pub fn capacitive_reactance(frequency: f64, capacitance: f64) -> f64 {
    1.0 / (2.0 * PI * frequency * capacitance)
}

/// Calculate inductive reactance: X_L = 2πfL
///
/// # Arguments
/// * `frequency` - Frequency in hertz (Hz)
/// * `inductance` - Inductance in henries (H)
///
/// # Returns
/// Inductive reactance in ohms (Ω)
///
/// # Example
/// ```
/// use computational_engine::electrical::ac_analysis::inductive_reactance;
/// let xl = inductive_reactance(60.0, 0.1); // 60Hz, 0.1H
/// assert!((xl - 37.699).abs() < 0.001);
/// ```
pub fn inductive_reactance(frequency: f64, inductance: f64) -> f64 {
    2.0 * PI * frequency * inductance
}

/// Calculate impedance of a resistor: Z = R
pub fn resistor_impedance(resistance: f64) -> Complex64 {
    Complex64::new(resistance, 0.0)
}

/// Calculate impedance of a capacitor: Z = -j / (2πfC)
pub fn capacitor_impedance(frequency: f64, capacitance: f64) -> Complex64 {
    let xc = capacitive_reactance(frequency, capacitance);
    Complex64::new(0.0, -xc)
}

/// Calculate impedance of an inductor: Z = j * 2πfL
pub fn inductor_impedance(frequency: f64, inductance: f64) -> Complex64 {
    let xl = inductive_reactance(frequency, inductance);
    Complex64::new(0.0, xl)
}

/// Calculate total impedance of series RLC circuit: Z = R + j(X_L - X_C)
///
/// # Arguments
/// * `frequency` - Frequency in hertz (Hz)
/// * `resistance` - Resistance in ohms (Ω)
/// * `inductance` - Inductance in henries (H)
/// * `capacitance` - Capacitance in farads (F)
///
/// # Returns
/// Complex impedance in ohms (Ω)
pub fn series_rlc_impedance(
    frequency: f64,
    resistance: f64,
    inductance: f64,
    capacitance: f64,
) -> Complex64 {
    let xl = inductive_reactance(frequency, inductance);
    let xc = capacitive_reactance(frequency, capacitance);
    Complex64::new(resistance, xl - xc)
}

/// Calculate parallel combination of two impedances
///
/// Z_parallel = (Z1 * Z2) / (Z1 + Z2)
pub fn parallel_impedance(z1: Complex64, z2: Complex64) -> Complex64 {
    (z1 * z2) / (z1 + z2)
}

/// Calculate magnitude of complex impedance
pub fn impedance_magnitude(z: Complex64) -> f64 {
    z.norm()
}

/// Calculate phase angle of complex impedance (in radians)
pub fn impedance_phase(z: Complex64) -> f64 {
    z.arg()
}

/// Calculate phase angle in degrees
pub fn impedance_phase_degrees(z: Complex64) -> f64 {
    z.arg() * 180.0 / PI
}

/// Calculate resonant frequency of series RLC circuit: f_0 = 1 / (2π√LC)
///
/// # Arguments
/// * `inductance` - Inductance in henries (H)
/// * `capacitance` - Capacitance in farads (F)
///
/// # Returns
/// Resonant frequency in hertz (Hz)
///
/// # Example
/// ```
/// use computational_engine::electrical::ac_analysis::resonant_frequency;
/// let f0 = resonant_frequency(0.1, 10e-6); // 0.1H, 10μF
/// assert!((f0 - 159.15).abs() < 0.01);
/// ```
pub fn resonant_frequency(inductance: f64, capacitance: f64) -> f64 {
    1.0 / (2.0 * PI * (inductance * capacitance).sqrt())
}

/// Calculate quality factor (Q) of series RLC circuit at resonance
///
/// Q = (ω₀L) / R = 1 / (ω₀RC) = (1/R)√(L/C)
///
/// # Arguments
/// * `resistance` - Resistance in ohms (Ω)
/// * `inductance` - Inductance in henries (H)
/// * `capacitance` - Capacitance in farads (F)
///
/// # Returns
/// Quality factor (dimensionless)
pub fn quality_factor(resistance: f64, inductance: f64, capacitance: f64) -> f64 {
    (1.0 / resistance) * (inductance / capacitance).sqrt()
}

/// Calculate bandwidth of RLC circuit: BW = f₀ / Q
///
/// # Arguments
/// * `resonant_freq` - Resonant frequency in Hz
/// * `q_factor` - Quality factor
///
/// # Returns
/// Bandwidth in Hz
pub fn bandwidth(resonant_freq: f64, q_factor: f64) -> f64 {
    resonant_freq / q_factor
}

/// Calculate current in series RLC circuit given voltage
///
/// I = V / Z
pub fn series_rlc_current(voltage: Complex64, impedance: Complex64) -> Complex64 {
    voltage / impedance
}

/// Calculate voltage across component given current and impedance
///
/// V = I * Z
pub fn component_voltage(current: Complex64, impedance: Complex64) -> Complex64 {
    current * impedance
}

/// Calculate time constant for RC circuit: τ = RC
///
/// # Arguments
/// * `resistance` - Resistance in ohms (Ω)
/// * `capacitance` - Capacitance in farads (F)
///
/// # Returns
/// Time constant in seconds (s)
pub fn rc_time_constant(resistance: f64, capacitance: f64) -> f64 {
    resistance * capacitance
}

/// Calculate time constant for RL circuit: τ = L/R
///
/// # Arguments
/// * `inductance` - Inductance in henries (H)
/// * `resistance` - Resistance in ohms (Ω)
///
/// # Returns
/// Time constant in seconds (s)
pub fn rl_time_constant(inductance: f64, resistance: f64) -> f64 {
    inductance / resistance
}

/// Calculate cutoff frequency for RC low-pass filter: f_c = 1 / (2πRC)
///
/// # Arguments
/// * `resistance` - Resistance in ohms (Ω)
/// * `capacitance` - Capacitance in farads (F)
///
/// # Returns
/// Cutoff frequency (-3dB point) in Hz
pub fn rc_cutoff_frequency(resistance: f64, capacitance: f64) -> f64 {
    1.0 / (2.0 * PI * resistance * capacitance)
}

/// Calculate cutoff frequency for RL circuit: f_c = R / (2πL)
///
/// # Arguments
/// * `resistance` - Resistance in ohms (Ω)
/// * `inductance` - Inductance in henries (H)
///
/// # Returns
/// Cutoff frequency in Hz
pub fn rl_cutoff_frequency(resistance: f64, inductance: f64) -> f64 {
    resistance / (2.0 * PI * inductance)
}

/// Calculate transfer function magnitude for RC low-pass filter
///
/// H(f) = 1 / √(1 + (f/f_c)²)
///
/// # Arguments
/// * `frequency` - Signal frequency in Hz
/// * `cutoff_freq` - Cutoff frequency in Hz
///
/// # Returns
/// Transfer function magnitude (gain, 0 to 1)
pub fn lowpass_transfer_function(frequency: f64, cutoff_freq: f64) -> f64 {
    1.0 / (1.0 + (frequency / cutoff_freq).powi(2)).sqrt()
}

/// Calculate transfer function magnitude for RC high-pass filter
///
/// H(f) = (f/f_c) / √(1 + (f/f_c)²)
pub fn highpass_transfer_function(frequency: f64, cutoff_freq: f64) -> f64 {
    let ratio = frequency / cutoff_freq;
    ratio / (1.0 + ratio.powi(2)).sqrt()
}

/// Convert transfer function magnitude to decibels
///
/// dB = 20 * log₁₀(|H|)
pub fn magnitude_to_db(magnitude: f64) -> f64 {
    20.0 * magnitude.log10()
}

/// Convert decibels to magnitude
///
/// |H| = 10^(dB/20)
pub fn db_to_magnitude(db: f64) -> f64 {
    10_f64.powf(db / 20.0)
}

/// Calculate damping ratio for RLC circuit: ζ = R/(2√(L/C))
///
/// ζ < 1: Underdamped (oscillatory)
/// ζ = 1: Critically damped
/// ζ > 1: Overdamped
///
/// # Arguments
/// * `resistance` - Resistance in ohms (Ω)
/// * `inductance` - Inductance in henries (H)
/// * `capacitance` - Capacitance in farads (F)
///
/// # Returns
/// Damping ratio (dimensionless)
pub fn damping_ratio(resistance: f64, inductance: f64, capacitance: f64) -> f64 {
    resistance / (2.0 * (inductance / capacitance).sqrt())
}

/// Calculate natural frequency: ω_n = 1/√(LC)
///
/// # Returns
/// Natural frequency in rad/s
pub fn natural_frequency(inductance: f64, capacitance: f64) -> f64 {
    1.0 / (inductance * capacitance).sqrt()
}

/// Calculate admittance (inverse of impedance): Y = 1/Z
///
/// # Returns
/// Admittance in siemens (S)
pub fn admittance(impedance: Complex64) -> Complex64 {
    Complex64::new(1.0, 0.0) / impedance
}

/// Calculate conductance (real part of admittance)
pub fn conductance(impedance: Complex64) -> f64 {
    admittance(impedance).re
}

/// Calculate susceptance (imaginary part of admittance)
pub fn susceptance(impedance: Complex64) -> f64 {
    admittance(impedance).im
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rms_conversions() {
        assert!((peak_to_rms(10.0) - 7.071).abs() < 0.001);
        assert!((rms_to_peak(7.071) - 10.0).abs() < 0.01);
    }

    #[test]
    fn test_reactance() {
        let xc = capacitive_reactance(60.0, 10e-6);
        assert!((xc - 265.26).abs() < 0.1);

        let xl = inductive_reactance(60.0, 0.1);
        assert!((xl - 37.699).abs() < 0.01);
    }

    #[test]
    fn test_impedance() {
        let z_r = resistor_impedance(100.0);
        assert_eq!(z_r.re, 100.0);
        assert_eq!(z_r.im, 0.0);

        let z_c = capacitor_impedance(60.0, 10e-6);
        assert!(z_c.re.abs() < 0.01);
        assert!(z_c.im < 0.0); // Capacitive reactance is negative
    }

    #[test]
    fn test_resonance() {
        let f0 = resonant_frequency(0.1, 10e-6);
        assert!((f0 - 159.15).abs() < 0.1);

        let q = quality_factor(10.0, 0.1, 10e-6);
        assert!(q > 0.0);
    }

    #[test]
    fn test_time_constants() {
        let tau_rc = rc_time_constant(1000.0, 10e-6);
        assert_eq!(tau_rc, 0.01); // 10ms

        let tau_rl = rl_time_constant(0.1, 10.0);
        assert_eq!(tau_rl, 0.01); // 10ms
    }

    #[test]
    fn test_transfer_functions() {
        let cutoff = 1000.0;
        let h_low = lowpass_transfer_function(cutoff, cutoff);
        assert!((h_low - 0.707).abs() < 0.01); // -3dB point

        let h_high = highpass_transfer_function(cutoff, cutoff);
        assert!((h_high - 0.707).abs() < 0.01);
    }

    #[test]
    fn test_db_conversions() {
        assert!((magnitude_to_db(1.0) - 0.0).abs() < 0.001);
        assert!((magnitude_to_db(0.707) - (-3.0)).abs() < 0.1);
        assert!((db_to_magnitude(0.0) - 1.0).abs() < 0.001);
    }
}
