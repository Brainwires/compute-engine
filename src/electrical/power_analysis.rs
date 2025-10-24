//! Power Analysis
//!
//! Power calculations for AC circuits including real, reactive, apparent power,
//! and power factor calculations.

use num_complex::Complex64;
use std::f64::consts::PI;

/// Calculate real power (active power): P = V * I * cos(θ)
///
/// # Arguments
/// * `voltage` - RMS voltage in volts
/// * `current` - RMS current in amperes
/// * `power_factor` - Power factor (cos θ)
///
/// # Returns
/// Real power in watts (W)
pub fn real_power(voltage: f64, current: f64, power_factor: f64) -> f64 {
    voltage * current * power_factor
}

/// Calculate reactive power: Q = V * I * sin(θ)
///
/// # Returns
/// Reactive power in volt-amperes reactive (VAR)
pub fn reactive_power(voltage: f64, current: f64, power_factor: f64) -> f64 {
    let theta = power_factor.acos();
    voltage * current * theta.sin()
}

/// Calculate apparent power: S = V * I
///
/// # Returns
/// Apparent power in volt-amperes (VA)
pub fn apparent_power(voltage: f64, current: f64) -> f64 {
    voltage * current
}

/// Calculate power factor: PF = P / S = cos(θ)
pub fn power_factor_from_powers(real_power: f64, apparent_power: f64) -> f64 {
    real_power / apparent_power
}

/// Calculate complex power: S = P + jQ
pub fn complex_power(real_power: f64, reactive_power: f64) -> Complex64 {
    Complex64::new(real_power, reactive_power)
}

/// Calculate three-phase real power: P = √3 * V_L * I_L * cos(θ)
///
/// # Arguments
/// * `line_voltage` - Line-to-line voltage in volts
/// * `line_current` - Line current in amperes
/// * `power_factor` - Power factor
pub fn three_phase_power(line_voltage: f64, line_current: f64, power_factor: f64) -> f64 {
    3_f64.sqrt() * line_voltage * line_current * power_factor
}

/// Calculate three-phase apparent power: S = √3 * V_L * I_L
pub fn three_phase_apparent_power(line_voltage: f64, line_current: f64) -> f64 {
    3_f64.sqrt() * line_voltage * line_current
}

/// Convert line voltage to phase voltage (wye connection): V_phase = V_line / √3
pub fn line_to_phase_voltage(line_voltage: f64) -> f64 {
    line_voltage / 3_f64.sqrt()
}

/// Convert phase voltage to line voltage (wye connection): V_line = √3 * V_phase
pub fn phase_to_line_voltage(phase_voltage: f64) -> f64 {
    3_f64.sqrt() * phase_voltage
}

/// Calculate capacitor value for power factor correction
///
/// # Arguments
/// * `power_kw` - Real power in kilowatts
/// * `voltage` - System voltage
/// * `frequency` - System frequency (Hz)
/// * `pf_current` - Current power factor
/// * `pf_target` - Target power factor
///
/// # Returns
/// Required capacitance in farads
pub fn power_factor_correction_capacitor(
    power_kw: f64,
    voltage: f64,
    frequency: f64,
    pf_current: f64,
    pf_target: f64,
) -> f64 {
    let theta_current = pf_current.acos();
    let theta_target = pf_target.acos();
    let q_reduction = power_kw * (theta_current.tan() - theta_target.tan());
    let q_var = q_reduction * 1000.0; // Convert to VAR
    q_var / (2.0 * PI * frequency * voltage * voltage)
}

/// Calculate efficiency: η = P_out / P_in
pub fn efficiency(power_out: f64, power_in: f64) -> f64 {
    power_out / power_in
}

/// Calculate power loss: P_loss = P_in - P_out
pub fn power_loss(power_in: f64, power_out: f64) -> f64 {
    power_in - power_out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_power_calculations() {
        let p = real_power(120.0, 10.0, 0.8);
        assert_eq!(p, 960.0); // 960W

        let s = apparent_power(120.0, 10.0);
        assert_eq!(s, 1200.0); // 1200VA

        let pf = power_factor_from_powers(960.0, 1200.0);
        assert_eq!(pf, 0.8);
    }

    #[test]
    fn test_three_phase() {
        let p = three_phase_power(480.0, 100.0, 1.0);
        assert!((p - 83138.0).abs() < 1.0);
    }

    #[test]
    fn test_voltage_conversions() {
        let v_phase = line_to_phase_voltage(480.0);
        assert!((v_phase - 277.0).abs() < 1.0);

        let v_line = phase_to_line_voltage(277.0);
        assert!((v_line - 480.0).abs() < 1.0);
    }
}
