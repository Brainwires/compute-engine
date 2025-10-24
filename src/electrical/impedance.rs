//! Impedance Calculations
//!
//! Advanced impedance manipulation and conversion utilities

use num_complex::Complex64;
use std::f64::consts::PI;

/// Convert rectangular form (R + jX) to polar form (|Z|∠θ)
///
/// # Returns
/// (magnitude, phase in radians)
pub fn rectangular_to_polar(z: Complex64) -> (f64, f64) {
    (z.norm(), z.arg())
}

/// Convert polar form (|Z|∠θ) to rectangular form (R + jX)
///
/// # Arguments
/// * `magnitude` - Impedance magnitude in ohms
/// * `phase` - Phase angle in radians
pub fn polar_to_rectangular(magnitude: f64, phase: f64) -> Complex64 {
    Complex64::from_polar(magnitude, phase)
}

/// Calculate series combination of multiple impedances
pub fn series_impedances(impedances: &[Complex64]) -> Complex64 {
    impedances.iter().sum()
}

/// Calculate parallel combination of multiple impedances
///
/// 1/Z_total = 1/Z1 + 1/Z2 + ... + 1/Zn
pub fn parallel_impedances(impedances: &[Complex64]) -> Complex64 {
    let sum_admittances: Complex64 = impedances
        .iter()
        .map(|&z| Complex64::new(1.0, 0.0) / z)
        .sum();
    Complex64::new(1.0, 0.0) / sum_admittances
}

/// Calculate reflection coefficient: Γ = (Z_L - Z_0) / (Z_L + Z_0)
///
/// Used in transmission line analysis
///
/// # Arguments
/// * `z_load` - Load impedance
/// * `z_characteristic` - Characteristic impedance of line
pub fn reflection_coefficient(z_load: Complex64, z_characteristic: Complex64) -> Complex64 {
    (z_load - z_characteristic) / (z_load + z_characteristic)
}

/// Calculate VSWR (Voltage Standing Wave Ratio) from reflection coefficient
///
/// VSWR = (1 + |Γ|) / (1 - |Γ|)
pub fn vswr_from_reflection(gamma: Complex64) -> f64 {
    let gamma_mag = gamma.norm();
    (1.0 + gamma_mag) / (1.0 - gamma_mag)
}

/// Calculate return loss from reflection coefficient in dB
///
/// RL = -20 * log₁₀(|Γ|)
pub fn return_loss(gamma: Complex64) -> f64 {
    -20.0 * gamma.norm().log10()
}

/// Calculate input impedance of transmission line
///
/// Z_in = Z_0 * [(Z_L + jZ_0*tan(βl)) / (Z_0 + jZ_L*tan(βl))]
///
/// # Arguments
/// * `z_0` - Characteristic impedance
/// * `z_load` - Load impedance
/// * `beta_l` - Electrical length (β*l in radians)
pub fn transmission_line_input_impedance(
    z_0: Complex64,
    z_load: Complex64,
    beta_l: f64,
) -> Complex64 {
    let tan_bl = Complex64::new(0.0, beta_l.tan());
    let numerator = z_load + z_0 * tan_bl;
    let denominator = z_0 + z_load * tan_bl;
    z_0 * (numerator / denominator)
}

/// Calculate Smith chart normalized impedance: z = Z / Z_0
pub fn normalize_impedance(z: Complex64, z_characteristic: f64) -> Complex64 {
    z / z_characteristic
}

/// Denormalize impedance: Z = z * Z_0
pub fn denormalize_impedance(z_normalized: Complex64, z_characteristic: f64) -> Complex64 {
    z_normalized * z_characteristic
}

/// Calculate impedance matching using L-network
///
/// Returns (series_component, shunt_component) impedances
/// Positive imaginary = inductive, Negative imaginary = capacitive
pub fn l_network_match(z_source: Complex64, z_load: Complex64) -> (Complex64, Complex64) {
    let r_s = z_source.re;
    let r_l = z_load.re;
    let x_l = z_load.im;

    // Calculate Q factor
    let q = ((r_l / r_s) - 1.0).sqrt();

    // Series reactance
    let x_series = Complex64::new(0.0, q * r_s - x_l);

    // Shunt susceptance (convert to impedance)
    let b_shunt = q / r_l;
    let x_shunt = Complex64::new(0.0, -1.0 / b_shunt);

    (x_series, x_shunt)
}

/// Calculate impedance of quarter-wave transformer
///
/// Z_0 = √(Z_in * Z_load)
pub fn quarter_wave_transformer_impedance(z_in: f64, z_load: f64) -> f64 {
    (z_in * z_load).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rectangular_polar_conversion() {
        let z = Complex64::new(3.0, 4.0);
        let (mag, phase) = rectangular_to_polar(z);
        assert!((mag - 5.0).abs() < 0.001);
        assert!((phase - 0.927).abs() < 0.01); // atan(4/3) ≈ 0.927 rad

        let z2 = polar_to_rectangular(5.0, 0.927);
        assert!((z2.re - 3.0).abs() < 0.01);
        assert!((z2.im - 4.0).abs() < 0.01);
    }

    #[test]
    fn test_series_parallel() {
        let z1 = Complex64::new(10.0, 5.0);
        let z2 = Complex64::new(20.0, -10.0);

        let z_series = series_impedances(&[z1, z2]);
        assert_eq!(z_series.re, 30.0);
        assert_eq!(z_series.im, -5.0);

        let z_parallel = parallel_impedances(&[
            Complex64::new(10.0, 0.0),
            Complex64::new(10.0, 0.0),
        ]);
        assert!((z_parallel.re - 5.0).abs() < 0.01);
    }

    #[test]
    fn test_reflection_coefficient() {
        let z_load = Complex64::new(100.0, 0.0);
        let z_0 = Complex64::new(50.0, 0.0);
        let gamma = reflection_coefficient(z_load, z_0);

        assert!((gamma.re - 0.333).abs() < 0.01);
        assert!(gamma.im.abs() < 0.01);

        let vswr = vswr_from_reflection(gamma);
        assert!((vswr - 2.0).abs() < 0.01);
    }

    #[test]
    fn test_quarter_wave_transformer() {
        let z_0 = quarter_wave_transformer_impedance(50.0, 75.0);
        assert!((z_0 - 61.237).abs() < 0.01);
    }
}
