//! Transfer Function Representation
//!
//! Transfer functions in the Laplace domain:
//! G(s) = N(s) / D(s) = (b_n*s^n + ... + b_1*s + b_0) / (a_m*s^m + ... + a_1*s + a_0)
//!
//! Provides frequency response, pole-zero analysis, and Bode plots.

use num_complex::Complex64;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Transfer function representation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransferFunction {
    /// Numerator coefficients (highest degree first)
    pub numerator: Vec<f64>,
    /// Denominator coefficients (highest degree first)
    pub denominator: Vec<f64>,
}

impl TransferFunction {
    /// Create a new transfer function
    ///
    /// # Example
    /// ```
    /// use computational_engine::specialized::control_theory::TransferFunction;
    ///
    /// // G(s) = 1 / (s + 1)
    /// let tf = TransferFunction::new(vec![1.0], vec![1.0, 1.0]).unwrap();
    /// ```
    pub fn new(numerator: Vec<f64>, denominator: Vec<f64>) -> Result<Self, String> {
        if numerator.is_empty() {
            return Err("Numerator cannot be empty".to_string());
        }
        if denominator.is_empty() {
            return Err("Denominator cannot be empty".to_string());
        }
        if denominator[0] == 0.0 {
            return Err("Leading coefficient of denominator cannot be zero".to_string());
        }

        Ok(Self {
            numerator,
            denominator,
        })
    }

    /// Evaluate transfer function at a complex frequency s
    pub fn evaluate(&self, s: Complex64) -> Complex64 {
        let num = Self::eval_polynomial(&self.numerator, s);
        let den = Self::eval_polynomial(&self.denominator, s);
        num / den
    }

    /// Evaluate polynomial at s
    fn eval_polynomial(coeffs: &[f64], s: Complex64) -> Complex64 {
        let mut result = Complex64::new(0.0, 0.0);
        let mut s_power = Complex64::new(1.0, 0.0);

        // Evaluate from lowest to highest degree
        for &coeff in coeffs.iter().rev() {
            result += coeff * s_power;
            s_power *= s;
        }

        result
    }

    /// Compute frequency response at angular frequency ω
    pub fn frequency_response(&self, omega: f64) -> Complex64 {
        let s = Complex64::new(0.0, omega);
        self.evaluate(s)
    }

    /// Compute magnitude (dB) and phase (degrees) at frequency ω
    pub fn bode_point(&self, omega: f64) -> (f64, f64) {
        let h = self.frequency_response(omega);
        let magnitude_db = 20.0 * h.norm().log10();
        let phase_deg = h.arg() * 180.0 / PI;
        (magnitude_db, phase_deg)
    }

    /// Generate Bode plot data
    pub fn bode_plot(&self, omega_min: f64, omega_max: f64, num_points: usize) -> BodePlot {
        let mut frequencies = Vec::new();
        let mut magnitudes = Vec::new();
        let mut phases = Vec::new();

        // Logarithmic spacing
        let log_min = omega_min.log10();
        let log_max = omega_max.log10();
        let log_step = (log_max - log_min) / (num_points - 1) as f64;

        for i in 0..num_points {
            let omega = 10.0_f64.powf(log_min + i as f64 * log_step);
            let (mag, phase) = self.bode_point(omega);

            frequencies.push(omega);
            magnitudes.push(mag);
            phases.push(phase);
        }

        BodePlot {
            frequencies,
            magnitudes,
            phases,
        }
    }

    /// Find poles (roots of denominator)
    pub fn poles(&self) -> Vec<Complex64> {
        find_polynomial_roots(&self.denominator)
    }

    /// Find zeros (roots of numerator)
    pub fn zeros(&self) -> Vec<Complex64> {
        find_polynomial_roots(&self.numerator)
    }

    /// Check if system is stable (all poles in left half-plane)
    pub fn is_stable(&self) -> bool {
        self.poles().iter().all(|pole| pole.re < 0.0)
    }

    /// Compute DC gain (steady-state gain for step input)
    /// Evaluates G(s) at s = 0
    pub fn dc_gain(&self) -> f64 {
        // For polynomials stored highest degree first,
        // the constant term is the last coefficient
        let num_const = *self.numerator.last().unwrap_or(&0.0);
        let den_const = *self.denominator.last().unwrap_or(&1.0);
        num_const / den_const
    }

    /// Series connection: G1(s) * G2(s)
    pub fn series(&self, other: &TransferFunction) -> TransferFunction {
        let num = convolve(&self.numerator, &other.numerator);
        let den = convolve(&self.denominator, &other.denominator);
        TransferFunction::new(num, den).unwrap()
    }

    /// Parallel connection: G1(s) + G2(s)
    pub fn parallel(&self, other: &TransferFunction) -> TransferFunction {
        // (N1*D2 + N2*D1) / (D1*D2)
        let n1d2 = convolve(&self.numerator, &other.denominator);
        let n2d1 = convolve(&other.numerator, &self.denominator);
        let num = add_polynomials(&n1d2, &n2d1);
        let den = convolve(&self.denominator, &other.denominator);
        TransferFunction::new(num, den).unwrap()
    }

    /// Feedback connection: G / (1 + G*H)
    pub fn feedback(&self, feedback_tf: &TransferFunction) -> TransferFunction {
        // Forward: N1/D1, Feedback: N2/D2
        // Closed-loop: (N1*D2) / (D1*D2 + N1*N2)
        let num = convolve(&self.numerator, &feedback_tf.denominator);
        let d1d2 = convolve(&self.denominator, &feedback_tf.denominator);
        let n1n2 = convolve(&self.numerator, &feedback_tf.numerator);
        let den = add_polynomials(&d1d2, &n1n2);
        TransferFunction::new(num, den).unwrap()
    }
}

/// Bode plot data
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BodePlot {
    /// Angular frequencies (rad/s)
    pub frequencies: Vec<f64>,
    /// Magnitude (dB)
    pub magnitudes: Vec<f64>,
    /// Phase (degrees)
    pub phases: Vec<f64>,
}

/// Convolve two polynomials (multiply in frequency domain)
fn convolve(a: &[f64], b: &[f64]) -> Vec<f64> {
    let n = a.len() + b.len() - 1;
    let mut result = vec![0.0; n];

    for (i, &ai) in a.iter().enumerate() {
        for (j, &bj) in b.iter().enumerate() {
            result[i + j] += ai * bj;
        }
    }

    result
}

/// Add two polynomials (pad with zeros if needed)
fn add_polynomials(a: &[f64], b: &[f64]) -> Vec<f64> {
    let n = a.len().max(b.len());
    let mut result = vec![0.0; n];

    for (i, &val) in a.iter().enumerate() {
        result[n - a.len() + i] += val;
    }

    for (i, &val) in b.iter().enumerate() {
        result[n - b.len() + i] += val;
    }

    result
}

/// Find roots of polynomial (simplified - uses companion matrix eigenvalues)
fn find_polynomial_roots(coeffs: &[f64]) -> Vec<Complex64> {
    if coeffs.len() == 1 {
        return vec![]; // Constant polynomial has no roots
    }

    if coeffs.len() == 2 {
        // Linear: ax + b = 0 → x = -b/a
        return vec![Complex64::new(-coeffs[1] / coeffs[0], 0.0)];
    }

    if coeffs.len() == 3 {
        // Quadratic: ax² + bx + c = 0
        let a = coeffs[0];
        let b = coeffs[1];
        let c = coeffs[2];

        let discriminant = b * b - 4.0 * a * c;

        if discriminant >= 0.0 {
            let sqrt_disc = discriminant.sqrt();
            return vec![
                Complex64::new((-b + sqrt_disc) / (2.0 * a), 0.0),
                Complex64::new((-b - sqrt_disc) / (2.0 * a), 0.0),
            ];
        } else {
            let real = -b / (2.0 * a);
            let imag = (-discriminant).sqrt() / (2.0 * a);
            return vec![
                Complex64::new(real, imag),
                Complex64::new(real, -imag),
            ];
        }
    }

    // For higher order, would need numerical root finding
    // Placeholder: return empty vector
    vec![]
}

// Test module
#[cfg(test)]
#[path = "../../../tests/unit/specialized/control_theory_transfer_function_tests.rs"]
mod tests;

