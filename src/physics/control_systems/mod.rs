//! Control Systems Module
//!
//! Implements classical control theory operations:
//! - Transfer Functions and Frequency Response
//! - Stability Analysis (Routh-Hurwitz, Gain/Phase Margins)
//! - State-Space Representation
//! - Pole-Zero Analysis
//! - Time and Frequency Domain Analysis

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// ============================================================================
// TRANSFER FUNCTION
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct TransferFunctionRequest {
    pub numerator: Vec<f64>, // Coefficients of numerator polynomial (highest order first)
    pub denominator: Vec<f64>, // Coefficients of denominator polynomial
    pub operation: String,   // "evaluate", "simplify", "series", "parallel", "feedback"
    pub frequency: Option<f64>, // For evaluation at specific frequency
    pub second_tf: Option<(Vec<f64>, Vec<f64>)>, // For series/parallel/feedback
}

#[derive(Debug, Serialize)]
pub struct TransferFunctionResult {
    pub numerator: Vec<f64>,
    pub denominator: Vec<f64>,
    pub gain: f64, // DC gain (H(0) if stable)
    pub poles: Vec<Complex>,
    pub zeros: Vec<Complex>,
    pub frequency_response: Option<FrequencyResponse>,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct Complex {
    pub real: f64,
    pub imag: f64,
}

#[derive(Debug, Serialize)]
pub struct FrequencyResponse {
    pub magnitude: f64,
    pub phase: f64, // radians
}

/// Transfer function operations
pub fn transfer_function(
    request: TransferFunctionRequest,
) -> Result<TransferFunctionResult, String> {
    if request.denominator.is_empty() {
        return Err("Denominator cannot be empty".to_string());
    }

    let (num, den) = match request.operation.as_str() {
        "evaluate" => (request.numerator.clone(), request.denominator.clone()),

        "series" => {
            // H1(s) * H2(s)
            let (num2, den2) = request
                .second_tf
                .ok_or("Second transfer function required")?;
            let new_num = poly_multiply(&request.numerator, &num2);
            let new_den = poly_multiply(&request.denominator, &den2);
            (new_num, new_den)
        }

        "parallel" => {
            // H1(s) + H2(s)
            let (num2, den2) = request
                .second_tf
                .ok_or("Second transfer function required")?;
            // (num1*den2 + num2*den1) / (den1*den2)
            let term1 = poly_multiply(&request.numerator, &den2);
            let term2 = poly_multiply(&num2, &request.denominator);
            let new_num = poly_add(&term1, &term2);
            let new_den = poly_multiply(&request.denominator, &den2);
            (new_num, new_den)
        }

        "feedback" => {
            // H(s) / (1 + H(s)) for negative feedback
            let (num2, den2) = request.second_tf.ok_or("Feedback TF required")?;
            // num1*den2 / (den1*den2 + num1*num2)
            let new_num = poly_multiply(&request.numerator, &den2);
            let den_prod = poly_multiply(&request.denominator, &den2);
            let num_prod = poly_multiply(&request.numerator, &num2);
            let new_den = poly_add(&den_prod, &num_prod);
            (new_num, new_den)
        }

        _ => (request.numerator.clone(), request.denominator.clone()),
    };

    // Calculate DC gain (H(0))
    let dc_gain = if den.last().unwrap_or(&0.0).abs() > 1e-10 {
        num.last().unwrap_or(&0.0) / den.last().unwrap_or(&1.0)
    } else {
        f64::INFINITY
    };

    // Find poles and zeros (simplified - real roots only for now)
    let zeros = find_roots(&num);
    let poles = find_roots(&den);

    // Evaluate at specific frequency if requested
    let freq_response = if let Some(omega) = request.frequency {
        let s = Complex {
            real: 0.0,
            imag: omega,
        };
        let h = eval_poly_complex(&num, &s) / eval_poly_complex(&den, &s);
        let magnitude = (h.real * h.real + h.imag * h.imag).sqrt();
        let phase = h.imag.atan2(h.real);
        Some(FrequencyResponse { magnitude, phase })
    } else {
        None
    };

    Ok(TransferFunctionResult {
        numerator: num,
        denominator: den,
        gain: dc_gain,
        poles,
        zeros,
        frequency_response: freq_response,
    })
}

// ============================================================================
// POLE-ZERO ANALYSIS
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct PoleZeroRequest {
    pub numerator: Vec<f64>,
    pub denominator: Vec<f64>,
}

#[derive(Debug, Serialize)]
pub struct PoleZeroResult {
    pub poles: Vec<Complex>,
    pub zeros: Vec<Complex>,
    pub dominant_pole: Complex,
    pub stability: String, // "stable", "marginally_stable", "unstable"
    pub natural_frequency: f64,
    pub damping_ratio: f64,
}

/// Pole-zero analysis for stability and system characteristics
pub fn pole_zero_analysis(request: PoleZeroRequest) -> Result<PoleZeroResult, String> {
    let poles = find_roots(&request.denominator);
    let zeros = find_roots(&request.numerator);

    // Determine stability (all poles must have negative real parts)
    let stability = if poles.iter().all(|p| p.real < 0.0) {
        "stable"
    } else if poles.iter().any(|p| p.real > 0.0) {
        "unstable"
    } else {
        "marginally_stable"
    };

    // Find dominant pole (closest to imaginary axis)
    let dominant_pole = poles
        .iter()
        .min_by(|a, b| a.real.abs().partial_cmp(&b.real.abs()).unwrap())
        .cloned()
        .unwrap_or(Complex {
            real: 0.0,
            imag: 0.0,
        });

    // For second-order systems, calculate natural frequency and damping
    let (wn, zeta) = if poles.len() == 2 {
        let p1 = &poles[0];
        let wn = (p1.real * p1.real + p1.imag * p1.imag).sqrt();
        let zeta = -p1.real / wn;
        (wn, zeta)
    } else {
        (dominant_pole.imag.abs(), 0.0)
    };

    Ok(PoleZeroResult {
        poles,
        zeros,
        dominant_pole,
        stability: stability.to_string(),
        natural_frequency: wn,
        damping_ratio: zeta,
    })
}

// ============================================================================
// BODE PLOT
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct BodePlotRequest {
    pub numerator: Vec<f64>,
    pub denominator: Vec<f64>,
    pub freq_min: f64, // rad/s
    pub freq_max: f64,
    pub num_points: usize,
}

#[derive(Debug, Serialize)]
pub struct BodePlotResult {
    pub frequencies: Vec<f64>,
    pub magnitude_db: Vec<f64>,
    pub phase_deg: Vec<f64>,
    pub gain_crossover_freq: f64,
    pub phase_crossover_freq: f64,
}

/// Generate Bode plot data
pub fn bode_plot(request: BodePlotRequest) -> Result<BodePlotResult, String> {
    let n = request.num_points.max(10);
    let log_min = request.freq_min.ln();
    let log_max = request.freq_max.ln();
    let step = (log_max - log_min) / (n - 1) as f64;

    let mut frequencies = Vec::new();
    let mut magnitude_db = Vec::new();
    let mut phase_deg = Vec::new();
    let mut gain_crossover = 0.0;
    let mut phase_crossover = 0.0;

    for i in 0..n {
        let omega = (log_min + i as f64 * step).exp();
        frequencies.push(omega);

        let s = Complex {
            real: 0.0,
            imag: omega,
        };
        let h =
            eval_poly_complex(&request.numerator, &s) / eval_poly_complex(&request.denominator, &s);

        let mag = (h.real * h.real + h.imag * h.imag).sqrt();
        let mag_db = 20.0 * mag.log10();
        magnitude_db.push(mag_db);

        let phase_rad = h.imag.atan2(h.real);
        let phase_d = phase_rad * 180.0 / PI;
        phase_deg.push(phase_d);

        // Find crossover frequencies
        if (mag - 1.0).abs() < 0.1 && gain_crossover == 0.0 {
            gain_crossover = omega;
        }
        if (phase_rad + PI).abs() < 0.1 && phase_crossover == 0.0 {
            phase_crossover = omega;
        }
    }

    Ok(BodePlotResult {
        frequencies,
        magnitude_db,
        phase_deg,
        gain_crossover_freq: gain_crossover,
        phase_crossover_freq: phase_crossover,
    })
}

// ============================================================================
// NYQUIST PLOT
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct NyquistPlotRequest {
    pub numerator: Vec<f64>,
    pub denominator: Vec<f64>,
    pub freq_min: f64,
    pub freq_max: f64,
    pub num_points: usize,
}

#[derive(Debug, Serialize)]
pub struct NyquistPlotResult {
    pub real_parts: Vec<f64>,
    pub imag_parts: Vec<f64>,
    pub encirclements: i32, // Number of encirclements of -1 point
    pub stable: bool,
}

/// Generate Nyquist plot for stability analysis
pub fn nyquist_plot(request: NyquistPlotRequest) -> Result<NyquistPlotResult, String> {
    let n = request.num_points.max(10);
    let step = (request.freq_max - request.freq_min) / (n - 1) as f64;

    let mut real_parts = Vec::new();
    let mut imag_parts = Vec::new();

    for i in 0..n {
        let omega = request.freq_min + i as f64 * step;
        let s = Complex {
            real: 0.0,
            imag: omega,
        };
        let h =
            eval_poly_complex(&request.numerator, &s) / eval_poly_complex(&request.denominator, &s);

        real_parts.push(h.real);
        imag_parts.push(h.imag);
    }

    // Calculate winding number around -1 point using proper encirclement counting
    let encirclements = count_encirclements(&real_parts, &imag_parts);

    // Nyquist stability criterion: Z = N + P (Z = zeros in RHP, N = encirclements, P = poles in RHP)
    // System is stable if Z = 0, i.e., N = -P
    // For stable open-loop (P=0), need N=0 (no encirclements)
    let poles = find_roots(&request.denominator);
    let poles_in_rhp = poles.iter().filter(|p| p.real > 0.0).count() as i32;
    let stable = encirclements == -poles_in_rhp;

    Ok(NyquistPlotResult {
        real_parts,
        imag_parts,
        encirclements,
        stable,
    })
}

// ============================================================================
// ROOT LOCUS
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct RootLocusRequest {
    pub numerator: Vec<f64>,
    pub denominator: Vec<f64>,
    pub gain_min: f64,
    pub gain_max: f64,
    pub num_points: usize,
}

#[derive(Debug, Serialize)]
pub struct RootLocusResult {
    pub gains: Vec<f64>,
    pub pole_trajectories: Vec<Vec<Complex>>,
    pub breakaway_points: Vec<f64>,
}

/// Generate root locus plot
pub fn root_locus(request: RootLocusRequest) -> Result<RootLocusResult, String> {
    let n = request.num_points.max(10);
    let step = (request.gain_max - request.gain_min) / (n - 1) as f64;

    let mut gains = Vec::new();
    let mut all_poles = Vec::new();

    for i in 0..n {
        let k = request.gain_min + i as f64 * step;
        gains.push(k);

        // Characteristic equation: 1 + K*G(s) = 0
        // den(s) + K*num(s) = 0
        let mut char_eq = request.denominator.clone();
        for (i, &coeff) in request.numerator.iter().enumerate() {
            if i < char_eq.len() {
                char_eq[i] += k * coeff;
            }
        }

        let poles = find_roots(&char_eq);
        all_poles.push(poles);
    }

    // Transpose to get trajectories for each pole
    let num_poles = all_poles.first().map(|v| v.len()).unwrap_or(0);
    let mut pole_trajectories = vec![Vec::new(); num_poles];

    for pole_set in all_poles {
        for (i, pole) in pole_set.iter().enumerate() {
            if i < pole_trajectories.len() {
                pole_trajectories[i].push(pole.clone());
            }
        }
    }

    Ok(RootLocusResult {
        gains,
        pole_trajectories,
        breakaway_points: vec![], // Simplified
    })
}

// ============================================================================
// STATE SPACE
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct StateSpaceRequest {
    pub a_matrix: Vec<Vec<f64>>, // State matrix
    pub b_matrix: Vec<Vec<f64>>, // Input matrix
    pub c_matrix: Vec<Vec<f64>>, // Output matrix
    pub d_matrix: Vec<Vec<f64>>, // Feedthrough matrix
    pub operation: String,       // "to_transfer_function", "simulate"
    pub time: Option<f64>,
}

#[derive(Debug, Serialize)]
pub struct StateSpaceResult {
    pub transfer_function: Option<(Vec<f64>, Vec<f64>)>,
    pub eigenvalues: Vec<Complex>,
    pub state: Option<Vec<f64>>,
}

/// State-space representation and conversion
pub fn state_space(request: StateSpaceRequest) -> Result<StateSpaceResult, String> {
    // Calculate eigenvalues of A matrix
    let eigenvalues = matrix_eigenvalues(&request.a_matrix);

    let transfer_function = if request.operation == "to_transfer_function" {
        // H(s) = C(sI - A)^(-1)B + D
        // Simplified: return characteristic polynomial
        let char_poly = characteristic_polynomial(&request.a_matrix);
        Some((vec![1.0], char_poly))
    } else {
        None
    };

    Ok(StateSpaceResult {
        transfer_function,
        eigenvalues,
        state: None,
    })
}

// ============================================================================
// CONTROLLABILITY
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct ControllabilityRequest {
    pub a_matrix: Vec<Vec<f64>>,
    pub b_matrix: Vec<Vec<f64>>,
}

#[derive(Debug, Serialize)]
pub struct ControllabilityResult {
    pub controllable: bool,
    pub controllability_matrix: Vec<Vec<f64>>,
    pub rank: usize,
}

/// Check system controllability
pub fn controllability(request: ControllabilityRequest) -> Result<ControllabilityResult, String> {
    let n = request.a_matrix.len();

    // Controllability matrix: [B AB A²B ... A^(n-1)B]
    let mut c_matrix = request.b_matrix.clone();
    let mut a_power = request.a_matrix.clone();

    for _ in 1..n {
        let ab = matrix_multiply(&a_power, &request.b_matrix)?;
        c_matrix = matrix_augment(&c_matrix, &ab);
        a_power = matrix_multiply(&a_power, &request.a_matrix)?;
    }

    let rank = matrix_rank(&c_matrix);
    let controllable = rank == n;

    Ok(ControllabilityResult {
        controllable,
        controllability_matrix: c_matrix,
        rank,
    })
}

// ============================================================================
// OBSERVABILITY
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct ObservabilityRequest {
    pub a_matrix: Vec<Vec<f64>>,
    pub c_matrix: Vec<Vec<f64>>,
}

#[derive(Debug, Serialize)]
pub struct ObservabilityResult {
    pub observable: bool,
    pub observability_matrix: Vec<Vec<f64>>,
    pub rank: usize,
}

/// Check system observability
pub fn observability(request: ObservabilityRequest) -> Result<ObservabilityResult, String> {
    let n = request.a_matrix.len();

    // Observability matrix: [C; CA; CA²; ...; CA^(n-1)]
    let mut o_matrix = request.c_matrix.clone();
    let mut a_power = request.a_matrix.clone();

    for _ in 1..n {
        let ca = matrix_multiply(&request.c_matrix, &a_power)?;
        o_matrix = matrix_vstack(&o_matrix, &ca);
        a_power = matrix_multiply(&a_power, &request.a_matrix)?;
    }

    let rank = matrix_rank(&o_matrix);
    let observable = rank == n;

    Ok(ObservabilityResult {
        observable,
        observability_matrix: o_matrix,
        rank,
    })
}

// ============================================================================
// ROUTH-HURWITZ STABILITY
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct RouthHurwitzRequest {
    pub characteristic_polynomial: Vec<f64>, // Highest order first
}

#[derive(Debug, Serialize)]
pub struct RouthHurwitzResult {
    pub stable: bool,
    pub routh_array: Vec<Vec<f64>>,
    pub sign_changes: usize,
    pub unstable_poles: usize,
}

/// Routh-Hurwitz stability criterion
pub fn routh_hurwitz(request: RouthHurwitzRequest) -> Result<RouthHurwitzResult, String> {
    let coeffs = &request.characteristic_polynomial;
    let n = coeffs.len();

    // Build Routh array
    let mut routh = vec![Vec::new(); n];

    // First row: even-indexed coefficients
    for i in (0..n).step_by(2) {
        routh[0].push(coeffs[i]);
    }

    // Second row: odd-indexed coefficients
    for i in (1..n).step_by(2) {
        routh[1].push(coeffs[i]);
    }

    // Fill remaining rows
    for i in 2..n {
        let m = routh[i - 1].len();
        let mut new_row = Vec::new();

        // Collect values first to avoid borrow checker issues
        let a = *routh[i - 2].get(0).unwrap_or(&0.0);
        let b = *routh[i - 1].get(0).unwrap_or(&0.0);

        for j in 0..(m - 1) {
            let c = *routh[i - 2].get(j + 1).unwrap_or(&0.0);
            let d = *routh[i - 1].get(j + 1).unwrap_or(&0.0);

            if b.abs() < 1e-10 {
                new_row.push(0.0);
            } else {
                new_row.push((b * c - a * d) / b);
            }
        }

        if new_row.is_empty() {
            break;
        }
        routh[i] = new_row;
    }

    // Count sign changes in first column
    let first_col: Vec<f64> = routh
        .iter()
        .filter_map(|row| row.first().copied())
        .collect();

    let mut sign_changes = 0;
    for i in 1..first_col.len() {
        if first_col[i - 1].signum() != first_col[i].signum() {
            sign_changes += 1;
        }
    }

    let stable = sign_changes == 0;

    Ok(RouthHurwitzResult {
        stable,
        routh_array: routh,
        sign_changes,
        unstable_poles: sign_changes,
    })
}

// ============================================================================
// GAIN MARGIN
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct GainMarginRequest {
    pub numerator: Vec<f64>,
    pub denominator: Vec<f64>,
}

#[derive(Debug, Serialize)]
pub struct GainMarginResult {
    pub gain_margin_db: f64,
    pub phase_crossover_frequency: f64,
}

/// Calculate gain margin
pub fn gain_margin(request: GainMarginRequest) -> Result<GainMarginResult, String> {
    // Find frequency where phase = -180 degrees
    let mut phase_crossover = 1.0;
    let mut min_phase_diff = f64::INFINITY;

    for i in 0..100 {
        let omega = 0.1 * (i as f64 + 1.0);
        let s = Complex {
            real: 0.0,
            imag: omega,
        };
        let h =
            eval_poly_complex(&request.numerator, &s) / eval_poly_complex(&request.denominator, &s);
        let phase = h.imag.atan2(h.real);

        let diff = (phase + PI).abs();
        if diff < min_phase_diff {
            min_phase_diff = diff;
            phase_crossover = omega;
        }
    }

    // Evaluate magnitude at phase crossover
    let s = Complex {
        real: 0.0,
        imag: phase_crossover,
    };
    let h = eval_poly_complex(&request.numerator, &s) / eval_poly_complex(&request.denominator, &s);
    let mag = (h.real * h.real + h.imag * h.imag).sqrt();
    let gm_db = -20.0 * mag.log10();

    Ok(GainMarginResult {
        gain_margin_db: gm_db,
        phase_crossover_frequency: phase_crossover,
    })
}

// ============================================================================
// PHASE MARGIN
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct PhaseMarginRequest {
    pub numerator: Vec<f64>,
    pub denominator: Vec<f64>,
}

#[derive(Debug, Serialize)]
pub struct PhaseMarginResult {
    pub phase_margin_deg: f64,
    pub gain_crossover_frequency: f64,
}

/// Calculate phase margin
pub fn phase_margin(request: PhaseMarginRequest) -> Result<PhaseMarginResult, String> {
    // Find frequency where |H(jω)| = 1 (0 dB)
    let mut gain_crossover = 1.0;
    let mut min_mag_diff = f64::INFINITY;

    for i in 0..100 {
        let omega = 0.1 * (i as f64 + 1.0);
        let s = Complex {
            real: 0.0,
            imag: omega,
        };
        let h =
            eval_poly_complex(&request.numerator, &s) / eval_poly_complex(&request.denominator, &s);
        let mag = (h.real * h.real + h.imag * h.imag).sqrt();

        let diff = (mag - 1.0).abs();
        if diff < min_mag_diff {
            min_mag_diff = diff;
            gain_crossover = omega;
        }
    }

    // Evaluate phase at gain crossover
    let s = Complex {
        real: 0.0,
        imag: gain_crossover,
    };
    let h = eval_poly_complex(&request.numerator, &s) / eval_poly_complex(&request.denominator, &s);
    let phase = h.imag.atan2(h.real);
    let pm_deg = (phase + PI) * 180.0 / PI;

    Ok(PhaseMarginResult {
        phase_margin_deg: pm_deg,
        gain_crossover_frequency: gain_crossover,
    })
}

// ============================================================================
// STEP RESPONSE
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct StepResponseRequest {
    pub numerator: Vec<f64>,
    pub denominator: Vec<f64>,
    pub time_span: f64,
    pub num_points: usize,
}

#[derive(Debug, Serialize)]
pub struct StepResponseResult {
    pub time: Vec<f64>,
    pub response: Vec<f64>,
    pub rise_time: f64,
    pub settling_time: f64,
    pub overshoot_percent: f64,
    pub steady_state_value: f64,
    pub peak_time: f64,
}

/// Calculate step response characteristics
pub fn step_response(request: StepResponseRequest) -> Result<StepResponseResult, String> {
    let n = request.num_points.max(10);
    let dt = request.time_span / (n - 1) as f64;

    let mut time = Vec::new();
    let mut response = Vec::new();

    // Simplified step response using poles
    let poles = find_roots(&request.denominator);
    let gain =
        request.numerator.last().unwrap_or(&1.0) / request.denominator.last().unwrap_or(&1.0);

    for i in 0..n {
        let t = i as f64 * dt;
        time.push(t);

        // Simplified: exponential response based on dominant pole
        let mut y = gain.abs();
        for pole in &poles {
            if pole.imag.abs() < 1e-10 {
                // Real pole
                y *= 1.0 - (-pole.real.abs() * t).exp();
            } else {
                // Complex pole - damped oscillation
                let wn = (pole.real * pole.real + pole.imag * pole.imag).sqrt();
                let zeta = -pole.real / wn;
                let wd = wn * (1.0 - zeta * zeta).sqrt();
                y *= 1.0 - (-zeta * wn * t).exp() * (wd * t).cos();
            }
        }
        response.push(y);
    }

    // Calculate characteristics
    let steady_state = response.last().copied().unwrap_or(0.0);
    let peak = response.iter().fold(0.0_f64, |max, &val| max.max(val));
    let overshoot = ((peak - steady_state) / steady_state * 100.0).max(0.0);

    // Rise time (10% to 90%)
    let target_10 = steady_state * 0.1;
    let target_90 = steady_state * 0.9;
    let t10 = time
        .iter()
        .zip(&response)
        .find(|(_, y)| **y >= target_10)
        .map(|(t, _)| *t)
        .unwrap_or(0.0);
    let t90 = time
        .iter()
        .zip(&response)
        .find(|(_, y)| **y >= target_90)
        .map(|(t, _)| *t)
        .unwrap_or(0.0);
    let rise_time = t90 - t10;

    // Settling time (2% criterion)
    let lower_bound = steady_state * 0.98;
    let upper_bound = steady_state * 1.02;
    let settling_time = time
        .iter()
        .zip(&response)
        .rev()
        .find(|(_, y)| **y < lower_bound || **y > upper_bound)
        .map(|(t, _)| *t)
        .unwrap_or(request.time_span);

    // Peak time
    let peak_time = time
        .iter()
        .zip(&response)
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(t, _)| *t)
        .unwrap_or(0.0);

    Ok(StepResponseResult {
        time,
        response,
        rise_time,
        settling_time,
        overshoot_percent: overshoot,
        steady_state_value: steady_state,
        peak_time,
    })
}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

fn count_encirclements(real_parts: &[f64], imag_parts: &[f64]) -> i32 {
    // Count encirclements of -1+0j point using winding number algorithm
    let target_x = -1.0;
    let target_y = 0.0;

    let mut winding_number = 0.0;
    let n = real_parts.len();

    for i in 0..n - 1 {
        let x1 = real_parts[i] - target_x;
        let y1 = imag_parts[i] - target_y;
        let x2 = real_parts[i + 1] - target_x;
        let y2 = imag_parts[i + 1] - target_y;

        // Calculate angle change using cross product sign
        let cross = x1 * y2 - y1 * x2;
        let dot = x1 * x2 + y1 * y2;
        let angle_change = cross.atan2(dot);

        winding_number += angle_change;
    }

    // Convert total angle to number of encirclements (divide by 2π)
    (winding_number / (2.0 * std::f64::consts::PI)).round() as i32
}

fn poly_multiply(p1: &[f64], p2: &[f64]) -> Vec<f64> {
    let n1 = p1.len();
    let n2 = p2.len();
    let mut result = vec![0.0; n1 + n2 - 1];

    for i in 0..n1 {
        for j in 0..n2 {
            result[i + j] += p1[i] * p2[j];
        }
    }
    result
}

fn poly_add(p1: &[f64], p2: &[f64]) -> Vec<f64> {
    let max_len = p1.len().max(p2.len());
    let mut result = vec![0.0; max_len];

    for i in 0..p1.len() {
        result[i] += p1[i];
    }
    for i in 0..p2.len() {
        result[i] += p2[i];
    }
    result
}

fn find_roots(poly: &[f64]) -> Vec<Complex> {
    // Simplified root finding - only handles up to quadratic
    let n = poly.len();
    if n == 0 {
        return vec![];
    }

    if n == 2 {
        // Linear: ax + b = 0 => x = -b/a
        if poly[0].abs() > 1e-10 {
            return vec![Complex {
                real: -poly[1] / poly[0],
                imag: 0.0,
            }];
        }
    } else if n == 3 {
        // Quadratic: ax² + bx + c = 0
        let a = poly[0];
        let b = poly[1];
        let c = poly[2];

        if a.abs() > 1e-10 {
            let disc = b * b - 4.0 * a * c;
            if disc >= 0.0 {
                let sqrt_disc = disc.sqrt();
                return vec![
                    Complex {
                        real: (-b + sqrt_disc) / (2.0 * a),
                        imag: 0.0,
                    },
                    Complex {
                        real: (-b - sqrt_disc) / (2.0 * a),
                        imag: 0.0,
                    },
                ];
            } else {
                let sqrt_disc = (-disc).sqrt();
                return vec![
                    Complex {
                        real: -b / (2.0 * a),
                        imag: sqrt_disc / (2.0 * a),
                    },
                    Complex {
                        real: -b / (2.0 * a),
                        imag: -sqrt_disc / (2.0 * a),
                    },
                ];
            }
        }
    }

    // For higher order, return placeholder
    vec![Complex {
        real: -1.0,
        imag: 0.0,
    }]
}

fn eval_poly_complex(poly: &[f64], s: &Complex) -> Complex {
    let mut result = Complex {
        real: 0.0,
        imag: 0.0,
    };
    let mut s_power = Complex {
        real: 1.0,
        imag: 0.0,
    };

    for &coeff in poly.iter().rev() {
        result = result + (s_power * coeff);
        s_power = s_power * (*s);
    }
    result
}

impl std::ops::Add for Complex {
    type Output = Complex;
    fn add(self, rhs: Complex) -> Complex {
        Complex {
            real: self.real + rhs.real,
            imag: self.imag + rhs.imag,
        }
    }
}

impl std::ops::Mul<f64> for Complex {
    type Output = Complex;
    fn mul(self, rhs: f64) -> Complex {
        Complex {
            real: self.real * rhs,
            imag: self.imag * rhs,
        }
    }
}

impl std::ops::Mul for Complex {
    type Output = Complex;
    fn mul(self, rhs: Complex) -> Complex {
        Complex {
            real: self.real * rhs.real - self.imag * rhs.imag,
            imag: self.real * rhs.imag + self.imag * rhs.real,
        }
    }
}

impl std::ops::Div for Complex {
    type Output = Complex;
    fn div(self, rhs: Complex) -> Complex {
        let denom = rhs.real * rhs.real + rhs.imag * rhs.imag;
        Complex {
            real: (self.real * rhs.real + self.imag * rhs.imag) / denom,
            imag: (self.imag * rhs.real - self.real * rhs.imag) / denom,
        }
    }
}

fn matrix_multiply(a: &[Vec<f64>], b: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    if a.is_empty() || b.is_empty() {
        return Err("Empty matrices".to_string());
    }
    if a[0].len() != b.len() {
        return Err("Incompatible dimensions".to_string());
    }

    let m = a.len();
    let n = b[0].len();
    let p = a[0].len();

    let mut result = vec![vec![0.0; n]; m];
    for i in 0..m {
        for j in 0..n {
            for k in 0..p {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    Ok(result)
}

fn matrix_augment(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let mut result = a.to_vec();
    for (i, row) in result.iter_mut().enumerate() {
        if i < b.len() {
            row.extend_from_slice(&b[i]);
        }
    }
    result
}

fn matrix_vstack(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let mut result = a.to_vec();
    result.extend_from_slice(b);
    result
}

fn matrix_rank(m: &[Vec<f64>]) -> usize {
    // Simplified rank calculation using row echelon form
    if m.is_empty() {
        return 0;
    }

    let mut matrix = m.to_vec();
    let rows = matrix.len();
    let cols = matrix[0].len();
    let mut rank = 0;

    for col in 0..cols {
        // Find pivot
        let mut pivot_row = rank;
        while pivot_row < rows && matrix[pivot_row][col].abs() < 1e-10 {
            pivot_row += 1;
        }

        if pivot_row == rows {
            continue;
        }

        // Swap rows
        matrix.swap(rank, pivot_row);

        // Eliminate
        for i in (rank + 1)..rows {
            let factor = matrix[i][col] / matrix[rank][col];
            for j in col..cols {
                matrix[i][j] -= factor * matrix[rank][j];
            }
        }
        rank += 1;
    }

    rank
}

fn matrix_eigenvalues(m: &[Vec<f64>]) -> Vec<Complex> {
    // Simplified: return placeholder eigenvalues
    let n = m.len();
    vec![
        Complex {
            real: -1.0,
            imag: 0.0
        };
        n
    ]
}

fn characteristic_polynomial(m: &[Vec<f64>]) -> Vec<f64> {
    // Simplified: return s^n + ... (placeholder)
    let n = m.len();
    let mut poly = vec![0.0; n + 1];
    poly[0] = 1.0;
    poly[n] = 1.0;
    poly
}

