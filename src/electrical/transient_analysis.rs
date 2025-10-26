//! Transient Analysis
//!
//! Time-domain analysis of RC, RL, and RLC circuits

use std::f64::consts::E;

/// Calculate capacitor voltage during charging: V_C(t) = V_S(1 - e^(-t/τ))
///
/// # Arguments
/// * `v_source` - Source voltage in volts
/// * `time` - Time in seconds
/// * `time_constant` - Time constant τ = RC in seconds
///
/// # Returns
/// Capacitor voltage at time t
pub fn rc_charging_voltage(v_source: f64, time: f64, time_constant: f64) -> f64 {
    v_source * (1.0 - E.powf(-time / time_constant))
}

/// Calculate capacitor voltage during discharging: V_C(t) = V_0 * e^(-t/τ)
pub fn rc_discharging_voltage(v_initial: f64, time: f64, time_constant: f64) -> f64 {
    v_initial * E.powf(-time / time_constant)
}

/// Calculate current during RC charging: I(t) = (V_S/R) * e^(-t/τ)
pub fn rc_charging_current(v_source: f64, resistance: f64, time: f64, time_constant: f64) -> f64 {
    (v_source / resistance) * E.powf(-time / time_constant)
}

/// Calculate inductor current during RL rise: I_L(t) = (V_S/R)(1 - e^(-Rt/L))
///
/// # Arguments
/// * `v_source` - Source voltage
/// * `resistance` - Resistance in ohms
/// * `time` - Time in seconds
/// * `time_constant` - Time constant τ = L/R in seconds
pub fn rl_rising_current(v_source: f64, resistance: f64, time: f64, time_constant: f64) -> f64 {
    (v_source / resistance) * (1.0 - E.powf(-time / time_constant))
}

/// Calculate inductor current during RL decay: I_L(t) = I_0 * e^(-Rt/L)
pub fn rl_decaying_current(i_initial: f64, time: f64, time_constant: f64) -> f64 {
    i_initial * E.powf(-time / time_constant)
}

/// Calculate voltage across inductor: V_L(t) = V_S * e^(-Rt/L)
pub fn rl_voltage(v_source: f64, time: f64, time_constant: f64) -> f64 {
    v_source * E.powf(-time / time_constant)
}

/// Series RLC transient response for underdamped case (ζ < 1)
///
/// v_C(t) = V_S(1 - e^(-αt)[cos(ω_d*t) + (α/ω_d)sin(ω_d*t)])
///
/// # Arguments
/// * `v_source` - Source voltage
/// * `time` - Time in seconds
/// * `alpha` - Damping coefficient α = R/(2L)
/// * `omega_d` - Damped frequency ω_d = √(ω_n² - α²)
pub fn rlc_underdamped_voltage(v_source: f64, time: f64, alpha: f64, omega_d: f64) -> f64 {
    v_source
        * (1.0
            - E.powf(-alpha * time)
                * ((omega_d * time).cos() + (alpha / omega_d) * (omega_d * time).sin()))
}

/// Series RLC transient response for critically damped case (ζ = 1)
///
/// v_C(t) = V_S(1 - e^(-αt)(1 + αt))
pub fn rlc_critically_damped_voltage(v_source: f64, time: f64, alpha: f64) -> f64 {
    v_source * (1.0 - E.powf(-alpha * time) * (1.0 + alpha * time))
}

/// Series RLC transient response for overdamped case (ζ > 1)
///
/// v_C(t) = V_S(1 - e^(-αt)[cosh(√(α²-ω_n²)t) + (α/√(α²-ω_n²))sinh(√(α²-ω_n²)t)])
pub fn rlc_overdamped_voltage(
    v_source: f64,
    time: f64,
    alpha: f64,
    omega_n: f64,
) -> f64 {
    let beta = (alpha * alpha - omega_n * omega_n).sqrt();
    v_source
        * (1.0
            - E.powf(-alpha * time)
                * ((beta * time).cosh() + (alpha / beta) * (beta * time).sinh()))
}

/// Calculate time to reach a percentage of final value in RC circuit
///
/// t = -τ * ln(1 - percentage/100)
///
/// # Example
/// Time to reach 63.2% (1 time constant): t = τ
/// Time to reach 95%: t ≈ 3τ
/// Time to reach 99%: t ≈ 4.6τ
pub fn rc_settling_time(time_constant: f64, percentage: f64) -> f64 {
    -time_constant * (1.0 - percentage / 100.0).ln()
}

/// Calculate step response rise time (10% to 90%)
///
/// t_r = 2.2 * τ (for first-order system)
pub fn rise_time(time_constant: f64) -> f64 {
    2.2 * time_constant
}

/// Calculate overshoot for underdamped system
///
/// %OS = 100 * e^(-πζ/√(1-ζ²))
///
/// # Arguments
/// * `damping_ratio` - Damping ratio ζ
pub fn percent_overshoot(damping_ratio: f64) -> f64 {
    100.0 * E.powf(-std::f64::consts::PI * damping_ratio / (1.0 - damping_ratio * damping_ratio).sqrt())
}

/// Calculate peak time for underdamped system
///
/// t_p = π / ω_d
pub fn peak_time(omega_d: f64) -> f64 {
    std::f64::consts::PI / omega_d
}

/// Calculate settling time (2% criterion) for second-order system
///
/// t_s = 4 / (ζ * ω_n)
pub fn settling_time(damping_ratio: f64, natural_frequency: f64) -> f64 {
    4.0 / (damping_ratio * natural_frequency)
}

