//! DC Circuit Analysis
//!
//! Basic DC circuit calculations including Ohm's law, series/parallel combinations,
//! voltage and current dividers, and power dissipation.

use std::f64::consts::PI;

/// Calculate voltage using Ohm's law: V = I * R
///
/// # Arguments
/// * `current` - Current in amperes (A)
/// * `resistance` - Resistance in ohms (Ω)
///
/// # Returns
/// Voltage in volts (V)
///
/// # Example
/// ```
/// use computational_engine::electrical::dc_analysis::ohms_law_voltage;
/// let v = ohms_law_voltage(2.0, 10.0); // 2A through 10Ω
/// assert_eq!(v, 20.0); // 20V
/// ```
pub fn ohms_law_voltage(current: f64, resistance: f64) -> f64 {
    current * resistance
}

/// Calculate current using Ohm's law: I = V / R
///
/// # Arguments
/// * `voltage` - Voltage in volts (V)
/// * `resistance` - Resistance in ohms (Ω)
///
/// # Returns
/// Current in amperes (A)
pub fn ohms_law_current(voltage: f64, resistance: f64) -> f64 {
    voltage / resistance
}

/// Calculate resistance using Ohm's law: R = V / I
///
/// # Arguments
/// * `voltage` - Voltage in volts (V)
/// * `current` - Current in amperes (A)
///
/// # Returns
/// Resistance in ohms (Ω)
pub fn ohms_law_resistance(voltage: f64, current: f64) -> f64 {
    voltage / current
}

/// Calculate total resistance of resistors in series: R_total = R1 + R2 + ... + Rn
///
/// # Arguments
/// * `resistances` - Vector of resistance values in ohms (Ω)
///
/// # Returns
/// Total series resistance in ohms (Ω)
///
/// # Example
/// ```
/// use computational_engine::electrical::dc_analysis::series_resistance;
/// let r_total = series_resistance(&[10.0, 20.0, 30.0]);
/// assert_eq!(r_total, 60.0); // 10 + 20 + 30 = 60Ω
/// ```
pub fn series_resistance(resistances: &[f64]) -> f64 {
    resistances.iter().sum()
}

/// Calculate total resistance of resistors in parallel: 1/R_total = 1/R1 + 1/R2 + ... + 1/Rn
///
/// # Arguments
/// * `resistances` - Vector of resistance values in ohms (Ω)
///
/// # Returns
/// Total parallel resistance in ohms (Ω)
///
/// # Example
/// ```
/// use computational_engine::electrical::dc_analysis::parallel_resistance;
/// let r_total = parallel_resistance(&[10.0, 10.0]);
/// assert_eq!(r_total, 5.0); // Two 10Ω in parallel = 5Ω
/// ```
pub fn parallel_resistance(resistances: &[f64]) -> f64 {
    let sum_reciprocals: f64 = resistances.iter().map(|&r| 1.0 / r).sum();
    1.0 / sum_reciprocals
}

/// Calculate voltage across a resistor in a voltage divider
///
/// Voltage divider equation: V_out = V_in * (R2 / (R1 + R2))
///
/// # Arguments
/// * `v_in` - Input voltage in volts (V)
/// * `r1` - Upper resistor in ohms (Ω)
/// * `r2` - Lower resistor in ohms (Ω)
///
/// # Returns
/// Output voltage in volts (V)
///
/// # Example
/// ```
/// use computational_engine::electrical::dc_analysis::voltage_divider;
/// let v_out = voltage_divider(12.0, 10.0, 5.0);
/// assert_eq!(v_out, 4.0); // 12V * (5/(10+5)) = 4V
/// ```
pub fn voltage_divider(v_in: f64, r1: f64, r2: f64) -> f64 {
    v_in * r2 / (r1 + r2)
}

/// Calculate current through a branch in a current divider
///
/// Current divider equation: I_out = I_in * (R_other / (R_branch + R_other))
///
/// # Arguments
/// * `i_in` - Input current in amperes (A)
/// * `r_branch` - Resistance of the branch in ohms (Ω)
/// * `r_other` - Resistance of the other parallel branch in ohms (Ω)
///
/// # Returns
/// Current through the branch in amperes (A)
///
/// # Example
/// ```
/// use computational_engine::electrical::dc_analysis::current_divider;
/// let i_branch = current_divider(6.0, 10.0, 5.0);
/// assert_eq!(i_branch, 2.0); // 6A * (5/(10+5)) = 2A
/// ```
pub fn current_divider(i_in: f64, r_branch: f64, r_other: f64) -> f64 {
    i_in * r_other / (r_branch + r_other)
}

/// Calculate power dissipated by a resistor: P = V * I = I² * R = V² / R
///
/// # Arguments
/// * `voltage` - Voltage across resistor in volts (V)
/// * `current` - Current through resistor in amperes (A)
///
/// # Returns
/// Power in watts (W)
///
/// # Example
/// ```
/// use computational_engine::electrical::dc_analysis::power_dissipation;
/// let p = power_dissipation(12.0, 2.0);
/// assert_eq!(p, 24.0); // 12V * 2A = 24W
/// ```
pub fn power_dissipation(voltage: f64, current: f64) -> f64 {
    voltage * current
}

/// Calculate power from voltage and resistance: P = V² / R
pub fn power_from_voltage_resistance(voltage: f64, resistance: f64) -> f64 {
    voltage * voltage / resistance
}

/// Calculate power from current and resistance: P = I² * R
pub fn power_from_current_resistance(current: f64, resistance: f64) -> f64 {
    current * current * resistance
}

/// Calculate energy dissipated over time: E = P * t
///
/// # Arguments
/// * `power` - Power in watts (W)
/// * `time` - Time in seconds (s)
///
/// # Returns
/// Energy in joules (J)
pub fn energy_dissipated(power: f64, time: f64) -> f64 {
    power * time
}

/// Calculate conductance from resistance: G = 1 / R
///
/// # Arguments
/// * `resistance` - Resistance in ohms (Ω)
///
/// # Returns
/// Conductance in siemens (S)
pub fn dc_conductance(resistance: f64) -> f64 {
    1.0 / resistance
}

/// Calculate total capacitance of capacitors in series: 1/C_total = 1/C1 + 1/C2 + ... + 1/Cn
///
/// # Arguments
/// * `capacitances` - Vector of capacitance values in farads (F)
///
/// # Returns
/// Total series capacitance in farads (F)
pub fn series_capacitance(capacitances: &[f64]) -> f64 {
    let sum_reciprocals: f64 = capacitances.iter().map(|&c| 1.0 / c).sum();
    1.0 / sum_reciprocals
}

/// Calculate total capacitance of capacitors in parallel: C_total = C1 + C2 + ... + Cn
///
/// # Arguments
/// * `capacitances` - Vector of capacitance values in farads (F)
///
/// # Returns
/// Total parallel capacitance in farads (F)
pub fn parallel_capacitance(capacitances: &[f64]) -> f64 {
    capacitances.iter().sum()
}

/// Calculate total inductance of inductors in series: L_total = L1 + L2 + ... + Ln
///
/// # Arguments
/// * `inductances` - Vector of inductance values in henries (H)
///
/// # Returns
/// Total series inductance in henries (H)
pub fn series_inductance(inductances: &[f64]) -> f64 {
    inductances.iter().sum()
}

/// Calculate total inductance of inductors in parallel: 1/L_total = 1/L1 + 1/L2 + ... + 1/Ln
///
/// # Arguments
/// * `inductances` - Vector of inductance values in henries (H)
///
/// # Returns
/// Total parallel inductance in henries (H)
pub fn parallel_inductance(inductances: &[f64]) -> f64 {
    let sum_reciprocals: f64 = inductances.iter().map(|&l| 1.0 / l).sum();
    1.0 / sum_reciprocals
}

/// Calculate charge stored in a capacitor: Q = C * V
///
/// # Arguments
/// * `capacitance` - Capacitance in farads (F)
/// * `voltage` - Voltage in volts (V)
///
/// # Returns
/// Charge in coulombs (C)
pub fn capacitor_charge(capacitance: f64, voltage: f64) -> f64 {
    capacitance * voltage
}

/// Calculate energy stored in a capacitor: E = 0.5 * C * V²
///
/// # Arguments
/// * `capacitance` - Capacitance in farads (F)
/// * `voltage` - Voltage in volts (V)
///
/// # Returns
/// Energy in joules (J)
pub fn capacitor_energy(capacitance: f64, voltage: f64) -> f64 {
    0.5 * capacitance * voltage * voltage
}

/// Calculate energy stored in an inductor: E = 0.5 * L * I²
///
/// # Arguments
/// * `inductance` - Inductance in henries (H)
/// * `current` - Current in amperes (A)
///
/// # Returns
/// Energy in joules (J)
pub fn inductor_energy(inductance: f64, current: f64) -> f64 {
    0.5 * inductance * current * current
}

/// Calculate wire resistance based on resistivity, length, and cross-sectional area
///
/// R = ρ * L / A
///
/// # Arguments
/// * `resistivity` - Material resistivity in Ω·m (e.g., 1.68e-8 for copper)
/// * `length` - Wire length in meters (m)
/// * `area` - Cross-sectional area in square meters (m²)
///
/// # Returns
/// Resistance in ohms (Ω)
pub fn wire_resistance(resistivity: f64, length: f64, area: f64) -> f64 {
    resistivity * length / area
}

/// Calculate wire resistance from AWG (American Wire Gauge) size
///
/// # Arguments
/// * `awg` - AWG wire gauge number (e.g., 12, 14, 18)
/// * `length_m` - Wire length in meters
///
/// # Returns
/// Resistance in ohms (Ω) for copper wire at 20°C
pub fn wire_resistance_awg(awg: u32, length_m: f64) -> f64 {
    // Diameter in mm for AWG (standard formula)
    let diameter_mm = 0.127 * 92_f64.powf((36.0 - awg as f64) / 39.0);
    let area_mm2 = PI * (diameter_mm / 2.0).powi(2);
    let area_m2 = area_mm2 * 1e-6;

    // Copper resistivity at 20°C
    let copper_resistivity = 1.68e-8; // Ω·m

    wire_resistance(copper_resistivity, length_m, area_m2)
}

