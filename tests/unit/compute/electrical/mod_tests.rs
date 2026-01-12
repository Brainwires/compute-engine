//! Comprehensive tests for electrical module
//!
//! Tests for AC/DC analysis, impedance, filters, network analysis,
//! power calculations, transient analysis, and NEC code calculations.

use computational_engine::electrical::*;
use num_complex::Complex64;
use std::f64::consts::PI;

use super::*;

// ============================================================================
// AC ANALYSIS TESTS
// ============================================================================

#[test]
fn test_peak_to_rms() {
    let peak = 100.0;
    let rms = ac_analysis::peak_to_rms(peak);
    assert!((rms - 70.71).abs() < 0.01, "RMS should be peak/√2");
}

#[test]
fn test_rms_to_peak() {
    let rms = 120.0;
    let peak = ac_analysis::rms_to_peak(rms);
    assert!((peak - 169.71).abs() < 0.01, "Peak should be RMS*√2");
}

#[test]
fn test_capacitive_reactance() {
    let xc = ac_analysis::capacitive_reactance(60.0, 10e-6); // 60Hz, 10μF
    assert!((xc - 265.26).abs() < 0.01, "Capacitive reactance at 60Hz");
}

#[test]
fn test_inductive_reactance() {
    let xl = ac_analysis::inductive_reactance(60.0, 0.1); // 60Hz, 100mH
    assert!((xl - 37.70).abs() < 0.01, "Inductive reactance at 60Hz");
}

#[test]
fn test_series_rlc_impedance() {
    let z = ac_analysis::series_rlc_impedance(10.0, 0.1, 100e-6, 60.0);
    assert!(z.re > 0.0 && z.im.abs() > 0.0, "Impedance should be complex");
}

#[test]
fn test_parallel_rlc_impedance() {
    let z = ac_analysis::parallel_rlc_impedance(100.0, 0.1, 10e-6, 1000.0);
    assert!(z.norm() > 0.0, "Parallel impedance should be non-zero");
}

#[test]
fn test_resonant_frequency() {
    let f = ac_analysis::resonant_frequency(0.1, 10e-6); // 100mH, 10μF
    assert!((f - 159.15).abs() < 0.1, "LC resonant frequency");
}

#[test]
fn test_quality_factor() {
    let q = ac_analysis::quality_factor(10.0, 0.1, 10e-6);
    assert!(q > 0.0, "Quality factor should be positive");
}

// ============================================================================
// DC ANALYSIS TESTS
// ============================================================================

#[test]
fn test_ohms_law_voltage() {
    let v = dc_analysis::voltage(10.0, 5.0); // I=10A, R=5Ω
    assert_eq!(v, 50.0, "V = IR");
}

#[test]
fn test_ohms_law_current() {
    let i = dc_analysis::current(120.0, 60.0); // V=120V, R=60Ω
    assert_eq!(i, 2.0, "I = V/R");
}

#[test]
fn test_ohms_law_resistance() {
    let r = dc_analysis::resistance(12.0, 3.0); // V=12V, I=3A
    assert_eq!(r, 4.0, "R = V/I");
}

#[test]
fn test_power() {
    let p = dc_analysis::power(120.0, 10.0); // V=120V, I=10A
    assert_eq!(p, 1200.0, "P = VI");
}

#[test]
fn test_power_from_resistance() {
    let p = dc_analysis::power_from_resistance(120.0, 60.0); // V=120V, R=60Ω
    assert_eq!(p, 240.0, "P = V²/R");
}

#[test]
fn test_conductance() {
    let g = dc_analysis::conductance(50.0); // R=50Ω
    assert_eq!(g, 0.02, "G = 1/R");
}

#[test]
fn test_series_resistance() {
    let r_total = dc_analysis::series_resistance(&[10.0, 20.0, 30.0]);
    assert_eq!(r_total, 60.0, "Series: R_total = R1 + R2 + R3");
}

#[test]
fn test_parallel_resistance() {
    let r_total = dc_analysis::parallel_resistance(&[60.0, 30.0, 20.0]);
    assert_eq!(r_total, 10.0, "Parallel resistors");
}

// ============================================================================
// IMPEDANCE TESTS
// ============================================================================

#[test]
fn test_impedance_from_resistance() {
    let z = impedance::from_resistance(100.0);
    assert_eq!(z, Complex64::new(100.0, 0.0), "Resistive impedance");
}

#[test]
fn test_impedance_from_capacitance() {
    let z = impedance::from_capacitance(10e-6, 1000.0); // 10μF, 1kHz
    assert!(z.im < 0.0, "Capacitive impedance has negative imaginary part");
}

#[test]
fn test_impedance_from_inductance() {
    let z = impedance::from_inductance(0.1, 1000.0); // 100mH, 1kHz
    assert!(z.im > 0.0, "Inductive impedance has positive imaginary part");
}

#[test]
fn test_impedance_magnitude() {
    let z = Complex64::new(3.0, 4.0);
    let mag = impedance::magnitude(z);
    assert_eq!(mag, 5.0, "|Z| for 3+4j should be 5");
}

#[test]
fn test_impedance_phase() {
    let z = Complex64::new(1.0, 1.0);
    let phase = impedance::phase(z);
    assert!((phase - PI / 4.0).abs() < 0.001, "Phase of 1+j should be 45°");
}

// ============================================================================
// FILTER DESIGN TESTS
// ============================================================================

#[test]
fn test_lowpass_rc_cutoff() {
    let fc = filter_design::rc_lowpass_cutoff(1000.0, 10e-6); // 1kΩ, 10μF
    assert!((fc - 15.92).abs() < 0.01, "RC lowpass cutoff frequency");
}

#[test]
fn test_highpass_rc_cutoff() {
    let fc = filter_design::rc_highpass_cutoff(1000.0, 100e-9); // 1kΩ, 100nF
    assert!((fc - 1591.55).abs() < 1.0, "RC highpass cutoff frequency");
}

#[test]
fn test_butterworth_poles() {
    let poles = filter_design::butterworth_poles(3);
    assert_eq!(poles.len(), 3, "3rd order Butterworth should have 3 poles");
}

#[test]
fn test_sallen_key_lowpass() {
    let (r1, r2, c1, c2) = filter_design::sallen_key_lowpass(1000.0, 0.707);
    assert!(r1 > 0.0 && r2 > 0.0 && c1 > 0.0 && c2 > 0.0, "Sallen-Key component values");
}

// ============================================================================
// NETWORK ANALYSIS TESTS
// ============================================================================

#[test]
fn test_voltage_divider() {
    let v_out = network_analysis::voltage_divider(12.0, 100.0, 100.0); // Equal resistors
    assert_eq!(v_out, 6.0, "Voltage divider with equal R should give V/2");
}

#[test]
fn test_current_divider() {
    let i_out = network_analysis::current_divider(10.0, 50.0, 50.0); // Equal resistors
    assert_eq!(i_out, 5.0, "Current divider with equal R should give I/2");
}

#[test]
fn test_thevenin_resistance() {
    // Simple example: two resistors in series seen from middle point
    let r_th = network_analysis::thevenin_resistance(vec![
        vec![10.0, -10.0],
        vec![-10.0, 20.0],
    ]);
    assert!(r_th > 0.0, "Thevenin resistance should be positive");
}

#[test]
fn test_norton_current() {
    let i_n = network_analysis::norton_current(12.0, 100.0);
    assert_eq!(i_n, 0.12, "Norton current I_N = V_oc / R_th");
}

#[test]
fn test_superposition_principle() {
    // Test with multiple sources
    let voltages = network_analysis::superposition(&[12.0, 5.0], &[100.0, 200.0]);
    assert!(voltages.len() > 0, "Superposition should return voltage array");
}

// ============================================================================
// POWER ANALYSIS TESTS
// ============================================================================

#[test]
fn test_apparent_power() {
    let s = power_analysis::apparent_power(120.0, 10.0); // V_rms, I_rms
    assert_eq!(s, 1200.0, "Apparent power S = V*I");
}

#[test]
fn test_real_power() {
    let p = power_analysis::real_power(120.0, 10.0, 0.8); // PF = 0.8
    assert_eq!(p, 960.0, "Real power P = V*I*cos(θ)");
}

#[test]
fn test_reactive_power() {
    let q = power_analysis::reactive_power(120.0, 10.0, 0.8);
    assert!((q - 720.0).abs() < 0.1, "Reactive power Q = V*I*sin(θ)");
}

#[test]
fn test_power_factor() {
    let pf = power_analysis::power_factor(800.0, 1000.0); // P=800W, S=1000VA
    assert_eq!(pf, 0.8, "Power factor PF = P/S");
}

#[test]
fn test_power_factor_correction() {
    let c = power_analysis::capacitor_for_pf_correction(10000.0, 0.7, 0.95, 60.0);
    assert!(c > 0.0, "Capacitor value for PF correction should be positive");
}

// ============================================================================
// TRANSIENT ANALYSIS TESTS
// ============================================================================

#[test]
fn test_rc_time_constant() {
    let tau = transient_analysis::rc_time_constant(1000.0, 100e-6); // 1kΩ, 100μF
    assert_eq!(tau, 0.1, "τ = RC");
}

#[test]
fn test_rl_time_constant() {
    let tau = transient_analysis::rl_time_constant(0.1, 10.0); // 100mH, 10Ω
    assert_eq!(tau, 0.01, "τ = L/R");
}

#[test]
fn test_rc_charging() {
    let v = transient_analysis::rc_charging_voltage(12.0, 0.1, 0.05); // V=12V, t=0.05s, τ=0.1s
    assert!(v > 0.0 && v < 12.0, "Charging voltage should be between 0 and V_final");
}

#[test]
fn test_rc_discharging() {
    let v = transient_analysis::rc_discharging_voltage(12.0, 0.1, 0.05);
    assert!(v > 0.0 && v < 12.0, "Discharging voltage should decrease");
}

#[test]
fn test_rl_current_rise() {
    let i = transient_analysis::rl_current_rise(10.0, 0.01, 0.005); // I_final=10A
    assert!(i > 0.0 && i < 10.0, "Current should rise toward final value");
}

#[test]
fn test_rl_current_decay() {
    let i = transient_analysis::rl_current_decay(10.0, 0.01, 0.005);
    assert!(i > 0.0 && i < 10.0, "Current should decay");
}

// ============================================================================
// NEC CALCULATIONS TESTS
// ============================================================================

#[test]
fn test_conductor_ampacity() {
    let ampacity = nec_calculations::conductor_ampacity("12 AWG", "THHN", 75.0);
    assert!(ampacity > 0.0, "Conductor ampacity should be positive");
}

#[test]
fn test_voltage_drop() {
    let drop = nec_calculations::voltage_drop(20.0, 100.0, 1.68, 2.0); // 20A, 100ft, copper
    assert!(drop > 0.0, "Voltage drop should be positive");
}

#[test]
fn test_conductor_resistance() {
    let r = nec_calculations::conductor_resistance("12 AWG", 100.0); // 100 feet
    assert!(r > 0.0, "Resistance should be positive");
}

#[test]
fn test_conduit_fill() {
    let fill = nec_calculations::conduit_fill(0.75, &[0.01, 0.01, 0.01]); // 3/4" conduit
    assert!(fill >= 0.0 && fill <= 100.0, "Conduit fill should be 0-100%");
}

#[test]
fn test_box_fill_calculation() {
    let volume = nec_calculations::box_fill_calculation(5, 2, 1, 3); // conductors, devices, clamps, grounds
    assert!(volume > 0.0, "Box fill volume should be positive");
}

#[test]
fn test_overcurrent_protection() {
    let ocpd = nec_calculations::overcurrent_protection_size(23.5); // 23.5A load
    assert!(ocpd >= 25.0, "OCPD should be next standard size up");
}

#[test]
fn test_grounding_electrode_conductor() {
    let size = nec_calculations::grounding_electrode_conductor_size("350 kcmil");
    assert!(!size.is_empty(), "GEC size should be returned");
}

#[test]
fn test_equipment_grounding_conductor() {
    let size = nec_calculations::equipment_grounding_conductor_size(100.0); // 100A OCPD
    assert!(!size.is_empty(), "EGC size should be returned");
}
