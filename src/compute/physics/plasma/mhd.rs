//! Magnetohydrodynamics (MHD) - Single fluid description of plasma

use super::PlasmaParams;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// MHD state (velocity, magnetic field, pressure, density)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MHDState {
    pub velocity: [f64; 3],      // Fluid velocity (m/s)
    pub b_field: [f64; 3],        // Magnetic field (T)
    pub pressure: f64,            // Pressure (Pa)
    pub density: f64,             // Mass density (kg/m³)
}

/// Ideal MHD assumptions
pub fn check_ideal_mhd(params: &PlasmaParams) -> bool {
    // Check if plasma is collisional enough
    let debye_len = super::debye_length(params);
    let system_size = 1.0; // Assume 1m system

    // Check Lundquist number (ratio of resistive to Alfvén timescales)
    let lundquist = lundquist_number(params, system_size);

    debye_len < system_size * 0.001 && lundquist > 100.0
}

/// Lundquist number: S = τ_R / τ_A
fn lundquist_number(params: &PlasmaParams, length: f64) -> f64 {
    #[allow(dead_code)]
    const MU_0: f64 = 4.0e-7 * PI;

    let rho = params.n_i * params.ion_mass;
    let v_a = super::alfven_velocity(params.b_field, rho);
    let tau_a = length / v_a;

    // Resistive time (simplified)
    let eta = 1e-6; // Typical resistivity
    let tau_r = length * length / eta;

    tau_r / tau_a
}

/// Magnetic pressure: p_B = B²/(2μ₀)
pub fn magnetic_pressure(b_field: f64) -> f64 {
    const MU_0: f64 = 4.0e-7 * PI;
    b_field * b_field / (2.0 * MU_0)
}

/// Magnetic tension force: (B·∇)B/μ₀
pub fn magnetic_tension(b: &[f64; 3], grad_b: &[[f64; 3]; 3]) -> [f64; 3] {
    const MU_0: f64 = 4.0e-7 * PI;

    let mut tension = [0.0; 3];
    for i in 0..3 {
        for j in 0..3 {
            tension[i] += b[j] * grad_b[j][i];
        }
        tension[i] /= MU_0;
    }
    tension
}

/// Lorentz force: J × B where J = (∇ × B)/μ₀
pub fn lorentz_force(b: &[f64; 3], curl_b: &[f64; 3]) -> [f64; 3] {
    const MU_0: f64 = 4.0e-7 * PI;

    // J = curl(B)/μ₀
    let j = [curl_b[0] / MU_0, curl_b[1] / MU_0, curl_b[2] / MU_0];

    // F = J × B
    [
        j[1] * b[2] - j[2] * b[1],
        j[2] * b[0] - j[0] * b[2],
        j[0] * b[1] - j[1] * b[0],
    ]
}

/// Frozen-in condition: B/ρ is conserved along fluid elements
pub fn frozen_in_parameter(state: &MHDState) -> f64 {
    let b_mag = (state.b_field[0].powi(2) + state.b_field[1].powi(2) + state.b_field[2].powi(2)).sqrt();
    b_mag / state.density
}

/// MHD equilibrium condition: ∇p = J × B
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MHDEquilibrium {
    pub grad_p: [f64; 3],        // Pressure gradient
    pub lorentz_force: [f64; 3], // J × B force
    pub is_equilibrium: bool,
}

pub fn check_mhd_equilibrium(grad_p: &[f64; 3], b: &[f64; 3], curl_b: &[f64; 3]) -> MHDEquilibrium {
    let lorentz = lorentz_force(b, curl_b);

    // Check if ∇p ≈ J × B
    let diff = [
        (grad_p[0] - lorentz[0]).abs(),
        (grad_p[1] - lorentz[1]).abs(),
        (grad_p[2] - lorentz[2]).abs(),
    ];

    let max_diff = diff.iter().cloned().fold(0.0, f64::max);
    let is_eq = max_diff < 0.01 * grad_p.iter().map(|x| x.abs()).fold(0.0, f64::max);

    MHDEquilibrium {
        grad_p: *grad_p,
        lorentz_force: lorentz,
        is_equilibrium: is_eq,
    }
}

/// Plasma current density from Ampère's law: J = ∇ × B / μ₀
pub fn current_density(curl_b: &[f64; 3]) -> [f64; 3] {
    const MU_0: f64 = 4.0e-7 * PI;
    [curl_b[0] / MU_0, curl_b[1] / MU_0, curl_b[2] / MU_0]
}

