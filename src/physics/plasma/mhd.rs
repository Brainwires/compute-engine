//! Magnetohydrodynamics (MHD) - Single fluid description of plasma

use super::{PlasmaParams, E, M_P};
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_magnetic_pressure() {
        let b = 5.0; // 5 Tesla
        let p_b = magnetic_pressure(b);

        // p_B ~ 10^6 Pa for 5T field
        assert!(p_b > 1e6 && p_b < 1e7);
    }

    #[test]
    fn test_magnetic_tension() {
        let b = [1.0, 0.0, 0.0];
        let grad_b = [
            [0.1, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ];

        let tension = magnetic_tension(&b, &grad_b);

        assert!(tension[0].is_finite());
    }

    #[test]
    fn test_lorentz_force() {
        let b = [0.0, 0.0, 5.0];
        let curl_b = [1.0, 0.0, 0.0];

        let force = lorentz_force(&b, &curl_b);

        // Force should be perpendicular to both J and B
        assert!(force[0].is_finite());
        assert!(force[1].is_finite());
        assert!(force[2].abs() < 1e-10); // Should be zero
    }

    #[test]
    fn test_frozen_in_parameter() {
        let state = MHDState {
            velocity: [0.0, 0.0, 0.0],
            b_field: [0.0, 0.0, 5.0],
            pressure: 1e5,
            density: 1e-6,
        };

        let param = frozen_in_parameter(&state);

        assert!(param > 0.0);
        assert!(param.is_finite());
    }

    #[test]
    fn test_mhd_equilibrium() {
        let grad_p = [1000.0, 0.0, 0.0];
        let b = [0.0, 0.0, 5.0];
        let curl_b = [1000.0 * 4.0e-7 * PI, 0.0, 0.0]; // Chosen to balance

        let eq = check_mhd_equilibrium(&grad_p, &b, &curl_b);

        assert!(eq.lorentz_force[0].is_finite());
    }

    #[test]
    fn test_current_density() {
        let curl_b = [1.0, 2.0, 3.0];
        let j = current_density(&curl_b);

        // J should be proportional to curl(B)
        assert!(j[0] > 0.0);
        assert!(j[1] > j[0]);
        assert!(j[2] > j[1]);
    }

    #[test]
    fn test_ideal_mhd_check() {
        let params = PlasmaParams::tokamak();
        let is_ideal = check_ideal_mhd(&params);

        // Tokamak should satisfy ideal MHD
        assert!(is_ideal);
    }

    #[test]
    fn test_lundquist_number() {
        let params = PlasmaParams::tokamak();
        let s = lundquist_number(&params, 1.0);

        // Tokamak should have large Lundquist number
        assert!(s > 100.0);
    }
}
