//! Plasma Physics and Magnetohydrodynamics (MHD)
//!
//! **Plasma:** The fourth state of matter - ionized gas with collective behavior
//!
//! **Key Physics:**
//! - Debye shielding and plasma frequency
//! - Magnetohydrodynamics (MHD) equations
//! - Plasma waves (Alfvén, Langmuir)
//! - Magnetic confinement (tokamaks)
//! - Plasma instabilities
//!
//! **Applications:**
//! - Fusion reactors (ITER, tokamaks)
//! - Space plasmas (solar wind, magnetosphere)
//! - Plasma propulsion

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

pub mod mhd;
pub mod waves;
pub mod confinement;

pub use mhd::*;
pub use waves::*;
pub use confinement::*;

/// Physical constants
pub const E: f64 = 1.602176634e-19;      // Elementary charge (C)
pub const M_E: f64 = 9.1093837015e-31;   // Electron mass (kg)
pub const M_P: f64 = 1.67262192369e-27;  // Proton mass (kg)
pub const EPSILON_0: f64 = 8.8541878128e-12; // Permittivity (F/m)
pub const K_B: f64 = 1.380649e-23;       // Boltzmann constant
pub const C: f64 = 299792458.0;          // Speed of light (m/s)

/// Plasma parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PlasmaParams {
    pub n_e: f64,           // Electron density (m⁻³)
    pub n_i: f64,           // Ion density (m⁻³)
    pub t_e: f64,           // Electron temperature (eV)
    pub t_i: f64,           // Ion temperature (eV)
    pub b_field: f64,       // Magnetic field strength (T)
    pub z_ion: f64,         // Ion charge number
    pub ion_mass: f64,      // Ion mass (kg)
}

impl PlasmaParams {
    /// Typical tokamak plasma parameters
    pub fn tokamak() -> Self {
        Self {
            n_e: 1e20,      // 10^20 m⁻³
            n_i: 1e20,
            t_e: 10000.0,   // 10 keV
            t_i: 10000.0,
            b_field: 5.0,   // 5 Tesla
            z_ion: 1.0,     // Deuterium
            ion_mass: 2.0 * M_P,
        }
    }

    /// Solar corona parameters
    pub fn solar_corona() -> Self {
        Self {
            n_e: 1e15,      // 10^15 m⁻³
            n_i: 1e15,
            t_e: 100.0,     // 100 eV ~ 1 MK
            t_i: 100.0,
            b_field: 0.01,  // 100 Gauss
            z_ion: 1.0,     // Hydrogen
            ion_mass: M_P,
        }
    }

    /// Convert temperature from eV to Kelvin
    pub fn electron_temp_kelvin(&self) -> f64 {
        self.t_e * E / K_B
    }

    pub fn ion_temp_kelvin(&self) -> f64 {
        self.t_i * E / K_B
    }
}

/// Debye length: λ_D = √(ε₀kT_e / n_e e²)
pub fn debye_length(params: &PlasmaParams) -> f64 {
    let t_e_joules = params.t_e * E;
    ((EPSILON_0 * t_e_joules) / (params.n_e * E * E)).sqrt()
}

/// Plasma frequency: ω_pe = √(n_e e² / ε₀m_e)
pub fn plasma_frequency(n_e: f64) -> f64 {
    ((n_e * E * E) / (EPSILON_0 * M_E)).sqrt()
}

/// Electron cyclotron frequency: ω_ce = eB/m_e
pub fn electron_cyclotron_frequency(b_field: f64) -> f64 {
    (E * b_field) / M_E
}

/// Ion cyclotron frequency: ω_ci = ZeB/m_i
pub fn ion_cyclotron_frequency(b_field: f64, z_ion: f64, ion_mass: f64) -> f64 {
    (z_ion * E * b_field) / ion_mass
}

/// Larmor radius (gyroradius): r_L = v⊥/ω_c
pub fn larmor_radius(v_perp: f64, omega_c: f64) -> f64 {
    v_perp / omega_c
}

/// Thermal velocity: v_th = √(kT/m)
pub fn thermal_velocity(temp_ev: f64, mass: f64) -> f64 {
    let t_joules = temp_ev * E;
    (t_joules / mass).sqrt()
}

/// Plasma beta: β = p / (B²/2μ₀)
pub fn plasma_beta(params: &PlasmaParams) -> f64 {
    const MU_0: f64 = 4.0e-7 * PI; // Permeability

    let p_e = params.n_e * params.t_e * E;
    let p_i = params.n_i * params.t_i * E;
    let p_total = p_e + p_i;

    let b_pressure = params.b_field * params.b_field / (2.0 * MU_0);

    p_total / b_pressure
}

/// Alfvén velocity: v_A = B/√(μ₀ρ)
pub fn alfven_velocity(b_field: f64, density: f64) -> f64 {
    const MU_0: f64 = 4.0e-7 * PI;
    b_field / (MU_0 * density).sqrt()
}

/// Number of particles in Debye sphere
pub fn debye_number(params: &PlasmaParams) -> f64 {
    let lambda_d = debye_length(params);
    (4.0 / 3.0) * PI * params.n_e * lambda_d.powi(3)
}

/// Check if system is a plasma (N_D >> 1)
pub fn is_plasma(params: &PlasmaParams) -> bool {
    debye_number(params) > 10.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_debye_length_tokamak() {
        let params = PlasmaParams::tokamak();
        let lambda_d = debye_length(&params);

        // Tokamak: λ_D ~ 10^-5 m
        assert!(lambda_d > 1e-6 && lambda_d < 1e-4);
    }

    #[test]
    fn test_plasma_frequency() {
        let n_e = 1e20; // Tokamak density
        let omega_pe = plasma_frequency(n_e);

        // ω_pe ~ 10^12 rad/s
        assert!(omega_pe > 1e11 && omega_pe < 1e13);
    }

    #[test]
    fn test_cyclotron_frequencies() {
        let b = 5.0; // 5 Tesla

        let omega_ce = electron_cyclotron_frequency(b);
        let omega_ci = ion_cyclotron_frequency(b, 1.0, M_P);

        // Electron cyclotron frequency much higher than ion
        assert!(omega_ce > omega_ci * 1000.0);
    }

    #[test]
    fn test_thermal_velocity() {
        let t_e = 1000.0; // 1 keV
        let v_th_e = thermal_velocity(t_e, M_E);

        // Electron thermal velocity ~ 10^7 m/s
        assert!(v_th_e > 1e6 && v_th_e < 1e8);
    }

    #[test]
    fn test_larmor_radius() {
        let v_perp = 1e6; // 1 Mm/s
        let omega_c = 1e8; // rad/s

        let r_l = larmor_radius(v_perp, omega_c);

        assert!(r_l > 0.0);
        assert!(r_l.is_finite());
    }

    #[test]
    fn test_plasma_beta() {
        let params = PlasmaParams::tokamak();
        let beta = plasma_beta(&params);

        // Tokamak typically has β < 0.1
        assert!(beta > 0.0);
        assert!(beta < 1.0);
    }

    #[test]
    fn test_alfven_velocity() {
        let b = 5.0; // Tesla
        let rho = 1e-6; // kg/m³ (low density plasma)

        let v_a = alfven_velocity(b, rho);

        // Alfvén velocity ~ 10^6 m/s for typical plasma
        assert!(v_a > 0.0);
        assert!(v_a.is_finite());
    }

    #[test]
    fn test_debye_number() {
        let params = PlasmaParams::tokamak();
        let n_d = debye_number(&params);

        // Should be >> 1 for good plasma
        assert!(n_d > 1000.0);
    }

    #[test]
    fn test_is_plasma() {
        let tokamak = PlasmaParams::tokamak();
        assert!(is_plasma(&tokamak));

        let corona = PlasmaParams::solar_corona();
        assert!(is_plasma(&corona));
    }

    #[test]
    fn test_temperature_conversion() {
        let params = PlasmaParams::tokamak();
        let t_k = params.electron_temp_kelvin();

        // 10 keV ~ 100 million K
        assert!(t_k > 1e8 && t_k < 1e9);
    }
}
