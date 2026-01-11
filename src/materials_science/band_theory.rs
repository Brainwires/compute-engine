//! Band Theory and Electronic Structure
//!
//! Calculations for band gaps, density of states, Fermi surfaces,
//! effective mass, and electronic properties of materials.

use super::constants::*;
use super::MaterialClass;
use std::f64::consts::PI;

/// Calculate Fermi energy for free electron gas
///
/// E_F = (ℏ²/2m)(3π²n)^(2/3)
///
/// # Arguments
/// * `electron_density` - Electron density (electrons/m³)
///
/// # Returns
/// Fermi energy in eV
pub fn fermi_energy_free_electron(electron_density: f64) -> f64 {
    let energy_joules = (HBAR * HBAR / (2.0 * M_E)) * (3.0 * PI * PI * electron_density).powf(2.0 / 3.0);
    energy_joules / E // Convert to eV
}

/// Calculate Fermi temperature
///
/// T_F = E_F / k_B
pub fn fermi_temperature(fermi_energy_ev: f64) -> f64 {
    (fermi_energy_ev * E) / K_B
}

/// Calculate Fermi velocity
///
/// v_F = √(2E_F/m)
pub fn fermi_velocity(fermi_energy_ev: f64) -> f64 {
    let energy_joules = fermi_energy_ev * E;
    (2.0 * energy_joules / M_E).sqrt()
}

/// Calculate density of states at Fermi level for 3D free electron gas
///
/// g(E_F) = (3n)/(2E_F)
///
/// # Returns
/// Density of states in states/(eV·atom)
pub fn density_of_states_3d(electron_density: f64, fermi_energy_ev: f64) -> f64 {
    (3.0 * electron_density) / (2.0 * fermi_energy_ev)
}

/// Calculate Thomas-Fermi screening length
///
/// λ_TF = √(ε₀/(e²·g(E_F)))
pub fn thomas_fermi_screening_length(fermi_energy_ev: f64, electron_density: f64) -> f64 {
    let dos = density_of_states_3d(electron_density, fermi_energy_ev);
    (EPSILON_0 / (E * E * dos)).sqrt()
}

/// Calculate band gap from conductivity temperature dependence
///
/// σ(T) = σ₀ exp(-E_g/(2k_B T))
///
/// This function calculates E_g from two (σ, T) measurements
pub fn band_gap_from_conductivity(
    sigma1: f64,
    temp1_k: f64,
    sigma2: f64,
    temp2_k: f64,
) -> f64 {
    let ln_ratio = (sigma2 / sigma1).ln();
    let t_ratio = 1.0 / temp1_k - 1.0 / temp2_k;
    let eg_joules = 2.0 * K_B * ln_ratio / t_ratio;
    eg_joules / E // Convert to eV
}

/// Classify material based on band gap
///
/// - Insulator: E_g > 3 eV
/// - Semiconductor: 0 < E_g < 3 eV
/// - Conductor: E_g = 0 (overlapping bands)
pub fn classify_material(band_gap_ev: f64) -> MaterialClass {
    if band_gap_ev < 0.001 {
        MaterialClass::Conductor
    } else if band_gap_ev < 3.0 {
        MaterialClass::Semiconductor
    } else {
        MaterialClass::Insulator
    }
}

/// Calculate intrinsic carrier concentration for semiconductor
///
/// n_i = √(N_c N_v) exp(-E_g/(2k_B T))
///
/// Simplified version assuming N_c ≈ N_v ≈ 10²⁵ m⁻³
pub fn intrinsic_carrier_concentration(band_gap_ev: f64, temperature_k: f64) -> f64 {
    let nc_nv: f64 = 1e25 * 1e25; // Typical values
    nc_nv.sqrt() * (-(band_gap_ev * E) / (2.0 * K_B * temperature_k)).exp()
}

/// Calculate effective mass from band curvature
///
/// 1/m* = (1/ℏ²) d²E/dk²
///
/// Approximate using parabolic band: E(k) = E_0 + (ℏ²k²)/(2m*)
/// Returns m*/m_e ratio
pub fn effective_mass_ratio(band_curvature: f64) -> f64 {
    // band_curvature in eV·Ų
    let curvature_joules_m2 = band_curvature * E * 1e-20;
    (HBAR * HBAR) / (curvature_joules_m2 * M_E)
}

/// Calculate cyclotron frequency
///
/// ω_c = eB/m*
///
/// # Arguments
/// * `magnetic_field` - Magnetic field in Tesla
/// * `effective_mass_ratio` - m*/m_e ratio
///
/// # Returns
/// Cyclotron frequency in rad/s
pub fn cyclotron_frequency(magnetic_field_t: f64, effective_mass_ratio: f64) -> f64 {
    (E * magnetic_field_t) / (effective_mass_ratio * M_E)
}

/// Calculate plasma frequency for metal
///
/// ω_p = √(ne²/(ε₀m))
///
/// # Returns
/// Plasma frequency in rad/s
pub fn plasma_frequency(electron_density: f64) -> f64 {
    ((electron_density * E * E) / (EPSILON_0 * M_E)).sqrt()
}

/// Calculate Debye frequency (phonons)
///
/// ω_D = v_s (6π²n)^(1/3)
///
/// # Arguments
/// * `sound_velocity` - Speed of sound in material (m/s)
/// * `atom_density` - Atomic density (atoms/m³)
pub fn debye_frequency(sound_velocity: f64, atom_density: f64) -> f64 {
    sound_velocity * (6.0 * PI * PI * atom_density).powf(1.0 / 3.0)
}

/// Calculate electron mobility from conductivity
///
/// μ = σ / (ne)
///
/// # Arguments
/// * `conductivity` - Electrical conductivity (S/m)
/// * `carrier_density` - Carrier density (carriers/m³)
///
/// # Returns
/// Mobility in m²/(V·s)
pub fn electron_mobility(conductivity: f64, carrier_density: f64) -> f64 {
    conductivity / (carrier_density * E)
}

/// Calculate Hall coefficient
///
/// R_H = 1/(ne) for single carrier type
///
/// # Returns
/// Hall coefficient in m³/C
pub fn hall_coefficient(carrier_density: f64) -> f64 {
    1.0 / (carrier_density * E)
}

/// Calculate drift velocity
///
/// v_d = μE
///
/// # Arguments
/// * `mobility` - Carrier mobility (m²/(V·s))
/// * `electric_field` - Electric field (V/m)
pub fn drift_velocity(mobility: f64, electric_field: f64) -> f64 {
    mobility * electric_field
}

/// Calculate mean free path
///
/// λ = v_F τ = v_F (μ m*) / e
pub fn mean_free_path(fermi_velocity: f64, mobility: f64, effective_mass_ratio: f64) -> f64 {
    let tau = (mobility * effective_mass_ratio * M_E) / E;
    fermi_velocity * tau
}

/// Calculate Richardson constant for thermionic emission
///
/// A = (4πmk²e)/(h³) ≈ 120 A/(cm²·K²)
pub const RICHARDSON_CONSTANT: f64 = 1.20173e6; // A/(m²·K²)

/// Calculate thermionic emission current density
///
/// J = A T² exp(-φ/(k_B T))
///
/// # Arguments
/// * `temperature_k` - Temperature in Kelvin
/// * `work_function_ev` - Work function in eV
///
/// # Returns
/// Current density in A/m²
pub fn thermionic_emission_current(temperature_k: f64, work_function_ev: f64) -> f64 {
    RICHARDSON_CONSTANT * temperature_k.powi(2)
        * (-(work_function_ev * E) / (K_B * temperature_k)).exp()
}

/// Calculate Seebeck coefficient (thermoelectric power)
///
/// S ≈ (π²k_B²T)/(3eE_F) for metals (simplified)
///
/// # Returns
/// Seebeck coefficient in V/K
pub fn seebeck_coefficient_metal(temperature_k: f64, fermi_energy_ev: f64) -> f64 {
    let fermi_energy_j = fermi_energy_ev * E;
    (PI * PI * K_B * K_B * temperature_k) / (3.0 * E * fermi_energy_j)
}

/// Calculate thermoelectric figure of merit
///
/// ZT = S²σT/κ
///
/// # Arguments
/// * `seebeck` - Seebeck coefficient (V/K)
/// * `conductivity` - Electrical conductivity (S/m)
/// * `thermal_conductivity` - Thermal conductivity (W/(m·K))
/// * `temperature` - Temperature (K)
pub fn figure_of_merit_zt(
    seebeck: f64,
    conductivity: f64,
    thermal_conductivity: f64,
    temperature: f64,
) -> f64 {
    (seebeck * seebeck * conductivity * temperature) / thermal_conductivity
}

/// Calculate Wiedemann-Franz ratio
///
/// L = κ/(σT) ≈ (π²/3)(k_B/e)² = 2.44×10⁻⁸ W·Ω/K²
///
/// (Lorenz number for metals)
pub const LORENZ_NUMBER: f64 = 2.44e-8; // W·Ω/K²

pub fn thermal_conductivity_from_lorenz(
    electrical_conductivity: f64,
    temperature_k: f64,
) -> f64 {
    LORENZ_NUMBER * electrical_conductivity * temperature_k
}

