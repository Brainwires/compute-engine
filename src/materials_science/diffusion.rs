//! Diffusion in Materials
//!
//! Fick's laws, diffusion coefficients, concentration profiles,
//! and mass transport in solids.

use super::constants::*;
use std::f64::consts::PI;

/// Calculate diffusion flux (Fick's first law)
///
/// J = -D(dC/dx)
///
/// # Arguments
/// * `diffusion_coeff` - Diffusion coefficient (m²/s)
/// * `concentration_gradient` - dC/dx (concentration/m)
///
/// # Returns
/// Flux (particles/(m²·s) or mol/(m²·s))
pub fn fick_first_law(diffusion_coeff: f64, concentration_gradient: f64) -> f64 {
    -diffusion_coeff * concentration_gradient
}

/// Calculate diffusion coefficient from Arrhenius equation
///
/// D = D₀ exp(-Q/(RT))
///
/// # Arguments
/// * `pre_exponential` - D₀ (m²/s)
/// * `activation_energy` - Q (J/mol)
/// * `temperature_k` - Temperature (K)
///
/// # Returns
/// Diffusion coefficient (m²/s)
pub fn arrhenius_diffusion_coefficient(
    pre_exponential: f64,
    activation_energy: f64,
    temperature_k: f64,
) -> f64 {
    pre_exponential * (-activation_energy / (R * temperature_k)).exp()
}

/// Calculate diffusion distance (root mean square displacement)
///
/// <x²>^(1/2) = √(2Dt)
///
/// # Arguments
/// * `diffusion_coeff` - Diffusion coefficient (m²/s)
/// * `time` - Diffusion time (s)
///
/// # Returns
/// RMS diffusion distance (m)
pub fn diffusion_distance(diffusion_coeff: f64, time: f64) -> f64 {
    (2.0 * diffusion_coeff * time).sqrt()
}

/// Calculate time required to diffuse a given distance
///
/// t = x²/(2D)
pub fn diffusion_time(diffusion_coeff: f64, distance: f64) -> f64 {
    (distance * distance) / (2.0 * diffusion_coeff)
}

/// Calculate concentration profile for semi-infinite diffusion
///
/// C(x,t) = C₀ erfc(x/(2√(Dt)))
///
/// where erfc is the complementary error function
///
/// # Arguments
/// * `initial_concentration` - C₀ (mol/m³ or similar)
/// * `position` - x (m)
/// * `diffusion_coeff` - D (m²/s)
/// * `time` - t (s)
pub fn concentration_semi_infinite(
    initial_concentration: f64,
    position: f64,
    diffusion_coeff: f64,
    time: f64,
) -> f64 {
    let arg = position / (2.0 * (diffusion_coeff * time).sqrt());
    initial_concentration * erfc_approx(arg)
}

/// Approximate complementary error function
///
/// erfc(x) = 1 - erf(x)
///
/// Using Abramowitz and Stegun approximation
fn erfc_approx(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc_approx(-x);
    }

    // Constants for approximation
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;
    let p = 0.3275911;

    let t = 1.0 / (1.0 + p * x);
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;

    let erf = 1.0 - (a1 * t + a2 * t2 + a3 * t3 + a4 * t4 + a5 * t5) * (-x * x).exp();
    1.0 - erf
}

/// Calculate concentration for diffusion couple (two semi-infinite solids)
///
/// C(x,t) = (C_A + C_B)/2 + (C_A - C_B)/2 · erf(x/(2√(Dt)))
///
/// # Arguments
/// * `conc_a` - Concentration in material A
/// * `conc_b` - Concentration in material B
/// * `position` - x (distance from interface, can be negative)
/// * `diffusion_coeff` - D (m²/s)
/// * `time` - t (s)
pub fn diffusion_couple_concentration(
    conc_a: f64,
    conc_b: f64,
    position: f64,
    diffusion_coeff: f64,
    time: f64,
) -> f64 {
    let arg = position / (2.0 * (diffusion_coeff * time).sqrt());
    let erf = 1.0 - erfc_approx(arg);
    (conc_a + conc_b) / 2.0 + (conc_a - conc_b) / 2.0 * erf
}

/// Calculate diffusion through a membrane (steady state)
///
/// J = D(C₁ - C₂)/L
///
/// # Arguments
/// * `diffusion_coeff` - D (m²/s)
/// * `conc_1` - Concentration on side 1
/// * `conc_2` - Concentration on side 2
/// * `thickness` - Membrane thickness L (m)
pub fn membrane_flux_steady_state(
    diffusion_coeff: f64,
    conc_1: f64,
    conc_2: f64,
    thickness: f64,
) -> f64 {
    diffusion_coeff * (conc_1 - conc_2) / thickness
}

/// Calculate permeability coefficient
///
/// P = DS (where S is solubility)
///
/// # Arguments
/// * `diffusion_coeff` - D (m²/s)
/// * `solubility` - S (dimensionless or concentration ratio)
pub fn permeability_coefficient(diffusion_coeff: f64, solubility: f64) -> f64 {
    diffusion_coeff * solubility
}

/// Calculate interdiffusion coefficient (Darken equation)
///
/// D̃ = x_B D_A + x_A D_B
///
/// # Arguments
/// * `mole_fraction_a` - x_A
/// * `mole_fraction_b` - x_B
/// * `diffusion_a` - Self-diffusion coefficient of A (m²/s)
/// * `diffusion_b` - Self-diffusion coefficient of B (m²/s)
pub fn interdiffusion_coefficient(
    mole_fraction_a: f64,
    mole_fraction_b: f64,
    diffusion_a: f64,
    diffusion_b: f64,
) -> f64 {
    mole_fraction_b * diffusion_a + mole_fraction_a * diffusion_b
}

/// Calculate Kirkendall velocity (marker velocity in diffusion couple)
///
/// v = (D_A - D_B)(dC_A/dx) / C_total
///
/// # Arguments
/// * `diffusion_a` - Diffusion coefficient of species A (m²/s)
/// * `diffusion_b` - Diffusion coefficient of species B (m²/s)
/// * `concentration_gradient_a` - dC_A/dx (concentration/m)
/// * `total_concentration` - Total concentration (atoms/m³)
pub fn kirkendall_velocity(
    diffusion_a: f64,
    diffusion_b: f64,
    concentration_gradient_a: f64,
    total_concentration: f64,
) -> f64 {
    (diffusion_a - diffusion_b) * concentration_gradient_a / total_concentration
}

/// Calculate vacancy diffusion coefficient
///
/// D_v = D₀ exp(-Q_v/(RT))
///
/// where Q_v = E_f + E_m (formation energy + migration energy)
///
/// # Arguments
/// * `pre_exponential` - D₀ (m²/s), typically a₀²·ν₀ where a₀ is lattice parameter
/// * `formation_energy` - Vacancy formation energy (J/mol)
/// * `migration_energy` - Vacancy migration energy (J/mol)
/// * `temperature_k` - Temperature (K)
pub fn vacancy_diffusion_coefficient(
    pre_exponential: f64,
    formation_energy: f64,
    migration_energy: f64,
    temperature_k: f64,
) -> f64 {
    let total_activation = formation_energy + migration_energy;
    pre_exponential * (-total_activation / (R * temperature_k)).exp()
}

/// Calculate grain boundary diffusion coefficient
///
/// Typically faster than bulk diffusion by several orders of magnitude
///
/// D_gb = D₀_gb exp(-Q_gb/(RT))
///
/// Generally: Q_gb ≈ 0.4 to 0.6 Q_bulk
pub fn grain_boundary_diffusion_coefficient(
    pre_exponential_gb: f64,
    activation_energy_gb: f64,
    temperature_k: f64,
) -> f64 {
    pre_exponential_gb * (-activation_energy_gb / (R * temperature_k)).exp()
}

/// Calculate effective diffusion coefficient (Harrison's model)
///
/// D_eff = D_bulk + (δ/d)D_gb
///
/// # Arguments
/// * `bulk_diffusion` - Bulk diffusion coefficient (m²/s)
/// * `grain_boundary_diffusion` - GB diffusion coefficient (m²/s)
/// * `grain_boundary_width` - δ (typically ~1 nm)
/// * `grain_size` - d (m)
pub fn effective_diffusion_coefficient(
    bulk_diffusion: f64,
    grain_boundary_diffusion: f64,
    grain_boundary_width: f64,
    grain_size: f64,
) -> f64 {
    bulk_diffusion + (grain_boundary_width / grain_size) * grain_boundary_diffusion
}

/// Calculate interstitial diffusion coefficient
///
/// Generally faster than vacancy diffusion
///
/// D_i = a² z ν exp(-E_m/(k_B T))
///
/// Simplified version using Arrhenius form
pub fn interstitial_diffusion_coefficient(
    lattice_param: f64,
    jump_frequency: f64,
    migration_energy_ev: f64,
    temperature_k: f64,
) -> f64 {
    let migration_energy_j = migration_energy_ev * E;
    lattice_param * lattice_param * jump_frequency * (-migration_energy_j / (K_B * temperature_k)).exp()
}

/// Calculate chemical diffusion coefficient (relates to chemical potential gradient)
///
/// D_chem = D (∂ln(γC)/∂ln(C))
///
/// where γ is activity coefficient. For ideal solution γ=1, D_chem = D
///
/// # Arguments
/// * `self_diffusion` - Self-diffusion coefficient (m²/s)
/// * `thermodynamic_factor` - ∂ln(a)/∂ln(C), where a is activity
pub fn chemical_diffusion_coefficient(self_diffusion: f64, thermodynamic_factor: f64) -> f64 {
    self_diffusion * thermodynamic_factor
}

/// Calculate Einstein relation for diffusion
///
/// D = (k_B T)/(6πηr)
///
/// (For spherical particles in a fluid)
///
/// # Arguments
/// * `temperature_k` - Temperature (K)
/// * `viscosity` - Dynamic viscosity (Pa·s)
/// * `particle_radius` - Particle radius (m)
pub fn einstein_diffusion_coefficient(
    temperature_k: f64,
    viscosity: f64,
    particle_radius: f64,
) -> f64 {
    (K_B * temperature_k) / (6.0 * PI * viscosity * particle_radius)
}

/// Calculate diffusion penetration depth (for case hardening, nitriding, etc.)
///
/// Depth at which concentration drops to a fraction f of surface concentration
///
/// x = 2√(Dt) · erfc⁻¹(f)
///
/// For f = 0.5 (half the surface concentration), erfc⁻¹(0.5) ≈ 0.477
pub fn penetration_depth(diffusion_coeff: f64, time: f64, concentration_fraction: f64) -> f64 {
    // Approximate inverse erfc for common values
    let inv_erfc = if concentration_fraction > 0.9 {
        0.089
    } else if concentration_fraction > 0.5 {
        0.477
    } else if concentration_fraction > 0.1 {
        1.163
    } else {
        2.0 // Very low concentration
    };

    2.0 * (diffusion_coeff * time).sqrt() * inv_erfc
}

/// Calculate ionic conductivity from diffusion (Nernst-Einstein relation)
///
/// σ = (n q² D)/(k_B T)
///
/// # Arguments
/// * `carrier_density` - n (carriers/m³)
/// * `charge` - q (C), typically e for electrons or ions
/// * `diffusion_coeff` - D (m²/s)
/// * `temperature_k` - T (K)
///
/// # Returns
/// Ionic conductivity (S/m)
pub fn ionic_conductivity_from_diffusion(
    carrier_density: f64,
    charge: f64,
    diffusion_coeff: f64,
    temperature_k: f64,
) -> f64 {
    (carrier_density * charge * charge * diffusion_coeff) / (K_B * temperature_k)
}

