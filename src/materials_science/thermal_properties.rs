//! Thermal Properties of Materials
//!
//! Heat capacity, thermal expansion, thermal conductivity,
//! and temperature-dependent properties.

use super::constants::*;
use std::f64::consts::PI;

/// Calculate Debye temperature from Debye frequency
///
/// Θ_D = ℏω_D / k_B
pub fn debye_temperature(debye_frequency: f64) -> f64 {
    (HBAR * debye_frequency) / K_B
}

/// Calculate Debye heat capacity (low temperature limit)
///
/// C_V = (12π⁴/5)Nk_B(T/Θ_D)³
///
/// # Arguments
/// * `temperature_k` - Temperature (K)
/// * `debye_temp` - Debye temperature (K)
/// * `num_atoms` - Number of atoms
///
/// # Returns
/// Heat capacity in J/K
pub fn debye_heat_capacity(temperature_k: f64, debye_temp: f64, num_atoms: f64) -> f64 {
    let ratio = temperature_k / debye_temp;
    (12.0 * PI.powi(4) / 5.0) * num_atoms * K_B * ratio.powi(3)
}

/// Calculate Einstein heat capacity
///
/// C_V = 3Nk_B(Θ_E/T)² exp(Θ_E/T) / [exp(Θ_E/T) - 1]²
///
/// # Arguments
/// * `temperature_k` - Temperature (K)
/// * `einstein_temp` - Einstein temperature (K)
/// * `num_atoms` - Number of atoms
pub fn einstein_heat_capacity(temperature_k: f64, einstein_temp: f64, num_atoms: f64) -> f64 {
    let x = einstein_temp / temperature_k;
    let exp_x = x.exp();
    3.0 * num_atoms * K_B * x.powi(2) * exp_x / (exp_x - 1.0).powi(2)
}

/// Dulong-Petit law (high temperature limit)
///
/// C_V = 3Nk_B (constant)
///
/// # Arguments
/// * `num_atoms` - Number of atoms
///
/// # Returns
/// Molar heat capacity in J/(mol·K)
pub fn dulong_petit_heat_capacity(num_atoms: f64) -> f64 {
    3.0 * num_atoms * K_B
}

/// Calculate molar heat capacity at constant pressure
///
/// C_p = C_V + TVα²B
///
/// where α is thermal expansion coefficient, B is bulk modulus
///
/// # Arguments
/// * `cv` - Heat capacity at constant volume (J/(mol·K))
/// * `temperature_k` - Temperature (K)
/// * `molar_volume` - Molar volume (m³/mol)
/// * `thermal_expansion` - Volumetric thermal expansion coefficient (1/K)
/// * `bulk_modulus` - Bulk modulus (Pa)
pub fn heat_capacity_cp(
    cv: f64,
    temperature_k: f64,
    molar_volume: f64,
    thermal_expansion: f64,
    bulk_modulus: f64,
) -> f64 {
    cv + temperature_k * molar_volume * thermal_expansion.powi(2) * bulk_modulus
}

/// Calculate linear thermal expansion coefficient
///
/// α_L = (1/L)(dL/dT)
///
/// # Arguments
/// * `length1` - Initial length (m)
/// * `length2` - Final length (m)
/// * `temp1_k` - Initial temperature (K)
/// * `temp2_k` - Final temperature (K)
pub fn linear_thermal_expansion_coefficient(
    length1: f64,
    length2: f64,
    temp1_k: f64,
    temp2_k: f64,
) -> f64 {
    (length2 - length1) / (length1 * (temp2_k - temp1_k))
}

/// Calculate volumetric thermal expansion coefficient
///
/// α_V = (1/V)(dV/dT) ≈ 3α_L
pub fn volumetric_thermal_expansion_coefficient(linear_coefficient: f64) -> f64 {
    3.0 * linear_coefficient
}

/// Calculate thermal expansion at a given temperature
///
/// ΔL = L₀ · α · ΔT
///
/// # Returns
/// Change in length (m)
pub fn thermal_expansion_length(
    original_length: f64,
    thermal_expansion_coeff: f64,
    delta_temp_k: f64,
) -> f64 {
    original_length * thermal_expansion_coeff * delta_temp_k
}

/// Calculate new length after thermal expansion
pub fn length_after_expansion(
    original_length: f64,
    thermal_expansion_coeff: f64,
    delta_temp_k: f64,
) -> f64 {
    original_length * (1.0 + thermal_expansion_coeff * delta_temp_k)
}

/// Calculate thermal stress due to constrained expansion
///
/// σ = E · α · ΔT
///
/// # Arguments
/// * `youngs_modulus` - Young's modulus (Pa)
/// * `thermal_expansion_coeff` - Linear thermal expansion coefficient (1/K)
/// * `delta_temp_k` - Temperature change (K)
///
/// # Returns
/// Thermal stress (Pa)
pub fn thermal_stress_constrained(
    youngs_modulus: f64,
    thermal_expansion_coeff: f64,
    delta_temp_k: f64,
) -> f64 {
    youngs_modulus * thermal_expansion_coeff * delta_temp_k
}

/// Calculate thermal conductivity from Wiedemann-Franz law
///
/// κ = L·σ·T
///
/// where L is Lorenz number (2.44×10⁻⁸ W·Ω/K²)
///
/// # Arguments
/// * `electrical_conductivity` - Electrical conductivity (S/m)
/// * `temperature_k` - Temperature (K)
///
/// # Returns
/// Thermal conductivity (W/(m·K))
pub fn thermal_conductivity_wiedemann_franz(
    electrical_conductivity: f64,
    temperature_k: f64,
) -> f64 {
    const LORENZ_NUMBER: f64 = 2.44e-8; // W·Ω/K²
    LORENZ_NUMBER * electrical_conductivity * temperature_k
}

/// Calculate phonon thermal conductivity (simplified kinetic theory)
///
/// κ_ph = (1/3) C v l
///
/// # Arguments
/// * `heat_capacity_volume` - Volumetric heat capacity (J/(m³·K))
/// * `sound_velocity` - Average sound velocity (m/s)
/// * `mean_free_path` - Phonon mean free path (m)
pub fn phonon_thermal_conductivity(
    heat_capacity_volume: f64,
    sound_velocity: f64,
    mean_free_path: f64,
) -> f64 {
    (1.0 / 3.0) * heat_capacity_volume * sound_velocity * mean_free_path
}

/// Calculate thermal diffusivity
///
/// α = κ/(ρ·C_p)
///
/// # Arguments
/// * `thermal_conductivity` - Thermal conductivity (W/(m·K))
/// * `density` - Density (kg/m³)
/// * `specific_heat` - Specific heat capacity (J/(kg·K))
///
/// # Returns
/// Thermal diffusivity (m²/s)
pub fn thermal_diffusivity(thermal_conductivity: f64, density: f64, specific_heat: f64) -> f64 {
    thermal_conductivity / (density * specific_heat)
}

/// Calculate Biot number (for heat transfer analysis)
///
/// Bi = hL_c/κ
///
/// # Arguments
/// * `heat_transfer_coeff` - Convective heat transfer coefficient (W/(m²·K))
/// * `characteristic_length` - Characteristic length (m)
/// * `thermal_conductivity` - Thermal conductivity (W/(m·K))
pub fn biot_number(
    heat_transfer_coeff: f64,
    characteristic_length: f64,
    thermal_conductivity: f64,
) -> f64 {
    (heat_transfer_coeff * characteristic_length) / thermal_conductivity
}

/// Calculate Fourier number (dimensionless time for heat conduction)
///
/// Fo = αt/L²
///
/// # Arguments
/// * `thermal_diffusivity` - Thermal diffusivity (m²/s)
/// * `time` - Time (s)
/// * `characteristic_length` - Characteristic length (m)
pub fn fourier_number(thermal_diffusivity: f64, time: f64, characteristic_length: f64) -> f64 {
    (thermal_diffusivity * time) / (characteristic_length * characteristic_length)
}

/// Calculate Stefan-Boltzmann radiation heat transfer
///
/// q = εσA(T⁴ - T_surr⁴)
///
/// # Arguments
/// * `emissivity` - Surface emissivity (0-1)
/// * `area` - Surface area (m²)
/// * `temperature_k` - Surface temperature (K)
/// * `surrounding_temp_k` - Surrounding temperature (K)
///
/// # Returns
/// Heat transfer rate (W)
pub fn stefan_boltzmann_radiation(
    emissivity: f64,
    area: f64,
    temperature_k: f64,
    surrounding_temp_k: f64,
) -> f64 {
    const STEFAN_BOLTZMANN: f64 = 5.670374419e-8; // W/(m²·K⁴)
    emissivity * STEFAN_BOLTZMANN * area * (temperature_k.powi(4) - surrounding_temp_k.powi(4))
}

/// Calculate heat transfer through a slab (1D steady-state conduction)
///
/// q = κA(T₁ - T₂)/L
///
/// # Arguments
/// * `thermal_conductivity` - Thermal conductivity (W/(m·K))
/// * `area` - Cross-sectional area (m²)
/// * `temp1_k` - Temperature on side 1 (K)
/// * `temp2_k` - Temperature on side 2 (K)
/// * `thickness` - Slab thickness (m)
///
/// # Returns
/// Heat transfer rate (W)
pub fn heat_transfer_slab(
    thermal_conductivity: f64,
    area: f64,
    temp1_k: f64,
    temp2_k: f64,
    thickness: f64,
) -> f64 {
    thermal_conductivity * area * (temp1_k - temp2_k) / thickness
}

/// Calculate thermal resistance
///
/// R = L/(κA)
///
/// # Arguments
/// * `thickness` - Material thickness (m)
/// * `thermal_conductivity` - Thermal conductivity (W/(m·K))
/// * `area` - Cross-sectional area (m²)
///
/// # Returns
/// Thermal resistance (K/W)
pub fn thermal_resistance(thickness: f64, thermal_conductivity: f64, area: f64) -> f64 {
    thickness / (thermal_conductivity * area)
}

/// Calculate thermal interface resistance (Kapitza resistance)
///
/// Often important at material interfaces, especially at low temperatures
///
/// # Returns
/// Interface thermal resistance (K·m²/W)
pub fn kapitza_resistance_estimate(temperature_k: f64) -> f64 {
    // Simplified empirical relation: R_K ∝ T⁻³
    // Typical values: 10⁻⁴ to 10⁻⁸ K·m²/W
    1e-6 / temperature_k.powi(3)
}

/// Calculate specific heat from temperature dependence (polynomial fit)
///
/// C_p(T) = a + bT + cT² + dT³
///
/// Common for experimental data fitting
pub fn specific_heat_polynomial(
    temperature_k: f64,
    a: f64,
    b: f64,
    c: f64,
    d: f64,
) -> f64 {
    a + b * temperature_k + c * temperature_k.powi(2) + d * temperature_k.powi(3)
}

/// Calculate Grüneisen parameter (relates thermal and mechanical properties)
///
/// γ = (3α_V B V)/(C_V)
///
/// # Arguments
/// * `volumetric_expansion` - Volumetric thermal expansion coefficient (1/K)
/// * `bulk_modulus` - Bulk modulus (Pa)
/// * `molar_volume` - Molar volume (m³/mol)
/// * `heat_capacity_cv` - Heat capacity at constant volume (J/(mol·K))
pub fn gruneisen_parameter(
    volumetric_expansion: f64,
    bulk_modulus: f64,
    molar_volume: f64,
    heat_capacity_cv: f64,
) -> f64 {
    (3.0 * volumetric_expansion * bulk_modulus * molar_volume) / heat_capacity_cv
}

/// Calculate temperature rise from adiabatic compression/expansion
///
/// ΔT = T(γ - 1)(ΔV/V)
///
/// # Arguments
/// * `initial_temp_k` - Initial temperature (K)
/// * `gruneisen` - Grüneisen parameter
/// * `volume_change_fraction` - ΔV/V (fractional volume change)
pub fn adiabatic_temperature_change(
    initial_temp_k: f64,
    gruneisen: f64,
    volume_change_fraction: f64,
) -> f64 {
    initial_temp_k * (gruneisen - 1.0) * volume_change_fraction
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dulong_petit() {
        // For 1 mole of atoms (N_A atoms)
        let cv = dulong_petit_heat_capacity(N_A);
        let expected = 3.0 * R; // 3R ≈ 24.9 J/(mol·K)
        assert!((cv - expected).abs() < 0.1);
    }

    #[test]
    fn test_thermal_expansion() {
        // Steel: α ≈ 12×10⁻⁶ /K
        let alpha = 12e-6;
        let original_length = 1.0; // 1 meter
        let delta_t = 100.0; // 100K temperature rise

        let delta_l = thermal_expansion_length(original_length, alpha, delta_t);
        assert!((delta_l - 0.0012).abs() < 1e-6); // 1.2 mm expansion

        let new_length = length_after_expansion(original_length, alpha, delta_t);
        assert!((new_length - 1.0012).abs() < 1e-6);
    }

    #[test]
    fn test_thermal_stress() {
        // Steel constrained expansion
        let e = 200e9; // 200 GPa
        let alpha = 12e-6; // /K
        let delta_t = 100.0; // K

        let stress = thermal_stress_constrained(e, alpha, delta_t);
        assert!((stress - 240e6).abs() < 1e6); // 240 MPa
    }

    #[test]
    fn test_thermal_diffusivity() {
        // Copper: κ ≈ 400 W/(m·K), ρ ≈ 8960 kg/m³, C_p ≈ 385 J/(kg·K)
        let kappa = 400.0;
        let rho = 8960.0;
        let cp = 385.0;

        let alpha = thermal_diffusivity(kappa, rho, cp);
        assert!(alpha > 1e-4 && alpha < 2e-4); // ~1.16×10⁻⁴ m²/s
    }

    #[test]
    fn test_wiedemann_franz() {
        // Copper at 300K: σ ≈ 6×10⁷ S/m
        let sigma = 6e7;
        let temp = 300.0;

        let kappa = thermal_conductivity_wiedemann_franz(sigma, temp);
        assert!(kappa > 400.0 && kappa < 500.0); // Should be ~440 W/(m·K)
    }

    #[test]
    fn test_heat_transfer_slab() {
        // 1 W/(m·K) material, 1 m² area, 1 m thick, 100K difference
        let q = heat_transfer_slab(1.0, 1.0, 400.0, 300.0, 1.0);
        assert!((q - 100.0).abs() < 0.01); // 100 W
    }

    #[test]
    fn test_biot_number() {
        // h = 100 W/(m²·K), L = 0.01 m, κ = 50 W/(m·K)
        let bi = biot_number(100.0, 0.01, 50.0);
        assert!((bi - 0.02).abs() < 0.001); // Bi << 1 (lumped capacitance valid)
    }

    #[test]
    fn test_stefan_boltzmann() {
        // Black body (ε=1) at 1000K, 1 m², surrounding at 300K
        let q = stefan_boltzmann_radiation(1.0, 1.0, 1000.0, 300.0);
        assert!(q > 50000.0 && q < 60000.0); // ~56 kW
    }
}
