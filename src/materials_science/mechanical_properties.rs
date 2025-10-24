//! Mechanical Properties of Materials
//!
//! Stress, strain, elasticity, plasticity, fracture mechanics,
//! and material failure criteria.

use super::Matrix3x3;

/// Calculate engineering stress
///
/// σ = F/A₀
///
/// # Arguments
/// * `force` - Applied force (N)
/// * `original_area` - Original cross-sectional area (m²)
///
/// # Returns
/// Engineering stress (Pa)
pub fn engineering_stress(force: f64, original_area: f64) -> f64 {
    force / original_area
}

/// Calculate engineering strain
///
/// ε = ΔL/L₀
pub fn engineering_strain(original_length: f64, final_length: f64) -> f64 {
    (final_length - original_length) / original_length
}

/// Calculate true stress
///
/// σ_true = F/A (using current area)
pub fn true_stress(force: f64, current_area: f64) -> f64 {
    force / current_area
}

/// Calculate true strain
///
/// ε_true = ln(L/L₀)
pub fn true_strain(original_length: f64, final_length: f64) -> f64 {
    (final_length / original_length).ln()
}

/// Calculate Young's modulus from stress-strain data
///
/// E = σ/ε (in elastic region)
pub fn youngs_modulus(stress: f64, strain: f64) -> f64 {
    stress / strain
}

/// Calculate shear modulus
///
/// G = τ/γ
pub fn shear_modulus(shear_stress: f64, shear_strain: f64) -> f64 {
    shear_stress / shear_strain
}

/// Calculate bulk modulus
///
/// K = -V(dP/dV) = stress / volumetric_strain
pub fn bulk_modulus(pressure_change: f64, volumetric_strain: f64) -> f64 {
    pressure_change / (-volumetric_strain)
}

/// Calculate Poisson's ratio
///
/// ν = -ε_lateral / ε_axial
pub fn poissons_ratio(lateral_strain: f64, axial_strain: f64) -> f64 {
    -lateral_strain / axial_strain
}

/// Relationship between elastic constants
///
/// E = 2G(1 + ν) = 3K(1 - 2ν)
///
/// Calculate Young's modulus from shear modulus and Poisson's ratio
pub fn youngs_from_shear_poisson(shear_modulus: f64, poissons_ratio: f64) -> f64 {
    2.0 * shear_modulus * (1.0 + poissons_ratio)
}

/// Calculate shear modulus from Young's modulus and Poisson's ratio
pub fn shear_from_youngs_poisson(youngs_modulus: f64, poissons_ratio: f64) -> f64 {
    youngs_modulus / (2.0 * (1.0 + poissons_ratio))
}

/// Calculate bulk modulus from Young's modulus and Poisson's ratio
pub fn bulk_from_youngs_poisson(youngs_modulus: f64, poissons_ratio: f64) -> f64 {
    youngs_modulus / (3.0 * (1.0 - 2.0 * poissons_ratio))
}

/// Calculate strain energy density
///
/// U = (1/2)σε = σ²/(2E)
pub fn strain_energy_density(stress: f64, strain: f64) -> f64 {
    0.5 * stress * strain
}

/// Calculate toughness (area under stress-strain curve up to fracture)
///
/// Approximated as: K = (1/2)(σ_f + σ_y)ε_f
pub fn toughness(yield_stress: f64, fracture_stress: f64, fracture_strain: f64) -> f64 {
    0.5 * (fracture_stress + yield_stress) * fracture_strain
}

/// Calculate resilience (energy absorbed in elastic region)
///
/// U_r = σ_y²/(2E)
pub fn resilience(yield_stress: f64, youngs_modulus: f64) -> f64 {
    (yield_stress * yield_stress) / (2.0 * youngs_modulus)
}

/// Calculate hardness from yield strength (approximation)
///
/// H ≈ 3σ_y (for metals)
pub fn hardness_from_yield(yield_stress: f64) -> f64 {
    3.0 * yield_stress
}

/// Calculate stress concentration factor for elliptical hole
///
/// K_t = 1 + 2(a/b)
///
/// where a is major axis, b is minor axis
pub fn stress_concentration_ellipse(major_axis: f64, minor_axis: f64) -> f64 {
    1.0 + 2.0 * (major_axis / minor_axis)
}

/// Calculate stress concentration factor for circular hole
pub fn stress_concentration_circular() -> f64 {
    3.0 // K_t = 3 for circular hole
}

/// Calculate stress intensity factor (Mode I)
///
/// K_I = Y σ √(πa)
///
/// # Arguments
/// * `geometry_factor` - Y (depends on geometry, typically 1.0-1.2)
/// * `applied_stress` - Applied stress (Pa)
/// * `crack_length` - Crack length (m)
pub fn stress_intensity_factor(
    geometry_factor: f64,
    applied_stress: f64,
    crack_length: f64,
) -> f64 {
    geometry_factor * applied_stress * (std::f64::consts::PI * crack_length).sqrt()
}

/// Check Griffith fracture criterion
///
/// Fracture occurs when K_I ≥ K_Ic (fracture toughness)
pub fn check_fracture_griffith(stress_intensity: f64, fracture_toughness: f64) -> bool {
    stress_intensity >= fracture_toughness
}

/// Calculate critical crack length
///
/// a_c = (1/π)(K_Ic/(Yσ))²
pub fn critical_crack_length(
    fracture_toughness: f64,
    applied_stress: f64,
    geometry_factor: f64,
) -> f64 {
    let ratio = fracture_toughness / (geometry_factor * applied_stress);
    (ratio * ratio) / std::f64::consts::PI
}

/// von Mises stress (equivalent stress for multiaxial loading)
///
/// σ_vm = √[(σ₁-σ₂)² + (σ₂-σ₃)² + (σ₃-σ₁)²]/√2
pub fn von_mises_stress(sigma1: f64, sigma2: f64, sigma3: f64) -> f64 {
    let term1 = (sigma1 - sigma2).powi(2);
    let term2 = (sigma2 - sigma3).powi(2);
    let term3 = (sigma3 - sigma1).powi(2);
    ((term1 + term2 + term3) / 2.0).sqrt()
}

/// Check von Mises yield criterion
///
/// Yields when σ_vm ≥ σ_y
pub fn check_yield_von_mises(von_mises_stress: f64, yield_stress: f64) -> bool {
    von_mises_stress >= yield_stress
}

/// Tresca stress (maximum shear stress)
///
/// τ_max = (σ_max - σ_min)/2
pub fn tresca_stress(sigma_max: f64, sigma_min: f64) -> f64 {
    (sigma_max - sigma_min) / 2.0
}

/// Calculate Brinell hardness number
///
/// BHN = 2F/(πD(D - √(D² - d²)))
///
/// # Arguments
/// * `force` - Applied force (kgf)
/// * `indenter_diameter` - Indenter ball diameter (mm)
/// * `indent_diameter` - Diameter of indent (mm)
pub fn brinell_hardness(force_kgf: f64, indenter_diameter_mm: f64, indent_diameter_mm: f64) -> f64 {
    let d = indenter_diameter_mm;
    let d_indent = indent_diameter_mm;
    (2.0 * force_kgf) / (std::f64::consts::PI * d * (d - (d * d - d_indent * d_indent).sqrt()))
}

/// Convert between hardness scales (approximation)
///
/// Rockwell C to Brinell: BHN ≈ 14.68 HRC + 215.6 (for steels)
pub fn rockwell_c_to_brinell(hrc: f64) -> f64 {
    14.68 * hrc + 215.6
}

/// Calculate creep rate (power law creep)
///
/// dε/dt = A σⁿ exp(-Q/(RT))
///
/// # Arguments
/// * `stress` - Applied stress (Pa)
/// * `temperature` - Temperature (K)
/// * `activation_energy` - Activation energy (J/mol)
/// * `stress_exponent` - n (typically 3-8)
/// * `prefactor` - A (material constant)
pub fn creep_rate_power_law(
    stress: f64,
    temperature_k: f64,
    activation_energy: f64,
    stress_exponent: f64,
    prefactor: f64,
) -> f64 {
    const R: f64 = 8.314; // Gas constant
    prefactor * stress.powf(stress_exponent) * (-activation_energy / (R * temperature_k)).exp()
}

/// Calculate fatigue life (Basquin equation)
///
/// Δσ = σ_f' (2N_f)^b
///
/// Solving for N_f: N_f = (1/2)(Δσ/σ_f')^(1/b)
pub fn fatigue_life_basquin(
    stress_amplitude: f64,
    fatigue_strength_coeff: f64,
    fatigue_strength_exp: f64,
) -> f64 {
    0.5 * (stress_amplitude / fatigue_strength_coeff).powf(1.0 / fatigue_strength_exp)
}

/// Calculate endurance limit for steels
///
/// S_e ≈ 0.5 σ_UTS (for σ_UTS < 1400 MPa)
pub fn endurance_limit_steel(ultimate_tensile_strength: f64) -> f64 {
    0.5 * ultimate_tensile_strength
}

/// Calculate percent elongation
///
/// %EL = (L_f - L_0)/L_0 × 100%
pub fn percent_elongation(original_length: f64, final_length: f64) -> f64 {
    ((final_length - original_length) / original_length) * 100.0
}

/// Calculate percent reduction in area
///
/// %RA = (A_0 - A_f)/A_0 × 100%
pub fn percent_reduction_area(original_area: f64, final_area: f64) -> f64 {
    ((original_area - final_area) / original_area) * 100.0
}

/// Calculate ductility from percent elongation and reduction in area
///
/// Higher values indicate more ductile material
pub fn ductility_index(percent_elongation: f64, percent_reduction_area: f64) -> f64 {
    (percent_elongation + percent_reduction_area) / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stress_strain() {
        let stress = engineering_stress(1000.0, 0.001); // 1000N on 1mm²
        assert_eq!(stress, 1e6); // 1 MPa

        let strain = engineering_strain(100.0, 101.0);
        assert_eq!(strain, 0.01); // 1% strain
    }

    #[test]
    fn test_elastic_constants() {
        let e = youngs_from_shear_poisson(80e9, 0.3);
        assert!((e - 208e9).abs() < 1e9);

        let g = shear_from_youngs_poisson(200e9, 0.3);
        assert!((g - 76.9e9).abs() < 1e9);
    }

    #[test]
    fn test_von_mises() {
        // Uniaxial tension: σ₁=σ, σ₂=σ₃=0
        let sigma_vm = von_mises_stress(100e6, 0.0, 0.0);
        assert!((sigma_vm - 100e6).abs() < 1e3);
    }

    #[test]
    fn test_stress_concentration() {
        let kt_circle = stress_concentration_circular();
        assert_eq!(kt_circle, 3.0);

        let kt_ellipse = stress_concentration_ellipse(10.0, 1.0);
        assert_eq!(kt_ellipse, 21.0); // 1 + 2*10
    }

    #[test]
    fn test_hardness_conversion() {
        let bhn = rockwell_c_to_brinell(50.0);
        assert!((bhn - 950.0).abs() < 50.0); // Approximate
    }
}
