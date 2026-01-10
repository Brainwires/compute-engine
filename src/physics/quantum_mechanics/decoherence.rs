//! Decoherence Calculations for Quantum-Classical Transition
//!
//! Implements decoherence length scales and rates relevant to quantum hydrodynamics:
//! - Thermal de Broglie wavelength (decoherence length scale)
//! - Joos-Zeh decoherence rate formula
//! - Classical validity regime estimates
//!
//! These calculations support the decoherence-corrected Navier-Stokes framework
//! for analyzing the quantum-classical transition in fluid dynamics.

use super::{H_BAR, K_B};
use serde::{Deserialize, Serialize};

/// Result of decoherence scale calculation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DecoherenceScaleResult {
    /// Thermal de Broglie wavelength (decoherence length scale) in meters
    pub decoherence_length: f64,

    /// Decoherence time scale in seconds (if scattering rate provided)
    pub decoherence_time: Option<f64>,

    /// Classical validity assessment
    pub classical_valid: bool,

    /// Ratio of decoherence length to typical molecular spacing
    pub length_ratio: f64,
}

/// Parameters for decoherence scale calculation
#[derive(Debug, Clone, Deserialize)]
pub struct DecoherenceParams {
    /// Particle/molecule mass in kg
    pub mass: f64,

    /// Temperature in Kelvin
    pub temperature: f64,

    /// Optional: scattering/collision rate in Hz (for decoherence time)
    pub scattering_rate: Option<f64>,

    /// Optional: typical length scale for comparison (default: 3e-10 m for molecular spacing)
    pub reference_length: Option<f64>,
}

/// Calculate thermal de Broglie wavelength (decoherence length scale)
///
/// Formula: L_D = ℏ / √(2 m k_B T)
///
/// This is the characteristic length scale below which quantum coherence effects
/// are significant. For scales L >> L_D, classical behavior dominates.
///
/// # Arguments
/// * `mass` - Particle mass in kg
/// * `temperature` - Temperature in Kelvin
///
/// # Returns
/// Decoherence length in meters
///
/// # Example
/// For water molecule (m ≈ 3×10⁻²⁶ kg) at T = 300K:
/// L_D ≈ 1.054×10⁻³⁴ / √(2 × 3×10⁻²⁶ × 1.38×10⁻²³ × 300)
///     ≈ 10⁻¹¹ m (0.1 Å)
pub fn decoherence_length(mass: f64, temperature: f64) -> f64 {
    if mass <= 0.0 || temperature <= 0.0 {
        return f64::NAN;
    }

    // L_D = ℏ / √(2 m k_B T)
    H_BAR / (2.0 * mass * K_B * temperature).sqrt()
}

/// Calculate full decoherence scale analysis
///
/// Provides decoherence length, time (if scattering rate known), and
/// assessment of classical validity regime.
pub fn calculate_decoherence_scale(params: DecoherenceParams) -> Result<DecoherenceScaleResult, String> {
    if params.mass <= 0.0 {
        return Err("Mass must be positive".to_string());
    }
    if params.temperature <= 0.0 {
        return Err("Temperature must be positive".to_string());
    }

    let l_d = decoherence_length(params.mass, params.temperature);

    // Calculate decoherence time if scattering rate provided
    // τ_D ~ τ_R × (λ_th / Δx)² where τ_R = 1/scattering_rate
    let decoherence_time = params.scattering_rate.map(|gamma| {
        // For localization decoherence: τ_D ≈ 1 / (γ × (Δx/λ_th)²)
        // Assuming Δx ~ reference_length
        let ref_length = params.reference_length.unwrap_or(3e-10);
        1.0 / (gamma * (ref_length / l_d).powi(2))
    });

    // Reference length for comparison (default: typical molecular spacing ~3 Å)
    let ref_length = params.reference_length.unwrap_or(3e-10);
    let length_ratio = l_d / ref_length;

    // Classical validity: L_D << molecular spacing means quantum effects negligible
    // Typically valid when L_D < 0.1 × molecular spacing
    let classical_valid = length_ratio < 0.1;

    Ok(DecoherenceScaleResult {
        decoherence_length: l_d,
        decoherence_time,
        classical_valid,
        length_ratio,
    })
}

/// Joos-Zeh decoherence rate for environmental scattering
///
/// Formula: Γ_D = Γ_scatt × (Δx / λ_th)²
///
/// This gives the rate at which spatial superpositions decohere due to
/// environmental scattering (photons, air molecules, etc.)
///
/// # Arguments
/// * `scattering_rate` - Rate of environmental scattering events (Hz)
/// * `separation` - Spatial separation of superposition (m)
/// * `thermal_wavelength` - Thermal de Broglie wavelength (m)
///
/// # Returns
/// Decoherence rate in Hz
pub fn joos_zeh_decoherence_rate(
    scattering_rate: f64,
    separation: f64,
    thermal_wavelength: f64,
) -> f64 {
    if thermal_wavelength <= 0.0 {
        return f64::NAN;
    }
    scattering_rate * (separation / thermal_wavelength).powi(2)
}

/// Calculate quantum correction magnitude ratio
///
/// Estimates |F_Q| / |viscous terms| as function of scale L
///
/// F_Q ~ ℏ²/(m) × ∇⁴ρ / ρ ~ ℏ²/(m L⁴) × Δρ/ρ
/// Viscous ~ μ ∇²u ~ μ U / L²
///
/// Ratio ~ (ℏ²/m μ) × (1/L²) × (Δρ/ρ U)
///
/// # Arguments
/// * `length_scale` - Characteristic length scale L (m)
/// * `mass` - Particle mass (kg)
/// * `viscosity` - Dynamic viscosity (Pa·s)
/// * `density_variation` - Δρ/ρ (dimensionless, typically ~0.01-0.1)
/// * `velocity_scale` - Characteristic velocity U (m/s)
///
/// # Returns
/// Ratio of quantum to viscous forces (dimensionless)
pub fn quantum_viscous_ratio(
    length_scale: f64,
    mass: f64,
    viscosity: f64,
    density_variation: f64,
    velocity_scale: f64,
) -> f64 {
    if length_scale <= 0.0 || mass <= 0.0 || viscosity <= 0.0 || velocity_scale <= 0.0 {
        return f64::NAN;
    }

    let hbar_sq = H_BAR * H_BAR;
    (hbar_sq / (mass * viscosity)) * (1.0 / (length_scale * length_scale)) * (density_variation / velocity_scale)
}

/// Kolmogorov dissipation scale
///
/// Formula: η = (ν³/ε)^(1/4)
///
/// The smallest scale of turbulent motion where viscous dissipation dominates.
///
/// # Arguments
/// * `kinematic_viscosity` - ν = μ/ρ (m²/s)
/// * `dissipation_rate` - ε, turbulent kinetic energy dissipation rate (m²/s³)
///
/// # Returns
/// Kolmogorov scale in meters
pub fn kolmogorov_scale(kinematic_viscosity: f64, dissipation_rate: f64) -> f64 {
    if kinematic_viscosity <= 0.0 || dissipation_rate <= 0.0 {
        return f64::NAN;
    }

    (kinematic_viscosity.powi(3) / dissipation_rate).powf(0.25)
}

/// Compare decoherence length to Kolmogorov scale
///
/// For the research paper's key prediction: For normal fluids, L_D << η,
/// so classical N-S is always valid before molecular effects dominate.
///
/// # Arguments
/// * `mass` - Particle mass (kg)
/// * `temperature` - Temperature (K)
/// * `kinematic_viscosity` - ν (m²/s)
/// * `dissipation_rate` - ε (m²/s³)
///
/// # Returns
/// Tuple of (L_D, η, L_D/η ratio)
pub fn decoherence_vs_kolmogorov(
    mass: f64,
    temperature: f64,
    kinematic_viscosity: f64,
    dissipation_rate: f64,
) -> (f64, f64, f64) {
    let l_d = decoherence_length(mass, temperature);
    let eta = kolmogorov_scale(kinematic_viscosity, dissipation_rate);
    let ratio = l_d / eta;

    (l_d, eta, ratio)
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER_MOLECULE_MASS: f64 = 2.99e-26; // kg (H2O)
    const ROOM_TEMP: f64 = 300.0; // K

    #[test]
    fn test_decoherence_length_water() {
        // For water at room temperature, L_D should be ~10^-11 m
        let l_d = decoherence_length(WATER_MOLECULE_MASS, ROOM_TEMP);

        assert!(l_d > 1e-12, "L_D should be > 10^-12 m");
        assert!(l_d < 1e-10, "L_D should be < 10^-10 m");

        // More precise check: should be approximately 1e-11 m
        let expected = H_BAR / (2.0 * WATER_MOLECULE_MASS * K_B * ROOM_TEMP).sqrt();
        assert!((l_d - expected).abs() < 1e-15);
    }

    #[test]
    fn test_decoherence_length_helium() {
        // He-4 atom: m = 6.646e-27 kg
        // At T = 2K (superfluid regime), should have larger L_D
        let helium_mass = 6.646e-27;
        let superfluid_temp = 2.0;

        let l_d = decoherence_length(helium_mass, superfluid_temp);

        // Should be much larger than water at room temp
        let l_d_water = decoherence_length(WATER_MOLECULE_MASS, ROOM_TEMP);
        assert!(l_d > l_d_water * 10.0, "He-4 at 2K should have L_D >> water at 300K");
    }

    #[test]
    fn test_decoherence_scale_full() {
        let params = DecoherenceParams {
            mass: WATER_MOLECULE_MASS,
            temperature: ROOM_TEMP,
            scattering_rate: Some(1e12), // typical collision rate in liquid
            reference_length: Some(3e-10), // molecular spacing
        };

        let result = calculate_decoherence_scale(params).unwrap();

        // Classical should be valid for water at room temp
        assert!(result.classical_valid, "Classical physics should be valid for water at 300K");
        assert!(result.length_ratio < 0.1, "L_D/molecular_spacing should be << 1");
    }

    #[test]
    fn test_joos_zeh_rate() {
        let scatt_rate = 1e10; // Hz
        let separation = 1e-6; // 1 micron
        let lambda_th = 1e-11; // m

        let gamma_d = joos_zeh_decoherence_rate(scatt_rate, separation, lambda_th);

        // (1e-6 / 1e-11)^2 = 1e10, so Γ_D ~ 1e10 × 1e10 = 1e20 Hz
        assert!(gamma_d > 1e19);
        assert!(gamma_d < 1e21);
    }

    #[test]
    fn test_kolmogorov_scale() {
        // Water: ν ≈ 1e-6 m²/s
        // Moderate turbulence: ε ≈ 0.1 m²/s³
        let nu = 1e-6;
        let epsilon = 0.1;

        let eta = kolmogorov_scale(nu, epsilon);

        // η = (10^-18 / 0.1)^0.25 = (10^-17)^0.25 ≈ 5.6e-5 m
        assert!(eta > 1e-5, "Kolmogorov scale should be > 10 μm");
        assert!(eta < 1e-3, "Kolmogorov scale should be < 1 mm");
    }

    #[test]
    fn test_decoherence_vs_kolmogorov() {
        // Water at room temp with moderate turbulence
        let (l_d, eta, ratio) = decoherence_vs_kolmogorov(
            WATER_MOLECULE_MASS,
            ROOM_TEMP,
            1e-6,  // ν
            0.1,   // ε
        );

        // L_D ~ 10^-11 m, η ~ 10^-5 m, ratio ~ 10^-6
        assert!(ratio < 1e-4, "L_D should be << η for water");

        // This confirms the paper's prediction: classical N-S valid at all observable scales
        eprintln!("L_D = {:.2e} m, η = {:.2e} m, ratio = {:.2e}", l_d, eta, ratio);
    }

    #[test]
    fn test_invalid_inputs() {
        assert!(decoherence_length(-1.0, 300.0).is_nan());
        assert!(decoherence_length(1e-26, -100.0).is_nan());
        assert!(decoherence_length(0.0, 300.0).is_nan());

        assert!(calculate_decoherence_scale(DecoherenceParams {
            mass: -1.0,
            temperature: 300.0,
            scattering_rate: None,
            reference_length: None,
        }).is_err());
    }
}
