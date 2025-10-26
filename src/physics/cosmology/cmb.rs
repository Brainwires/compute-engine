//! Cosmic Microwave Background (CMB) Physics

use super::{CosmologyParams, Redshift, C, K_B};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Planck blackbody spectrum
/// B_ν(T) = (2hν³/c²) / (exp(hν/kT) - 1)
pub fn planck_spectrum(frequency: f64, temperature: f64) -> f64 {
    const H: f64 = 6.62607015e-34; // Planck constant

    let x = H * frequency / (K_B * temperature);

    if x > 100.0 {
        return 0.0; // Exponential cutoff
    }

    let numerator = 2.0 * H * frequency.powi(3) / (C * C);
    let denominator = x.exp() - 1.0;

    numerator / denominator
}

/// Peak frequency of CMB spectrum (Wien's law)
/// ν_peak = 2.821 kT/h
pub fn cmb_peak_frequency(temperature: f64) -> f64 {
    const H: f64 = 6.62607015e-34;
    2.821 * K_B * temperature / H
}

/// Number density of CMB photons
/// n_γ = (2ζ(3)/π²) (kT/ℏc)³
pub fn photon_number_density(temperature: f64) -> f64 {
    const HBAR: f64 = 1.054571817e-34;
    const ZETA_3: f64 = 1.202056903; // Riemann zeta(3)

    let coeff = 2.0 * ZETA_3 / (PI * PI);
    let kt_over_hbar_c = K_B * temperature / (HBAR * C);

    coeff * kt_over_hbar_c.powi(3)
}

/// CMB energy density
/// ρ_γ = (π²/15) (kT)⁴/(ℏc)³
pub fn cmb_energy_density(temperature: f64) -> f64 {
    const HBAR: f64 = 1.054571817e-34;

    let coeff = PI * PI / 15.0;
    let kt_over_hbar_c = K_B * temperature / (HBAR * C);

    coeff * (K_B * temperature) * kt_over_hbar_c.powi(3)
}

/// Recombination epoch calculations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RecombinationEpoch {
    pub redshift: f64,
    pub temperature: f64,
    pub age: f64,
    pub scale_factor: f64,
}

/// Calculate recombination parameters
pub fn recombination(params: &CosmologyParams) -> RecombinationEpoch {
    // Recombination occurs at z ~ 1100, T ~ 3000 K
    let z_recomb = 1100.0;
    let z = Redshift(z_recomb);

    RecombinationEpoch {
        redshift: z_recomb,
        temperature: z.temperature(params),
        age: z.age(params),
        scale_factor: z.to_scale_factor(),
    }
}

/// Matter-radiation equality
pub fn matter_radiation_equality(params: &CosmologyParams) -> f64 {
    // z_eq = Ω_m / Ω_r - 1
    params.omega_m / params.omega_r - 1.0
}

/// Sound horizon at recombination (characteristic scale for CMB fluctuations)
pub fn sound_horizon_recombination(params: &CosmologyParams) -> f64 {
    // Simplified: r_s ~ c / (H₀ √(Ω_m z_recomb))
    let z_recomb = 1100.0;
    let h0 = params.hubble_constant_si();

    C / (h0 * (params.omega_m * z_recomb).sqrt())
}

/// Angular scale of first acoustic peak in CMB
pub fn first_acoustic_peak_angle(params: &CosmologyParams) -> f64 {
    let z_recomb = Redshift(1100.0);
    let r_s = sound_horizon_recombination(params);
    let d_a = super::friedmann::angular_diameter_distance(z_recomb, params);

    // θ = r_s / D_A (in radians)
    r_s / d_a
}

/// Sachs-Wolfe effect (temperature fluctuation from gravitational potential)
pub fn sachs_wolfe_amplitude(potential: f64) -> f64 {
    // ΔT/T ~ Φ/c²
    potential / (C * C)
}

