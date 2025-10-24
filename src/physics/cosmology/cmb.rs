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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_planck_spectrum() {
        let t_cmb = 2.7255;

        // Peak should be at microwave frequencies (~160 GHz)
        let freq_peak = cmb_peak_frequency(t_cmb);
        assert!(freq_peak > 1e11 && freq_peak < 2e11);

        // Spectrum at peak should be non-zero
        let intensity = planck_spectrum(freq_peak, t_cmb);
        assert!(intensity > 0.0);
        assert!(intensity.is_finite());
    }

    #[test]
    fn test_cmb_peak_frequency() {
        let params = CosmologyParams::planck_2018();
        let freq = cmb_peak_frequency(params.t_cmb);

        // Should be ~160 GHz
        assert!(freq > 1.4e11 && freq < 1.8e11);
    }

    #[test]
    fn test_photon_number_density() {
        let params = CosmologyParams::planck_2018();
        let n_gamma = photon_number_density(params.t_cmb);

        // CMB photon density should be ~400 photons/cm³ = 4e8 photons/m³
        assert!(n_gamma > 1e8 && n_gamma < 1e9);
    }

    #[test]
    fn test_cmb_energy_density() {
        let params = CosmologyParams::planck_2018();
        let rho_gamma = cmb_energy_density(params.t_cmb);

        // Should be finite and positive
        assert!(rho_gamma > 0.0);
        assert!(rho_gamma.is_finite());
    }

    #[test]
    fn test_recombination() {
        let params = CosmologyParams::planck_2018();
        let recomb = recombination(&params);

        // Recombination at z ~ 1100
        assert!((recomb.redshift - 1100.0).abs() < 1.0);

        // Temperature ~ 3000 K
        assert!(recomb.temperature > 2900.0 && recomb.temperature < 3100.0);

        // Scale factor ~ 1/1101
        assert!(recomb.scale_factor > 0.0009 && recomb.scale_factor < 0.001);
    }

    #[test]
    fn test_matter_radiation_equality() {
        let params = CosmologyParams::planck_2018();
        let z_eq = matter_radiation_equality(&params);

        // Matter-radiation equality at z ~ 3400
        assert!(z_eq > 3000.0 && z_eq < 4000.0);
    }

    #[test]
    fn test_sound_horizon() {
        let params = CosmologyParams::planck_2018();
        let r_s = sound_horizon_recombination(&params);

        // Sound horizon should be positive and finite
        assert!(r_s > 0.0);
        assert!(r_s.is_finite());
    }

    #[test]
    fn test_first_acoustic_peak() {
        let params = CosmologyParams::planck_2018();
        let theta = first_acoustic_peak_angle(&params);

        // First peak should be a positive, finite angle
        assert!(theta > 0.0);
        assert!(theta.is_finite());
    }

    #[test]
    fn test_sachs_wolfe() {
        let potential = 1e-5; // Typical gravitational potential
        let dt_t = sachs_wolfe_amplitude(potential);

        // Temperature fluctuation should be tiny
        assert!(dt_t.abs() < 1e-10);
    }
}
