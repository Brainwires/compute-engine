//! Magnetic Confinement - Tokamaks and Stellarators

use super::{PlasmaParams, E, K_B};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Tokamak configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TokamakConfig {
    pub major_radius: f64,  // R (m)
    pub minor_radius: f64,  // a (m)
    pub toroidal_field: f64, // B_T (T)
    pub plasma_current: f64, // I_p (A)
    pub elongation: f64,    // κ
    pub triangularity: f64, // δ
}

/// Safety factor q(r): measures field line winding
/// q = (r B_T) / (R B_p)
pub fn safety_factor(config: &TokamakConfig, r: f64) -> f64 {
    const MU_0: f64 = 4.0e-7 * PI;

    // Poloidal field from plasma current
    let b_p = (MU_0 * config.plasma_current) / (2.0 * PI * config.major_radius);

    if b_p.abs() < 1e-10 {
        return f64::INFINITY;
    }

    (r * config.toroidal_field) / (config.major_radius * b_p)
}

/// Aspect ratio: A = R/a
pub fn aspect_ratio(config: &TokamakConfig) -> f64 {
    config.major_radius / config.minor_radius
}

/// Plasma volume for tokamak: V = 2π²Raκ
pub fn plasma_volume(config: &TokamakConfig) -> f64 {
    2.0 * PI * PI * config.major_radius * config.minor_radius * config.minor_radius
        * config.elongation
}

/// Plasma surface area: S = 4π²Ra
pub fn plasma_surface_area(config: &TokamakConfig) -> f64 {
    4.0 * PI * PI * config.major_radius * config.minor_radius
}

/// Toroidal magnetic field energy
pub fn magnetic_field_energy(config: &TokamakConfig) -> f64 {
    const MU_0: f64 = 4.0e-7 * PI;
    let volume = plasma_volume(config);
    let b_sq = config.toroidal_field * config.toroidal_field;

    (b_sq / (2.0 * MU_0)) * volume
}

/// Kruskal-Shafranov limit for stability: q(a) > 2
pub fn kruskal_shafranov_stable(config: &TokamakConfig) -> bool {
    let q_edge = safety_factor(config, config.minor_radius);
    q_edge > 2.0
}

/// Troyon beta limit: β < β_N I_p / (a B_T)
/// where β_N ~ 2.8-3.5 for tokamaks
pub fn troyon_beta_limit(config: &TokamakConfig, beta_n: f64) -> f64 {
    let ip_ma = config.plasma_current / 1e6; // Convert to MA
    beta_n * ip_ma / (config.minor_radius * config.toroidal_field)
}

/// Greenwald density limit: n_G = I_p / (πa²)
pub fn greenwald_density_limit(config: &TokamakConfig) -> f64 {
    let ip_ma = config.plasma_current / 1e6; // MA
    let area = PI * config.minor_radius * config.minor_radius; // m²

    ip_ma * 1e20 / area // Return in m⁻³
}

/// Energy confinement time scaling (ITER98y2)
/// τ_E ∝ I_p^0.93 B_T^0.15 P^-0.69 n^0.41 M^0.19 R^1.97 ε^0.58 κ^0.78
pub fn iter98_confinement_time(
    config: &TokamakConfig,
    params: &PlasmaParams,
    heating_power: f64,
) -> f64 {
    let ip_ma = config.plasma_current / 1e6;
    let p_mw = heating_power / 1e6;
    let n_20 = params.n_e / 1e20;
    let m = params.ion_mass / super::M_P; // Ion mass in proton masses
    let epsilon = config.minor_radius / config.major_radius;

    // ITER98y2 scaling
    0.0562 * ip_ma.powf(0.93)
        * config.toroidal_field.powf(0.15)
        * p_mw.powf(-0.69)
        * n_20.powf(0.41)
        * m.powf(0.19)
        * config.major_radius.powf(1.97)
        * epsilon.powf(0.58)
        * config.elongation.powf(0.78)
}

/// Triple product (Lawson criterion): n T τ_E
/// Need n T τ_E > 3×10²¹ keV·s/m³ for fusion
pub fn fusion_triple_product(
    params: &PlasmaParams,
    tau_e: f64,
) -> f64 {
    let n_20 = params.n_e / 1e20;
    let t_kev = params.t_e / 1000.0;

    n_20 * t_kev * tau_e * 1e20 // Return in m⁻³ keV s
}

/// Fusion power density (D-T reaction)
/// P_fusion ∝ n² <σv> Q
pub fn dt_fusion_power_density(params: &PlasmaParams) -> f64 {
    let t_kev = params.t_e / 1000.0;

    // D-T reactivity <σv> (approximate, valid for 10-100 keV)
    let sigma_v = if t_kev > 5.0 && t_kev < 100.0 {
        1e-22 * t_kev * t_kev / (1.0 + 0.01 * t_kev) // Simplified fit
    } else {
        0.0
    };

    let n_i = params.n_i; // Assume equal D and T densities
    let q_fusion = 17.6e6 * E; // 17.6 MeV per reaction

    0.25 * n_i * n_i * sigma_v * q_fusion // 0.25 for n_D = n_T = n/2
}

/// H-mode confinement enhancement factor
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConfinementMode {
    pub h_factor: f64, // H98(y,2) enhancement factor
    pub is_hmode: bool,
}

pub fn estimate_confinement_mode(
    config: &TokamakConfig,
    params: &PlasmaParams,
    heating_power: f64,
) -> ConfinementMode {
    let tau_scaling = iter98_confinement_time(config, params, heating_power);

    // Simple H-mode threshold: P > P_threshold
    let p_threshold_mw = 2.84; // Simplified threshold
    let p_mw = heating_power / 1e6;

    let is_hmode = p_mw > p_threshold_mw;
    let h_factor = if is_hmode { 1.0 } else { 0.85 }; // H-mode vs L-mode

    ConfinementMode { h_factor, is_hmode }
}

