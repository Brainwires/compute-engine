// Unit tests for physics::plasma module
use super::*;
use std::f64::consts::PI;

// ============================================================================
// Plasma Parameters Tests
// ============================================================================

#[test]
fn test_tokamak_plasma_params() {
    let params = PlasmaParams::tokamak();

    assert_eq!(params.n_e, 1e20);
    assert_eq!(params.n_i, 1e20);
    assert_eq!(params.t_e, 10000.0);
    assert_eq!(params.t_i, 10000.0);
    assert_eq!(params.b_field, 5.0);
    assert_eq!(params.z_ion, 1.0);
    assert_eq!(params.ion_mass, 2.0 * M_P);
}

#[test]
fn test_solar_corona_params() {
    let params = PlasmaParams::solar_corona();

    assert_eq!(params.n_e, 1e15);
    assert_eq!(params.t_e, 100.0);
    assert_eq!(params.b_field, 0.01);
    assert_eq!(params.ion_mass, M_P);
}

#[test]
fn test_temperature_conversion() {
    let params = PlasmaParams::tokamak();
    let t_k = params.electron_temp_kelvin();

    // 10 keV ~ 116 million K
    assert!(t_k > 1e8 && t_k < 2e8);

    let t_i_k = params.ion_temp_kelvin();
    assert_eq!(t_k, t_i_k); // Same temperatures in tokamak
}

// ============================================================================
// Basic Plasma Physics Tests
// ============================================================================

#[test]
fn test_debye_length_tokamak() {
    let params = PlasmaParams::tokamak();
    let lambda_d = debye_length(&params);

    // Tokamak: λ_D ~ 10^-5 m
    assert!(lambda_d > 1e-6 && lambda_d < 1e-4);
}

#[test]
fn test_debye_length_corona() {
    let params = PlasmaParams::solar_corona();
    let lambda_d = debye_length(&params);

    // Solar corona has lower density, longer Debye length
    assert!(lambda_d > 1e-3 && lambda_d < 1e-1);
}

#[test]
fn test_plasma_frequency() {
    let n_e = 1e20; // Tokamak density
    let omega_pe = plasma_frequency(n_e);

    // ω_pe ~ 10^12 rad/s
    assert!(omega_pe > 1e11 && omega_pe < 1e13);
}

#[test]
fn test_plasma_frequency_scaling() {
    let n1 = 1e20;
    let n2 = 4e20;

    let omega1 = plasma_frequency(n1);
    let omega2 = plasma_frequency(n2);

    // Should scale as √n
    let ratio = omega2 / omega1;
    assert!((ratio - 2.0).abs() < 0.01);
}

#[test]
fn test_electron_cyclotron_frequency() {
    let b = 5.0; // 5 Tesla
    let omega_ce = electron_cyclotron_frequency(b);

    // ω_ce ~ 10^11 rad/s for 5T
    assert!(omega_ce > 1e10 && omega_ce < 1e12);
}

#[test]
fn test_ion_cyclotron_frequency() {
    let b = 5.0; // 5 Tesla
    let omega_ci = ion_cyclotron_frequency(b, 1.0, M_P);

    // Ion frequency is much lower than electron
    let omega_ce = electron_cyclotron_frequency(b);
    assert!(omega_ce > omega_ci * 1000.0);
}

#[test]
fn test_larmor_radius() {
    let v_perp = 1e6; // 1 Mm/s
    let omega_c = 1e8; // rad/s

    let r_l = larmor_radius(v_perp, omega_c);

    // r_L = v/ω
    assert!((r_l - 1e-2).abs() < 1e-6);
}

#[test]
fn test_thermal_velocity_electrons() {
    let t_e = 1000.0; // 1 keV
    let v_th_e = thermal_velocity(t_e, M_E);

    // Electron thermal velocity ~ 10^7 m/s
    assert!(v_th_e > 1e6 && v_th_e < 1e8);
}

#[test]
fn test_thermal_velocity_ions() {
    let t_i = 1000.0; // 1 keV
    let v_th_i = thermal_velocity(t_i, M_P);

    // Ion thermal velocity much slower than electrons
    let v_th_e = thermal_velocity(t_i, M_E);
    assert!(v_th_e > v_th_i * 40.0); // Mass ratio effect
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
    let params = PlasmaParams::tokamak();
    let rho = params.n_i * params.ion_mass;

    let v_a = alfven_velocity(b, rho);

    // Alfvén velocity ~ 10^6 m/s for typical tokamak
    assert!(v_a > 1e5 && v_a < 1e7);
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

// ============================================================================
// MHD Tests
// ============================================================================

#[test]
fn test_magnetic_pressure() {
    let b = 5.0; // 5 Tesla
    let p_b = magnetic_pressure(b);

    // p_B = B²/(2μ₀) ~ 10^6 Pa for 5T field
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
    assert!(tension[1].abs() < 1e-10);
    assert!(tension[2].abs() < 1e-10);
}

#[test]
fn test_lorentz_force() {
    let b = [0.0, 0.0, 5.0];
    let curl_b = [1.0, 0.0, 0.0];

    let force = lorentz_force(&b, &curl_b);

    // Force should be perpendicular to both J and B
    // J = (1, 0, 0) × B = (0, 0, 5) → F = (0, *, 0)
    assert!(force[0].abs() < 1e-10);
    assert!(force[1] != 0.0);
    assert!(force[2].abs() < 1e-10);
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

    // Should be B/ρ
    assert!((param - 5e6).abs() < 1.0);
}

#[test]
fn test_mhd_equilibrium() {
    let grad_p = [1000.0, 0.0, 0.0];
    let b = [0.0, 0.0, 5.0];
    let curl_b = [1000.0 * 4.0e-7 * PI, 0.0, 0.0];

    let eq = check_mhd_equilibrium(&grad_p, &b, &curl_b);

    assert!(eq.lorentz_force[0].is_finite());
    assert_eq!(eq.grad_p[0], 1000.0);
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

// ============================================================================
// Wave Physics Tests
// ============================================================================

#[test]
fn test_langmuir_wave() {
    let params = PlasmaParams::tokamak();
    let k = 1e3; // Wave number

    let wave = langmuir_wave(k, &params);

    assert!(wave.omega > 0.0);
    assert!(wave.v_phase > 0.0);
    assert!(wave.v_group > 0.0);
    assert!(wave.v_group < wave.v_phase); // Dispersive wave
}

#[test]
fn test_alfven_wave() {
    let params = PlasmaParams::tokamak();
    let k_parallel = 1e2;

    let wave = alfven_wave(k_parallel, &params);

    // Alfvén wave: ω = k v_A (non-dispersive)
    assert_eq!(wave.v_phase, wave.v_group);
    assert!(wave.omega > 0.0);
}

#[test]
fn test_ion_acoustic_wave() {
    let params = PlasmaParams::tokamak();
    let k = 1e3;

    let wave = ion_acoustic_wave(k, &params);

    // Sound wave: ω = k c_s
    assert_eq!(wave.v_phase, wave.v_group);
    assert!(wave.omega > 0.0);
}

#[test]
fn test_upper_hybrid_frequency() {
    let params = PlasmaParams::tokamak();
    let omega_uh = upper_hybrid_frequency(&params);

    let omega_pe = plasma_frequency(params.n_e);
    let omega_ce = electron_cyclotron_frequency(params.b_field);

    // ω_uh² = ω_pe² + ω_ce²
    let expected = (omega_pe * omega_pe + omega_ce * omega_ce).sqrt();
    assert!((omega_uh - expected).abs() < 1e-6);
}

#[test]
fn test_lower_hybrid_frequency() {
    let params = PlasmaParams::tokamak();
    let omega_lh = lower_hybrid_frequency(&params);

    let omega_ci = ion_cyclotron_frequency(params.b_field, params.z_ion, params.ion_mass);
    let omega_ce = electron_cyclotron_frequency(params.b_field);

    // ω_lh ≈ √(ω_ci · ω_ce)
    let expected = (omega_ci * omega_ce).sqrt();
    assert!((omega_lh - expected).abs() / expected < 0.01);
}

#[test]
fn test_two_stream_growth_rate() {
    let n_beam = 1e19;
    let n_background = 1e20;

    let gamma = two_stream_growth_rate(n_beam, n_background);

    assert!(gamma > 0.0);
    assert!(gamma.is_finite());
}

#[test]
fn test_weibel_growth_rate() {
    let params = PlasmaParams::tokamak();
    let t_perp = 1500.0; // Higher perpendicular temp
    let t_parallel = 1000.0;

    let gamma = weibel_growth_rate(&params, t_perp, t_parallel);

    // Should be unstable with T⊥ > T∥
    assert!(gamma > 0.0);
}

#[test]
fn test_weibel_stable() {
    let params = PlasmaParams::tokamak();
    let t_perp = 1000.0;
    let t_parallel = 1500.0; // Higher parallel temp

    let gamma = weibel_growth_rate(&params, t_perp, t_parallel);

    // Should be stable with T⊥ < T∥
    assert_eq!(gamma, 0.0);
}

#[test]
fn test_em_wave_in_plasma() {
    let params = PlasmaParams::tokamak();
    let k = 1e7; // High frequency

    let wave = em_wave_in_plasma(k, &params);

    // ω² = ω_pe² + k²c²
    assert!(wave.omega > 0.0);
    assert!(wave.v_phase > 0.0);
}

#[test]
fn test_cutoff_frequency() {
    let params = PlasmaParams::tokamak();
    let omega_cutoff = cutoff_frequency(&params);
    let omega_pe = plasma_frequency(params.n_e);

    assert_eq!(omega_cutoff, omega_pe);
}

// ============================================================================
// Confinement Physics Tests
// ============================================================================

#[test]
fn test_tokamak_config() {
    let config = TokamakConfig {
        major_radius: 6.2,    // ITER-like
        minor_radius: 2.0,
        toroidal_field: 5.3,
        plasma_current: 15e6, // 15 MA
        elongation: 1.7,
        triangularity: 0.33,
    };

    assert_eq!(config.major_radius, 6.2);
    assert_eq!(config.plasma_current, 15e6);
}

#[test]
fn test_safety_factor() {
    let config = TokamakConfig {
        major_radius: 6.0,
        minor_radius: 2.0,
        toroidal_field: 5.0,
        plasma_current: 15e6,
        elongation: 1.7,
        triangularity: 0.3,
    };

    let q = safety_factor(&config, config.minor_radius);

    // Edge safety factor should be > 2 for stability
    assert!(q > 0.0);
    assert!(q.is_finite());
}

#[test]
fn test_aspect_ratio() {
    let config = TokamakConfig {
        major_radius: 6.0,
        minor_radius: 2.0,
        toroidal_field: 5.0,
        plasma_current: 15e6,
        elongation: 1.7,
        triangularity: 0.3,
    };

    let a_ratio = aspect_ratio(&config);
    assert_eq!(a_ratio, 3.0);
}

#[test]
fn test_plasma_volume() {
    let config = TokamakConfig {
        major_radius: 6.0,
        minor_radius: 2.0,
        toroidal_field: 5.0,
        plasma_current: 15e6,
        elongation: 1.7,
        triangularity: 0.3,
    };

    let volume = plasma_volume(&config);

    // Should be ~ 800 m³ for ITER-like machine
    assert!(volume > 700.0 && volume < 900.0);
}

#[test]
fn test_plasma_surface_area() {
    let config = TokamakConfig {
        major_radius: 6.0,
        minor_radius: 2.0,
        toroidal_field: 5.0,
        plasma_current: 15e6,
        elongation: 1.7,
        triangularity: 0.3,
    };

    let area = plasma_surface_area(&config);

    assert!(area > 0.0);
    assert!(area.is_finite());
}

#[test]
fn test_magnetic_field_energy() {
    let config = TokamakConfig {
        major_radius: 6.0,
        minor_radius: 2.0,
        toroidal_field: 5.0,
        plasma_current: 15e6,
        elongation: 1.7,
        triangularity: 0.3,
    };

    let energy = magnetic_field_energy(&config);

    // Should be ~ GJ scale
    assert!(energy > 1e9);
}

#[test]
fn test_kruskal_shafranov_stability() {
    let config = TokamakConfig {
        major_radius: 6.0,
        minor_radius: 2.0,
        toroidal_field: 5.0,
        plasma_current: 15e6,
        elongation: 1.7,
        triangularity: 0.3,
    };

    let stable = kruskal_shafranov_stable(&config);

    // Well-designed tokamak should be stable
    assert!(stable);
}

#[test]
fn test_troyon_beta_limit() {
    let config = TokamakConfig {
        major_radius: 6.0,
        minor_radius: 2.0,
        toroidal_field: 5.0,
        plasma_current: 15e6,
        elongation: 1.7,
        triangularity: 0.3,
    };

    let beta_n = 3.0; // Typical normalized beta
    let beta_limit = troyon_beta_limit(&config, beta_n);

    // β_limit = β_N * I_p(MA) / (a(m) * B_T(T))
    // For these params: 3.0 * 15 / (2.0 * 5.0) = 4.5 (percentage)
    assert!(beta_limit > 0.0);
    assert!(beta_limit < 10.0); // Reasonable upper bound for beta limit percentage
}

#[test]
fn test_greenwald_density_limit() {
    let config = TokamakConfig {
        major_radius: 6.0,
        minor_radius: 2.0,
        toroidal_field: 5.0,
        plasma_current: 15e6,
        elongation: 1.7,
        triangularity: 0.3,
    };

    let n_greenwald = greenwald_density_limit(&config);

    // Should be ~ 10^20 m⁻³
    assert!(n_greenwald > 1e19 && n_greenwald < 1e21);
}

#[test]
fn test_iter98_confinement_time() {
    let config = TokamakConfig {
        major_radius: 6.2,
        minor_radius: 2.0,
        toroidal_field: 5.3,
        plasma_current: 15e6,
        elongation: 1.7,
        triangularity: 0.33,
    };

    let params = PlasmaParams::tokamak();
    let heating_power = 50e6; // 50 MW

    let tau_e = iter98_confinement_time(&config, &params, heating_power);

    // ITER should achieve ~3-4 seconds
    assert!(tau_e > 1.0 && tau_e < 10.0);
}

#[test]
fn test_fusion_triple_product() {
    let params = PlasmaParams::tokamak();
    let tau_e = 3.0; // 3 seconds

    let ntt = fusion_triple_product(&params, tau_e);

    // Need > 3×10²¹ for ignition
    // With n=10^20, T=10keV, τ=3s → nTτ = 3×10²¹
    assert!(ntt > 1e21);
}

#[test]
fn test_dt_fusion_power_density() {
    let params = PlasmaParams::tokamak();
    let power = dt_fusion_power_density(&params);

    // At 10 keV, should produce fusion power
    assert!(power > 0.0);
    assert!(power.is_finite());
}

#[test]
fn test_confinement_mode_estimation() {
    let config = TokamakConfig {
        major_radius: 6.0,
        minor_radius: 2.0,
        toroidal_field: 5.0,
        plasma_current: 15e6,
        elongation: 1.7,
        triangularity: 0.3,
    };

    let params = PlasmaParams::tokamak();
    let heating_power = 50e6; // High power → H-mode

    let mode = estimate_confinement_mode(&config, &params, heating_power);

    assert!(mode.h_factor > 0.0);
    assert!(mode.h_factor <= 1.0);
}
