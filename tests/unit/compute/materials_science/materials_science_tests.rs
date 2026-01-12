//! Comprehensive unit tests for materials_science module
//!
//! Tests cover:
//! - Crystal structures and crystallography
//! - Band theory and electronic structure
//! - Mechanical properties
//! - Thermal properties
//! - Diffusion
//! - XRD analysis

use super::*;

// =============================================================================
// Crystal Structures Tests
// =============================================================================

mod crystal_structures {
    use super::super::*;

    #[test]
    fn test_cubic_unit_cell() {
        let cell = UnitCell::cubic(4.0);
        assert_eq!(cell.a, 4.0);
        assert_eq!(cell.b, 4.0);
        assert_eq!(cell.c, 4.0);
        assert_eq!(cell.alpha, 90.0);
        assert_eq!(cell.beta, 90.0);
        assert_eq!(cell.gamma, 90.0);
    }

    #[test]
    fn test_tetragonal_unit_cell() {
        let cell = UnitCell::tetragonal(3.0, 5.0);
        assert_eq!(cell.a, 3.0);
        assert_eq!(cell.b, 3.0);
        assert_eq!(cell.c, 5.0);
        assert_eq!(cell.gamma, 90.0);
    }

    #[test]
    fn test_hexagonal_unit_cell() {
        let cell = UnitCell::hexagonal(2.5, 4.0);
        assert_eq!(cell.a, 2.5);
        assert_eq!(cell.b, 2.5);
        assert_eq!(cell.c, 4.0);
        assert_eq!(cell.gamma, 120.0);
    }

    #[test]
    fn test_unit_cell_volume_cubic() {
        let cell = UnitCell::cubic(4.0);
        let volume = unit_cell_volume(&cell);
        assert!((volume - 64.0).abs() < 0.001);
    }

    #[test]
    fn test_unit_cell_volume_hexagonal() {
        let cell = UnitCell::hexagonal(3.0, 5.0);
        let volume = unit_cell_volume(&cell);
        // V = a²c√3/2 for hexagonal
        let expected = 3.0 * 3.0 * 5.0 * (3.0_f64.sqrt() / 2.0);
        assert!((volume - expected).abs() < 0.01);
    }

    #[test]
    fn test_interplanar_spacing_cubic() {
        // (100) plane
        let d = interplanar_spacing_cubic(4.0, 1, 0, 0);
        assert!((d - 4.0).abs() < 0.001);

        // (110) plane
        let d = interplanar_spacing_cubic(4.0, 1, 1, 0);
        assert!((d - 2.828).abs() < 0.01);

        // (111) plane
        let d = interplanar_spacing_cubic(4.0, 1, 1, 1);
        assert!((d - 2.309).abs() < 0.01);
    }

    #[test]
    fn test_bragg_angle() {
        // Cu Kα wavelength = 1.5418 Å, d = 2.0 Å
        let theta = bragg_angle(2.0, 1.5418, 1);
        assert!(theta.is_some());
        let angle = theta.unwrap();
        assert!(angle > 20.0 && angle < 50.0);
    }

    #[test]
    fn test_bragg_angle_impossible() {
        // λ > 2d, Bragg condition cannot be satisfied
        let theta = bragg_angle(0.5, 1.5418, 1);
        assert!(theta.is_none());
    }

    #[test]
    fn test_atomic_packing_factor() {
        let apf_sc = atomic_packing_factor(BravaisLattice::SimpleCubic);
        assert!((apf_sc - 0.524).abs() < 0.001);

        let apf_bcc = atomic_packing_factor(BravaisLattice::BodyCenteredCubic);
        assert!((apf_bcc - 0.680).abs() < 0.001);

        let apf_fcc = atomic_packing_factor(BravaisLattice::FaceCenteredCubic);
        assert!((apf_fcc - 0.740).abs() < 0.001);
    }

    #[test]
    fn test_coordination_number() {
        assert_eq!(coordination_number(BravaisLattice::SimpleCubic), 6);
        assert_eq!(coordination_number(BravaisLattice::BodyCenteredCubic), 8);
        assert_eq!(coordination_number(BravaisLattice::FaceCenteredCubic), 12);
    }

    #[test]
    fn test_atoms_per_unit_cell() {
        assert!((atoms_per_unit_cell(BravaisLattice::SimpleCubic) - 1.0).abs() < 0.001);
        assert!((atoms_per_unit_cell(BravaisLattice::BodyCenteredCubic) - 2.0).abs() < 0.001);
        assert!((atoms_per_unit_cell(BravaisLattice::FaceCenteredCubic) - 4.0).abs() < 0.001);
    }

    #[test]
    fn test_systematic_absences_bcc() {
        // BCC: only h+k+l even allowed
        assert!(is_reflection_allowed_bcc(1, 1, 0)); // 2 - even
        assert!(is_reflection_allowed_bcc(2, 0, 0)); // 2 - even
        assert!(!is_reflection_allowed_bcc(1, 0, 0)); // 1 - odd
        assert!(!is_reflection_allowed_bcc(1, 1, 1)); // 3 - odd
    }

    #[test]
    fn test_systematic_absences_fcc() {
        // FCC: all even or all odd
        assert!(is_reflection_allowed_fcc(2, 0, 0)); // all even
        assert!(is_reflection_allowed_fcc(1, 1, 1)); // all odd
        assert!(!is_reflection_allowed_fcc(1, 0, 0)); // mixed
        assert!(!is_reflection_allowed_fcc(2, 1, 1)); // mixed
    }

    #[test]
    fn test_crystal_density() {
        // Iron BCC: a=2.87 Å, M=55.85 g/mol, 2 atoms/cell
        let vol = 2.87_f64.powi(3);
        let density = crystal_density(2.0, 55.85, vol);
        assert!((density - 7.87).abs() < 0.2); // Iron density ~7.87 g/cm³
    }

    #[test]
    fn test_angle_between_planes() {
        // (100) and (010) should be 90°
        let angle = angle_between_planes_cubic(1, 0, 0, 0, 1, 0);
        assert!((angle - 90.0).abs() < 0.001);

        // (100) and (110) should be 45°
        let angle = angle_between_planes_cubic(1, 0, 0, 1, 1, 0);
        assert!((angle - 45.0).abs() < 0.001);
    }
}

// =============================================================================
// Band Theory Tests
// =============================================================================

mod band_theory {
    use super::super::*;

    #[test]
    fn test_fermi_energy() {
        // Typical metal: n ~ 10^28 electrons/m³
        let n = 1e28;
        let ef = fermi_energy_free_electron(n);
        assert!(ef > 0.0 && ef < 20.0); // Typical range 2-10 eV
    }

    #[test]
    fn test_fermi_temperature() {
        let ef = 5.0; // 5 eV
        let tf = fermi_temperature(ef);
        assert!(tf > 10000.0); // Should be tens of thousands of K
    }

    #[test]
    fn test_fermi_velocity() {
        let ef = 5.0; // 5 eV
        let vf = fermi_velocity(ef);
        assert!(vf > 1e5 && vf < 1e7); // Typical ~10^6 m/s
    }

    #[test]
    fn test_material_classification() {
        assert_eq!(classify_material(0.0), MaterialClass::Conductor);
        assert_eq!(classify_material(1.1), MaterialClass::Semiconductor); // Si
        assert_eq!(classify_material(5.0), MaterialClass::Insulator);
    }

    #[test]
    fn test_band_gap_from_conductivity() {
        // Simulate semiconductor measurements
        let sigma1 = 1.0;
        let temp1 = 300.0;
        let sigma2 = 10.0;
        let temp2 = 400.0;

        let eg = band_gap_from_conductivity(sigma1, temp1, sigma2, temp2);
        assert!(eg > 0.0 && eg < 3.0); // Reasonable for semiconductor
    }

    #[test]
    fn test_intrinsic_carrier_concentration() {
        // Silicon at room temperature
        let ni = intrinsic_carrier_concentration(1.1, 300.0);
        assert!(ni > 1e15 && ni < 1e17); // Si: ~10^16 m^-3 at 300K
    }

    #[test]
    fn test_plasma_frequency() {
        let n = 1e28; // Typical metal
        let omega_p = plasma_frequency(n);
        assert!(omega_p > 1e15); // Typical plasma frequency
    }

    #[test]
    fn test_electron_mobility() {
        let sigma = 1e7; // S/m
        let n = 1e28;
        let mu = electron_mobility(sigma, n);
        assert!(mu > 0.0);
    }

    #[test]
    fn test_drift_velocity() {
        let mobility = 0.1; // m²/(V·s)
        let field = 1000.0; // V/m
        let vd = drift_velocity(mobility, field);
        assert!((vd - 100.0).abs() < 0.001);
    }

    #[test]
    fn test_thermionic_emission() {
        let temp = 1000.0; // K
        let work_function = 4.5; // eV
        let j = thermionic_emission_current(temp, work_function);
        assert!(j > 0.0);
    }
}

// =============================================================================
// Mechanical Properties Tests
// =============================================================================

mod mechanical_properties {
    use super::super::*;

    #[test]
    fn test_engineering_stress() {
        let stress = engineering_stress(1000.0, 0.001); // 1000 N, 1 mm²
        assert!((stress - 1e6).abs() < 0.1); // 1 MPa
    }

    #[test]
    fn test_engineering_strain() {
        let strain = engineering_strain(100.0, 102.0);
        assert!((strain - 0.02).abs() < 0.0001);
    }

    #[test]
    fn test_true_strain() {
        let strain = true_strain(100.0, 110.0);
        let expected = (1.1_f64).ln();
        assert!((strain - expected).abs() < 0.0001);
    }

    #[test]
    fn test_youngs_modulus() {
        let e = youngs_modulus(200e6, 0.001); // 200 MPa, 0.1% strain
        assert!((e - 200e9).abs() < 1e6); // 200 GPa
    }

    #[test]
    fn test_poissons_ratio() {
        let nu = poissons_ratio(-0.001, 0.003);
        assert!((nu - 0.333).abs() < 0.01); // Typical ~0.3
    }

    #[test]
    fn test_elastic_constant_relationships() {
        let g = 80e9; // Shear modulus 80 GPa
        let nu = 0.3;
        let e = youngs_from_shear_poisson(g, nu);

        // Check round-trip
        let g_calc = shear_from_youngs_poisson(e, nu);
        assert!((g - g_calc).abs() / g < 0.01);
    }

    #[test]
    fn test_bulk_modulus() {
        let k = bulk_modulus(1e6, -0.001); // 1 MPa, -0.1% volume change
        assert!((k - 1e9).abs() < 1e6);
    }

    #[test]
    fn test_von_mises_stress() {
        // Uniaxial tension
        let sigma_vm = von_mises_stress(100e6, 0.0, 0.0);
        assert!((sigma_vm - 100e6).abs() < 1.0);

        // Pure shear
        let sigma_vm = von_mises_stress(50e6, -50e6, 0.0);
        assert!(sigma_vm > 80e6 && sigma_vm < 90e6);
    }

    #[test]
    fn test_von_mises_yield_criterion() {
        let sigma_vm = 250e6; // Pa
        let sigma_y = 200e6; // Pa
        assert!(check_yield_von_mises(sigma_vm, sigma_y));
        assert!(!check_yield_von_mises(150e6, sigma_y));
    }

    #[test]
    fn test_tresca_stress() {
        let tau_max = tresca_stress(100e6, 20e6);
        assert!((tau_max - 40e6).abs() < 1.0);
    }

    #[test]
    fn test_stress_concentration() {
        let kt_circle = stress_concentration_circular();
        assert!((kt_circle - 3.0).abs() < 0.001);

        let kt_ellipse = stress_concentration_ellipse(10.0, 2.0);
        assert!((kt_ellipse - 11.0).abs() < 0.001); // 1 + 2*(10/2) = 1 + 10 = 11
    }

    #[test]
    fn test_stress_intensity_factor() {
        let ki = stress_intensity_factor(1.0, 100e6, 0.001); // 1 mm crack
        assert!(ki > 0.0);
    }

    #[test]
    fn test_critical_crack_length() {
        let kic = 50e6; // Pa·√m
        let sigma = 100e6; // Pa
        let ac = critical_crack_length(kic, sigma, 1.0);
        assert!(ac > 0.0);
    }

    #[test]
    fn test_hardness_conversions() {
        let hrc = 40.0; // Rockwell C
        let bhn = rockwell_c_to_brinell(hrc);
        // BHN ≈ 14.68 * 40 + 215.6 = 587.2 + 215.6 = 802.8
        assert!(bhn > 800.0 && bhn < 810.0);
    }

    #[test]
    fn test_ductility_measures() {
        let pct_el = percent_elongation(100.0, 115.0);
        assert!((pct_el - 15.0).abs() < 0.001);

        let pct_ra = percent_reduction_area(100.0, 70.0);
        assert!((pct_ra - 30.0).abs() < 0.001);
    }
}

// =============================================================================
// Thermal Properties Tests
// =============================================================================

mod thermal_properties {
    use super::super::*;

    #[test]
    fn test_dulong_petit() {
        let cv = dulong_petit_heat_capacity(6.022e23); // 1 mole
        let expected = 3.0 * 6.022e23 * 1.380649e-23; // 3Nk_B
        assert!((cv - expected).abs() / cv < 0.01);
    }

    #[test]
    fn test_debye_heat_capacity() {
        // Low temperature: C ∝ T³
        let cv1 = debye_heat_capacity(10.0, 300.0, 1e23);
        let cv2 = debye_heat_capacity(20.0, 300.0, 1e23);
        let ratio = cv2 / cv1;
        assert!((ratio - 8.0).abs() < 1.0); // (20/10)³ = 8
    }

    #[test]
    fn test_einstein_heat_capacity() {
        // High temperature limit should approach Dulong-Petit
        let num_atoms = 6.022e23;
        let cv_high = einstein_heat_capacity(1000.0, 100.0, num_atoms);
        let cv_dp = dulong_petit_heat_capacity(num_atoms);
        assert!((cv_high - cv_dp).abs() / cv_dp < 0.1);
    }

    #[test]
    fn test_thermal_expansion_coefficient() {
        let alpha = linear_thermal_expansion_coefficient(1.0, 1.00001, 300.0, 301.0);
        assert!(alpha > 0.0 && alpha < 1e-4); // Typical 10^-5 to 10^-6 /K
    }

    #[test]
    fn test_volumetric_expansion() {
        let alpha_l = 1e-5;
        let alpha_v = volumetric_thermal_expansion_coefficient(alpha_l);
        assert!((alpha_v - 3e-5).abs() < 1e-7);
    }

    #[test]
    fn test_thermal_expansion_length() {
        let delta_l = thermal_expansion_length(1.0, 1e-5, 100.0);
        assert!((delta_l - 0.001).abs() < 1e-6); // 1 mm expansion
    }

    #[test]
    fn test_thermal_stress() {
        let sigma = thermal_stress_constrained(200e9, 1e-5, 100.0);
        assert!(sigma > 1e8); // Significant stress
    }

    #[test]
    fn test_thermal_conductivity_wiedemann_franz() {
        let sigma = 6e7; // Copper electrical conductivity
        let kappa = thermal_conductivity_wiedemann_franz(sigma, 300.0);
        assert!(kappa > 300.0 && kappa < 500.0); // Copper ~400 W/(m·K)
    }

    #[test]
    fn test_thermal_diffusivity() {
        let alpha = thermal_diffusivity(400.0, 8960.0, 385.0); // Copper
        assert!(alpha > 1e-4 && alpha < 2e-4);
    }

    #[test]
    fn test_biot_number() {
        let bi = biot_number(10.0, 0.01, 50.0);
        assert!(bi > 0.0);
    }

    #[test]
    fn test_fourier_number() {
        let fo = fourier_number(1e-5, 100.0, 0.01);
        assert!(fo > 0.0);
    }

    #[test]
    fn test_stefan_boltzmann_radiation() {
        let q = stefan_boltzmann_radiation(0.9, 1.0, 1000.0, 300.0);
        assert!(q > 0.0); // Net radiation from hot to cold
    }

    #[test]
    fn test_heat_transfer_slab() {
        let q = heat_transfer_slab(50.0, 1.0, 400.0, 300.0, 0.1);
        assert!((q - 50000.0).abs() < 1.0); // 50 kW
    }

    #[test]
    fn test_thermal_resistance() {
        let r = thermal_resistance(0.1, 50.0, 1.0);
        assert!((r - 0.002).abs() < 1e-6); // 0.002 K/W
    }
}

// =============================================================================
// Diffusion Tests
// =============================================================================

mod diffusion {
    use super::super::*;

    #[test]
    fn test_fick_first_law() {
        let j = fick_first_law(1e-14, 1e20); // D=10^-14 m²/s, dC/dx=10^20 /m⁴
        assert!(j < 0.0); // Flux opposite to gradient
    }

    #[test]
    fn test_arrhenius_diffusion() {
        let d1 = arrhenius_diffusion_coefficient(1e-5, 100e3, 300.0);
        let d2 = arrhenius_diffusion_coefficient(1e-5, 100e3, 600.0);
        assert!(d2 > d1); // Higher temperature = faster diffusion
    }

    #[test]
    fn test_diffusion_distance() {
        let d = diffusion_distance(1e-14, 1000.0); // D=10^-14 m²/s, t=1000s
        assert!(d > 0.0 && d < 1e-5); // Nanometer scale
    }

    #[test]
    fn test_diffusion_time() {
        let t = diffusion_time(1e-14, 1e-6); // 1 μm diffusion distance
        assert!(t > 0.0);
    }

    #[test]
    fn test_concentration_semi_infinite() {
        let c0 = 1000.0;
        let c_surface = concentration_semi_infinite(c0, 0.0, 1e-14, 1000.0);
        assert!((c_surface - c0).abs() < 0.1); // At surface, C ≈ C₀
    }

    #[test]
    fn test_membrane_flux() {
        let j = membrane_flux_steady_state(1e-14, 1000.0, 100.0, 0.001);
        assert!(j > 0.0);
    }

    #[test]
    fn test_permeability() {
        let p = permeability_coefficient(1e-14, 0.5);
        assert!(p > 0.0);
    }

    #[test]
    fn test_interdiffusion_coefficient() {
        let d_tilde = interdiffusion_coefficient(0.6, 0.4, 1e-14, 2e-14);
        assert!(d_tilde > 1e-14 && d_tilde < 2e-14);
    }

    #[test]
    fn test_einstein_diffusion() {
        // Water at room temperature
        let d = einstein_diffusion_coefficient(300.0, 1e-3, 1e-9); // 1 nm particle
        // D = kT/(6πηr) ≈ 2.2e-10 m²/s for 1 nm particle in water
        assert!(d > 2e-10 && d < 2.3e-10);
    }

    #[test]
    fn test_ionic_conductivity() {
        let sigma = ionic_conductivity_from_diffusion(1e27, 1.6e-19, 1e-14, 300.0);
        assert!(sigma > 0.0);
    }
}

// =============================================================================
// XRD Analysis Tests
// =============================================================================

mod xrd_analysis {
    use super::super::*;

    #[test]
    fn test_two_theta_calculation() {
        let two_theta = calculate_two_theta(2.0, 1.5418);
        assert!(two_theta.is_some());
        let angle = two_theta.unwrap();
        assert!(angle > 0.0 && angle < 90.0);
    }

    #[test]
    fn test_d_spacing_from_two_theta() {
        let d = calculate_d_spacing(45.0, 1.5418);
        assert!(d > 0.0 && d < 5.0);
    }

    #[test]
    fn test_round_trip_d_and_two_theta() {
        let d_original = 2.5;
        let two_theta = calculate_two_theta(d_original, 1.5418).unwrap();
        let d_calculated = calculate_d_spacing(two_theta, 1.5418);
        assert!((d_original - d_calculated).abs() < 0.001);
    }

    #[test]
    fn test_scherrer_crystallite_size() {
        // Convert 0.5° FWHM to radians
        let fwhm_rad = 0.5_f64.to_radians();
        let theta_rad = 22.5_f64.to_radians();
        let size = scherrer_crystallite_size(1.5418, fwhm_rad, theta_rad, 0.9);
        assert!(size > 0.0 && size < 1000.0); // Reasonable crystallite size
    }

    #[test]
    fn test_instrumental_broadening_correction() {
        let beta_corrected = correct_instrumental_broadening(0.5, 0.3);
        assert!(beta_corrected > 0.0 && beta_corrected < 0.5);
    }

    #[test]
    fn test_microstrain_calculation() {
        let fwhm_rad = 0.5_f64.to_radians();
        let theta_rad = 22.5_f64.to_radians();
        let strain = calculate_microstrain(fwhm_rad, theta_rad);
        assert!(strain > 0.0 && strain < 0.01); // Typical microstrain
    }

    #[test]
    fn test_lattice_parameter_from_d() {
        let d = 2.0;
        let a = lattice_parameter_from_d_cubic(d, 1, 1, 0);
        let expected = 2.0 * 2.0_f64.sqrt();
        assert!((a - expected).abs() < 0.001);
    }

    #[test]
    fn test_lorentz_polarization_factor() {
        let lp = lorentz_polarization_factor(45.0_f64.to_radians());
        assert!(lp > 0.0);
    }

    #[test]
    fn test_temperature_factor() {
        let t = temperature_factor(1.0, 22.5_f64.to_radians(), 1.5418);
        assert!(t > 0.0 && t <= 1.0); // Reduces intensity
    }

    #[test]
    fn test_texture_coefficient() {
        let tc = texture_coefficient(100.0, 80.0, 1.2);
        assert!(tc > 0.0);
    }

    #[test]
    fn test_degree_of_crystallinity() {
        let xc = degree_of_crystallinity(70.0, 30.0);
        assert!((xc - 0.7).abs() < 0.001);
    }

    #[test]
    fn test_generate_cubic_pattern() {
        let peaks = generate_cubic_pattern(
            3.5,
            BravaisLattice::FaceCenteredCubic,
            1.5418,
            90.0,
        );
        assert!(!peaks.is_empty());

        // Check peaks are sorted by 2θ
        for i in 1..peaks.len() {
            assert!(peaks[i].two_theta >= peaks[i - 1].two_theta);
        }

        // Check intensities are normalized
        let max_intensity = peaks.iter().map(|p| p.intensity).fold(0.0, f64::max);
        assert!((max_intensity - 100.0).abs() < 0.001);
    }
}
