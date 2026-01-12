//! Comprehensive tests for Electromagnetism module
//!
//! Tests Maxwell's equations, EM waves, antennas, transmission lines,
//! waveguides, scattering, and related electromagnetic phenomena.

use computational_engine::compute::physics::electromagnetism::*;

// ============================================================================
// MAXWELL'S EQUATIONS TESTS
// ============================================================================

#[test]
fn test_gauss_law_electric() {
    let result = maxwell_equations(MaxwellRequest {
        equation: "gauss_electric".to_string(),
        electric_field: Some(vec![1.0, 2.0, 3.0]),
        magnetic_field: None,
        charge_density: Some(1.0e-9),
        current_density: None,
        position: vec![0.0, 0.0, 0.0],
    })
    .unwrap();

    assert_eq!(result.equation, "Gauss's Law (Electric)");
    assert!(result.divergence.is_some());
    assert!(result.field_strength > 0.0);
}

#[test]
fn test_gauss_law_magnetic() {
    let result = maxwell_equations(MaxwellRequest {
        equation: "gauss_magnetic".to_string(),
        electric_field: None,
        magnetic_field: Some(vec![0.5, 0.5, 0.5]),
        charge_density: None,
        current_density: None,
        position: vec![0.0, 0.0, 0.0],
    })
    .unwrap();

    assert_eq!(result.equation, "Gauss's Law (Magnetic)");
    assert!(result.divergence.is_some());
}

#[test]
fn test_faraday_law() {
    let result = maxwell_equations(MaxwellRequest {
        equation: "faraday".to_string(),
        electric_field: Some(vec![1.0, 0.0, 0.0]),
        magnetic_field: None,
        charge_density: None,
        current_density: None,
        position: vec![0.0, 0.0, 0.0],
    })
    .unwrap();

    assert_eq!(result.equation, "Faraday's Law");
    assert!(result.curl.is_some());
    assert_eq!(result.curl.as_ref().unwrap().len(), 3);
}

#[test]
fn test_ampere_maxwell_law() {
    let result = maxwell_equations(MaxwellRequest {
        equation: "ampere".to_string(),
        electric_field: None,
        magnetic_field: Some(vec![0.0, 1.0, 0.0]),
        charge_density: None,
        current_density: Some(vec![0.0, 0.0, 1.0]),
        position: vec![0.0, 0.0, 0.0],
    })
    .unwrap();

    assert_eq!(result.equation, "Ampere-Maxwell Law");
    assert!(result.curl.is_some());
}

// ============================================================================
// ELECTROMAGNETIC WAVE TESTS
// ============================================================================

#[test]
fn test_em_wave_vacuum() {
    let result = em_wave(WaveRequest {
        frequency: 1e9, // 1 GHz
        wavelength: None,
        medium: "vacuum".to_string(),
        permittivity: None,
        permeability: None,
    })
    .unwrap();

    // λ = c/f = 3e8/1e9 = 0.3 m
    assert!((result.wavelength - 0.3).abs() < 0.01);
    assert!((result.frequency - 1e9).abs() < 1.0);
    assert!(result.phase_velocity > 2.9e8);
}

#[test]
fn test_em_wave_microwave() {
    let result = em_wave(WaveRequest {
        frequency: 2.45e9, // Microwave oven frequency
        wavelength: None,
        medium: "air".to_string(),
        permittivity: None,
        permeability: None,
    })
    .unwrap();

    // λ ≈ 12.2 cm
    assert!((result.wavelength - 0.122).abs() < 0.01);
    assert!(result.wave_number > 0.0);
    assert!(result.angular_frequency > 0.0);
}

#[test]
fn test_em_wave_radio() {
    let result = em_wave(WaveRequest {
        frequency: 100e6, // FM radio (100 MHz)
        wavelength: None,
        medium: "vacuum".to_string(),
        permittivity: None,
        permeability: None,
    })
    .unwrap();

    // λ = 3 m
    assert!((result.wavelength - 3.0).abs() < 0.1);
}

#[test]
fn test_em_wave_dielectric() {
    let eps = 4.0 * 8.854187817e-12; // εr = 4
    let result = em_wave(WaveRequest {
        frequency: 1e9,
        wavelength: None,
        medium: "dielectric".to_string(),
        permittivity: Some(eps),
        permeability: None,
    })
    .unwrap();

    // In dielectric, v = c/√εr, so λ = λ0/√εr
    assert!(result.wavelength < 0.3);
    assert!(result.phase_velocity < 3e8);
}

// ============================================================================
// ANTENNA TESTS
// ============================================================================

#[test]
fn test_half_wave_dipole() {
    let result = antenna_analysis(AntennaRequest {
        antenna_type: "dipole".to_string(),
        frequency: 100e6, // 100 MHz
        length: None,
        distance: 1000.0, // 1 km
        theta: std::f64::consts::PI / 2.0,
        power: 100.0, // 100 W
    })
    .unwrap();

    assert!(result.gain > 0.0);
    assert!((result.directivity - 1.64).abs() < 0.1);
    assert!((result.radiation_resistance - 73.0).abs() < 1.0);
    assert!(result.efficiency > 0.9);
}

#[test]
fn test_monopole_antenna() {
    let result = antenna_analysis(AntennaRequest {
        antenna_type: "monopole".to_string(),
        frequency: 900e6, // 900 MHz (cellular)
        length: None,
        distance: 100.0,
        theta: 0.0,
        power: 1.0, // 1 W
    })
    .unwrap();

    assert!(result.gain > 0.0);
    assert!(result.directivity > 1.64); // Higher than dipole
    assert!(result.beam_width > 0.0);
}

#[test]
fn test_patch_antenna() {
    let result = antenna_analysis(AntennaRequest {
        antenna_type: "patch".to_string(),
        frequency: 2.4e9, // WiFi
        length: None,
        distance: 10.0,
        theta: 0.0,
        power: 0.1, // 100 mW
    })
    .unwrap();

    assert!(result.gain > 3.0);
    assert!(result.efficiency > 0.8);
    assert!(result.electric_field_strength > 0.0);
}

#[test]
fn test_antenna_far_field() {
    let result = antenna_analysis(AntennaRequest {
        antenna_type: "dipole".to_string(),
        frequency: 1e9,
        length: None,
        distance: 10000.0, // 10 km
        theta: 0.0,
        power: 1000.0, // 1 kW
    })
    .unwrap();

    // Field should decrease with distance
    assert!(result.electric_field_strength > 0.0);
    assert!(result.electric_field_strength < 1.0); // Not too strong
}

// ============================================================================
// TRANSMISSION LINE TESTS
// ============================================================================

#[test]
fn test_transmission_line_matched() {
    let result = transmission_line(TransmissionLineRequest {
        line_type: "coax".to_string(),
        frequency: 1e9,
        length: 10.0,
        z0: 50.0,
        load_impedance: Some(50.0), // Matched
        attenuation: Some(0.1),
    })
    .unwrap();

    // Perfect match: VSWR = 1, reflection = 0
    assert!((result.vswr - 1.0).abs() < 0.01);
    assert!(result.reflection_coefficient < 0.01);
    assert!(result.return_loss > 30.0); // High return loss = good
}

#[test]
fn test_transmission_line_mismatched() {
    let result = transmission_line(TransmissionLineRequest {
        line_type: "coax".to_string(),
        frequency: 1e9,
        length: 1.0,
        z0: 50.0,
        load_impedance: Some(75.0), // Mismatched
        attenuation: None,
    })
    .unwrap();

    // VSWR > 1 for mismatch
    assert!(result.vswr > 1.0);
    assert!(result.reflection_coefficient > 0.0);
    assert!(result.return_loss < 20.0);
}

#[test]
fn test_transmission_line_open() {
    let result = transmission_line(TransmissionLineRequest {
        line_type: "coax".to_string(),
        frequency: 100e6,
        length: 1.0,
        z0: 50.0,
        load_impedance: Some(1e6), // Very high (open)
        attenuation: None,
    })
    .unwrap();

    // High VSWR for open circuit
    assert!(result.vswr > 10.0);
}

#[test]
fn test_transmission_line_attenuation() {
    let result = transmission_line(TransmissionLineRequest {
        line_type: "coax".to_string(),
        frequency: 10e9, // 10 GHz
        length: 100.0,   // 100 m
        z0: 50.0,
        load_impedance: Some(50.0),
        attenuation: Some(0.5), // 0.5 dB/m
    })
    .unwrap();

    // 50 dB insertion loss
    assert!((result.insertion_loss - 50.0).abs() < 1.0);
}

// ============================================================================
// WAVEGUIDE TESTS
// ============================================================================

#[test]
fn test_rectangular_waveguide_te10() {
    let result = waveguide(WaveguideRequest {
        guide_type: "rectangular".to_string(),
        frequency: 10e9,    // 10 GHz (X-band)
        width: 0.023,       // 23 mm (standard WR-90)
        height: Some(0.01), // 10 mm
        mode: "TE10".to_string(),
    })
    .unwrap();

    // Cutoff frequency for TE10: fc = c/(2a)
    let expected_fc = 3e8 / (2.0 * 0.023);
    assert!((result.cutoff_frequency - expected_fc).abs() < 1e8);
    assert!(result.guide_wavelength > 0.0);
    assert!(result.phase_velocity > 3e8); // vp > c in waveguide
    assert!(result.group_velocity < 3e8); // vg < c
}

#[test]
fn test_rectangular_waveguide_below_cutoff() {
    let result = waveguide(WaveguideRequest {
        guide_type: "rectangular".to_string(),
        frequency: 1e9, // Too low for X-band
        width: 0.023,
        height: Some(0.01),
        mode: "TE10".to_string(),
    });

    // Should fail below cutoff
    assert!(result.is_err());
}

#[test]
fn test_circular_waveguide() {
    let result = waveguide(WaveguideRequest {
        guide_type: "circular".to_string(),
        frequency: 10e9,
        width: 0.01, // 10 mm radius
        height: None,
        mode: "TE11".to_string(),
    })
    .unwrap();

    assert!(result.cutoff_frequency > 0.0);
    assert!(result.guide_wavelength > 0.0);
    assert!(result.attenuation < 0.01); // Lower than rectangular
}

#[test]
fn test_coaxial_waveguide() {
    let result = waveguide(WaveguideRequest {
        guide_type: "coaxial".to_string(),
        frequency: 1e9,
        width: 0.001,        // 1 mm inner radius
        height: Some(0.003), // 3 mm outer radius
        mode: "TEM".to_string(),
    })
    .unwrap();

    // TEM mode has no cutoff
    assert_eq!(result.cutoff_frequency, 0.0);
    assert!((result.phase_velocity - 3e8).abs() < 1e6);
}

#[test]
fn test_dielectric_waveguide() {
    let result = waveguide(WaveguideRequest {
        guide_type: "dielectric".to_string(),
        frequency: 200e12, // 200 THz (optical)
        width: 10e-6,      // 10 μm slab
        height: Some(1.5), // Refractive index
        mode: "TE0".to_string(),
    })
    .unwrap();

    assert!(result.cutoff_frequency > 0.0);
    assert!(result.attenuation > 0.0); // Has some loss
    assert!(result.attenuation < 100.0); // Reasonable loss value
}

// ============================================================================
// SCATTERING TESTS
// ============================================================================

#[test]
fn test_rayleigh_scattering() {
    let result = scattering(ScatteringRequest {
        scatterer_type: "sphere".to_string(),
        radius: 0.001,   // 1 mm
        wavelength: 0.1, // 10 cm (ka << 1)
        incident_angle: 0.0,
        polarization: "TE".to_string(),
    })
    .unwrap();

    // Rayleigh regime: small ka
    let ka = 2.0 * std::f64::consts::PI * 0.001 / 0.1;
    assert!(ka < 1.0);
    assert!(result.radar_cross_section > 0.0);
    assert_eq!(result.scattering_pattern.len(), 36);
}

#[test]
fn test_geometric_scattering() {
    let result = scattering(ScatteringRequest {
        scatterer_type: "sphere".to_string(),
        radius: 1.0,     // 1 m
        wavelength: 0.1, // 10 cm (ka >> 1)
        incident_angle: 0.0,
        polarization: "TM".to_string(),
    })
    .unwrap();

    // Geometric optics: RCS ≈ πa²
    let expected_rcs = std::f64::consts::PI;
    assert!((result.radar_cross_section - expected_rcs).abs() < 0.5);
}

#[test]
fn test_scattering_pattern() {
    let result = scattering(ScatteringRequest {
        scatterer_type: "sphere".to_string(),
        radius: 0.1,
        wavelength: 0.05,
        incident_angle: 0.0,
        polarization: "TE".to_string(),
    })
    .unwrap();

    // Check pattern has correct structure
    assert_eq!(result.scattering_pattern.len(), 36);
    for (angle, intensity) in &result.scattering_pattern {
        assert!(*angle >= 0.0);
        assert!(*angle <= 2.0 * std::f64::consts::PI); // Full circle in radians
        assert!(*intensity >= 0.0);
        assert!(*intensity <= 1.0);
    }
}

// ============================================================================
// POYNTING VECTOR TESTS
// ============================================================================

#[test]
fn test_poynting_vector_plane_wave() {
    let result = poynting_vector(PoyntingRequest {
        electric_field: vec![1.0, 0.0, 0.0],
        magnetic_field: vec![0.0, 1.0 / 377.0, 0.0], // H = E/Z0
    })
    .unwrap();

    assert_eq!(result.poynting_vector.len(), 3);
    assert!(result.power_density > 0.0);
}

#[test]
fn test_poynting_vector_perpendicular() {
    let result = poynting_vector(PoyntingRequest {
        electric_field: vec![10.0, 0.0, 0.0],
        magnetic_field: vec![0.0, 0.1, 0.0],
    })
    .unwrap();

    // S = E × H, should point in z direction
    assert!((result.poynting_vector[0]).abs() < 1e-10);
    assert!((result.poynting_vector[1]).abs() < 1e-10);
    assert!(result.poynting_vector[2].abs() > 0.0);
}

#[test]
fn test_poynting_vector_parallel() {
    let result = poynting_vector(PoyntingRequest {
        electric_field: vec![1.0, 0.0, 0.0],
        magnetic_field: vec![1.0, 0.0, 0.0],
    })
    .unwrap();

    // Parallel E and H: S = 0
    assert!(result.power_density < 1e-10);
}

#[test]
fn test_poynting_power_density() {
    let e_mag = 100.0; // V/m
    let h_mag = e_mag / 377.0;

    let result = poynting_vector(PoyntingRequest {
        electric_field: vec![e_mag, 0.0, 0.0],
        magnetic_field: vec![0.0, h_mag, 0.0],
    })
    .unwrap();

    // Power density = E²/Z0
    let expected = e_mag * e_mag / 377.0;
    assert!((result.power_density - expected).abs() < 1.0);
}

// ============================================================================
// SKIN EFFECT TESTS
// ============================================================================

#[test]
fn test_skin_effect_copper_1mhz() {
    let result = skin_effect(SkinEffectRequest {
        frequency: 1e6,      // 1 MHz
        conductivity: 5.8e7, // Copper
        permeability: None,
    })
    .unwrap();

    // At 1 MHz, copper skin depth ≈ 66 μm
    assert!((result.skin_depth - 66e-6).abs() < 10e-6);
    assert!(result.surface_resistance > 0.0);
}

#[test]
fn test_skin_effect_copper_1ghz() {
    let result = skin_effect(SkinEffectRequest {
        frequency: 1e9, // 1 GHz
        conductivity: 5.8e7,
        permeability: None,
    })
    .unwrap();

    // At 1 GHz, skin depth is much smaller
    assert!(result.skin_depth < 10e-6);
    assert!(result.skin_depth > 1e-6);
}

#[test]
fn test_skin_effect_aluminum() {
    let result = skin_effect(SkinEffectRequest {
        frequency: 1e6,
        conductivity: 3.5e7, // Aluminum
        permeability: None,
    })
    .unwrap();

    // Aluminum has larger skin depth than copper
    assert!(result.skin_depth > 66e-6);
}

#[test]
fn test_skin_effect_frequency_scaling() {
    let result1 = skin_effect(SkinEffectRequest {
        frequency: 1e6,
        conductivity: 5.8e7,
        permeability: None,
    })
    .unwrap();

    let result2 = skin_effect(SkinEffectRequest {
        frequency: 4e6, // 4x frequency
        conductivity: 5.8e7,
        permeability: None,
    })
    .unwrap();

    // Skin depth scales as 1/√f
    let ratio = result1.skin_depth / result2.skin_depth;
    assert!((ratio - 2.0).abs() < 0.1);
}

// ============================================================================
// PLASMA FREQUENCY TESTS
// ============================================================================

#[test]
fn test_plasma_frequency_ionosphere() {
    let n_e = 1e11; // electrons/m³ (ionosphere)
    let fp = plasma_frequency(n_e);

    // Typical ionosphere: 1-10 MHz
    assert!(fp > 1e6);
    assert!(fp < 10e6);
}

#[test]
fn test_plasma_frequency_metal() {
    let n_e = 8.5e28; // electrons/m³ (copper)
    let fp = plasma_frequency(n_e);

    // Metal plasma frequency in UV range
    assert!(fp > 1e15);
}

#[test]
fn test_plasma_frequency_scaling() {
    let n1 = 1e10;
    let n2 = 4e10; // 4x density

    let fp1 = plasma_frequency(n1);
    let fp2 = plasma_frequency(n2);

    // fp scales as √n
    let ratio = fp2 / fp1;
    assert!((ratio - 2.0).abs() < 0.1);
}
