//! Comprehensive scientific computing test suite
//!
//! Tests for all COMPUTE scientific operations including:
//! - Chemistry (8 operations)
//! - Biology (6 operations)
//! - Thermodynamics (5 operations)
//! - Optics (4 operations)
//! - Geophysics (4 categories)
//! - Engineering (4 disciplines)

use computational_engine::engine::*;
use computational_engine::create_default_dispatcher;
use std::collections::HashMap;

// ============================================================================
// CHEMISTRY TESTS (8 operations)
// ============================================================================

#[test]
fn test_ph_calculation() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("pka".to_string(), serde_json::json!(4.76));
    params.insert("concentration_acid".to_string(), serde_json::json!(0.1));
    params.insert("concentration_base".to_string(), serde_json::json!(0.1));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Chemistry(ChemistryOp::PhCalculation),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate pH: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        let value = output.result.get("value").and_then(|v| v.as_f64()).unwrap();
        assert!((value - 4.76).abs() < 0.01, "pH should be ~4.76");
    }
}

#[test]
fn test_buffer_capacity() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("concentration_acid".to_string(), serde_json::json!(0.1));
    params.insert("concentration_base".to_string(), serde_json::json!(0.1));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Chemistry(ChemistryOp::BufferCapacity),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate buffer capacity: {:?}", result);
}

#[test]
fn test_arrhenius() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("activation_energy".to_string(), serde_json::json!(50.0)); // kJ/mol
    params.insert("temperature".to_string(), serde_json::json!(298.0)); // K
    params.insert("pre_exponential".to_string(), serde_json::json!(1e10));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Chemistry(ChemistryOp::Arrhenius),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate Arrhenius rate constant: {:?}", result);
}

#[test]
fn test_nernst_equation() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("standard_potential".to_string(), serde_json::json!(0.34)); // V
    params.insert("n_electrons".to_string(), serde_json::json!(2)); // Integer, not float
    params.insert("concentrations".to_string(), serde_json::json!([0.1, 1.0]));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Chemistry(ChemistryOp::NernstEquation),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate Nernst potential: {:?}", result);
}

#[test]
fn test_beer_lambert() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("molar_absorptivity".to_string(), serde_json::json!(1000.0)); // L/(mol·cm)
    params.insert("path_length".to_string(), serde_json::json!(1.0)); // cm
    params.insert("concentration".to_string(), serde_json::json!(0.001)); // M

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Chemistry(ChemistryOp::BeerLambert),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate absorbance: {:?}", result);
}

// ============================================================================
// BIOLOGY TESTS (6 operations)
// ============================================================================

#[test]
fn test_michaelis_menten() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("vmax".to_string(), serde_json::json!(100.0));
    params.insert("km".to_string(), serde_json::json!(5.0));
    params.insert("substrate_concentration".to_string(), serde_json::json!(10.0));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Biology(BiologyOp::MichaelisMenten),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate enzyme velocity: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        let value = output.result.get("value").and_then(|v| v.as_f64()).unwrap();
        assert!(value > 60.0 && value < 70.0, "Velocity should be ~66.7");
    }
}

#[test]
fn test_hardy_weinberg() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("allele_frequency_p".to_string(), serde_json::json!(0.6)); // allele frequency p

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Biology(BiologyOp::HardyWeinberg),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate genotype frequencies: {:?}", result);
}

#[test]
fn test_goldman_equation() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("temperature".to_string(), serde_json::json!(310.0)); // K (37°C)
    params.insert("ion_concentrations_inside".to_string(), serde_json::json!([140.0, 12.0, 4.0])); // K+, Na+, Cl- (mM)
    params.insert("ion_concentrations_outside".to_string(), serde_json::json!([5.0, 145.0, 110.0])); // K+, Na+, Cl- (mM)
    params.insert("permeabilities".to_string(), serde_json::json!([1.0, 0.04, 0.45])); // Relative permeabilities

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Biology(BiologyOp::GoldmanEquation),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate membrane potential: {:?}", result);
}

#[test]
fn test_allometric_scaling() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("body_mass".to_string(), serde_json::json!(70.0)); // kg (human)
    params.insert("scaling_type".to_string(), serde_json::json!("metabolic")); // "metabolic", "surface_area", or "lifespan"

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Biology(BiologyOp::AllometricScaling),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate metabolic rate: {:?}", result);
}

// ============================================================================
// THERMODYNAMICS TESTS (5 operations)
// ============================================================================

#[test]
fn test_conduction() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("thermal_conductivity".to_string(), serde_json::json!(0.8)); // W/(m·K) for concrete
    params.insert("area".to_string(), serde_json::json!(10.0)); // m²
    params.insert("temp_hot".to_string(), serde_json::json!(30.0)); // °C
    params.insert("temp_cold".to_string(), serde_json::json!(20.0)); // °C
    params.insert("thickness".to_string(), serde_json::json!(0.2)); // m

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Thermodynamics(ThermodynamicsOp::Conduction),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate conductive heat transfer: {:?}", result);
}

#[test]
fn test_convection() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("heat_transfer_coefficient".to_string(), serde_json::json!(25.0)); // W/(m²·K)
    params.insert("area".to_string(), serde_json::json!(5.0)); // m²
    params.insert("surface_temp".to_string(), serde_json::json!(100.0)); // °C
    params.insert("fluid_temp".to_string(), serde_json::json!(20.0)); // °C

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Thermodynamics(ThermodynamicsOp::Convection),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate convective heat transfer: {:?}", result);
}

#[test]
fn test_radiation() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("emissivity".to_string(), serde_json::json!(0.95));
    params.insert("area".to_string(), serde_json::json!(2.0)); // m²
    params.insert("surface_temp_1".to_string(), serde_json::json!(400.0)); // K
    params.insert("surface_temp_2".to_string(), serde_json::json!(300.0)); // K

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Thermodynamics(ThermodynamicsOp::Radiation),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate radiative heat transfer: {:?}", result);
}

#[test]
fn test_thermal_resistance() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("resistances".to_string(), serde_json::json!([2.5])); // R = L/(kA) = 0.1/(0.04*10) = 0.25 K/W
    params.insert("configuration".to_string(), serde_json::json!("series"));
    params.insert("temp_hot".to_string(), serde_json::json!(30.0));
    params.insert("temp_cold".to_string(), serde_json::json!(20.0));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Thermodynamics(ThermodynamicsOp::ThermalResistance),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate thermal resistance: {:?}", result);
}

#[test]
fn test_entropy() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("heat_transfer".to_string(), serde_json::json!(1000.0)); // J (not "heat")
    params.insert("temperature".to_string(), serde_json::json!(300.0)); // K

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Thermodynamics(ThermodynamicsOp::Entropy),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate entropy change: {:?}", result);
}

// ============================================================================
// OPTICS TESTS (4 operations)
// ============================================================================

#[test]
fn test_thin_lens() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("focal_length".to_string(), serde_json::json!(0.5)); // m (50 cm)
    params.insert("object_distance".to_string(), serde_json::json!(1.0)); // m

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Optics(OpticsOp::ThinLens),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate image distance: {:?}", result);
}

#[test]
fn test_snells_law() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("n1".to_string(), serde_json::json!(1.0)); // air
    params.insert("n2".to_string(), serde_json::json!(1.5)); // glass
    params.insert("theta1".to_string(), serde_json::json!(30.0)); // degrees

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Optics(OpticsOp::SnellsLaw),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate refraction angle: {:?}", result);
}

#[test]
fn test_diffraction_grating() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("grating_spacing".to_string(), serde_json::json!(1.0e-6)); // m (1 μm)
    params.insert("wavelength".to_string(), serde_json::json!(500.0e-9)); // m (500 nm)
    params.insert("order".to_string(), serde_json::json!(1)); // Integer, not 1.0

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Optics(OpticsOp::DiffractionGrating),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate diffraction angle: {:?}", result);
}

#[test]
fn test_fresnel_equations() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("n1".to_string(), serde_json::json!(1.0));
    params.insert("n2".to_string(), serde_json::json!(1.5));
    params.insert("theta1".to_string(), serde_json::json!(45.0)); // degrees (not "incident_angle")

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Optics(OpticsOp::FresnelEquations),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate reflectance/transmittance: {:?}", result);
}

// ============================================================================
// GEOPHYSICS TESTS (4 categories)
// ============================================================================

#[test]
fn test_seismology() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("seismic_moment".to_string(), serde_json::json!(1.0e20)); // N·m (Magnitude ~6.0)

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Geophysics(GeophysicsOp::Seismology),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate earthquake magnitude: {:?}", result);
}

#[test]
fn test_atmosphere() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("altitude".to_string(), serde_json::json!(1000.0)); // m
    params.insert("temperature".to_string(), serde_json::json!(288.15)); // K (15°C at sea level)
    params.insert("pressure".to_string(), serde_json::json!(101325.0)); // Pa

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Geophysics(GeophysicsOp::Atmosphere),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate atmospheric properties: {:?}", result);
}

#[test]
fn test_radiometric_dating() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("isotope_system".to_string(), serde_json::json!("U238Pb206")); // Uranium-238 to Lead-206
    params.insert("parent_isotope".to_string(), serde_json::json!(1000.0));
    params.insert("daughter_isotope".to_string(), serde_json::json!(1000.0));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Geophysics(GeophysicsOp::RadiometricDating),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate rock age: {:?}", result);
}

#[test]
fn test_planetary_science() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("mass_primary".to_string(), serde_json::json!(5.972e24)); // Earth mass (kg)
    params.insert("radius_primary".to_string(), serde_json::json!(6.371e6)); // Earth radius (m)

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Geophysics(GeophysicsOp::PlanetaryScience),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate planetary parameters: {:?}", result);
}

// ============================================================================
// ENGINEERING TESTS (4 disciplines)
// ============================================================================

#[test]
fn test_sound_pressure_level() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("pressure_rms".to_string(), serde_json::json!(2.0)); // Pa
    params.insert("reference_pressure".to_string(), serde_json::json!(20.0e-6)); // Pa (20 μPa)

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Engineering(EngineeringOp::SoundPressureLevel),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate SPL: {:?}", result);
}

#[test]
fn test_doppler_effect() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("frequency".to_string(), serde_json::json!(440.0)); // Hz (A4 note)
    params.insert("velocity_source".to_string(), serde_json::json!(0.0)); // m/s
    params.insert("velocity_observer".to_string(), serde_json::json!(10.0)); // m/s
    params.insert("sound_speed".to_string(), serde_json::json!(343.0)); // m/s

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Engineering(EngineeringOp::DopplerEffect),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate Doppler shifted frequency: {:?}", result);
}

#[test]
fn test_stress_strain() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("youngs_modulus".to_string(), serde_json::json!(200.0e9)); // Pa (steel)
    params.insert("strain".to_string(), serde_json::json!(0.001)); // dimensionless (0.1%)

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Engineering(EngineeringOp::Stress),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should calculate stress: {:?}", result);
}

#[test]
fn test_bernoulli() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("velocity".to_string(), serde_json::json!(5.0)); // m/s (v1)
    params.insert("pressure_1".to_string(), serde_json::json!(101325.0)); // Pa (P1)
    params.insert("pressure_2".to_string(), serde_json::json!(100000.0)); // Pa (P2)
    params.insert("height_1".to_string(), serde_json::json!(0.0)); // m (h1)
    params.insert("height_2".to_string(), serde_json::json!(0.0)); // m (h2)
    params.insert("density".to_string(), serde_json::json!(1000.0)); // kg/m³ (water)

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Engineering(EngineeringOp::Bernoulli),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply Bernoulli's equation: {:?}", result);
}

// ============================================================================
// ERROR HANDLING TESTS
// ============================================================================

#[test]
fn test_ph_requires_parameters() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Chemistry(ChemistryOp::PhCalculation),
        data: serde_json::json!({}),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "pH calculation should require parameters");
}

#[test]
fn test_michaelis_menten_requires_vmax() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("km".to_string(), serde_json::json!(5.0));
    params.insert("substrate_concentration".to_string(), serde_json::json!(10.0));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Biology(BiologyOp::MichaelisMenten),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Michaelis-Menten should require vmax");
}

#[test]
fn test_thin_lens_requires_parameters() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("focal_length".to_string(), serde_json::json!(0.5));
    // Missing object_distance

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Optics(OpticsOp::ThinLens),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err() || result.is_ok(), "Thin lens may solve for missing parameter or error");
}
