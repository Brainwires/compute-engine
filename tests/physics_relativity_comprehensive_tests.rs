//! Comprehensive relativity test suite
//!
//! Tests for all relativity operations including:
//! - Special Relativity (5 operations)
//! - General Relativity (4 operations)
//! - Black Hole Physics (1 operation)

use computational_engine::create_default_dispatcher;
use computational_engine::engine::*;
use std::collections::HashMap;

// Physical constants for reference
const C: f64 = 299792458.0; // Speed of light (m/s)
const SOLAR_MASS: f64 = 1.989e30; // kg

// ============================================================================
// SPECIAL RELATIVITY TESTS (5 operations)
// ============================================================================

#[test]
fn test_lorentz_transform() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("velocity".to_string(), serde_json::json!(0.6 * C)); // 0.6c
    params.insert("position".to_string(), serde_json::json!([100.0, 0.0, 0.0])); // 100m in x
    params.insert("time".to_string(), serde_json::json!(1.0e-6)); // 1 microsecond

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::Relativity(RelativityOp::LorentzTransform)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Lorentz transform should succeed: {:?}",
        result
    );

    if let Ok(ToolResponse::Compute(output)) = result {
        // Check that gamma is correct (should be 1.25 for 0.6c)
        let gamma = output.result.get("gamma").and_then(|v| v.as_f64()).unwrap();
        assert!(
            (gamma - 1.25).abs() < 0.01,
            "Gamma should be ~1.25 for 0.6c"
        );

        // Check that position transformed
        let position_prime = output.result.get("position_prime").unwrap();
        assert!(
            position_prime.is_array(),
            "Transformed position should be an array"
        );
    }
}

#[test]
fn test_time_dilation() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("proper_time".to_string(), serde_json::json!(1.0)); // 1 second
    params.insert("velocity".to_string(), serde_json::json!(0.8 * C)); // 0.8c

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::Relativity(RelativityOp::TimeDilation)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Time dilation should succeed: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        let dilated_time = output
            .result
            .get("dilated_time")
            .and_then(|v| v.as_f64())
            .unwrap();
        // At 0.8c, gamma ≈ 1.667, so dilated time should be ~1.667s
        assert!(
            dilated_time > 1.6 && dilated_time < 1.7,
            "Dilated time should be ~1.667s"
        );

        let time_difference = output
            .result
            .get("time_difference")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(time_difference > 0.0, "Time difference should be positive");
    }
}

#[test]
fn test_length_contraction() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("proper_length".to_string(), serde_json::json!(10.0)); // 10 meters
    params.insert("velocity".to_string(), serde_json::json!(0.6 * C)); // 0.6c

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::Relativity(RelativityOp::LengthContraction)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Length contraction should succeed: {:?}",
        result
    );

    if let Ok(ToolResponse::Compute(output)) = result {
        let contracted_length = output
            .result
            .get("contracted_length")
            .and_then(|v| v.as_f64())
            .unwrap();
        // At 0.6c, gamma = 1.25, so contracted length = 10/1.25 = 8m
        assert!(
            (contracted_length - 8.0).abs() < 0.1,
            "Contracted length should be ~8m"
        );

        let contraction_factor = output
            .result
            .get("contraction_factor")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(
            contraction_factor < 1.0,
            "Contraction factor should be less than 1"
        );
    }
}

#[test]
fn test_relativistic_energy() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("mass".to_string(), serde_json::json!(1.0)); // 1 kg
    params.insert("velocity".to_string(), serde_json::json!(0.5 * C)); // 0.5c

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::Relativity(RelativityOp::RelativisticEnergy)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Relativistic energy should succeed: {:?}",
        result
    );

    if let Ok(ToolResponse::Compute(output)) = result {
        let rest_energy = output
            .result
            .get("rest_energy")
            .and_then(|v| v.as_f64())
            .unwrap();
        let expected_rest = C * C;
        assert!(
            (rest_energy - expected_rest).abs() / expected_rest < 0.01,
            "Rest energy should be E=mc²"
        );

        let total_energy = output
            .result
            .get("total_energy")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(
            total_energy > rest_energy,
            "Total energy should exceed rest energy"
        );

        let kinetic_energy = output
            .result
            .get("kinetic_energy")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(kinetic_energy > 0.0, "Kinetic energy should be positive");

        let momentum = output
            .result
            .get("momentum")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(momentum > 0.0, "Momentum should be positive");
    }
}

#[test]
fn test_velocity_addition() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("velocity1".to_string(), serde_json::json!(0.6 * C)); // 0.6c
    params.insert("velocity2".to_string(), serde_json::json!(0.6 * C)); // 0.6c

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::Relativity(RelativityOp::VelocityAddition)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Velocity addition should succeed: {:?}",
        result
    );

    if let Ok(ToolResponse::Compute(output)) = result {
        let relativistic_sum = output
            .result
            .get("relativistic_sum")
            .and_then(|v| v.as_f64())
            .unwrap();
        // Relativistic: (0.6c + 0.6c) / (1 + 0.6*0.6) = 1.2c / 1.36 ≈ 0.882c
        assert!(relativistic_sum < C, "Relativistic sum must be less than c");
        assert!(
            (relativistic_sum / C - 0.882).abs() < 0.01,
            "Should be ~0.882c"
        );

        let classical_sum = output
            .result
            .get("classical_sum")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(
            classical_sum > relativistic_sum,
            "Classical sum should exceed relativistic"
        );
    }
}

// ============================================================================
// GENERAL RELATIVITY TESTS (4 operations)
// ============================================================================

#[test]
fn test_schwarzschild_metric() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("mass".to_string(), serde_json::json!(SOLAR_MASS)); // Solar mass
    params.insert("radius".to_string(), serde_json::json!(1.0e10)); // 10,000 km from center

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::Relativity(RelativityOp::SchwarzschildMetric)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Schwarzschild metric should succeed: {:?}",
        result
    );

    if let Ok(ToolResponse::Compute(output)) = result {
        let rs = output
            .result
            .get("schwarzschild_radius")
            .and_then(|v| v.as_f64())
            .unwrap();
        // Sun's Schwarzschild radius ~2.95 km
        assert!(
            (rs - 2950.0).abs() < 100.0,
            "Sun's Schwarzschild radius should be ~2.95 km"
        );

        let g00 = output
            .result
            .get("metric_g00")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(g00 < 0.0, "Time component should be negative");

        let grr = output
            .result
            .get("metric_grr")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(grr > 0.0, "Radial component should be positive");

        let escape_velocity = output
            .result
            .get("escape_velocity")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(
            escape_velocity > 0.0 && escape_velocity < C,
            "Escape velocity should be less than c"
        );
    }
}

#[test]
fn test_gravitational_time_dilation() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("mass".to_string(), serde_json::json!(SOLAR_MASS));
    params.insert("radius".to_string(), serde_json::json!(7.0e8)); // Sun's radius
    params.insert("proper_time".to_string(), serde_json::json!(1.0)); // 1 second

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::Relativity(
            RelativityOp::GravitationalTimeDilation,
        )),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Gravitational time dilation should succeed: {:?}",
        result
    );

    if let Ok(ToolResponse::Compute(output)) = result {
        let time_dilation_factor = output
            .result
            .get("time_dilation_factor")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(
            time_dilation_factor > 0.0 && time_dilation_factor < 1.0,
            "Dilation factor should be < 1"
        );

        let coordinate_time = output
            .result
            .get("coordinate_time")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(
            coordinate_time > 1.0,
            "Coordinate time should be greater than proper time"
        );

        let time_difference = output
            .result
            .get("time_difference")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(time_difference > 0.0, "Time difference should be positive");
    }
}

#[test]
fn test_orbital_precession() {
    let dispatcher = create_default_dispatcher();

    // Mercury's orbit parameters
    let mut params = HashMap::new();
    params.insert("mass".to_string(), serde_json::json!(SOLAR_MASS));
    params.insert("semi_major_axis".to_string(), serde_json::json!(5.79e10)); // ~0.387 AU
    params.insert("eccentricity".to_string(), serde_json::json!(0.206)); // Mercury's eccentricity

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::Relativity(RelativityOp::OrbitalPrecession)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Orbital precession should succeed: {:?}",
        result
    );

    if let Ok(ToolResponse::Compute(output)) = result {
        let precession_per_orbit = output
            .result
            .get("precession_per_orbit")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(precession_per_orbit > 0.0, "Precession should be positive");

        let precession_per_century = output
            .result
            .get("precession_per_century")
            .and_then(|v| v.as_f64())
            .unwrap();
        // Mercury's precession is ~43 arcseconds/century
        assert!(
            precession_per_century > 30.0 && precession_per_century < 60.0,
            "Mercury's precession should be ~43 arcsec/century, got {}",
            precession_per_century
        );

        let orbital_period = output
            .result
            .get("orbital_period")
            .and_then(|v| v.as_f64())
            .unwrap();
        // Mercury's period ~88 days = ~7.6 million seconds
        assert!(
            orbital_period > 7.0e6 && orbital_period < 8.0e6,
            "Mercury's period should be ~88 days"
        );
    }
}

#[test]
fn test_gravitational_lensing() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("lens_mass".to_string(), serde_json::json!(SOLAR_MASS));
    params.insert("impact_parameter".to_string(), serde_json::json!(7.0e8)); // Sun's radius

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::Relativity(RelativityOp::GravitationalLensing)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Gravitational lensing should succeed: {:?}",
        result
    );

    if let Ok(ToolResponse::Compute(output)) = result {
        let deflection_angle = output
            .result
            .get("deflection_angle")
            .and_then(|v| v.as_f64())
            .unwrap();
        // Einstein's prediction for light grazing the sun is ~1.75 arcseconds = 8.5e-6 radians
        assert!(
            deflection_angle > 0.0,
            "Deflection angle should be positive"
        );
        assert!(
            deflection_angle < 1.0e-4,
            "Deflection should be small for solar mass"
        );

        let einstein_radius = output
            .result
            .get("einstein_radius")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(einstein_radius > 0.0, "Einstein radius should be positive");
    }
}

// ============================================================================
// BLACK HOLE PHYSICS TESTS (1 operation)
// ============================================================================

#[test]
fn test_black_hole_properties() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("mass".to_string(), serde_json::json!(SOLAR_MASS)); // Solar mass black hole

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::Relativity(RelativityOp::BlackHoleProperties)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(
        result.is_ok(),
        "Black hole properties should succeed: {:?}",
        result
    );

    if let Ok(ToolResponse::Compute(output)) = result {
        let schwarzschild_radius = output
            .result
            .get("schwarzschild_radius")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(
            (schwarzschild_radius - 2950.0).abs() < 100.0,
            "Event horizon should be ~2.95 km"
        );

        let surface_gravity = output
            .result
            .get("surface_gravity")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(surface_gravity > 0.0, "Surface gravity should be positive");

        let hawking_temperature = output
            .result
            .get("hawking_temperature")
            .and_then(|v| v.as_f64())
            .unwrap();
        // Solar mass black hole has incredibly low temperature (~60 nanokelvin)
        assert!(
            hawking_temperature > 0.0 && hawking_temperature < 1.0e-6,
            "Hawking temperature should be extremely low for solar mass"
        );

        let hawking_luminosity = output
            .result
            .get("hawking_luminosity")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(
            hawking_luminosity > 0.0,
            "Hawking radiation should be positive"
        );

        let evaporation_time = output
            .result
            .get("evaporation_time")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(
            evaporation_time > 1.0e50,
            "Evaporation time should be astronomical"
        );

        let entropy = output
            .result
            .get("entropy")
            .and_then(|v| v.as_f64())
            .unwrap();
        assert!(
            entropy > 0.0,
            "Bekenstein-Hawking entropy should be positive"
        );
    }
}

// ============================================================================
// ERROR HANDLING TESTS
// ============================================================================

#[test]
fn test_time_dilation_requires_velocity() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("proper_time".to_string(), serde_json::json!(1.0));
    // Missing velocity parameter

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::Relativity(RelativityOp::TimeDilation)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Should require velocity parameter");
}

#[test]
fn test_velocity_must_be_less_than_c() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("proper_time".to_string(), serde_json::json!(1.0));
    params.insert("velocity".to_string(), serde_json::json!(1.5 * C)); // Faster than light!

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::Relativity(RelativityOp::TimeDilation)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Should reject superluminal velocities");
}

#[test]
fn test_schwarzschild_radius_must_be_outside_horizon() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("mass".to_string(), serde_json::json!(SOLAR_MASS));
    params.insert("radius".to_string(), serde_json::json!(1000.0)); // Inside event horizon!

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::Relativity(RelativityOp::SchwarzschildMetric)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Should reject radius inside event horizon");
}
