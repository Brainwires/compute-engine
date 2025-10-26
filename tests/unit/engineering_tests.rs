use super::{
    calculate_engineering, EngineeringDiscipline, EngineeringInput, EngineeringParams,
};

// ============================================================================
// ACOUSTICS TESTS
// ============================================================================

#[test]
fn test_sound_pressure_level_quiet() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::Acoustics,
        parameters: EngineeringParams {
            pressure_rms: Some(0.002), // 2 mPa
            reference_pressure: Some(20e-6), // 20 μPa standard
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "dB SPL");
    assert!(result.value > 30.0 && result.value < 50.0); // Should be in "Quiet" range
    assert!(result.formula_used.contains("SPL"));
    assert!(result.classification.is_some());
}

#[test]
fn test_sound_pressure_level_loud() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::Acoustics,
        parameters: EngineeringParams {
            pressure_rms: Some(0.2), // 200 mPa - loud
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "dB SPL");
    assert!(result.value > 70.0); // Should be loud
    assert!(result.additional.is_some());
}

#[test]
fn test_doppler_effect_approaching() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::Acoustics,
        parameters: EngineeringParams {
            frequency: Some(1000.0),    // 1 kHz
            sound_speed: Some(343.0),   // m/s (20°C)
            velocity_source: Some(20.0), // Source approaching (positive in formula)
            velocity_observer: Some(0.0),
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "Hz");
    assert!(result.value > 1000.0); // Frequency should increase
    assert!(result.formula_used.contains("Doppler"));
    assert!(result.additional.is_some());
}

#[test]
fn test_doppler_effect_receding() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::Acoustics,
        parameters: EngineeringParams {
            frequency: Some(1000.0),   // 1 kHz
            sound_speed: Some(343.0),  // m/s
            velocity_source: Some(-20.0), // Source receding (negative)
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert!(result.value < 1000.0); // Frequency should decrease
}

#[test]
fn test_reverberation_time_dry_room() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::Acoustics,
        parameters: EngineeringParams {
            room_volume: Some(200.0),           // 200 m³
            absorption_coefficient: Some(0.4),  // High absorption
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "seconds");
    assert!(result.value < 1.0); // Should be dry
    assert!(result.formula_used.contains("Sabine"));
}

#[test]
fn test_reverberation_time_reverberant_room() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::Acoustics,
        parameters: EngineeringParams {
            room_volume: Some(1000.0),          // Large volume
            absorption_coefficient: Some(0.05), // Low absorption
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert!(result.value > 1.0); // Should be reverberant
    assert!(result.classification.is_some());
}

// ============================================================================
// MATERIALS SCIENCE TESTS
// ============================================================================

#[test]
fn test_hookes_law_stress_from_strain() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::Materials,
        parameters: EngineeringParams {
            youngs_modulus: Some(200e9),  // 200 GPa (steel)
            strain: Some(0.001),           // 0.1% strain
            yield_strength: Some(250e6),  // 250 MPa
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "MPa");
    assert!((result.value - 200.0).abs() < 1.0); // 200 MPa expected
    assert!(result.formula_used.contains("Hooke"));
    assert!(result.classification.unwrap().contains("Elastic"));
}

#[test]
fn test_hookes_law_strain_from_stress() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::Materials,
        parameters: EngineeringParams {
            youngs_modulus: Some(70e9), // 70 GPa (aluminum)
            stress: Some(100e6),         // 100 MPa
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "% strain");
    assert!(result.value > 0.0);
    assert!(result.formula_used.contains("Hooke"));
}

#[test]
fn test_fracture_mechanics_safe_crack() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::Materials,
        parameters: EngineeringParams {
            stress: Some(100e6),         // 100 MPa
            crack_length: Some(0.001),   // 1 mm crack
            geometry_factor: Some(1.0),  // Simple geometry
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "MPa·√m");
    assert!(result.formula_used.contains("Fracture"));
    assert!(result.classification.is_some());
}

#[test]
fn test_fracture_mechanics_critical_crack() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::Materials,
        parameters: EngineeringParams {
            stress: Some(300e6),         // 300 MPa - high stress
            crack_length: Some(0.01),    // 10 mm crack - large
            geometry_factor: Some(1.12), // Edge crack
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert!(result.value > 10.0); // Should show significant stress intensity
    assert!(result.additional.is_some());
}

// ============================================================================
// FLUID MECHANICS TESTS
// ============================================================================

#[test]
fn test_bernoulli_equation() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::FluidMechanics,
        parameters: EngineeringParams {
            density: Some(1000.0),    // Water
            pressure_1: Some(200000.0), // 200 kPa
            velocity: Some(2.0),      // 2 m/s
            height_1: Some(0.0),      // Ground level
            pressure_2: Some(100000.0), // 100 kPa
            height_2: Some(0.0),
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "m/s");
    assert!(result.value > 2.0); // Velocity should increase with pressure drop
    assert!(result.formula_used.contains("Bernoulli"));
}

#[test]
fn test_poiseuille_flow_laminar() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::FluidMechanics,
        parameters: EngineeringParams {
            pressure_1: Some(1000.0),    // 1 kPa pressure difference
            radius: Some(0.005),         // 5 mm pipe
            length: Some(1.0),           // 1 m long
            viscosity: Some(0.001),      // Water viscosity (1 cP)
            density: Some(1000.0),       // Water density
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "L/s");
    assert!(result.value > 0.0);
    assert!(result.formula_used.contains("Poiseuille"));
    assert!(result.classification.is_some());
    assert!(result.additional.is_some());
}

#[test]
fn test_poiseuille_flow_turbulent_detection() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::FluidMechanics,
        parameters: EngineeringParams {
            pressure_1: Some(100000.0),  // High pressure - will cause high flow
            radius: Some(0.05),          // 50 mm pipe
            length: Some(1.0),
            viscosity: Some(0.001),
            density: Some(1000.0),
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    // Should detect turbulent flow and warn
    let additional = result.additional.unwrap();
    assert!(additional.contains_key("reynolds_number"));
}

#[test]
fn test_drag_force_low_speed() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::FluidMechanics,
        parameters: EngineeringParams {
            density: Some(1.225),              // Air at sea level
            velocity: Some(10.0),              // 10 m/s
            drag_coefficient: Some(0.47),      // Sphere
            cross_sectional_area: Some(1.0),   // 1 m²
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "N");
    assert!(result.value > 0.0);
    assert!(result.formula_used.contains("Drag"));
    assert!(result.additional.is_some());
}

#[test]
fn test_drag_force_high_speed() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::FluidMechanics,
        parameters: EngineeringParams {
            density: Some(1.225),
            velocity: Some(30.0),             // 30 m/s (108 km/h)
            drag_coefficient: Some(0.3),      // Streamlined car
            cross_sectional_area: Some(2.2),  // Typical car
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert!(result.value > 100.0); // Should have significant drag
    let additional = result.additional.unwrap();
    assert!(additional.contains_key("drag_power_w"));
}

// ============================================================================
// CONTROL THEORY TESTS
// ============================================================================

#[test]
fn test_pid_error_calculation() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::ControlTheory,
        parameters: EngineeringParams {
            setpoint: Some(100.0),
            process_variable: Some(80.0),
            kp: Some(1.5),
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "control output");
    assert!(result.formula_used.contains("PID"));
    assert!(result.additional.is_some());
    let additional = result.additional.unwrap();
    assert_eq!(additional.get("error").unwrap(), &20.0);
}

#[test]
fn test_ziegler_nichols_tuning() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::ControlTheory,
        parameters: EngineeringParams {
            setpoint: Some(50.0),
            process_variable: Some(45.0),
            time_constant: Some(10.0), // τ = 10s
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "gain");
    assert!(result.formula_used.contains("Ziegler-Nichols"));
    assert!(result.classification.is_some());
    assert!(result.additional.is_some());
    let additional = result.additional.unwrap();
    assert!(additional.contains_key("zn_kp_suggested"));
    assert!(additional.contains_key("zn_ti_suggested"));
}

#[test]
fn test_pid_missing_setpoint_error() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::ControlTheory,
        parameters: EngineeringParams {
            process_variable: Some(80.0),
            // Missing setpoint
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("Setpoint required"));
}

#[test]
fn test_pid_missing_process_variable_error() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::ControlTheory,
        parameters: EngineeringParams {
            setpoint: Some(100.0),
            // Missing process variable
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("Process variable required"));
}

// ============================================================================
// ERROR HANDLING TESTS
// ============================================================================

#[test]
fn test_acoustics_insufficient_parameters() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::Acoustics,
        parameters: EngineeringParams {
            // No parameters provided
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("Insufficient acoustic parameters"));
}

#[test]
fn test_materials_insufficient_parameters() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::Materials,
        parameters: EngineeringParams {
            // Only one parameter, need two for any calculation
            stress: Some(100e6),
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("Insufficient materials parameters"));
}

#[test]
fn test_fluid_mechanics_insufficient_parameters() {
    let input = EngineeringInput {
        discipline: EngineeringDiscipline::FluidMechanics,
        parameters: EngineeringParams {
            density: Some(1000.0),
            // Not enough for any calculation
            ..Default::default()
        },
    };

    let result = calculate_engineering(input);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("Insufficient fluid mechanics parameters"));
}
