use super::{calculate_thermodynamics, ThermodynamicsInput, ThermodynamicsOperation, ThermodynamicsParams};

// ============================================================================
// CONDUCTION TESTS
// ============================================================================

#[test]
fn test_conduction_with_gradient() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Conduction,
        parameters: ThermodynamicsParams {
            thermal_conductivity: Some(200.0), // W/(m·K) - Aluminum
            area: Some(0.01),                  // m² - 10cm × 10cm
            temperature_gradient: Some(-100.0), // K/m - negative gradient
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "W");
    assert!(result.value > 0.0);
    assert!(result.formula_used.contains("Fourier"));
}

#[test]
fn test_conduction_with_temperatures() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Conduction,
        parameters: ThermodynamicsParams {
            thermal_conductivity: Some(0.6), // W/(m·K) - Wood
            area: Some(2.0),                 // m²
            thickness: Some(0.05),           // m - 5cm thick
            temp_hot: Some(293.15),          // K - 20°C
            temp_cold: Some(273.15),         // K - 0°C
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "W");
    assert!(result.value > 0.0);
    assert!(result.interpretation.contains("Conduction"));
}

#[test]
fn test_conduction_zero_thickness_error() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Conduction,
        parameters: ThermodynamicsParams {
            thermal_conductivity: Some(1.0),
            area: Some(1.0),
            thickness: Some(0.0), // Invalid
            temp_hot: Some(373.0),
            temp_cold: Some(293.0),
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("positive"));
}

#[test]
fn test_conduction_missing_parameters() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Conduction,
        parameters: ThermodynamicsParams {
            thermal_conductivity: Some(1.0),
            // Missing area
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_err());
}

// ============================================================================
// CONVECTION TESTS
// ============================================================================

#[test]
fn test_convection_natural() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Convection,
        parameters: ThermodynamicsParams {
            heat_transfer_coefficient: Some(5.0), // W/(m²·K) - Natural convection
            area: Some(1.0),                      // m²
            surface_temp: Some(323.15),           // K - 50°C
            fluid_temp: Some(293.15),             // K - 20°C
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "W");
    assert!(result.value > 0.0);
    assert!(result.interpretation.contains("Natural convection"));
}

#[test]
fn test_convection_forced() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Convection,
        parameters: ThermodynamicsParams {
            heat_transfer_coefficient: Some(50.0), // W/(m²·K) - Forced convection
            area: Some(2.5),                       // m²
            surface_temp: Some(373.15),            // K - 100°C
            fluid_temp: Some(293.15),              // K - 20°C
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "W");
    assert!(result.interpretation.contains("Forced convection"));
}

#[test]
fn test_convection_phase_change() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Convection,
        parameters: ThermodynamicsParams {
            heat_transfer_coefficient: Some(1000.0), // W/(m²·K) - Condensation
            area: Some(0.5),                         // m²
            surface_temp: Some(373.15),              // K
            fluid_temp: Some(373.15),                // K - no temp difference
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert!(result.interpretation.contains("Phase change") || result.interpretation.contains("high velocity"));
}

#[test]
fn test_convection_missing_temperature() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Convection,
        parameters: ThermodynamicsParams {
            heat_transfer_coefficient: Some(10.0),
            area: Some(1.0),
            surface_temp: Some(300.0),
            // Missing fluid_temp
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_err());
}

// ============================================================================
// RADIATION TESTS
// ============================================================================

#[test]
fn test_radiation_black_body() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Radiation,
        parameters: ThermodynamicsParams {
            emissivity: Some(1.0),      // Black body
            area: Some(1.0),            // m²
            surface_temp_1: Some(373.15), // K - 100°C
            surface_temp_2: Some(293.15), // K - 20°C
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "W");
    assert!(result.value > 0.0);
    assert!(result.interpretation.contains("Black body"));
}

#[test]
fn test_radiation_gray_body() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Radiation,
        parameters: ThermodynamicsParams {
            emissivity: Some(0.7),        // Gray body
            area: Some(2.0),              // m²
            surface_temp_1: Some(500.0),  // K
            surface_temp_2: Some(300.0),  // K
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert!(result.interpretation.contains("Gray body"));
}

#[test]
fn test_radiation_to_space() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Radiation,
        parameters: ThermodynamicsParams {
            emissivity: Some(0.9),
            area: Some(1.0),
            surface_temp_1: Some(300.0), // K
            // surface_temp_2 defaults to 0.0 (space)
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert!(result.value > 0.0);
}

#[test]
fn test_radiation_invalid_emissivity() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Radiation,
        parameters: ThermodynamicsParams {
            emissivity: Some(1.5), // Invalid: > 1
            area: Some(1.0),
            surface_temp_1: Some(300.0),
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("between 0 and 1"));
}

#[test]
fn test_radiation_negative_temperature() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Radiation,
        parameters: ThermodynamicsParams {
            emissivity: Some(0.8),
            area: Some(1.0),
            surface_temp_1: Some(-100.0), // Invalid: must be Kelvin
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("Kelvin"));
}

// ============================================================================
// THERMAL RESISTANCE TESTS
// ============================================================================

#[test]
fn test_resistance_series() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::ThermalResistance,
        parameters: ThermodynamicsParams {
            resistances: Some(vec![0.1, 0.2, 0.3]), // K/W
            configuration: Some("series".to_string()),
            temp_hot: Some(373.0),  // K
            temp_cold: Some(293.0), // K
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "W");
    assert!(result.formula_used.contains("series"));
    // Total R = 0.1 + 0.2 + 0.3 = 0.6 K/W
    // Q = (373 - 293) / 0.6 = 133.33 W
    assert!((result.value - 133.33).abs() < 0.1);
}

#[test]
fn test_resistance_parallel() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::ThermalResistance,
        parameters: ThermodynamicsParams {
            resistances: Some(vec![0.3, 0.3, 0.3]), // K/W
            configuration: Some("parallel".to_string()),
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert!(result.formula_used.contains("parallel"));
    // Total R = 1 / (1/0.3 + 1/0.3 + 1/0.3) = 0.1 K/W
}

#[test]
fn test_resistance_default_series() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::ThermalResistance,
        parameters: ThermodynamicsParams {
            resistances: Some(vec![0.5, 0.5]),
            // No configuration specified - defaults to series
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_ok());
}

#[test]
fn test_resistance_invalid_configuration() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::ThermalResistance,
        parameters: ThermodynamicsParams {
            resistances: Some(vec![0.1, 0.2]),
            configuration: Some("invalid".to_string()),
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_err());
    let err_msg = result.unwrap_err();
    assert!(err_msg.contains("series") || err_msg.contains("parallel"));
}

#[test]
fn test_resistance_empty_array() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::ThermalResistance,
        parameters: ThermodynamicsParams {
            resistances: Some(vec![]),
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("At least one"));
}

// ============================================================================
// ENTROPY TESTS
// ============================================================================

#[test]
fn test_entropy_clausius_heat_added() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Entropy,
        parameters: ThermodynamicsParams {
            heat_transfer: Some(1000.0), // J - heat added
            temperature: Some(300.0),    // K
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "J/K");
    assert!(result.value > 0.0);
    assert!(result.formula_used.contains("Clausius"));
    assert!(result.interpretation.contains("increases"));
}

#[test]
fn test_entropy_clausius_heat_removed() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Entropy,
        parameters: ThermodynamicsParams {
            heat_transfer: Some(-500.0), // J - heat removed
            temperature: Some(298.15),   // K
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert!(result.value < 0.0);
    assert!(result.interpretation.contains("decreases"));
}

#[test]
fn test_entropy_boltzmann() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Entropy,
        parameters: ThermodynamicsParams {
            num_microstates: Some(1e23), // Ω
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "J/K");
    assert!(result.formula_used.contains("Boltzmann"));
    assert!(result.interpretation.contains("microstates"));
}

#[test]
fn test_entropy_thermal_heating() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Entropy,
        parameters: ThermodynamicsParams {
            mass: Some(1.0),            // kg
            specific_heat: Some(4186.0), // J/(kg·K) - water
            initial_temp: Some(293.15),  // K - 20°C
            final_temp: Some(373.15),    // K - 100°C
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "J/K");
    assert!(result.value > 0.0);
    assert!(result.formula_used.contains("Thermal"));
    assert!(result.interpretation.contains("Heating"));
}

#[test]
fn test_entropy_thermal_cooling() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Entropy,
        parameters: ThermodynamicsParams {
            mass: Some(2.0),             // kg
            specific_heat: Some(900.0),  // J/(kg·K) - aluminum
            initial_temp: Some(400.0),   // K
            final_temp: Some(300.0),     // K
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert!(result.value < 0.0);
    assert!(result.interpretation.contains("Cooling"));
}

#[test]
fn test_entropy_invalid_temperature() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Entropy,
        parameters: ThermodynamicsParams {
            heat_transfer: Some(100.0),
            temperature: Some(-10.0), // Invalid: negative temperature
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("positive"));
}

#[test]
fn test_entropy_insufficient_parameters() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Entropy,
        parameters: ThermodynamicsParams {
            mass: Some(1.0),
            // Missing other required parameters
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("Insufficient parameters"));
}

#[test]
fn test_entropy_boltzmann_invalid_microstates() {
    let input = ThermodynamicsInput {
        operation: ThermodynamicsOperation::Entropy,
        parameters: ThermodynamicsParams {
            num_microstates: Some(-1.0), // Invalid: negative
            ..Default::default()
        },
    };

    let result = calculate_thermodynamics(input);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("positive"));
}
