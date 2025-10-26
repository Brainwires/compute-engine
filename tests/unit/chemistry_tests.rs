use super::{calculate_chemistry, ChemistryInput, ChemistryOperation, ChemistryParams};

#[test]
fn test_ph_calculation() {
    let input = ChemistryInput {
        operation: ChemistryOperation::PhCalculation,
        parameters: ChemistryParams {
            pka: Some(4.76),
            concentration_acid: Some(0.1),
            concentration_base: Some(0.1),
            ..Default::default()
        },
    };

    let result = calculate_chemistry(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "pH units");
}

#[test]
fn test_buffer_capacity() {
    let input = ChemistryInput {
        operation: ChemistryOperation::BufferCapacity,
        parameters: ChemistryParams {
            concentration_acid: Some(0.1),
            concentration_base: Some(0.1),
            ..Default::default()
        },
    };

    let result = calculate_chemistry(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "mol/(L·pH)");
    assert!(result.value > 0.0);
}

#[test]
fn test_arrhenius() {
    let input = ChemistryInput {
        operation: ChemistryOperation::Arrhenius,
        parameters: ChemistryParams {
            activation_energy: Some(50.0), // kJ/mol
            temperature: Some(298.15),     // K
            pre_exponential: Some(1e10),   // A
            ..Default::default()
        },
    };

    let result = calculate_chemistry(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert!(result.value > 0.0);
    assert!(result.formula_used.contains("Arrhenius"));
}

#[test]
fn test_rate_law() {
    let input = ChemistryInput {
        operation: ChemistryOperation::RateLaw,
        parameters: ChemistryParams {
            rate_constant: Some(0.5),
            concentration: Some(0.1),
            order: Some(1),
            ..Default::default()
        },
    };

    let result = calculate_chemistry(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "M/s");
    assert!(result.value > 0.0);
}

#[test]
fn test_gibbs_free_energy() {
    let input = ChemistryInput {
        operation: ChemistryOperation::GibbsFreeEnergy,
        parameters: ChemistryParams {
            enthalpy: Some(-50.0),    // kJ/mol
            entropy: Some(100.0),     // J/(mol·K)
            temperature: Some(298.15), // K
            ..Default::default()
        },
    };

    let result = calculate_chemistry(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "kJ/mol");
    assert!(result.interpretation.len() > 0);
}

#[test]
fn test_nernst_equation() {
    let input = ChemistryInput {
        operation: ChemistryOperation::NernstEquation,
        parameters: ChemistryParams {
            standard_potential: Some(1.1), // V
            n_electrons: Some(2),
            temperature: Some(298.15),
            concentrations: Some(vec![0.1, 0.01]),
            ..Default::default()
        },
    };

    let result = calculate_chemistry(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "V");
    assert!(result.formula_used.contains("Nernst"));
}

#[test]
fn test_beer_lambert_absorbance() {
    let input = ChemistryInput {
        operation: ChemistryOperation::BeerLambert,
        parameters: ChemistryParams {
            molar_absorptivity: Some(1000.0), // L/(mol·cm)
            path_length: Some(1.0),            // cm
            concentration: Some(0.001),        // M
            ..Default::default()
        },
    };

    let result = calculate_chemistry(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "AU (absorbance units)");
    assert!(result.value > 0.0);
}

#[test]
fn test_beer_lambert_concentration() {
    let input = ChemistryInput {
        operation: ChemistryOperation::BeerLambert,
        parameters: ChemistryParams {
            absorbance: Some(0.5),
            molar_absorptivity: Some(1000.0),
            path_length: Some(1.0),
            ..Default::default()
        },
    };

    let result = calculate_chemistry(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "M");
    assert!(result.value > 0.0);
}

#[test]
fn test_van_der_waals() {
    let input = ChemistryInput {
        operation: ChemistryOperation::VanDerWaals,
        parameters: ChemistryParams {
            pressure: Some(1.0),    // atm
            volume: Some(22.4),     // L/mol
            temperature: Some(273.15), // K
            a_constant: Some(1.36), // CO2
            b_constant: Some(0.0318), // CO2
            ..Default::default()
        },
    };

    let result = calculate_chemistry(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.unit, "% deviation from ideal");
}

#[test]
fn test_rate_law_second_order() {
    let input = ChemistryInput {
        operation: ChemistryOperation::RateLaw,
        parameters: ChemistryParams {
            rate_constant: Some(0.5),
            concentration: Some(0.2),
            order: Some(2),
            ..Default::default()
        },
    };

    let result = calculate_chemistry(input);
    assert!(result.is_ok());
    let result = result.unwrap();
    assert!(result.interpretation.contains("Second order"));
}
