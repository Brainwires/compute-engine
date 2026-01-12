// Unit tests for chemistry::mod
use computational_engine::chemistry::mod::*;

use super::*;

    // pH Calculation Tests
    #[test]
    fn test_ph_calculation_equal_concentrations() {
        let params = ChemistryParams {
            pka: Some(4.76),
            concentration_acid: Some(0.1),
            concentration_base: Some(0.1),
            ..Default::default()
        };
        let result = calculate_ph(&params).unwrap();
        assert!((result.value - 4.76).abs() < 0.01);
    }

    #[test]
    fn test_ph_calculation_acidic() {
        let params = ChemistryParams {
            pka: Some(4.76),
            concentration_acid: Some(0.2),
            concentration_base: Some(0.02),
            ..Default::default()
        };
        let result = calculate_ph(&params).unwrap();
        assert!(result.value < 4.76); // More acid = lower pH
        assert!(result.interpretation.contains("Acidic"));
    }

    #[test]
    fn test_ph_calculation_basic() {
        let params = ChemistryParams {
            pka: Some(4.76),
            concentration_acid: Some(0.01),
            concentration_base: Some(0.1),
            ..Default::default()
        };
        let result = calculate_ph(&params).unwrap();
        assert!(result.value > 4.76); // More base = higher pH
    }

    // Buffer Capacity Tests
    #[test]
    fn test_buffer_capacity_equal_concentrations() {
        let params = ChemistryParams {
            concentration_acid: Some(0.1),
            concentration_base: Some(0.1),
            ..Default::default()
        };
        let result = calculate_buffer_capacity(&params).unwrap();
        assert!(result.value > 0.0);
        assert_eq!(result.unit, "mol/(L·pH)");
    }

    #[test]
    fn test_buffer_capacity_unequal_concentrations() {
        let params = ChemistryParams {
            concentration_acid: Some(0.2),
            concentration_base: Some(0.05),
            ..Default::default()
        };
        let result = calculate_buffer_capacity(&params).unwrap();
        assert!(result.value > 0.0);
    }

    // Arrhenius Equation Tests
    #[test]
    fn test_arrhenius_room_temperature() {
        let params = ChemistryParams {
            activation_energy: Some(50.0),
            temperature: Some(298.0),
            pre_exponential: Some(1e10),
            ..Default::default()
        };
        let result = calculate_arrhenius(&params).unwrap();
        assert!(result.value > 0.0);
    }

    #[test]
    fn test_arrhenius_high_temperature() {
        let params_low = ChemistryParams {
            activation_energy: Some(50.0),
            temperature: Some(298.0),
            pre_exponential: Some(1e10),
            ..Default::default()
        };
        let result_low = calculate_arrhenius(&params_low).unwrap();

        let params_high = ChemistryParams {
            activation_energy: Some(50.0),
            temperature: Some(500.0), // Higher temp = faster reaction
            pre_exponential: Some(1e10),
            ..Default::default()
        };
        let result_high = calculate_arrhenius(&params_high).unwrap();

        // At higher temperature, rate constant should be significantly higher
        assert!(result_high.value > result_low.value * 10.0);
    }

    #[test]
    fn test_arrhenius_low_activation_energy() {
        let params = ChemistryParams {
            activation_energy: Some(10.0), // Low Ea = fast reaction
            temperature: Some(298.0),
            pre_exponential: Some(1e10),
            ..Default::default()
        };
        let result = calculate_arrhenius(&params).unwrap();
        assert!(result.value > 0.0);
    }

    // Rate Law Tests
    #[test]
    fn test_rate_law_first_order() {
        let params = ChemistryParams {
            rate_constant: Some(0.1),
            concentration: Some(2.0),
            order: Some(1),
            ..Default::default()
        };
        let result = calculate_rate_law(&params).unwrap();
        assert!((result.value - 0.2).abs() < 1e-10); // rate = 0.1 * 2^1 = 0.2
        assert!(result.interpretation.contains("First order"));
    }

    #[test]
    fn test_rate_law_second_order() {
        let params = ChemistryParams {
            rate_constant: Some(0.5),
            concentration: Some(3.0),
            order: Some(2),
            ..Default::default()
        };
        let result = calculate_rate_law(&params).unwrap();
        assert!((result.value - 4.5).abs() < 1e-10); // rate = 0.5 * 3^2 = 4.5
        assert!(result.interpretation.contains("Second order"));
    }

    #[test]
    fn test_rate_law_zero_order() {
        let params = ChemistryParams {
            rate_constant: Some(0.05),
            concentration: Some(10.0),
            order: Some(0),
            ..Default::default()
        };
        let result = calculate_rate_law(&params).unwrap();
        assert!((result.value - 0.05).abs() < 1e-10); // rate = 0.05 * 10^0 = 0.05
        assert!(result.interpretation.contains("Zero order"));
    }

    // Gibbs Free Energy Tests
    #[test]
    fn test_gibbs_spontaneous() {
        let params = ChemistryParams {
            enthalpy: Some(-100.0), // Exothermic
            entropy: Some(200.0),   // Entropy increase
            temperature: Some(298.0),
            ..Default::default()
        };
        let result = calculate_gibbs(&params).unwrap();
        assert!(result.value < 0.0); // Spontaneous
        assert!(result.interpretation.contains("Spontaneous"));
    }

    #[test]
    fn test_gibbs_nonspontaneous() {
        let params = ChemistryParams {
            enthalpy: Some(100.0), // Endothermic
            entropy: Some(-50.0),  // Entropy decrease
            temperature: Some(298.0),
            ..Default::default()
        };
        let result = calculate_gibbs(&params).unwrap();
        assert!(result.value > 0.0); // Non-spontaneous
        assert!(result.interpretation.contains("Non-spontaneous"));
    }

    #[test]
    fn test_gibbs_temperature_dependent() {
        // At low temp, enthalpy dominates; at high temp, entropy dominates
        let params_low = ChemistryParams {
            enthalpy: Some(-50.0),
            entropy: Some(-100.0),
            temperature: Some(100.0),
            ..Default::default()
        };
        let result_low = calculate_gibbs(&params_low).unwrap();
        assert!(result_low.value < 0.0); // Spontaneous at low T

        let params_high = ChemistryParams {
            enthalpy: Some(-50.0),
            entropy: Some(-100.0),
            temperature: Some(1000.0),
            ..Default::default()
        };
        let result_high = calculate_gibbs(&params_high).unwrap();
        assert!(result_high.value > 0.0); // Non-spontaneous at high T
    }

    // Nernst Equation Tests
    #[test]
    fn test_nernst_standard_conditions() {
        let params = ChemistryParams {
            standard_potential: Some(0.34), // Cu2+ reduction
            n_electrons: Some(2),
            temperature: Some(298.15),
            concentrations: Some(vec![1.0, 1.0]), // [reactants], [products]
            ..Default::default()
        };
        let result = calculate_nernst(&params).unwrap();
        // At Q=1, E = E° (ln(1) = 0)
        assert!((result.value - 0.34).abs() < 0.01);
    }

    #[test]
    fn test_nernst_shifted_potential() {
        let params = ChemistryParams {
            standard_potential: Some(1.0),
            n_electrons: Some(1),
            temperature: Some(298.15),
            concentrations: Some(vec![0.1, 1.0]), // Low reactant, high product
            ..Default::default()
        };
        let result = calculate_nernst(&params).unwrap();
        // Q > 1 should give E < E°
        assert!(result.value < 1.0);
    }

    // Beer-Lambert Law Tests
    #[test]
    fn test_beer_lambert_calculate_absorbance() {
        let params = ChemistryParams {
            molar_absorptivity: Some(1000.0),
            path_length: Some(1.0),
            concentration: Some(0.001),
            ..Default::default()
        };
        let result = calculate_beer_lambert(&params).unwrap();
        assert!((result.value - 1.0).abs() < 1e-10); // A = 1000 * 1 * 0.001 = 1.0
    }

    #[test]
    fn test_beer_lambert_calculate_concentration() {
        let params = ChemistryParams {
            absorbance: Some(0.5),
            molar_absorptivity: Some(500.0),
            path_length: Some(1.0),
            ..Default::default()
        };
        let result = calculate_beer_lambert(&params).unwrap();
        assert!((result.value - 0.001).abs() < 1e-10); // c = 0.5 / (500 * 1) = 0.001 M
    }

    #[test]
    fn test_beer_lambert_proportional() {
        // Double concentration = double absorbance
        let params1 = ChemistryParams {
            molar_absorptivity: Some(1000.0),
            path_length: Some(1.0),
            concentration: Some(0.001),
            ..Default::default()
        };
        let result1 = calculate_beer_lambert(&params1).unwrap();

        let params2 = ChemistryParams {
            molar_absorptivity: Some(1000.0),
            path_length: Some(1.0),
            concentration: Some(0.002),
            ..Default::default()
        };
        let result2 = calculate_beer_lambert(&params2).unwrap();

        assert!((result2.value / result1.value - 2.0).abs() < 1e-10);
    }

    // Van der Waals Tests
    #[test]
    fn test_van_der_waals_real_gas() {
        let params = ChemistryParams {
            pressure: Some(1.0),
            volume: Some(22.4), // Near ideal at STP
            temperature: Some(273.15),
            a_constant: Some(1.36), // CO2
            b_constant: Some(0.0318),
            ..Default::default()
        };
        let result = calculate_van_der_waals(&params).unwrap();
        assert!(result.value < 10.0); // Small deviation from ideal
    }

    #[test]
    fn test_van_der_waals_high_pressure() {
        let params = ChemistryParams {
            pressure: Some(100.0), // High pressure
            volume: Some(0.5),
            temperature: Some(300.0),
            a_constant: Some(1.36),
            b_constant: Some(0.0318),
            ..Default::default()
        };
        let result = calculate_van_der_waals(&params).unwrap();
        assert!(result.value > 0.0); // Larger deviation at high P
    }
