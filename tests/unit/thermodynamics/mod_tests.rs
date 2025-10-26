// Unit tests for thermodynamics::mod
use computational_engine::thermodynamics::mod::*;

use super::*;

    // Conduction Tests (3 tests)
    #[test]
    fn test_conduction_basic() {
        let params = ThermodynamicsParams {
            thermal_conductivity: Some(0.6), // W/(m·K) - insulation
            area: Some(10.0),
            temp_hot: Some(293.0),
            temp_cold: Some(273.0),
            thickness: Some(0.1),
            ..Default::default()
        };

        let result = calculate_conduction(&params).unwrap();
        assert!(result.value > 0.0);
        // q = k·A·ΔT/L = 0.6·10·20/0.1 = 1200 W
        assert!((result.value - 1200.0).abs() < 1.0);
        assert_eq!(result.unit, "W");
    }

    #[test]
    fn test_conduction_high_conductivity() {
        let params = ThermodynamicsParams {
            thermal_conductivity: Some(400.0), // W/(m·K) - copper
            area: Some(0.01),                  // 10 cm²
            temp_hot: Some(373.0),
            temp_cold: Some(293.0),
            thickness: Some(0.05), // 5 cm
            ..Default::default()
        };

        let result = calculate_conduction(&params).unwrap();
        // q = 400·0.01·80/0.05 = 6400 W
        assert!((result.value - 6400.0).abs() < 10.0);
    }

    #[test]
    fn test_conduction_with_gradient() {
        let params = ThermodynamicsParams {
            thermal_conductivity: Some(1.0),
            area: Some(5.0),
            temperature_gradient: Some(-100.0), // K/m (negative gradient)
            ..Default::default()
        };

        let result = calculate_conduction(&params).unwrap();
        // q = -k·A·grad = -1·5·(-100) = 500 W
        assert!((result.value - 500.0).abs() < 1.0);
    }

    // Convection Tests (3 tests)
    #[test]
    fn test_convection_natural() {
        let params = ThermodynamicsParams {
            heat_transfer_coefficient: Some(5.0), // W/(m²·K) - natural convection
            area: Some(2.0),
            surface_temp: Some(323.0), // 50°C
            fluid_temp: Some(293.0),   // 20°C
            ..Default::default()
        };

        let result = calculate_convection(&params).unwrap();
        // q = h·A·ΔT = 5·2·30 = 300 W
        assert!((result.value - 300.0).abs() < 1.0);
        assert!(result.interpretation.contains("Natural convection"));
    }

    #[test]
    fn test_convection_forced() {
        let params = ThermodynamicsParams {
            heat_transfer_coefficient: Some(50.0), // W/(m²·K) - forced convection
            area: Some(1.0),
            surface_temp: Some(400.0),
            fluid_temp: Some(300.0),
            ..Default::default()
        };

        let result = calculate_convection(&params).unwrap();
        // q = 50·1·100 = 5000 W
        assert!((result.value - 5000.0).abs() < 10.0);
        assert!(result.interpretation.contains("Forced convection"));
    }

    #[test]
    fn test_convection_high_h() {
        let params = ThermodynamicsParams {
            heat_transfer_coefficient: Some(500.0), // W/(m²·K) - phase change
            area: Some(0.5),
            surface_temp: Some(373.0), // Boiling water
            fluid_temp: Some(373.0),   // Same temp (minimal transfer)
            ..Default::default()
        };

        let result = calculate_convection(&params).unwrap();
        assert!(result.value.abs() < 1.0); // Nearly zero transfer
        assert!(result.interpretation.contains("Phase change"));
    }

    // Radiation Tests (3 tests)
    #[test]
    fn test_radiation_black_body() {
        let params = ThermodynamicsParams {
            emissivity: Some(1.0), // Black body
            area: Some(1.0),
            surface_temp_1: Some(373.0), // 100°C in K
            surface_temp_2: Some(293.0), // 20°C in K
            ..Default::default()
        };

        let result = calculate_radiation(&params).unwrap();
        assert!(result.value > 0.0);
        assert!(result.interpretation.contains("Black body"));
        // q = σ·ε·A·(T₁⁴ - T₂⁴) ≈ 5.67e-8·1·1·(373⁴ - 293⁴) ≈ 681 W
        assert!(result.value > 650.0 && result.value < 720.0);
    }

    #[test]
    fn test_radiation_to_space() {
        let params = ThermodynamicsParams {
            emissivity: Some(0.9),
            area: Some(1.0),
            surface_temp_1: Some(300.0), // Room temp
            surface_temp_2: None,        // Defaults to 0K (space)
            ..Default::default()
        };

        let result = calculate_radiation(&params).unwrap();
        // q = σ·ε·A·T₁⁴ ≈ 5.67e-8·0.9·1·300⁴ ≈ 413 W
        assert!(result.value > 400.0 && result.value < 450.0);
    }

    #[test]
    fn test_radiation_reflective_surface() {
        let params = ThermodynamicsParams {
            emissivity: Some(0.1), // Polished metal
            area: Some(1.0),
            surface_temp_1: Some(400.0),
            surface_temp_2: Some(300.0),
            ..Default::default()
        };

        let result = calculate_radiation(&params).unwrap();
        assert!(result.interpretation.contains("Reflective surface"));
        // Much lower than black body
        assert!(result.value > 0.0 && result.value < 100.0);
    }

    // Thermal Resistance Tests (3 tests)
    #[test]
    fn test_thermal_resistance_series() {
        let params = ThermodynamicsParams {
            resistances: Some(vec![0.1, 0.2, 0.3]), // K/W
            configuration: Some("series".to_string()),
            temp_hot: Some(373.0),
            temp_cold: Some(293.0),
            ..Default::default()
        };

        let result = calculate_resistance(&params).unwrap();
        // R_total = 0.1 + 0.2 + 0.3 = 0.6 K/W
        // q = ΔT/R = 80/0.6 = 133.33 W
        assert!((result.value - 133.33).abs() < 1.0);
        let info = result.additional_info.unwrap();
        assert!((info["thermal_resistance"].as_f64().unwrap() - 0.6).abs() < 0.01);
    }

    #[test]
    fn test_thermal_resistance_parallel() {
        let params = ThermodynamicsParams {
            resistances: Some(vec![0.3, 0.6]), // K/W
            configuration: Some("parallel".to_string()),
            temp_hot: Some(350.0),
            temp_cold: Some(300.0),
            ..Default::default()
        };

        let result = calculate_resistance(&params).unwrap();
        // 1/R_total = 1/0.3 + 1/0.6 = 3.333 + 1.667 = 5
        // R_total = 0.2 K/W
        // q = 50/0.2 = 250 W
        assert!((result.value - 250.0).abs() < 1.0);
        let info = result.additional_info.unwrap();
        assert!((info["thermal_resistance"].as_f64().unwrap() - 0.2).abs() < 0.01);
    }

    #[test]
    fn test_thermal_resistance_no_temps() {
        let params = ThermodynamicsParams {
            resistances: Some(vec![1.0, 2.0, 3.0]),
            configuration: Some("series".to_string()),
            ..Default::default()
        };

        let result = calculate_resistance(&params).unwrap();
        // Only R_total calculated, no heat rate
        assert_eq!(result.value, 0.0);
        let info = result.additional_info.unwrap();
        assert!((info["thermal_resistance"].as_f64().unwrap() - 6.0).abs() < 0.01);
    }

    // Entropy Tests (4 tests)
    #[test]
    fn test_entropy_clausius() {
        let params = ThermodynamicsParams {
            heat_transfer: Some(1000.0), // 1 kJ
            temperature: Some(300.0),    // 300 K
            ..Default::default()
        };

        let result = calculate_entropy(&params).unwrap();
        // ΔS = Q/T = 1000/300 = 3.333 J/K
        assert!((result.value - 3.333).abs() < 0.01);
        assert_eq!(result.unit, "J/K");
        assert!(result.interpretation.contains("increases"));
    }

    #[test]
    fn test_entropy_boltzmann() {
        let params = ThermodynamicsParams {
            num_microstates: Some(1e23), // Large number of states
            ..Default::default()
        };

        let result = calculate_entropy(&params).unwrap();
        // S = k·ln(Ω) = 1.381e-23·ln(1e23) ≈ 7.31e-22 J/K
        assert!(result.value > 0.0);
        assert!(result.formula_used.contains("Boltzmann"));
    }

    #[test]
    fn test_entropy_thermal_heating() {
        let params = ThermodynamicsParams {
            mass: Some(1.0),             // 1 kg water
            specific_heat: Some(4186.0), // J/(kg·K) for water
            initial_temp: Some(293.0),   // 20°C
            final_temp: Some(373.0),     // 100°C
            ..Default::default()
        };

        let result = calculate_entropy(&params).unwrap();
        // ΔS = m·c·ln(T₂/T₁) = 1·4186·ln(373/293) ≈ 1016 J/K
        assert!(result.value > 1000.0 && result.value < 1100.0);
        assert!(result.interpretation.contains("Heating"));
    }

    #[test]
    fn test_entropy_thermal_cooling() {
        let params = ThermodynamicsParams {
            mass: Some(2.0),
            specific_heat: Some(900.0), // J/(kg·K) for aluminum
            initial_temp: Some(400.0),
            final_temp: Some(300.0),
            ..Default::default()
        };

        let result = calculate_entropy(&params).unwrap();
        // ΔS = 2·900·ln(300/400) = 1800·ln(0.75) < 0 (cooling)
        assert!(result.value < 0.0);
        assert!(result.interpretation.contains("Cooling"));
    }
