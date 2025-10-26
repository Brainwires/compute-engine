// Unit tests for geophysics::mod
use computational_engine::geophysics::mod::*;

use super::*;

    // ===== SEISMOLOGY TESTS =====

    #[test]
    fn test_seismology_moment_magnitude_large_quake() {
        let params = GeophysicsParams {
            seismic_moment: Some(1e20), // Large earthquake
            ..Default::default()
        };

        let result = calculate_seismology(&params).unwrap();
        assert!(result.value > 6.0 && result.value < 8.0);
        assert_eq!(result.unit, "Mw");
        assert!(result.additional_data.is_some());
        let additional = result.additional_data.unwrap();
        assert!(additional.contains_key("moment_magnitude_mw"));
        assert!(additional.contains_key("energy_joules"));
        assert!(additional.contains_key("energy_tnt_equivalent_tons"));
    }

    #[test]
    fn test_seismology_moment_magnitude_micro() {
        let params = GeophysicsParams {
            seismic_moment: Some(1e12), // Micro earthquake
            ..Default::default()
        };

        let result = calculate_seismology(&params).unwrap();
        assert!(result.value < 3.0);
        assert!(result.interpretation.contains("Micro"));
    }

    #[test]
    fn test_seismology_moment_magnitude_moderate() {
        let params = GeophysicsParams {
            seismic_moment: Some(1e18), // Moderate earthquake
            ..Default::default()
        };

        let result = calculate_seismology(&params).unwrap();
        assert!(result.value >= 5.0 && result.value < 6.0);
        assert!(result.interpretation.contains("Moderate"));
    }

    #[test]
    fn test_seismology_moment_magnitude_great() {
        let params = GeophysicsParams {
            seismic_moment: Some(1e23), // Great earthquake (like 2004 Indian Ocean)
            ..Default::default()
        };

        let result = calculate_seismology(&params).unwrap();
        assert!(result.value >= 8.0);
        assert!(result.interpretation.contains("Great"));
    }

    #[test]
    fn test_seismology_richter_from_energy() {
        let params = GeophysicsParams {
            energy: Some(1e15), // ~7.1 magnitude
            ..Default::default()
        };

        let result = calculate_seismology(&params).unwrap();
        assert!(result.value > 6.5 && result.value < 7.5);
        assert_eq!(result.unit, "ML (Richter)");
        let additional = result.additional_data.unwrap();
        assert!(additional.contains_key("richter_magnitude"));
        assert!(additional.contains_key("tnt_equivalent_tons"));
    }

    #[test]
    fn test_seismology_energy_conversion() {
        let params = GeophysicsParams {
            seismic_moment: Some(1e20),
            ..Default::default()
        };

        let result = calculate_seismology(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let energy_joules = additional.get("energy_joules").unwrap();
        let tnt_tons = additional.get("energy_tnt_equivalent_tons").unwrap();

        // Verify conversion: 1 ton TNT = 4.184e9 J
        assert!((*energy_joules / 4.184e9 - tnt_tons).abs() < 1.0);
    }

    #[test]
    fn test_seismology_missing_parameters() {
        let params = GeophysicsParams::default();
        let result = calculate_seismology(&params);
        assert!(result.is_err());
    }

    // ===== ATMOSPHERE TESTS =====

    #[test]
    fn test_atmosphere_pressure_at_altitude() {
        let params = GeophysicsParams {
            pressure: Some(101325.0),  // Sea level pressure (Pa)
            temperature: Some(288.15), // 15°C in Kelvin
            altitude: Some(1000.0),    // 1 km
            ..Default::default()
        };

        let result = calculate_atmosphere(&params).unwrap();
        // Pressure should decrease with altitude
        assert!(result.value < 101325.0);
        assert!(result.value > 85000.0); // Reasonable range
        assert_eq!(result.unit, "Pa");
        assert!(result.interpretation.contains("1000m altitude"));
    }

    #[test]
    fn test_atmosphere_scale_height() {
        let params = GeophysicsParams {
            pressure: Some(101325.0),
            temperature: Some(288.15),
            altitude: Some(8500.0), // Approximately one scale height
            ..Default::default()
        };

        let result = calculate_atmosphere(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let scale_height = additional.get("scale_height_m").unwrap();

        // Scale height for Earth's atmosphere ~8400m
        assert!((*scale_height - 8400.0).abs() < 500.0);

        // After one scale height, pressure should drop to ~1/e of initial
        let pressure_drop = additional.get("pressure_drop_percent").unwrap();
        assert!((*pressure_drop - 63.2).abs() < 5.0); // 1-1/e ≈ 63.2%
    }

    #[test]
    fn test_atmosphere_high_altitude() {
        let params = GeophysicsParams {
            pressure: Some(101325.0),
            temperature: Some(288.15),
            altitude: Some(10000.0), // 10 km (cruising altitude)
            ..Default::default()
        };

        let result = calculate_atmosphere(&params).unwrap();
        // Pressure at 10 km should be ~30% of sea level (~31000 Pa)
        assert!(result.value < 32000.0);
        assert!(result.value > 29000.0);
    }

    #[test]
    fn test_atmosphere_temperature_celsius_conversion() {
        let params_celsius = GeophysicsParams {
            pressure: Some(101325.0),
            temperature: Some(15.0), // Celsius
            altitude: Some(1000.0),
            ..Default::default()
        };

        let params_kelvin = GeophysicsParams {
            pressure: Some(101325.0),
            temperature: Some(288.15), // Kelvin equivalent
            altitude: Some(1000.0),
            ..Default::default()
        };

        let result_c = calculate_atmosphere(&params_celsius).unwrap();
        let result_k = calculate_atmosphere(&params_kelvin).unwrap();

        // Should give same result
        assert!((result_c.value - result_k.value).abs() < 10.0);
    }

    #[test]
    fn test_atmosphere_dew_point() {
        let params = GeophysicsParams {
            temperature: Some(20.0),       // 20°C
            relative_humidity: Some(60.0), // 60% RH
            ..Default::default()
        };

        let result = calculate_atmosphere(&params).unwrap();
        // Dew point should be below air temperature
        assert!(result.value < 20.0);
        assert!(result.value > 10.0); // Reasonable range for 60% RH
        assert_eq!(result.unit, "°C");
        assert!(result.interpretation.contains("Dew point"));
        assert!(result.interpretation.contains("60% RH"));
    }

    #[test]
    fn test_atmosphere_saturation_vapor_pressure() {
        let params = GeophysicsParams {
            temperature: Some(25.0),        // 25°C
            relative_humidity: Some(100.0), // Saturated
            ..Default::default()
        };

        let result = calculate_atmosphere(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let e_sat = additional.get("saturation_vapor_pressure_pa").unwrap();
        let e_actual = additional.get("actual_vapor_pressure_pa").unwrap();

        // At 100% RH, actual = saturation
        assert!((e_sat - e_actual).abs() < 1.0);

        // Saturation vapor pressure at 25°C ~3167 Pa
        assert!((*e_sat - 3167.0).abs() < 100.0);
    }

    #[test]
    fn test_atmosphere_dew_point_low_humidity() {
        let params = GeophysicsParams {
            temperature: Some(30.0),       // 30°C
            relative_humidity: Some(20.0), // Low humidity
            ..Default::default()
        };

        let result = calculate_atmosphere(&params).unwrap();
        // Low humidity means large difference between temp and dew point
        assert!(30.0 - result.value > 15.0);
    }

    #[test]
    fn test_atmosphere_missing_parameters() {
        let params = GeophysicsParams {
            pressure: Some(101325.0),
            // Missing temperature and altitude
            ..Default::default()
        };

        let result = calculate_atmosphere(&params);
        assert!(result.is_err());
    }

    // ===== RADIOMETRIC DATING TESTS =====

    #[test]
    fn test_carbon_dating_one_half_life() {
        let params = GeophysicsParams {
            parent_isotope: Some(50.0),
            daughter_isotope: Some(50.0),
            isotope_system: Some("C14".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        // When P=D, age should be ~1 half-life
        assert!((result.value - 5730.0).abs() < 100.0);
        assert_eq!(result.unit, "years");
        assert!(result.interpretation.contains("years"));
    }

    #[test]
    fn test_carbon_dating_recent() {
        let params = GeophysicsParams {
            parent_isotope: Some(90.0),
            daughter_isotope: Some(10.0),
            isotope_system: Some("C14".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        // Young sample (less than 1 half-life)
        assert!(result.value < 5730.0);
        assert!(result.value > 0.0);
    }

    #[test]
    fn test_carbon_dating_old() {
        let params = GeophysicsParams {
            parent_isotope: Some(10.0),
            daughter_isotope: Some(90.0),
            isotope_system: Some("C14".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        // Old sample (multiple half-lives)
        assert!(result.value > 10000.0);
        assert!(result.value < 40000.0); // C14 limit
    }

    #[test]
    fn test_uranium_dating() {
        let params = GeophysicsParams {
            parent_isotope: Some(50.0),
            daughter_isotope: Some(50.0),
            isotope_system: Some("U238Pb206".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        // One half-life of U-238 is 4.468 billion years
        assert!((result.value - 4.468e9).abs() < 1e8);
        assert!(result.interpretation.contains("billion years"));
    }

    #[test]
    fn test_potassium_dating() {
        let params = GeophysicsParams {
            parent_isotope: Some(75.0),
            daughter_isotope: Some(25.0),
            isotope_system: Some("K40Ar40".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        // Younger than 1 half-life (1.25 billion years)
        assert!(result.value < 1.25e9);
        assert!(result.value > 0.0);
    }

    #[test]
    fn test_rubidium_dating() {
        let params = GeophysicsParams {
            parent_isotope: Some(50.0),
            daughter_isotope: Some(50.0),
            isotope_system: Some("Rb87Sr87".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        // Rb-87 has very long half-life (48.8 billion years)
        assert!(result.value > 4e10);
        assert!(result.value < 5e10);
        assert!(result.interpretation.contains("billion years"));
    }

    #[test]
    fn test_dating_custom_half_life() {
        let params = GeophysicsParams {
            parent_isotope: Some(50.0),
            daughter_isotope: Some(50.0),
            isotope_system: Some("Custom".to_string()),
            half_life: Some(10000.0),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        // Should use provided half-life
        assert!((result.value - 10000.0).abs() < 200.0);
    }

    #[test]
    fn test_dating_decay_constant() {
        let params = GeophysicsParams {
            parent_isotope: Some(50.0),
            daughter_isotope: Some(50.0),
            isotope_system: Some("C14".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let decay_const = additional.get("decay_constant").unwrap();
        let half_life = additional.get("half_life_years").unwrap();

        // λ = ln(2) / t½
        assert!((decay_const * half_life - 0.693147).abs() < 1e-5);
    }

    #[test]
    fn test_dating_daughter_parent_ratio() {
        let params = GeophysicsParams {
            parent_isotope: Some(25.0),
            daughter_isotope: Some(75.0), // 3:1 D/P ratio
            isotope_system: Some("C14".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let ratio = additional.get("daughter_parent_ratio").unwrap();

        assert!((*ratio - 3.0).abs() < 0.01);
        // High D/P ratio means old sample
        assert!(result.value > 11000.0); // ~2 half-lives
    }

    #[test]
    fn test_dating_missing_isotope_system() {
        let params = GeophysicsParams {
            parent_isotope: Some(50.0),
            daughter_isotope: Some(50.0),
            ..Default::default()
        };

        let result = calculate_dating(&params);
        assert!(result.is_err());
    }

    #[test]
    fn test_dating_missing_isotope_amounts() {
        let params = GeophysicsParams {
            isotope_system: Some("C14".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params);
        assert!(result.is_err());
    }

    // ===== PLANETARY SCIENCE TESTS =====

    #[test]
    fn test_escape_velocity_earth() {
        let params = GeophysicsParams {
            mass_primary: Some(5.972e24),  // Earth mass (kg)
            radius_primary: Some(6.371e6), // Earth radius (m)
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        // Earth's escape velocity ~11.2 km/s
        assert!((result.value - 11.2).abs() < 0.2);
        assert_eq!(result.unit, "km/s");
        assert!(result.interpretation.contains("Escape velocity"));
    }

    #[test]
    fn test_escape_velocity_mars() {
        let params = GeophysicsParams {
            mass_primary: Some(6.4171e23),  // Mars mass
            radius_primary: Some(3.3895e6), // Mars radius
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        // Mars escape velocity ~5.0 km/s
        assert!((result.value - 5.0).abs() < 0.3);
    }

    #[test]
    fn test_escape_velocity_jupiter() {
        let params = GeophysicsParams {
            mass_primary: Some(1.898e27),   // Jupiter mass
            radius_primary: Some(6.9911e7), // Jupiter radius
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        // Jupiter escape velocity ~59.5 km/s
        assert!((result.value - 59.5).abs() < 2.0);
    }

    #[test]
    fn test_escape_velocity_moon() {
        let params = GeophysicsParams {
            mass_primary: Some(7.342e22),  // Moon mass
            radius_primary: Some(1.737e6), // Moon radius
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        // Moon escape velocity ~2.38 km/s
        assert!((result.value - 2.38).abs() < 0.1);
    }

    #[test]
    fn test_orbital_velocity_surface() {
        let params = GeophysicsParams {
            mass_primary: Some(5.972e24),
            radius_primary: Some(6.371e6),
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let v_orbital = additional.get("orbital_velocity_surface").unwrap();

        // Earth surface orbital velocity ~7.9 km/s
        assert!((v_orbital / 1000.0 - 7.9).abs() < 0.2);

        // Escape velocity should be √2 times orbital velocity
        let v_escape = result.value * 1000.0; // Convert back to m/s
        assert!((v_escape / v_orbital - 2_f64.sqrt()).abs() < 0.01);
    }

    #[test]
    fn test_roche_limit_earth_moon() {
        let params = GeophysicsParams {
            radius_primary: Some(6.371e6),   // Earth radius
            density_primary: Some(5514.0),   // Earth density (kg/m³)
            density_secondary: Some(3344.0), // Moon density
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        // Earth-Moon Roche limit ~18,365 km (for fluid body)
        assert!((result.value - 18365.0).abs() < 500.0);
        assert_eq!(result.unit, "km");
    }

    #[test]
    fn test_roche_limit_saturn_rings() {
        let params = GeophysicsParams {
            radius_primary: Some(5.8232e7),  // Saturn radius
            density_primary: Some(687.0),    // Saturn density
            density_secondary: Some(1000.0), // Ice density
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let roche_radii = additional.get("roche_limit_primary_radii").unwrap();

        // Saturn's Roche limit is ~2.15 Saturn radii (density ratio < 1)
        assert!((*roche_radii - 2.15).abs() < 0.1);

        // Saturn's rings are inside Roche limit
        assert!(result.value > 100000.0); // >100,000 km
    }

    #[test]
    fn test_roche_limit_rigid_vs_fluid() {
        let params = GeophysicsParams {
            radius_primary: Some(6.371e6),
            density_primary: Some(5514.0),
            density_secondary: Some(3344.0),
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let rigid = additional.get("roche_limit_rigid").unwrap();
        let fluid = additional.get("roche_limit_fluid").unwrap();

        // Rigid body Roche limit (2.46) is slightly larger than fluid (2.44)
        assert!(rigid > fluid);
        assert!((rigid / fluid - 1.008).abs() < 0.01);
    }

    #[test]
    fn test_roche_limit_equal_densities() {
        let params = GeophysicsParams {
            radius_primary: Some(6.371e6),
            density_primary: Some(5000.0),
            density_secondary: Some(5000.0), // Same density
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let roche_radii = additional.get("roche_limit_primary_radii").unwrap();

        // When densities are equal, Roche limit = 2.44 radii
        assert!((*roche_radii - 2.44).abs() < 0.01);
    }

    #[test]
    fn test_planetary_missing_parameters() {
        let params = GeophysicsParams {
            mass_primary: Some(5.972e24),
            // Missing radius
            ..Default::default()
        };

        let result = calculate_planetary(&params);
        assert!(result.is_err());
    }

    // ===== INTEGRATION TESTS =====

    #[test]
    fn test_full_workflow_seismology() {
        let input = GeophysicsInput {
            category: GeophysicsCategory::Seismology,
            parameters: GeophysicsParams {
                seismic_moment: Some(5e19),
                ..Default::default()
            },
        };

        let result = calculate_geophysics(input).unwrap();
        assert!(result.value > 0.0);
        assert!(result.uncertainty.is_some());
        assert!(!result.interpretation.is_empty());
    }

    #[test]
    fn test_full_workflow_atmosphere() {
        let input = GeophysicsInput {
            category: GeophysicsCategory::Atmosphere,
            parameters: GeophysicsParams {
                pressure: Some(101325.0),
                temperature: Some(288.15),
                altitude: Some(5000.0),
                ..Default::default()
            },
        };

        let result = calculate_geophysics(input).unwrap();
        assert!(result.value > 0.0);
        assert_eq!(result.unit, "Pa");
    }

    #[test]
    fn test_full_workflow_dating() {
        let input = GeophysicsInput {
            category: GeophysicsCategory::RadiometricDating,
            parameters: GeophysicsParams {
                parent_isotope: Some(60.0),
                daughter_isotope: Some(40.0),
                isotope_system: Some("U238Pb206".to_string()),
                ..Default::default()
            },
        };

        let result = calculate_geophysics(input).unwrap();
        assert!(result.value > 0.0);
        assert_eq!(result.unit, "years");
    }

    #[test]
    fn test_full_workflow_planetary() {
        let input = GeophysicsInput {
            category: GeophysicsCategory::PlanetaryScience,
            parameters: GeophysicsParams {
                mass_primary: Some(5.972e24),
                radius_primary: Some(6.371e6),
                ..Default::default()
            },
        };

        let result = calculate_geophysics(input).unwrap();
        assert!(result.value > 0.0);
        assert_eq!(result.unit, "km/s");
    }
