// Unit tests for engineering::mod
use computational_engine::engineering::mod::*;

use super::*;

    // Acoustics Tests (6 tests)
    #[test]
    fn test_acoustics_spl_quiet() {
        let params = EngineeringParams {
            pressure_rms: Some(0.002), // 2 mPa
            ..Default::default()
        };
        let result = calculate_acoustics(&params).unwrap();
        // SPL = 20·log₁₀(0.002/20e-6) = 20·log₁₀(100) = 40 dB
        assert!((result.value - 40.0).abs() < 1.0);
        assert!(result.classification.unwrap().contains("Quiet"));
    }

    #[test]
    fn test_acoustics_spl_loud() {
        let params = EngineeringParams {
            pressure_rms: Some(0.2), // 200 mPa
            ..Default::default()
        };
        let result = calculate_acoustics(&params).unwrap();
        // SPL = 20·log₁₀(0.2/20e-6) = 20·log₁₀(10000) = 80 dB
        assert!((result.value - 80.0).abs() < 1.0);
        assert!(result.classification.unwrap().contains("Loud"));
    }

    #[test]
    fn test_acoustics_doppler_approaching() {
        let params = EngineeringParams {
            frequency: Some(1000.0),     // 1 kHz
            sound_speed: Some(343.0),    // Speed of sound in air
            velocity_source: Some(20.0), // Source moving toward observer (positive)
            ..Default::default()
        };
        let result = calculate_acoustics(&params).unwrap();
        // f' = f·(v + v_obs)/(v - v_src) = 1000·343/(343-20) = 1000·343/323 ≈ 1062 Hz
        assert!(result.value > 1000.0); // Frequency should increase
    }

    #[test]
    fn test_acoustics_doppler_receding() {
        let params = EngineeringParams {
            frequency: Some(1000.0),
            sound_speed: Some(343.0),
            velocity_source: Some(-20.0), // Source moving away (negative)
            ..Default::default()
        };
        let result = calculate_acoustics(&params).unwrap();
        // f' = 1000·343/(343-(-20)) = 1000·343/363 ≈ 945 Hz
        assert!(result.value < 1000.0); // Frequency should decrease
    }

    #[test]
    fn test_acoustics_reverberation_dry() {
        let params = EngineeringParams {
            room_volume: Some(100.0),          // 100 m³ room
            absorption_coefficient: Some(0.8), // High absorption (carpet, curtains)
            ..Default::default()
        };
        let result = calculate_acoustics(&params).unwrap();
        // RT60 should be short for dry room
        assert!(result.value < 1.0);
        assert!(result.classification.is_some());
    }

    #[test]
    fn test_acoustics_reverberation_reverberant() {
        let params = EngineeringParams {
            room_volume: Some(10000.0),        // Large cathedral
            absorption_coefficient: Some(0.1), // Low absorption (stone walls)
            ..Default::default()
        };
        let result = calculate_acoustics(&params).unwrap();
        // RT60 should be long for reverberant space
        assert!(result.value > 2.0);
    }

    // Materials Tests (5 tests)
    #[test]
    fn test_materials_hookes_law_elastic() {
        let params = EngineeringParams {
            youngs_modulus: Some(200e9), // Steel: 200 GPa
            strain: Some(0.0005),        // 0.05% strain (lower to stay in safe zone)
            ..Default::default()
        };
        let result = calculate_materials(&params).unwrap();
        // σ = E·ε = 200e9 · 0.0005 = 100e6 Pa = 100 MPa
        assert!((result.value - 100.0).abs() < 1.0);
        assert!(result.classification.unwrap().contains("safe"));
    }

    #[test]
    fn test_materials_hookes_law_strain_from_stress() {
        let params = EngineeringParams {
            youngs_modulus: Some(70e9), // Aluminum: 70 GPa
            stress: Some(140e6),        // 140 MPa
            ..Default::default()
        };
        let result = calculate_materials(&params).unwrap();
        // ε = σ/E = 140e6/70e9 = 0.002 = 0.2%
        assert!((result.value - 0.2).abs() < 0.01);
    }

    #[test]
    fn test_materials_fracture_mechanics_safe() {
        let params = EngineeringParams {
            stress: Some(100e6),       // 100 MPa
            crack_length: Some(0.001), // 1 mm crack
            geometry_factor: Some(1.0),
            ..Default::default()
        };
        let result = calculate_materials(&params).unwrap();
        // K = Y·σ·√(π·a) = 1.0 · 100e6 · √(π·0.001) ≈ 5.6 MPa·√m
        assert!(result.value < 20.0); // Safe range
        assert!(result.classification.unwrap().contains("Safe"));
    }

    #[test]
    fn test_materials_fracture_mechanics_critical() {
        let params = EngineeringParams {
            stress: Some(500e6),         // 500 MPa (very high)
            crack_length: Some(0.01),    // 10 mm crack
            geometry_factor: Some(1.12), // Edge crack
            ..Default::default()
        };
        let result = calculate_materials(&params).unwrap();
        // K = 1.12 · 500e6 · √(π·0.01) ≈ 99 MPa·√m (very high!)
        assert!(result.value > 40.0);
    }

    #[test]
    fn test_materials_yield_strength_check() {
        let params = EngineeringParams {
            youngs_modulus: Some(200e9),
            strain: Some(0.0025),        // 0.25% strain
            yield_strength: Some(250e6), // 250 MPa yield
            ..Default::default()
        };
        let result = calculate_materials(&params).unwrap();
        // σ = 200e9 · 0.0025 = 500e6 Pa = 500 MPa (exceeds yield!)
        assert!(result.value > 400.0);
        assert!(result.classification.unwrap().contains("Plastic"));
    }

    // Fluid Mechanics Tests (5 tests)
    #[test]
    fn test_fluid_bernoulli_exit_velocity() {
        let params = EngineeringParams {
            density: Some(1000.0),      // Water
            pressure_1: Some(200000.0), // 200 kPa
            pressure_2: Some(101325.0), // Atmospheric
            velocity: Some(1.0),        // 1 m/s inlet
            height_1: Some(0.0),
            height_2: Some(0.0), // Same height
            ..Default::default()
        };
        let result = calculate_fluid_mechanics(&params).unwrap();
        // P₁ + ½ρv₁² = P₂ + ½ρv₂²
        // 200000 + 500 = 101325 + 500·v₂²
        // v₂² = (200000 + 500 - 101325) / 500 ≈ 198.35
        // v₂ ≈ 14.08 m/s
        assert!(result.value > 10.0 && result.value < 20.0);
    }

    #[test]
    fn test_fluid_poiseuille_laminar() {
        let params = EngineeringParams {
            pressure_1: Some(1000.0), // 1 kPa pressure drop
            radius: Some(0.001),      // 1 mm radius pipe (smaller for laminar flow)
            length: Some(1.0),        // 1 meter long
            viscosity: Some(0.001),   // Water viscosity
            density: Some(1000.0),
            ..Default::default()
        };
        let result = calculate_fluid_mechanics(&params).unwrap();
        // Q = (π·ΔP·r⁴)/(8·η·L), Re = 250 (laminar)
        assert!(result.value > 0.0); // Positive flow rate
        assert!(result.classification.unwrap().contains("Laminar"));
    }

    #[test]
    fn test_fluid_drag_force() {
        let params = EngineeringParams {
            density: Some(1.225),            // Air at sea level
            velocity: Some(30.0),            // 30 m/s ≈ 108 km/h
            drag_coefficient: Some(0.3),     // Streamlined car
            cross_sectional_area: Some(2.0), // 2 m² frontal area
            ..Default::default()
        };
        let result = calculate_fluid_mechanics(&params).unwrap();
        // F_D = ½·ρ·v²·C_D·A = 0.5·1.225·900·0.3·2.0 = 330.75 N
        assert!((result.value - 330.75).abs() < 10.0);
    }

    #[test]
    fn test_fluid_drag_high_speed() {
        let params = EngineeringParams {
            density: Some(1.225),
            velocity: Some(100.0),       // 100 m/s ≈ 360 km/h
            drag_coefficient: Some(0.5), // Less streamlined
            cross_sectional_area: Some(1.5),
            ..Default::default()
        };
        let result = calculate_fluid_mechanics(&params).unwrap();
        // F_D = 0.5·1.225·10000·0.5·1.5 = 4593.75 N
        assert!(result.value > 4000.0);
    }

    #[test]
    fn test_fluid_poiseuille_reynolds() {
        let params = EngineeringParams {
            pressure_1: Some(1000.0),
            radius: Some(0.01), // 10 mm radius
            length: Some(2.0),
            viscosity: Some(0.001),
            density: Some(1000.0),
            ..Default::default()
        };
        let result = calculate_fluid_mechanics(&params).unwrap();
        let additional = result.additional.unwrap();
        // Should calculate Reynolds number
        assert!(additional.contains_key("reynolds_number"));
    }

    // Control Theory Tests (4 tests)
    #[test]
    fn test_control_pid_proportional() {
        let params = EngineeringParams {
            setpoint: Some(100.0),
            process_variable: Some(80.0),
            kp: Some(2.0),
            ..Default::default()
        };
        let result = calculate_control(&params).unwrap();
        // error = 100 - 80 = 20
        // P-term = 2.0 · 20 = 40
        assert!((result.value - 40.0).abs() < 0.1);
    }

    #[test]
    fn test_control_pid_negative_error() {
        let params = EngineeringParams {
            setpoint: Some(50.0),
            process_variable: Some(70.0), // Overshoot
            kp: Some(1.5),
            ..Default::default()
        };
        let result = calculate_control(&params).unwrap();
        // error = 50 - 70 = -20
        // P-term = 1.5 · (-20) = -30
        assert!((result.value - (-30.0)).abs() < 0.1);
    }

    #[test]
    fn test_control_ziegler_nichols_tuning() {
        let params = EngineeringParams {
            setpoint: Some(100.0),
            process_variable: Some(100.0), // At setpoint
            time_constant: Some(10.0),     // τ = 10 seconds
            ..Default::default()
        };
        let result = calculate_control(&params).unwrap();
        // Kp = 0.9·τ = 9.0
        assert!((result.value - 9.0).abs() < 0.1);
        assert!(result.interpretation.contains("Kp"));
    }

    #[test]
    fn test_control_pid_zero_error() {
        let params = EngineeringParams {
            setpoint: Some(75.0),
            process_variable: Some(75.0), // Perfect match
            kp: Some(3.0),
            ..Default::default()
        };
        let result = calculate_control(&params).unwrap();
        // error = 0, P-term = 0
        assert!((result.value).abs() < 0.001);
    }
