// Unit tests for physics::electromagnetism::mod
use computational_engine::physics::electromagnetism::mod::*;

use super::*;

    #[test]
    fn test_em_wave() {
        let result = em_wave(WaveRequest {
            frequency: 1e9, // 1 GHz
            wavelength: None,
            medium: "vacuum".to_string(),
            permittivity: None,
            permeability: None,
        })
        .unwrap();

        assert!((result.wavelength - 0.3).abs() < 0.01);
    }

    #[test]
    fn test_skin_effect() {
        let result = skin_effect(SkinEffectRequest {
            frequency: 1e6,      // 1 MHz
            conductivity: 5.8e7, // Copper
            permeability: None,
        })
        .unwrap();

        assert!(result.skin_depth > 0.0);
    }
