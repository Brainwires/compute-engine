// Unit tests for specialized::chemistry::mod
use computational_engine::specialized::chemistry::mod::*;

use super::*;

    #[test]
    fn test_molar_mass() {
        let result = molar_mass(MolarMassRequest {
            formula: "H2O".to_string(),
        })
        .unwrap();

        assert!((result.molar_mass - 18.015).abs() < 0.1);
    }

    #[test]
    fn test_gas_law() {
        let result = gas_law(GasLawRequest {
            pressure: Some(1.0),
            volume: Some(22.4),
            temperature: Some(273.15),
            moles: None,
            gas_type: None,
        })
        .unwrap();

        assert!((result.moles - 1.0).abs() < 0.01);
    }
