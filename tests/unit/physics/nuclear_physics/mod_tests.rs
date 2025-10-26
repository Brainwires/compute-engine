// Unit tests for physics::nuclear_physics::mod
use computational_engine::physics::nuclear_physics::mod::*;

use super::*;

    #[test]
    fn test_radioactive_decay() {
        let result = radioactive_decay(RadioactiveDecayRequest {
            initial_quantity: 1000.0,
            decay_constant: 0.693 / 5730.0, // C-14 (5730 years)
            time: 5730.0,                   // One half-life
        })
        .unwrap();

        assert!((result.fraction_remaining - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_half_life_calculation() {
        let result = half_life(HalfLifeRequest {
            decay_constant: Some(0.1),
            half_life: None,
            mean_lifetime: None,
        })
        .unwrap();

        assert!((result.half_life - 6.9315).abs() < 0.01);
    }

    #[test]
    fn test_mass_defect() {
        // Helium-4 nucleus
        let result = mass_defect(MassDefectRequest {
            protons: 2,
            neutrons: 2,
            nuclear_mass: 4.001506, // amu
        })
        .unwrap();

        assert!(result.mass_defect_amu > 0.0);
        assert!(result.binding_energy_mev > 28.0); // He-4 has ~28.3 MeV
    }
