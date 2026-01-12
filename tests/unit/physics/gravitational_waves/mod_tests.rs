// Unit tests for physics::gravitational_waves::mod
use computational_engine::compute::physics::gravitational_waves::mod::*;

use super::*;

    #[test]
    fn test_binary_mass_parameters() {
        let binary = BinarySystem {
            mass1: 36.0,
            mass2: 29.0,
            distance: 410.0, // Mpc (GW150914-like)
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryBlackHole,
        };

        let total = binary.total_mass();
        assert_eq!(total, 65.0);

        let chirp = binary.chirp_mass();
        assert!(chirp > 25.0 && chirp < 35.0); // ~30 M_sun for GW150914

        let eta = binary.symmetric_mass_ratio();
        assert!(eta > 0.2 && eta < 0.26); // Should be ~0.247
    }

    #[test]
    fn test_isco_frequency() {
        let binary = BinarySystem {
            mass1: 30.0,
            mass2: 30.0,
            distance: 100.0,
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryBlackHole,
        };

        let f_isco = binary.isco_frequency();

        // ISCO frequency should be positive and finite
        assert!(f_isco > 0.0);
        assert!(f_isco.is_finite());
    }

    #[test]
    fn test_time_to_coalescence() {
        let binary = BinarySystem {
            mass1: 1.4,
            mass2: 1.4,
            distance: 100.0,
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryNeutronStar,
        };

        // Time to coalescence from 10 Hz
        let t_coal = binary.time_to_coalescence(10.0);

        // For BNS, should be minutes to hours at 10 Hz
        assert!(t_coal > 60.0); // More than 1 minute
    }

    #[test]
    fn test_binary_types() {
        let bbh = BinarySystem {
            mass1: 30.0,
            mass2: 30.0,
            distance: 410.0,
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryBlackHole,
        };

        let bns = BinarySystem {
            mass1: 1.4,
            mass2: 1.4,
            distance: 40.0,
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryNeutronStar,
        };

        assert_eq!(bbh.system_type, BinaryType::BinaryBlackHole);
        assert_eq!(bns.system_type, BinaryType::BinaryNeutronStar);
    }
