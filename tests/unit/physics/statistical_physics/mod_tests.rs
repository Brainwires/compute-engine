// Unit tests for physics::statistical_physics::mod
use computational_engine::compute::physics::statistical_physics::mod::*;

use super::*;

    #[test]
    fn test_maxwell_boltzmann() {
        let result = maxwell_boltzmann(MaxwellBoltzmannRequest {
            temperature: 300.0,
            mass: 4.65e-26, // Nitrogen molecule
            velocity: Some(500.0),
        })
        .unwrap();

        assert!(result.probability_density.is_some());
        assert!(result.average_speed > 0.0);
        assert!(result.rms_speed > result.average_speed);
    }

    #[test]
    fn test_fermi_dirac() {
        let result = fermi_dirac(FermiDiracRequest {
            energy: 1.0e-19, // 1 eV in Joules
            temperature: 300.0,
            chemical_potential: 0.8e-19,
        })
        .unwrap();

        assert!(result.occupation_probability > 0.0 && result.occupation_probability < 1.0);
    }

    #[test]
    fn test_canonical_partition() {
        let result = canonical_partition(CanonicalPartitionRequest {
            temperature: 300.0,
            energy_levels: vec![0.0, 1.0e-20, 2.0e-20],
            degeneracies: vec![1, 2, 1],
        })
        .unwrap();

        assert!(result.partition_function > 0.0);
        assert!(result.internal_energy >= 0.0);
    }
