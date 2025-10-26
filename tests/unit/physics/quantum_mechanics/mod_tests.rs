// Unit tests for physics::quantum_mechanics::mod
use computational_engine::physics::quantum_mechanics::mod::*;

use super::*;

    #[test]
    fn test_schrodinger_infinite_well() {
        let mut params = std::collections::HashMap::new();
        params.insert("width".to_string(), 1.0e-9); // 1 nm
        params.insert("quantum_number".to_string(), 1.0);

        let result = schrodinger_equation(SchrodingerRequest {
            potential: "infinite_well".to_string(),
            energy: 0.0,
            position: 0.5e-9,
            parameters: params,
        })
        .unwrap();

        assert!(result.wavefunction.abs() > 0.0);
        assert!(result.energy_eigenvalue > 0.0);
    }

    #[test]
    fn test_harmonic_oscillator() {
        let result = harmonic_oscillator(HarmonicOscillatorRequest {
            quantum_number: 0,
            omega: 1.0e15,
            mass: M_E,
            position: Some(0.0),
        })
        .unwrap();

        assert_eq!(result.energy, H_BAR * 1.0e15 / 2.0);
        assert!(result.wavefunction.unwrap() > 0.0);
    }

    #[test]
    fn test_hydrogen_atom() {
        let result = hydrogen_atom(HydrogenAtomRequest {
            n: 1,
            l: 0,
            m: 0,
            r: Some(A_0),
        })
        .unwrap();

        assert!(result.energy < 0.0); // Bound state
        assert_eq!(result.bohr_radius, A_0);
    }
