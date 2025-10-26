// Unit tests for specialized::control_theory::state_space
use computational_engine::specialized::control_theory::state_space::*;

use super::*;

    #[test]
    fn test_state_space_creation() {
        // Simple integrator
        let sys = StateSpace::new(
            vec![vec![0.0]],
            vec![vec![1.0]],
            vec![vec![1.0]],
            vec![vec![0.0]],
        );

        assert!(sys.is_ok());
        let sys = sys.unwrap();
        assert_eq!(sys.num_states(), 1);
        assert_eq!(sys.num_inputs(), 1);
        assert_eq!(sys.num_outputs(), 1);
    }

    #[test]
    fn test_output_computation() {
        // y = 2x (C = [2], D = [0])
        let sys = StateSpace::new(
            vec![vec![0.0]],
            vec![vec![1.0]],
            vec![vec![2.0]],
            vec![vec![0.0]],
        )
        .unwrap();

        let mut sys_mut = sys.clone();
        sys_mut.set_state(vec![5.0]).unwrap();

        let output = sys_mut.output(&[0.0]).unwrap();
        assert_eq!(output[0], 10.0);
    }

    #[test]
    fn test_integrator_simulation() {
        // Integrator: ẋ = u
        let mut sys = StateSpace::new(
            vec![vec![0.0]],
            vec![vec![1.0]],
            vec![vec![1.0]],
            vec![vec![0.0]],
        )
        .unwrap();

        let dt = 0.1;
        let input = vec![1.0]; // Constant input

        // Step 10 times
        for _ in 0..10 {
            sys.step(&input, dt).unwrap();
        }

        // After 10 steps with dt=0.1 and u=1, x should be ~1.0
        assert!((sys.state[0] - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_siso_creation() {
        let sys = siso_state_space(vec![vec![-1.0]], vec![1.0], vec![1.0], 0.0);

        assert!(sys.is_ok());
        let sys = sys.unwrap();
        assert_eq!(sys.num_states(), 1);
        assert_eq!(sys.num_inputs(), 1);
        assert_eq!(sys.num_outputs(), 1);
    }

    #[test]
    fn test_dimension_validation() {
        // Mismatched B dimensions
        let sys = StateSpace::new(
            vec![vec![0.0, 1.0], vec![-1.0, 0.0]],
            vec![vec![1.0]],        // Should be 2×1, but only 1×1
            vec![vec![1.0, 0.0]],
            vec![vec![0.0]],
        );

        assert!(sys.is_err());
    }
