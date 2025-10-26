// Unit tests for specialized::machine_learning::optimization
use computational_engine::specialized::machine_learning::*;

use super::*;

    #[test]
    fn test_optimizer_state() {
        let mut optimizer = OptimizerState::new(OptimizerType::Adam, 0.01, 2);
        optimizer.initialize(&[(3, 2), (2, 1)]);

        assert!(optimizer.velocity.is_some());
        assert!(optimizer.second_moment.is_some());
    }

    #[test]
    fn test_sgd_update() {
        let optimizer = OptimizerState::new(OptimizerType::SGD, 0.1, 1);
        let gradients = vec![vec![vec![1.0, 2.0], vec![3.0, 4.0]]];

        let updates = optimizer.update_sgd(&gradients);

        assert_eq!(updates[0][0][0], 0.1);
        assert_eq!(updates[0][0][1], 0.2);
    }
