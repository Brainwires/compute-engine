// Unit tests for specialized::machine_learning::neural_network
use crate::ml::*;

    #[test]
    fn test_activation_functions() {
        // Test sigmoid
        let sig = activate_scalar(0.0, ActivationFunction::Sigmoid);
        assert!((sig - 0.5).abs() < 1e-10);

        // Test ReLU
        assert_eq!(activate_scalar(5.0, ActivationFunction::ReLU), 5.0);
        assert_eq!(activate_scalar(-5.0, ActivationFunction::ReLU), 0.0);

        // Test tanh
        let tanh_val = activate_scalar(0.0, ActivationFunction::Tanh);
        assert!((tanh_val - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_softmax() {
        let x = vec![1.0, 2.0, 3.0];
        let result = activate_softmax(&x);

        // Sum should be 1.0
        let sum: f64 = result.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);

        // Values should be in (0, 1)
        for &val in &result {
            assert!(val > 0.0 && val < 1.0);
        }
    }

    #[test]
    fn test_dense_layer() {
        let layer = DenseLayer::new(3, 2, ActivationFunction::ReLU);
        let input = vec![1.0, 2.0, 3.0];
        let output = layer.forward(&input);

        assert_eq!(output.len(), 2);
    }

    #[test]
    fn test_neural_network() {
        let network = NeuralNetwork::new(&[
            (2, ActivationFunction::Linear),
            (4, ActivationFunction::ReLU),
            (1, ActivationFunction::Sigmoid),
        ]);

        let input = vec![0.5, 0.3];
        let output = network.predict(&input);

        assert_eq!(output.len(), 1);
        assert!(output[0] >= 0.0 && output[0] <= 1.0); // Sigmoid output
    }
