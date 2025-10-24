//! Neural Network Primitives
//!
//! Basic building blocks for neural networks including:
//! - Dense (fully connected) layers
//! - Activation functions and their derivatives
//! - Forward and backward propagation

use super::ActivationFunction;
use serde::{Deserialize, Serialize};

/// Apply activation function element-wise to a vector
pub fn activate(x: &[f64], activation: ActivationFunction) -> Vec<f64> {
    x.iter().map(|&val| activate_scalar(val, activation)).collect()
}

/// Apply activation function to a single scalar
pub fn activate_scalar(x: f64, activation: ActivationFunction) -> f64 {
    match activation {
        ActivationFunction::Linear => x,
        ActivationFunction::Sigmoid => 1.0 / (1.0 + (-x).exp()),
        ActivationFunction::Tanh => x.tanh(),
        ActivationFunction::ReLU => x.max(0.0),
        ActivationFunction::LeakyReLU => {
            if x > 0.0 {
                x
            } else {
                0.01 * x
            }
        }
        ActivationFunction::Softmax => {
            // For single scalar, softmax doesn't make sense
            // This should be called on vectors via activate_softmax
            panic!("Use activate_softmax for softmax activation")
        }
    }
}

/// Apply softmax activation to a vector
pub fn activate_softmax(x: &[f64]) -> Vec<f64> {
    // Subtract max for numerical stability
    let max_val = x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let exp_vals: Vec<f64> = x.iter().map(|&val| (val - max_val).exp()).collect();
    let sum: f64 = exp_vals.iter().sum();
    exp_vals.iter().map(|&val| val / sum).collect()
}

/// Compute derivative of activation function
pub fn activate_derivative(x: &[f64], activation: ActivationFunction) -> Vec<f64> {
    x.iter()
        .map(|&val| activate_derivative_scalar(val, activation))
        .collect()
}

/// Compute derivative of activation function for a single scalar
pub fn activate_derivative_scalar(x: f64, activation: ActivationFunction) -> f64 {
    match activation {
        ActivationFunction::Linear => 1.0,
        ActivationFunction::Sigmoid => {
            let sig = 1.0 / (1.0 + (-x).exp());
            sig * (1.0 - sig)
        }
        ActivationFunction::Tanh => {
            let tanh_val = x.tanh();
            1.0 - tanh_val * tanh_val
        }
        ActivationFunction::ReLU => {
            if x > 0.0 {
                1.0
            } else {
                0.0
            }
        }
        ActivationFunction::LeakyReLU => {
            if x > 0.0 {
                1.0
            } else {
                0.01
            }
        }
        ActivationFunction::Softmax => {
            panic!("Use softmax_derivative for softmax activation")
        }
    }
}

/// Dense (fully connected) layer
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DenseLayer {
    /// Weight matrix (rows = output size, cols = input size)
    pub weights: Vec<Vec<f64>>,
    /// Bias vector (length = output size)
    pub biases: Vec<f64>,
    /// Activation function
    pub activation: ActivationFunction,
}

impl DenseLayer {
    /// Create a new dense layer with random initialization
    pub fn new(input_size: usize, output_size: usize, activation: ActivationFunction) -> Self {
        use rand::Rng;
        let mut rng = rand::thread_rng();

        // He initialization for ReLU, Xavier for others
        let scale = match activation {
            ActivationFunction::ReLU | ActivationFunction::LeakyReLU => {
                (2.0 / input_size as f64).sqrt()
            }
            _ => (1.0 / input_size as f64).sqrt(),
        };

        let weights: Vec<Vec<f64>> = (0..output_size)
            .map(|_| {
                (0..input_size)
                    .map(|_| rng.gen_range(-scale..scale))
                    .collect()
            })
            .collect();

        let biases = vec![0.0; output_size];

        Self {
            weights,
            biases,
            activation,
        }
    }

    /// Forward pass through the layer
    pub fn forward(&self, input: &[f64]) -> Vec<f64> {
        let z = self.linear_forward(input);

        if self.activation == ActivationFunction::Softmax {
            activate_softmax(&z)
        } else {
            activate(&z, self.activation)
        }
    }

    /// Linear forward pass (before activation)
    pub fn linear_forward(&self, input: &[f64]) -> Vec<f64> {
        self.weights
            .iter()
            .zip(&self.biases)
            .map(|(weights, &bias)| {
                let weighted_sum: f64 = weights.iter().zip(input).map(|(w, x)| w * x).sum();
                weighted_sum + bias
            })
            .collect()
    }

    /// Backward pass through the layer
    /// Returns (gradient w.r.t. input, gradient w.r.t. weights, gradient w.r.t. biases)
    pub fn backward(
        &self,
        input: &[f64],
        output_gradient: &[f64],
    ) -> (Vec<f64>, Vec<Vec<f64>>, Vec<f64>) {
        let z = self.linear_forward(input);

        // Compute activation gradient
        let activation_grad = if self.activation == ActivationFunction::Softmax {
            // For softmax with cross-entropy, the gradient simplifies
            output_gradient.to_vec()
        } else {
            let act_deriv = activate_derivative(&z, self.activation);
            output_gradient
                .iter()
                .zip(&act_deriv)
                .map(|(g, d)| g * d)
                .collect()
        };

        // Gradient w.r.t. weights: outer product of activation_grad and input
        let weight_gradient: Vec<Vec<f64>> = activation_grad
            .iter()
            .map(|&g| input.iter().map(|&x| g * x).collect())
            .collect();

        // Gradient w.r.t. biases is just the activation gradient
        let bias_gradient = activation_grad.clone();

        // Gradient w.r.t. input
        let input_gradient: Vec<f64> = (0..input.len())
            .map(|i| {
                self.weights
                    .iter()
                    .zip(&activation_grad)
                    .map(|(weights, &g)| weights[i] * g)
                    .sum()
            })
            .collect();

        (input_gradient, weight_gradient, bias_gradient)
    }

    /// Update weights and biases using gradients
    pub fn update(
        &mut self,
        weight_gradient: &[Vec<f64>],
        bias_gradient: &[f64],
        learning_rate: f64,
    ) {
        for (i, weight_row) in self.weights.iter_mut().enumerate() {
            for (j, weight) in weight_row.iter_mut().enumerate() {
                *weight -= learning_rate * weight_gradient[i][j];
            }
        }

        for (i, bias) in self.biases.iter_mut().enumerate() {
            *bias -= learning_rate * bias_gradient[i];
        }
    }
}

/// Simple feedforward neural network
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NeuralNetwork {
    /// Layers in the network
    pub layers: Vec<DenseLayer>,
}

impl NeuralNetwork {
    /// Create a new neural network with specified layer sizes
    pub fn new(layer_sizes: &[(usize, ActivationFunction)]) -> Self {
        let mut layers = Vec::new();

        for i in 0..layer_sizes.len() - 1 {
            let (input_size, _) = layer_sizes[i];
            let (output_size, activation) = layer_sizes[i + 1];

            layers.push(DenseLayer::new(input_size, output_size, activation));
        }

        Self { layers }
    }

    /// Forward pass through the entire network
    pub fn predict(&self, input: &[f64]) -> Vec<f64> {
        let mut activation = input.to_vec();

        for layer in &self.layers {
            activation = layer.forward(&activation);
        }

        activation
    }

    /// Forward pass storing intermediate activations for backpropagation
    pub fn forward_with_cache(&self, input: &[f64]) -> Vec<Vec<f64>> {
        let mut activations = vec![input.to_vec()];

        for layer in &self.layers {
            let next_activation = layer.forward(activations.last().unwrap());
            activations.push(next_activation);
        }

        activations
    }

    /// Backward pass through the network
    /// Returns gradients for each layer (weights, biases)
    pub fn backward(
        &self,
        activations: &[Vec<f64>],
        target: &[f64],
        loss_function: super::LossFunction,
    ) -> Vec<(Vec<Vec<f64>>, Vec<f64>)> {
        let output = activations.last().unwrap();

        // Compute output gradient based on loss function
        let mut gradient = compute_loss_gradient(output, target, loss_function);

        let mut gradients = Vec::new();

        // Backpropagate through layers in reverse
        for (i, layer) in self.layers.iter().enumerate().rev() {
            let layer_input = &activations[i];
            let (input_grad, weight_grad, bias_grad) = layer.backward(layer_input, &gradient);

            gradients.push((weight_grad, bias_grad));
            gradient = input_grad;
        }

        gradients.reverse();
        gradients
    }

    /// Update network weights using gradients
    pub fn update(&mut self, gradients: &[(Vec<Vec<f64>>, Vec<f64>)], learning_rate: f64) {
        for (layer, (weight_grad, bias_grad)) in self.layers.iter_mut().zip(gradients) {
            layer.update(weight_grad, bias_grad, learning_rate);
        }
    }
}

/// Compute gradient of loss function w.r.t. predictions
fn compute_loss_gradient(
    predictions: &[f64],
    targets: &[f64],
    loss_function: super::LossFunction,
) -> Vec<f64> {
    match loss_function {
        super::LossFunction::MSE => predictions
            .iter()
            .zip(targets)
            .map(|(&pred, &target)| 2.0 * (pred - target) / predictions.len() as f64)
            .collect(),

        super::LossFunction::MAE => predictions
            .iter()
            .zip(targets)
            .map(|(&pred, &target)| {
                if pred > target {
                    1.0 / predictions.len() as f64
                } else if pred < target {
                    -1.0 / predictions.len() as f64
                } else {
                    0.0
                }
            })
            .collect(),

        super::LossFunction::CrossEntropy => {
            // For softmax + cross-entropy, gradient simplifies to: pred - target
            predictions
                .iter()
                .zip(targets)
                .map(|(&pred, &target)| pred - target)
                .collect()
        }

        super::LossFunction::BinaryCrossEntropy => predictions
            .iter()
            .zip(targets)
            .map(|(&pred, &target)| {
                // Gradient: (pred - target) / (pred * (1 - pred))
                // Clamp predictions to avoid division by zero
                let pred_clamped = pred.max(1e-7).min(1.0 - 1e-7);
                (pred_clamped - target) / (pred_clamped * (1.0 - pred_clamped))
                    / predictions.len() as f64
            })
            .collect(),
    }
}

/// Compute loss value
pub fn compute_loss(
    predictions: &[f64],
    targets: &[f64],
    loss_function: super::LossFunction,
) -> f64 {
    match loss_function {
        super::LossFunction::MSE => {
            let sum: f64 = predictions
                .iter()
                .zip(targets)
                .map(|(&pred, &target)| (pred - target).powi(2))
                .sum();
            sum / predictions.len() as f64
        }

        super::LossFunction::MAE => {
            let sum: f64 = predictions
                .iter()
                .zip(targets)
                .map(|(&pred, &target)| (pred - target).abs())
                .sum();
            sum / predictions.len() as f64
        }

        super::LossFunction::CrossEntropy => {
            let sum: f64 = predictions
                .iter()
                .zip(targets)
                .map(|(&pred, &target)| {
                    if target > 0.0 {
                        -target * pred.max(1e-7).ln()
                    } else {
                        0.0
                    }
                })
                .sum();
            sum
        }

        super::LossFunction::BinaryCrossEntropy => {
            let sum: f64 = predictions
                .iter()
                .zip(targets)
                .map(|(&pred, &target)| {
                    let pred_clamped = pred.max(1e-7).min(1.0 - 1e-7);
                    -target * pred_clamped.ln() - (1.0 - target) * (1.0 - pred_clamped).ln()
                })
                .sum();
            sum / predictions.len() as f64
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
