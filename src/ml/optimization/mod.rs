//! Optimization Algorithms for Machine Learning
//!
//! Gradient descent variants:
//! - SGD (Stochastic Gradient Descent)
//! - SGD with Momentum
//! - RMSprop
//! - Adam (Adaptive Moment Estimation)

use crate::ml::{OptimizerType, TrainingConfig, TrainingHistory};
use crate::ml::neural_network::{NeuralNetwork, compute_loss};
use serde::{Deserialize, Serialize};

/// Optimizer state for maintaining momentum, velocity, etc.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptimizerState {
    /// Optimizer type
    pub optimizer_type: OptimizerType,
    /// Learning rate
    pub learning_rate: f64,
    /// Momentum velocity (for Momentum and Adam)
    pub velocity: Option<Vec<Vec<Vec<f64>>>>,
    /// Second moment estimate (for Adam and RMSprop)
    pub second_moment: Option<Vec<Vec<Vec<f64>>>>,
    /// Time step (for Adam bias correction)
    pub time_step: usize,
    /// Beta1 parameter (for Momentum and Adam)
    pub beta1: f64,
    /// Beta2 parameter (for Adam and RMSprop)
    pub beta2: f64,
    /// Epsilon for numerical stability
    pub epsilon: f64,
}

impl OptimizerState {
    /// Create a new optimizer state
    pub fn new(optimizer_type: OptimizerType, learning_rate: f64, _num_layers: usize) -> Self {
        let beta1 = match optimizer_type {
            OptimizerType::Momentum => 0.9,
            OptimizerType::Adam => 0.9,
            _ => 0.0,
        };

        let beta2 = match optimizer_type {
            OptimizerType::RMSprop => 0.999,
            OptimizerType::Adam => 0.999,
            _ => 0.0,
        };

        Self {
            optimizer_type,
            learning_rate,
            velocity: None,
            second_moment: None,
            time_step: 0,
            beta1,
            beta2,
            epsilon: 1e-8,
        }
    }

    /// Initialize optimizer state with layer dimensions
    pub fn initialize(&mut self, layer_dims: &[(usize, usize)]) {
        match self.optimizer_type {
            OptimizerType::Momentum | OptimizerType::Adam => {
                // Initialize velocity
                self.velocity = Some(
                    layer_dims
                        .iter()
                        .map(|&(rows, cols)| vec![vec![0.0; cols]; rows])
                        .collect(),
                );
            }
            _ => {}
        }

        match self.optimizer_type {
            OptimizerType::RMSprop | OptimizerType::Adam => {
                // Initialize second moment
                self.second_moment = Some(
                    layer_dims
                        .iter()
                        .map(|&(rows, cols)| vec![vec![0.0; cols]; rows])
                        .collect(),
                );
            }
            _ => {}
        }
    }

    /// Update parameters using the optimizer
    pub fn update(&mut self, gradients: &[Vec<Vec<f64>>]) -> Vec<Vec<Vec<f64>>> {
        self.time_step += 1;

        match self.optimizer_type {
            OptimizerType::SGD => self.update_sgd(gradients),
            OptimizerType::Momentum => self.update_momentum(gradients),
            OptimizerType::RMSprop => self.update_rmsprop(gradients),
            OptimizerType::Adam => self.update_adam(gradients),
        }
    }

    /// Standard SGD update
    pub(crate) fn update_sgd(&self, gradients: &[Vec<Vec<f64>>]) -> Vec<Vec<Vec<f64>>> {
        gradients
            .iter()
            .map(|layer_grad| {
                layer_grad
                    .iter()
                    .map(|row| row.iter().map(|&g| self.learning_rate * g).collect())
                    .collect()
            })
            .collect()
    }

    /// SGD with Momentum update
    fn update_momentum(&mut self, gradients: &[Vec<Vec<f64>>]) -> Vec<Vec<Vec<f64>>> {
        let velocity = self.velocity.as_mut().unwrap();

        let updates: Vec<Vec<Vec<f64>>> = gradients
            .iter()
            .zip(velocity.iter_mut())
            .map(|(layer_grad, layer_velocity)| {
                layer_grad
                    .iter()
                    .zip(layer_velocity.iter_mut())
                    .map(|(row_grad, row_velocity)| {
                        row_grad
                            .iter()
                            .zip(row_velocity.iter_mut())
                            .map(|(&g, v)| {
                                *v = self.beta1 * *v + (1.0 - self.beta1) * g;
                                self.learning_rate * *v
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect();

        updates
    }

    /// RMSprop update
    fn update_rmsprop(&mut self, gradients: &[Vec<Vec<f64>>]) -> Vec<Vec<Vec<f64>>> {
        let second_moment = self.second_moment.as_mut().unwrap();

        let updates: Vec<Vec<Vec<f64>>> = gradients
            .iter()
            .zip(second_moment.iter_mut())
            .map(|(layer_grad, layer_sm)| {
                layer_grad
                    .iter()
                    .zip(layer_sm.iter_mut())
                    .map(|(row_grad, row_sm)| {
                        row_grad
                            .iter()
                            .zip(row_sm.iter_mut())
                            .map(|(&g, sm)| {
                                *sm = self.beta2 * *sm + (1.0 - self.beta2) * g * g;
                                self.learning_rate * g / (sm.sqrt() + self.epsilon)
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect();

        updates
    }

    /// Adam update (combines Momentum and RMSprop)
    fn update_adam(&mut self, gradients: &[Vec<Vec<f64>>]) -> Vec<Vec<Vec<f64>>> {
        let velocity = self.velocity.as_mut().unwrap();
        let second_moment = self.second_moment.as_mut().unwrap();

        // Bias correction factors
        let bias_correction1 = 1.0 - self.beta1.powi(self.time_step as i32);
        let bias_correction2 = 1.0 - self.beta2.powi(self.time_step as i32);

        let updates: Vec<Vec<Vec<f64>>> = gradients
            .iter()
            .zip(velocity.iter_mut())
            .zip(second_moment.iter_mut())
            .map(|((layer_grad, layer_velocity), layer_sm)| {
                layer_grad
                    .iter()
                    .zip(layer_velocity.iter_mut())
                    .zip(layer_sm.iter_mut())
                    .map(|((row_grad, row_velocity), row_sm)| {
                        row_grad
                            .iter()
                            .zip(row_velocity.iter_mut())
                            .zip(row_sm.iter_mut())
                            .map(|((&g, v), sm)| {
                                // Update biased first moment estimate
                                *v = self.beta1 * *v + (1.0 - self.beta1) * g;

                                // Update biased second moment estimate
                                *sm = self.beta2 * *sm + (1.0 - self.beta2) * g * g;

                                // Bias-corrected estimates
                                let v_corrected = *v / bias_correction1;
                                let sm_corrected = *sm / bias_correction2;

                                // Parameter update
                                self.learning_rate * v_corrected / (sm_corrected.sqrt() + self.epsilon)
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect();

        updates
    }
}

/// Train a neural network using the specified optimizer
pub fn train_network(
    network: &mut NeuralNetwork,
    train_data: &[(Vec<f64>, Vec<f64>)],
    config: &TrainingConfig,
) -> TrainingHistory {
    // Get layer dimensions for optimizer initialization
    let layer_dims: Vec<(usize, usize)> = network
        .layers
        .iter()
        .map(|layer| (layer.weights.len(), layer.weights[0].len()))
        .collect();

    let mut optimizer = OptimizerState::new(config.optimizer, config.learning_rate, network.layers.len());
    optimizer.initialize(&layer_dims);

    let mut losses = Vec::new();
    let mut best_loss = f64::INFINITY;
    let mut patience_counter = 0;

    for epoch in 0..config.epochs {
        let mut epoch_loss = 0.0;
        // batch_count could be used for progress tracking or metrics
        let mut _batch_count = 0;

        // Mini-batch training
        for batch in train_data.chunks(config.batch_size) {
            let mut batch_weight_grads: Vec<Vec<Vec<f64>>> = layer_dims
                .iter()
                .map(|&(rows, cols)| vec![vec![0.0; cols]; rows])
                .collect();
            let mut batch_bias_grads: Vec<Vec<f64>> = layer_dims
                .iter()
                .map(|&(rows, _)| vec![0.0; rows])
                .collect();

            // Accumulate gradients over batch
            for (input, target) in batch {
                let activations = network.forward_with_cache(input);
                let output = activations.last().unwrap();

                epoch_loss += compute_loss(output, target, config.loss_function);

                let gradients = network.backward(&activations, target, config.loss_function);

                for (i, (weight_grad, bias_grad)) in gradients.iter().enumerate() {
                    for (j, row) in weight_grad.iter().enumerate() {
                        for (k, &val) in row.iter().enumerate() {
                            batch_weight_grads[i][j][k] += val;
                        }
                    }

                    for (j, &val) in bias_grad.iter().enumerate() {
                        batch_bias_grads[i][j] += val;
                    }
                }
            }

            // Average gradients
            let batch_size_f = batch.len() as f64;
            for layer_grad in &mut batch_weight_grads {
                for row in layer_grad {
                    for val in row {
                        *val /= batch_size_f;
                    }
                }
            }

            for bias_grad in &mut batch_bias_grads {
                for val in bias_grad {
                    *val /= batch_size_f;
                }
            }

            // Update weights using optimizer
            let weight_updates = optimizer.update(&batch_weight_grads);

            for (i, layer) in network.layers.iter_mut().enumerate() {
                for (j, row) in layer.weights.iter_mut().enumerate() {
                    for (k, weight) in row.iter_mut().enumerate() {
                        *weight -= weight_updates[i][j][k];
                    }
                }

                for (j, bias) in layer.biases.iter_mut().enumerate() {
                    *bias -= config.learning_rate * batch_bias_grads[i][j];
                }
            }

            _batch_count += 1;
        }

        epoch_loss /= train_data.len() as f64;
        losses.push(epoch_loss);

        // Early stopping check
        if config.early_stopping_patience > 0 {
            if epoch_loss < best_loss {
                best_loss = epoch_loss;
                patience_counter = 0;
            } else {
                patience_counter += 1;
                if patience_counter >= config.early_stopping_patience {
                    return TrainingHistory {
                        losses,
                        val_losses: None,
                        epochs_trained: epoch + 1,
                        early_stopped: true,
                    };
                }
            }
        }
    }

    TrainingHistory {
        losses,
        val_losses: None,
        epochs_trained: config.epochs,
        early_stopped: false,
    }
}
