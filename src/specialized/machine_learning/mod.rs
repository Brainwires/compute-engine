//! Machine Learning Foundations Module
//!
//! Provides fundamental machine learning algorithms and neural network primitives.
//!
//! # Modules
//!
//! - `neural_network`: Layer operations, activations, backpropagation
//! - `optimization`: Gradient descent variants (SGD, Adam, RMSprop)
//! - `dimensionality_reduction`: PCA and related techniques
//! - `clustering`: K-means, hierarchical clustering
//! - `regression`: Linear and logistic regression with gradient descent

pub mod neural_network;
pub mod optimization;
pub mod dimensionality_reduction;
pub mod clustering;
pub mod regression;

pub use neural_network::*;
pub use optimization::*;
pub use dimensionality_reduction::*;
pub use clustering::*;
pub use regression::*;

use serde::{Deserialize, Serialize};

/// Activation function types for neural networks
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum ActivationFunction {
    /// Linear activation: f(x) = x
    Linear,
    /// Sigmoid activation: f(x) = 1 / (1 + e^(-x))
    Sigmoid,
    /// Hyperbolic tangent: f(x) = tanh(x)
    Tanh,
    /// Rectified Linear Unit: f(x) = max(0, x)
    ReLU,
    /// Leaky ReLU: f(x) = max(αx, x) where α = 0.01
    LeakyReLU,
    /// Softmax activation (for output layers)
    Softmax,
}

/// Loss function types
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum LossFunction {
    /// Mean Squared Error: (1/n) Σ(y - ŷ)²
    MSE,
    /// Mean Absolute Error: (1/n) Σ|y - ŷ|
    MAE,
    /// Cross-Entropy Loss (for classification)
    CrossEntropy,
    /// Binary Cross-Entropy
    BinaryCrossEntropy,
}

/// Optimization algorithm types
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum OptimizerType {
    /// Stochastic Gradient Descent
    SGD,
    /// SGD with momentum
    Momentum,
    /// RMSprop
    RMSprop,
    /// Adam (Adaptive Moment Estimation)
    Adam,
}

/// Training configuration for ML models
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrainingConfig {
    /// Learning rate
    pub learning_rate: f64,
    /// Number of epochs
    pub epochs: usize,
    /// Batch size
    pub batch_size: usize,
    /// Optimizer type
    pub optimizer: OptimizerType,
    /// Loss function
    pub loss_function: LossFunction,
    /// Early stopping patience (0 = disabled)
    pub early_stopping_patience: usize,
}

impl Default for TrainingConfig {
    fn default() -> Self {
        Self {
            learning_rate: 0.01,
            epochs: 100,
            batch_size: 32,
            optimizer: OptimizerType::Adam,
            loss_function: LossFunction::MSE,
            early_stopping_patience: 0,
        }
    }
}

/// Training history and metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrainingHistory {
    /// Loss values per epoch
    pub losses: Vec<f64>,
    /// Validation losses per epoch (if validation data provided)
    pub val_losses: Option<Vec<f64>>,
    /// Number of epochs trained
    pub epochs_trained: usize,
    /// Whether early stopping was triggered
    pub early_stopped: bool,
}

// Test modules - part of this module, can access private functions
#[cfg(test)]
#[path = "../../../tests/unit/specialized/machine_learning/regression_tests.rs"]
mod regression_tests;

#[cfg(test)]
#[path = "../../../tests/unit/specialized/machine_learning/clustering_tests.rs"]
mod clustering_tests;

#[cfg(test)]
#[path = "../../../tests/unit/specialized/machine_learning/neural_network_tests.rs"]
mod neural_network_tests;

#[cfg(test)]
#[path = "../../../tests/unit/specialized/machine_learning/optimization_tests.rs"]
mod optimization_tests;

#[cfg(test)]
#[path = "../../../tests/unit/specialized/machine_learning/dimensionality_reduction_tests.rs"]
mod dimensionality_reduction_tests;
