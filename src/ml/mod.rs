//! Machine Learning Foundations Module
//!
//! Provides fundamental machine learning algorithms and neural network primitives.
//!
//! # Modules
//!
//! - `neural_network`: Layer operations, activations, backpropagation
//! - `optimization`: Gradient descent variants (SGD, Adam, RMSprop)
//! - `dim_reduction`: PCA and related techniques
//! - `clustering`: K-means, hierarchical clustering
//! - `regression`: Linear and logistic regression with gradient descent
//! - `classification`: Classification algorithms (placeholder)

use serde::{Deserialize, Serialize};

use crate::engine::*;

// Submodules
pub mod classification;
pub mod clustering;
pub mod dim_reduction;
pub mod neural_network;
pub mod optimization;
pub mod regression;

// Re-exports
pub use clustering::*;
pub use dim_reduction::*;
pub use neural_network::*;
pub use optimization::*;
pub use regression::*;

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

// ============================================================================
// UnifiedML - Tool Implementation
// ============================================================================

pub struct UnifiedML;

impl UnifiedML {
    pub fn new() -> Self {
        Self
    }

    /// Handle clustering operations
    fn handle_clustering(&self, method: &ClusteringMethod, input: &MLInput) -> ToolResult<MLOutput> {
        match method {
            ClusteringMethod::KMeans => {
                // Extract data from input
                let data = self.extract_matrix_data(&input.data)?;

                let k = input.parameters
                    .get("k")
                    .and_then(|v| v.as_u64())
                    .map(|v| v as usize)
                    .unwrap_or(3);

                let max_iterations = input.parameters
                    .get("max_iterations")
                    .and_then(|v| v.as_u64())
                    .map(|v| v as usize);

                let tolerance = input.parameters
                    .get("tolerance")
                    .and_then(|v| v.as_f64());

                let init_method = input.parameters
                    .get("init_method")
                    .and_then(|v| v.as_str());

                let result = clustering::kmeans(&data, k, max_iterations, tolerance, init_method)
                    .map_err(|e| format!("K-means clustering failed: {}", e))?;

                let mut metrics = std::collections::HashMap::new();
                metrics.insert("inertia".to_string(), result.inertia);
                metrics.insert("n_iterations".to_string(), result.n_iterations as f64);
                metrics.insert("converged".to_string(), if result.converged { 1.0 } else { 0.0 });

                Ok(MLOutput {
                    result: serde_json::json!({
                        "labels": result.labels,
                        "centroids": result.centroids
                    }),
                    model: Some(serde_json::json!({
                        "centroids": result.centroids,
                        "k": k
                    })),
                    metrics: Some(metrics),
                    metadata: Some(serde_json::json!({
                        "method": "kmeans",
                        "n_samples": data.len(),
                        "n_features": data.first().map(|v| v.len()).unwrap_or(0)
                    })),
                })
            }

            ClusteringMethod::SilhouetteScore => {
                let data = self.extract_matrix_data(&input.data)?;

                let labels: Vec<usize> = input.parameters
                    .get("labels")
                    .and_then(|v| v.as_array())
                    .map(|arr| arr.iter().filter_map(|v| v.as_u64().map(|n| n as usize)).collect())
                    .unwrap_or_default();

                if labels.is_empty() {
                    return Err("labels parameter required for silhouette score".to_string());
                }

                let score = clustering::silhouette_score(&data, &labels)
                    .map_err(|e| format!("Silhouette score failed: {}", e))?;

                let mut metrics = std::collections::HashMap::new();
                metrics.insert("silhouette_score".to_string(), score);

                Ok(MLOutput {
                    result: serde_json::json!({ "silhouette_score": score }),
                    model: None,
                    metrics: Some(metrics),
                    metadata: Some(serde_json::json!({
                        "method": "silhouette_score",
                        "n_samples": data.len()
                    })),
                })
            }

            ClusteringMethod::DBSCAN => {
                Err("DBSCAN clustering not yet implemented".to_string())
            }

            ClusteringMethod::Hierarchical => {
                Err("Hierarchical clustering not yet implemented".to_string())
            }

            ClusteringMethod::GaussianMixture => {
                Err("Gaussian Mixture Models not yet implemented".to_string())
            }
        }
    }

    /// Handle neural network operations
    fn handle_neural_network(&self, op: &NeuralNetworkOp, input: &MLInput) -> ToolResult<MLOutput> {
        match op {
            NeuralNetworkOp::CreateLayer => {
                let input_size = input.parameters
                    .get("input_size")
                    .and_then(|v| v.as_u64())
                    .map(|v| v as usize)
                    .unwrap_or(1);

                let output_size = input.parameters
                    .get("output_size")
                    .and_then(|v| v.as_u64())
                    .map(|v| v as usize)
                    .unwrap_or(1);

                let activation = input.parameters
                    .get("activation")
                    .and_then(|v| v.as_str())
                    .map(|s| match s.to_lowercase().as_str() {
                        "sigmoid" => ActivationFunction::Sigmoid,
                        "tanh" => ActivationFunction::Tanh,
                        "relu" => ActivationFunction::ReLU,
                        "leaky_relu" => ActivationFunction::LeakyReLU,
                        "softmax" => ActivationFunction::Softmax,
                        _ => ActivationFunction::Linear,
                    })
                    .unwrap_or(ActivationFunction::ReLU);

                let layer = neural_network::DenseLayer::new(input_size, output_size, activation);

                Ok(MLOutput {
                    result: serde_json::json!({
                        "input_size": input_size,
                        "output_size": output_size,
                        "activation": format!("{:?}", activation),
                        "weights_shape": [layer.weights.len(), layer.weights.first().map(|v| v.len()).unwrap_or(0)]
                    }),
                    model: Some(serde_json::json!({
                        "weights": layer.weights,
                        "biases": layer.biases,
                        "activation": format!("{:?}", activation)
                    })),
                    metrics: None,
                    metadata: Some(serde_json::json!({
                        "operation": "create_layer",
                        "n_parameters": layer.weights.len() * layer.weights.first().map(|v| v.len()).unwrap_or(0) + layer.biases.len()
                    })),
                })
            }

            NeuralNetworkOp::Train => {
                // Extract network architecture
                let layer_sizes: Vec<usize> = input.parameters
                    .get("layer_sizes")
                    .and_then(|v| v.as_array())
                    .map(|arr| arr.iter().filter_map(|v| v.as_u64().map(|n| n as usize)).collect())
                    .unwrap_or_else(|| vec![2, 4, 1]);

                // Extract training data
                let train_data = self.extract_training_data(&input.data)?;

                // Build network
                let activations: Vec<ActivationFunction> = input.parameters
                    .get("activations")
                    .and_then(|v| v.as_array())
                    .map(|arr| arr.iter().map(|v| {
                        match v.as_str().unwrap_or("relu").to_lowercase().as_str() {
                            "sigmoid" => ActivationFunction::Sigmoid,
                            "tanh" => ActivationFunction::Tanh,
                            "relu" => ActivationFunction::ReLU,
                            "leaky_relu" => ActivationFunction::LeakyReLU,
                            "softmax" => ActivationFunction::Softmax,
                            _ => ActivationFunction::ReLU,
                        }
                    }).collect())
                    .unwrap_or_else(|| vec![ActivationFunction::ReLU; layer_sizes.len()]);

                // Build layer specification for NeuralNetwork::new()
                let layer_spec: Vec<(usize, ActivationFunction)> = layer_sizes.iter()
                    .enumerate()
                    .map(|(i, &size)| {
                        let activation = activations.get(i).copied().unwrap_or(ActivationFunction::ReLU);
                        (size, activation)
                    })
                    .collect();

                let mut network = neural_network::NeuralNetwork::new(&layer_spec);

                // Training config
                let config = TrainingConfig {
                    learning_rate: input.parameters.get("learning_rate").and_then(|v| v.as_f64()).unwrap_or(0.01),
                    epochs: input.parameters.get("epochs").and_then(|v| v.as_u64()).map(|v| v as usize).unwrap_or(100),
                    batch_size: input.parameters.get("batch_size").and_then(|v| v.as_u64()).map(|v| v as usize).unwrap_or(32),
                    optimizer: match input.parameters.get("optimizer").and_then(|v| v.as_str()).unwrap_or("adam") {
                        "sgd" => OptimizerType::SGD,
                        "momentum" => OptimizerType::Momentum,
                        "rmsprop" => OptimizerType::RMSprop,
                        _ => OptimizerType::Adam,
                    },
                    loss_function: match input.parameters.get("loss").and_then(|v| v.as_str()).unwrap_or("mse") {
                        "mae" => LossFunction::MAE,
                        "cross_entropy" => LossFunction::CrossEntropy,
                        "binary_cross_entropy" => LossFunction::BinaryCrossEntropy,
                        _ => LossFunction::MSE,
                    },
                    early_stopping_patience: input.parameters.get("early_stopping_patience").and_then(|v| v.as_u64()).map(|v| v as usize).unwrap_or(0),
                };

                let history = optimization::train_network(&mut network, &train_data, &config);

                let mut metrics = std::collections::HashMap::new();
                metrics.insert("final_loss".to_string(), *history.losses.last().unwrap_or(&f64::NAN));
                metrics.insert("epochs_trained".to_string(), history.epochs_trained as f64);
                metrics.insert("early_stopped".to_string(), if history.early_stopped { 1.0 } else { 0.0 });

                Ok(MLOutput {
                    result: serde_json::json!({
                        "losses": history.losses,
                        "epochs_trained": history.epochs_trained,
                        "early_stopped": history.early_stopped
                    }),
                    model: Some(serde_json::json!({
                        "layers": network.layers.iter().map(|l| serde_json::json!({
                            "weights": l.weights,
                            "biases": l.biases
                        })).collect::<Vec<_>>()
                    })),
                    metrics: Some(metrics),
                    metadata: Some(serde_json::json!({
                        "operation": "train",
                        "architecture": layer_sizes
                    })),
                })
            }

            NeuralNetworkOp::Forward => {
                Err("Forward pass requires a trained model - use Train first".to_string())
            }

            NeuralNetworkOp::Backward => {
                Err("Backward pass is handled internally during training".to_string())
            }

            NeuralNetworkOp::Predict => {
                Err("Predict requires a trained model - use Train first".to_string())
            }
        }
    }

    /// Handle regression operations
    fn handle_regression(&self, method: &RegressionMethod, input: &MLInput) -> ToolResult<MLOutput> {
        // Extract X and y from input
        let (x, y) = self.extract_regression_data(&input.data)?;

        match method {
            RegressionMethod::Linear => {
                let model = regression::linear_regression(&x, &y)
                    .map_err(|e| format!("Linear regression failed: {}", e))?;

                let mut metrics = std::collections::HashMap::new();
                if let Some(r2) = model.r_squared {
                    metrics.insert("r_squared".to_string(), r2);
                }

                Ok(MLOutput {
                    result: serde_json::json!({
                        "coefficients": model.coefficients,
                        "intercept": model.intercept,
                        "equation": format!("y = {:.4} + {}", model.intercept,
                            model.coefficients.iter().enumerate()
                                .map(|(i, c)| format!("{:.4}*x{}", c, i))
                                .collect::<Vec<_>>().join(" + "))
                    }),
                    model: Some(serde_json::json!({
                        "type": "linear_regression",
                        "coefficients": model.coefficients,
                        "intercept": model.intercept
                    })),
                    metrics: Some(metrics),
                    metadata: Some(serde_json::json!({
                        "method": "linear_regression",
                        "n_samples": x.len(),
                        "n_features": x.first().map(|v| v.len()).unwrap_or(0)
                    })),
                })
            }

            RegressionMethod::Ridge => {
                let alpha = input.parameters
                    .get("alpha")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(1.0);

                let model = regression::ridge_regression(&x, &y, alpha)
                    .map_err(|e| format!("Ridge regression failed: {}", e))?;

                let mut metrics = std::collections::HashMap::new();
                if let Some(r2) = model.r_squared {
                    metrics.insert("r_squared".to_string(), r2);
                }

                Ok(MLOutput {
                    result: serde_json::json!({
                        "coefficients": model.coefficients,
                        "intercept": model.intercept
                    }),
                    model: Some(serde_json::json!({
                        "type": "ridge_regression",
                        "coefficients": model.coefficients,
                        "intercept": model.intercept,
                        "alpha": alpha
                    })),
                    metrics: Some(metrics),
                    metadata: Some(serde_json::json!({
                        "method": "ridge_regression",
                        "alpha": alpha
                    })),
                })
            }

            RegressionMethod::Logistic => {
                let learning_rate = input.parameters
                    .get("learning_rate")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.1);
                let max_iterations = input.parameters
                    .get("max_iterations")
                    .and_then(|v| v.as_u64())
                    .map(|v| v as usize)
                    .unwrap_or(1000);
                let tolerance = input.parameters
                    .get("tolerance")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(1e-4);

                let model = regression::logistic_regression(&x, &y, learning_rate, max_iterations, tolerance)
                    .map_err(|e| format!("Logistic regression failed: {}", e))?;

                let mut metrics = std::collections::HashMap::new();
                if let Some(acc) = model.accuracy {
                    metrics.insert("accuracy".to_string(), acc);
                }

                Ok(MLOutput {
                    result: serde_json::json!({
                        "coefficients": model.coefficients,
                        "intercept": model.intercept,
                        "accuracy": model.accuracy
                    }),
                    model: Some(serde_json::json!({
                        "type": "logistic_regression",
                        "coefficients": model.coefficients,
                        "intercept": model.intercept
                    })),
                    metrics: Some(metrics),
                    metadata: Some(serde_json::json!({
                        "method": "logistic_regression"
                    })),
                })
            }

            RegressionMethod::Lasso => {
                Err("Lasso regression not yet implemented".to_string())
            }

            RegressionMethod::ElasticNet => {
                Err("Elastic Net not yet implemented".to_string())
            }

            RegressionMethod::PolynomialRegression => {
                Err("Polynomial regression not yet implemented".to_string())
            }
        }
    }

    /// Handle dimensionality reduction operations
    fn handle_dim_reduction(&self, method: &DimReductionMethod, input: &MLInput) -> ToolResult<MLOutput> {
        let data = self.extract_matrix_data(&input.data)?;

        match method {
            DimReductionMethod::PCA => {
                let n_components = input.parameters
                    .get("n_components")
                    .and_then(|v| v.as_u64())
                    .map(|v| v as usize);

                let result = dim_reduction::pca(&data, n_components)
                    .map_err(|e| format!("PCA failed: {}", e))?;

                let mut metrics = std::collections::HashMap::new();
                for (i, &var) in result.explained_variance_ratio.iter().enumerate() {
                    metrics.insert(format!("explained_variance_ratio_{}", i), var);
                }
                metrics.insert("total_explained_variance".to_string(),
                    result.explained_variance_ratio.iter().sum::<f64>());

                Ok(MLOutput {
                    result: serde_json::json!({
                        "components": result.components,
                        "explained_variance": result.explained_variance,
                        "explained_variance_ratio": result.explained_variance_ratio,
                        "singular_values": result.singular_values
                    }),
                    model: Some(serde_json::json!({
                        "type": "pca",
                        "components": result.components,
                        "mean": result.mean
                    })),
                    metrics: Some(metrics),
                    metadata: Some(serde_json::json!({
                        "method": "pca",
                        "n_components": result.components.len(),
                        "n_features": result.mean.len()
                    })),
                })
            }

            DimReductionMethod::Transform => {
                // Transform data using existing PCA model
                let pca_model = input.parameters
                    .get("model")
                    .ok_or("model parameter required for transform")?;

                let components: Vec<Vec<f64>> = serde_json::from_value(
                    pca_model.get("components").cloned().unwrap_or_default()
                ).map_err(|e| format!("Invalid components: {}", e))?;

                let mean: Vec<f64> = serde_json::from_value(
                    pca_model.get("mean").cloned().unwrap_or_default()
                ).map_err(|e| format!("Invalid mean: {}", e))?;

                let pca_result = dim_reduction::PCAResult {
                    components,
                    explained_variance: vec![],
                    explained_variance_ratio: vec![],
                    mean,
                    singular_values: vec![],
                };

                let transformed = dim_reduction::transform(&data, &pca_result)
                    .map_err(|e| format!("Transform failed: {}", e))?;

                Ok(MLOutput {
                    result: serde_json::json!({ "transformed": transformed }),
                    model: None,
                    metrics: None,
                    metadata: Some(serde_json::json!({
                        "method": "pca_transform",
                        "n_samples": transformed.len()
                    })),
                })
            }

            DimReductionMethod::TSNE => {
                Err("t-SNE not yet implemented".to_string())
            }

            DimReductionMethod::UMAP => {
                Err("UMAP not yet implemented".to_string())
            }

            DimReductionMethod::LDA => {
                Err("LDA not yet implemented".to_string())
            }
        }
    }

    /// Handle classification operations
    fn handle_classification(&self, method: &ClassificationMethod, _input: &MLInput) -> ToolResult<MLOutput> {
        match method {
            ClassificationMethod::DecisionTree => {
                Err("Decision Tree not yet implemented".to_string())
            }
            ClassificationMethod::RandomForest => {
                Err("Random Forest not yet implemented".to_string())
            }
            ClassificationMethod::SVM => {
                Err("SVM not yet implemented".to_string())
            }
            ClassificationMethod::NaiveBayes => {
                Err("Naive Bayes not yet implemented".to_string())
            }
            ClassificationMethod::KNN => {
                Err("K-Nearest Neighbors not yet implemented".to_string())
            }
        }
    }

    // ============================================================================
    // Helper functions
    // ============================================================================

    /// Extract matrix data from JSON Value
    fn extract_matrix_data(&self, data: &serde_json::Value) -> ToolResult<Vec<Vec<f64>>> {
        if let Some(arr) = data.as_array() {
            let mut result = Vec::new();
            for row in arr {
                if let Some(row_arr) = row.as_array() {
                    let row_data: Vec<f64> = row_arr
                        .iter()
                        .filter_map(|v| v.as_f64())
                        .collect();
                    if !row_data.is_empty() {
                        result.push(row_data);
                    }
                }
            }
            if result.is_empty() {
                Err("Data must be a non-empty 2D array".to_string())
            } else {
                Ok(result)
            }
        } else {
            Err("Data must be a 2D array [[f64, ...], ...]".to_string())
        }
    }

    /// Extract regression data (X, y) from JSON Value
    fn extract_regression_data(&self, data: &serde_json::Value) -> ToolResult<(Vec<Vec<f64>>, Vec<f64>)> {
        // Try to get X and y from object format
        if let Some(obj) = data.as_object() {
            if let (Some(x_val), Some(y_val)) = (obj.get("X").or(obj.get("x")), obj.get("y").or(obj.get("Y"))) {
                let x = self.extract_matrix_data(x_val)?;
                let y: Vec<f64> = y_val.as_array()
                    .ok_or("y must be an array")?
                    .iter()
                    .filter_map(|v| v.as_f64())
                    .collect();
                return Ok((x, y));
            }
        }

        // Fallback: assume data is matrix where last column is y
        let matrix = self.extract_matrix_data(data)?;
        let x: Vec<Vec<f64>> = matrix.iter()
            .map(|row| row[..row.len()-1].to_vec())
            .collect();
        let y: Vec<f64> = matrix.iter()
            .map(|row| *row.last().unwrap_or(&0.0))
            .collect();

        Ok((x, y))
    }

    /// Extract training data for neural networks
    fn extract_training_data(&self, data: &serde_json::Value) -> ToolResult<Vec<(Vec<f64>, Vec<f64>)>> {
        if let Some(obj) = data.as_object() {
            if let (Some(x_val), Some(y_val)) = (obj.get("X").or(obj.get("x")), obj.get("y").or(obj.get("Y"))) {
                let x = self.extract_matrix_data(x_val)?;
                let y = self.extract_matrix_data(y_val)?;

                if x.len() != y.len() {
                    return Err("X and y must have same number of samples".to_string());
                }

                return Ok(x.into_iter().zip(y.into_iter()).collect());
            }
        }

        Err("Training data must be {X: [[...]], y: [[...]]}".to_string())
    }
}

impl Default for UnifiedML {
    fn default() -> Self {
        Self::new()
    }
}

impl MachineLearning for UnifiedML {
    fn ml(&self, input: &MLInput) -> ToolResult<MLOutput> {
        match &input.operation {
            MLOp::Clustering(method) => self.handle_clustering(method, input),
            MLOp::NeuralNetwork(op) => self.handle_neural_network(op, input),
            MLOp::Regression(method) => self.handle_regression(method, input),
            MLOp::DimReduction(method) => self.handle_dim_reduction(method, input),
            MLOp::Classification(method) => self.handle_classification(method, input),
        }
    }
}

// Test modules
#[cfg(test)]
#[path = "../../tests/unit/ml/regression/regression_tests.rs"]
mod regression_tests;

#[cfg(test)]
#[path = "../../tests/unit/ml/clustering/clustering_tests.rs"]
mod clustering_tests;

#[cfg(test)]
#[path = "../../tests/unit/ml/neural_network/neural_network_tests.rs"]
mod neural_network_tests;

#[cfg(test)]
#[path = "../../tests/unit/ml/optimization/optimization_tests.rs"]
mod optimization_tests;

#[cfg(test)]
#[path = "../../tests/unit/ml/dim_reduction/dim_reduction_tests.rs"]
mod dimensionality_reduction_tests;
