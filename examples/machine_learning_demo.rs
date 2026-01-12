//! Machine Learning Demo
//!
//! This example demonstrates the machine learning capabilities:
//! - Neural networks with backpropagation
//! - K-means clustering
//! - PCA (Principal Component Analysis)
//! - Linear and logistic regression
//!
//! Run with: cargo run --release --example machine_learning_demo

use computational_engine::ml::*;

fn main() {
    println!("=== Machine Learning Demo ===\n");

    // Example 1: Neural Network for XOR
    neural_network_xor();

    // Example 2: Linear Regression
    linear_regression_demo();

    // Example 3: Logistic Regression
    logistic_regression_demo();

    // Example 4: K-means Clustering
    kmeans_clustering_demo();

    // Example 5: PCA
    pca_demo();
}

/// Example 1: Neural Network learning XOR function
///
/// XOR is not linearly separable, so we need a hidden layer.
/// Network: 2 inputs → 4 hidden (ReLU) → 1 output (Sigmoid)
fn neural_network_xor() {
    println!("--- Example 1: Neural Network (XOR) ---");
    println!("Training a neural network to learn the XOR function");
    println!();

    // XOR dataset: [x1, x2] → y
    let train_data = vec![
        (vec![0.0, 0.0], vec![0.0]),
        (vec![0.0, 1.0], vec![1.0]),
        (vec![1.0, 0.0], vec![1.0]),
        (vec![1.0, 1.0], vec![0.0]),
    ];

    // Create network: 2 inputs → 4 hidden → 1 output
    let mut network = NeuralNetwork::new(&[
        (2, ActivationFunction::Linear), // Input layer size
        (4, ActivationFunction::ReLU),   // Hidden layer
        (1, ActivationFunction::Sigmoid), // Output layer
    ]);

    // Training configuration
    let config = TrainingConfig {
        learning_rate: 0.1,
        epochs: 1000,
        batch_size: 4,
        optimizer: OptimizerType::Adam,
        loss_function: LossFunction::BinaryCrossEntropy,
        early_stopping_patience: 100,
    };

    println!("Configuration:");
    println!("  Epochs: {}", config.epochs);
    println!("  Learning rate: {}", config.learning_rate);
    println!("  Optimizer: {:?}", config.optimizer);
    println!();

    // Train the network
    let history = train_network(&mut network, &train_data, &config);

    println!("Training complete:");
    println!("  Epochs trained: {}", history.epochs_trained);
    println!("  Final loss: {:.6}", history.losses.last().unwrap());
    println!("  Early stopped: {}", history.early_stopped);
    println!();

    // Test the network
    println!("Predictions:");
    for (input, expected) in &train_data {
        let output = network.predict(input);
        let predicted = if output[0] > 0.5 { 1.0 } else { 0.0 };
        println!(
            "  Input: [{}, {}] → Predicted: {:.4} ({}) | Expected: {}",
            input[0], input[1], output[0], predicted, expected[0]
        );
    }
    println!();
}

/// Example 2: Linear Regression
///
/// Fit a line to data: y = 2.5x + 3.0
fn linear_regression_demo() {
    println!("--- Example 2: Linear Regression ---");
    println!("Fitting a line to data: y = 2.5x + 3.0");
    println!();

    // Generate training data with some noise
    let x: Vec<Vec<f64>> = (1..=10).map(|i| vec![i as f64]).collect();
    let y: Vec<f64> = x
        .iter()
        .map(|xi| 2.5 * xi[0] + 3.0 + (rand::random::<f64>() - 0.5) * 0.5)
        .collect();

    // Fit model
    let model = linear_regression(&x, &y).unwrap();

    println!("Model:");
    println!("  Coefficient: {:.4}", model.coefficients[0]);
    println!("  Intercept: {:.4}", model.intercept);
    println!("  R² score: {:.4}", model.r_squared.unwrap());
    println!();

    // Make predictions
    let x_test = vec![vec![11.0], vec![12.0], vec![15.0]];
    let predictions = linear_regression_predict(&model, &x_test);

    println!("Predictions:");
    for (x_val, pred) in x_test.iter().zip(predictions.iter()) {
        let expected = 2.5 * x_val[0] + 3.0;
        println!(
            "  x = {:.0} → Predicted: {:.2} | Expected: {:.2}",
            x_val[0], pred, expected
        );
    }
    println!();
}

/// Example 3: Logistic Regression
///
/// Binary classification: classify points based on a threshold
fn logistic_regression_demo() {
    println!("--- Example 3: Logistic Regression ---");
    println!("Binary classification: x < 5 → class 0, x >= 5 → class 1");
    println!();

    // Generate training data
    let x: Vec<Vec<f64>> = (1..=10).map(|i| vec![i as f64]).collect();
    let y: Vec<f64> = x.iter().map(|xi| if xi[0] < 5.0 { 0.0 } else { 1.0 }).collect();

    // Fit model
    let model = logistic_regression(&x, &y, 0.1, 1000, 1e-6).unwrap();

    println!("Model:");
    println!("  Coefficient: {:.4}", model.coefficients[0]);
    println!("  Intercept: {:.4}", model.intercept);
    println!("  Training accuracy: {:.2}%", model.accuracy.unwrap() * 100.0);
    println!();

    // Make predictions
    let x_test = vec![vec![2.0], vec![5.0], vec![8.0]];
    let probs = logistic_regression_predict_proba(&model, &x_test);
    let predictions = logistic_regression_predict(&model, &x_test);

    println!("Predictions:");
    for (x_val, (prob, pred)) in x_test.iter().zip(probs.iter().zip(predictions.iter())) {
        println!(
            "  x = {:.0} → P(class=1): {:.4} | Predicted class: {}",
            x_val[0], prob, pred
        );
    }
    println!();
}

/// Example 4: K-means Clustering
///
/// Cluster 2D points into groups
fn kmeans_clustering_demo() {
    println!("--- Example 4: K-means Clustering ---");
    println!("Clustering 2D points into 3 groups");
    println!();

    // Generate data with 3 clusters
    let data = vec![
        // Cluster 1 (around [1, 1])
        vec![1.0, 1.0],
        vec![1.5, 2.0],
        vec![1.8, 1.5],
        vec![2.0, 2.0],
        // Cluster 2 (around [5, 5])
        vec![5.0, 5.0],
        vec![5.5, 5.5],
        vec![5.2, 4.8],
        vec![4.8, 5.2],
        // Cluster 3 (around [9, 1])
        vec![9.0, 1.0],
        vec![9.5, 1.5],
        vec![9.2, 0.8],
        vec![8.8, 1.2],
    ];

    // Perform clustering
    let result = kmeans(&data, 3, Some(100), Some(1e-4), Some("kmeans++")).unwrap();

    println!("Clustering result:");
    println!("  Converged: {}", result.converged);
    println!("  Iterations: {}", result.n_iterations);
    println!("  Inertia: {:.4}", result.inertia);
    println!();

    println!("Cluster centers:");
    for (i, centroid) in result.centroids.iter().enumerate() {
        println!("  Cluster {}: [{:.2}, {:.2}]", i, centroid[0], centroid[1]);
    }
    println!();

    println!("Point assignments:");
    for (i, (point, label)) in data.iter().zip(result.labels.iter()).enumerate() {
        println!(
            "  Point {}: [{:.1}, {:.1}] → Cluster {}",
            i, point[0], point[1], label
        );
    }
    println!();

    // Compute silhouette score
    let score = silhouette_score(&data, &result.labels).unwrap();
    println!("Silhouette score: {:.4} (closer to 1.0 is better)", score);
    println!();
}

/// Example 5: PCA (Principal Component Analysis)
///
/// Reduce dimensionality while preserving variance
fn pca_demo() {
    println!("--- Example 5: PCA (Dimensionality Reduction) ---");
    println!("Reducing 3D data to 2D");
    println!();

    // Generate 3D data (correlated features)
    let data = vec![
        vec![2.5, 2.4, 3.0],
        vec![0.5, 0.7, 1.0],
        vec![2.2, 2.9, 3.5],
        vec![1.9, 2.2, 2.8],
        vec![3.1, 3.0, 4.0],
        vec![2.3, 2.7, 3.2],
        vec![2.0, 1.6, 2.5],
        vec![1.0, 1.1, 1.5],
        vec![1.5, 1.6, 2.0],
        vec![1.1, 0.9, 1.3],
    ];

    // Perform PCA (keep all components initially)
    let result = pca(&data, None).unwrap();

    println!("PCA Analysis:");
    println!("  Original dimensions: 3");
    println!("  Components computed: {}", result.components.len());
    println!();

    println!("Explained variance:");
    for (i, (variance, ratio)) in result
        .explained_variance
        .iter()
        .zip(result.explained_variance_ratio.iter())
        .enumerate()
    {
        println!(
            "  PC{}: Variance = {:.4}, Ratio = {:.2}%",
            i + 1,
            variance,
            ratio * 100.0
        );
    }
    println!();

    // Cumulative variance
    let cumulative = cumulative_explained_variance(&result);
    println!("Cumulative explained variance:");
    for (i, cum) in cumulative.iter().enumerate() {
        println!("  First {} components: {:.2}%", i + 1, cum * 100.0);
    }
    println!();

    // Determine components needed for 95% variance
    let n_components_95 = components_for_variance(&result, 0.95);
    println!(
        "Components needed for 95% variance: {}",
        n_components_95
    );
    println!();

    // Reduce to 2D
    let result_2d = pca(&data, Some(2)).unwrap();
    let transformed = transform(&data, &result_2d).unwrap();

    println!("Transformed data (first 5 samples):");
    for (i, sample) in transformed.iter().take(5).enumerate() {
        println!(
            "  Sample {}: [{:.3}, {:.3}, {:.3}] → [{:.3}, {:.3}]",
            i, data[i][0], data[i][1], data[i][2], sample[0], sample[1]
        );
    }
    println!();

    // Reconstruct data
    let reconstructed = inverse_transform(&transformed, &result_2d).unwrap();

    // Compute reconstruction error
    let mut total_error = 0.0;
    for (original, recon) in data.iter().zip(reconstructed.iter()) {
        for (o, r) in original.iter().zip(recon.iter()) {
            total_error += (o - r).powi(2);
        }
    }
    let mse = total_error / (data.len() * data[0].len()) as f64;

    println!("Reconstruction error (MSE): {:.6}", mse);
    println!();
}
