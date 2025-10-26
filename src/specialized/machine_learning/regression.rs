//! Regression Models
//!
//! Supervised learning for continuous and binary outcomes:
//! - Linear regression (ordinary least squares)
//! - Logistic regression (binary classification)
//! - Ridge regression (L2 regularization)
//! - Gradient descent training

use nalgebra::{DMatrix, DVector};

/// Linear regression model
#[derive(Debug, Clone)]
pub struct LinearRegression {
    /// Coefficients (weights) for each feature
    pub coefficients: Vec<f64>,
    /// Intercept (bias) term
    pub intercept: f64,
    /// R² score (coefficient of determination)
    pub r_squared: Option<f64>,
}

/// Logistic regression model
#[derive(Debug, Clone)]
pub struct LogisticRegression {
    /// Coefficients (weights) for each feature
    pub coefficients: Vec<f64>,
    /// Intercept (bias) term
    pub intercept: f64,
    /// Training accuracy
    pub accuracy: Option<f64>,
}

/// Fit linear regression using ordinary least squares
///
/// Solves: β = (X^T X)^(-1) X^T y
///
/// # Arguments
/// * `x` - Feature matrix (samples as rows, features as columns)
/// * `y` - Target values
///
/// # Example
/// ```
/// use computational_engine::specialized::machine_learning::regression::linear_regression;
///
/// let x = vec![
///     vec![1.0],
///     vec![2.0],
///     vec![3.0],
///     vec![4.0],
/// ];
/// let y = vec![2.0, 4.0, 6.0, 8.0];
///
/// let model = linear_regression(&x, &y).unwrap();
/// assert!((model.coefficients[0] - 2.0).abs() < 0.01);
/// assert!(model.intercept.abs() < 0.01);
/// ```
pub fn linear_regression(x: &[Vec<f64>], y: &[f64]) -> Result<LinearRegression, String> {
    if x.is_empty() || y.is_empty() {
        return Err("Data cannot be empty".to_string());
    }

    if x.len() != y.len() {
        return Err("X and y must have same number of samples".to_string());
    }

    let n_samples = x.len();
    let n_features = x[0].len();

    // Add intercept column (column of 1s)
    let mut x_with_intercept = Vec::new();
    for row in x {
        let mut new_row = vec![1.0];
        new_row.extend(row);
        x_with_intercept.push(new_row);
    }

    // Convert to matrices
    let mut x_flat = Vec::new();
    for row in &x_with_intercept {
        x_flat.extend(row);
    }
    let x_matrix = DMatrix::from_row_slice(n_samples, n_features + 1, &x_flat);

    let y_vector = DVector::from_vec(y.to_vec());

    // Solve normal equation: (X^T X) β = X^T y
    let xt_x = x_matrix.transpose() * &x_matrix;
    let xt_y = x_matrix.transpose() * y_vector;

    // Use SVD for numerical stability
    let beta = match xt_x.svd(true, true).solve(&xt_y, 1e-10) {
        Ok(b) => b,
        Err(_) => return Err("Failed to solve normal equation (singular matrix)".to_string()),
    };

    let intercept = beta[0];
    let coefficients: Vec<f64> = beta.iter().skip(1).copied().collect();

    // Compute R² score
    let y_pred: Vec<f64> = x
        .iter()
        .map(|sample| {
            intercept
                + sample
                    .iter()
                    .zip(&coefficients)
                    .map(|(xi, &coef)| xi * coef)
                    .sum::<f64>()
        })
        .collect();

    let r_squared = compute_r_squared(y, &y_pred);

    Ok(LinearRegression {
        coefficients,
        intercept,
        r_squared: Some(r_squared),
    })
}

/// Predict using linear regression model
pub fn linear_regression_predict(model: &LinearRegression, x: &[Vec<f64>]) -> Vec<f64> {
    x.iter()
        .map(|sample| {
            model.intercept
                + sample
                    .iter()
                    .zip(&model.coefficients)
                    .map(|(xi, &coef)| xi * coef)
                    .sum::<f64>()
        })
        .collect()
}

/// Fit ridge regression (linear regression with L2 regularization)
///
/// Solves: β = (X^T X + λI)^(-1) X^T y
///
/// # Arguments
/// * `x` - Feature matrix
/// * `y` - Target values
/// * `alpha` - Regularization strength (λ)
pub fn ridge_regression(
    x: &[Vec<f64>],
    y: &[f64],
    alpha: f64,
) -> Result<LinearRegression, String> {
    if x.is_empty() || y.is_empty() {
        return Err("Data cannot be empty".to_string());
    }

    if x.len() != y.len() {
        return Err("X and y must have same number of samples".to_string());
    }

    let n_samples = x.len();
    let n_features = x[0].len();

    // Add intercept column
    let mut x_with_intercept = Vec::new();
    for row in x {
        let mut new_row = vec![1.0];
        new_row.extend(row);
        x_with_intercept.push(new_row);
    }

    // Convert to matrices
    let mut x_flat = Vec::new();
    for row in &x_with_intercept {
        x_flat.extend(row);
    }
    let x_matrix = DMatrix::from_row_slice(n_samples, n_features + 1, &x_flat);

    let y_vector = DVector::from_vec(y.to_vec());

    // Add regularization: (X^T X + αI) β = X^T y
    let xt_x = x_matrix.transpose() * &x_matrix;
    let regularization = DMatrix::identity(n_features + 1, n_features + 1) * alpha;
    let regularized = xt_x + regularization;

    let xt_y = x_matrix.transpose() * y_vector;

    let beta = match regularized.svd(true, true).solve(&xt_y, 1e-10) {
        Ok(b) => b,
        Err(_) => return Err("Failed to solve ridge regression".to_string()),
    };

    let intercept = beta[0];
    let coefficients: Vec<f64> = beta.iter().skip(1).copied().collect();

    // Compute R²
    let y_pred = linear_regression_predict(
        &LinearRegression {
            coefficients: coefficients.clone(),
            intercept,
            r_squared: None,
        },
        x,
    );
    let r_squared = compute_r_squared(y, &y_pred);

    Ok(LinearRegression {
        coefficients,
        intercept,
        r_squared: Some(r_squared),
    })
}

/// Fit logistic regression using gradient descent
///
/// Predicts probability using sigmoid: P(y=1) = 1 / (1 + exp(-(β^T x + β₀)))
///
/// # Arguments
/// * `x` - Feature matrix
/// * `y` - Binary labels (0 or 1)
/// * `learning_rate` - Learning rate for gradient descent
/// * `max_iterations` - Maximum iterations
/// * `tolerance` - Convergence tolerance
pub fn logistic_regression(
    x: &[Vec<f64>],
    y: &[f64],
    learning_rate: f64,
    max_iterations: usize,
    tolerance: f64,
) -> Result<LogisticRegression, String> {
    if x.is_empty() || y.is_empty() {
        return Err("Data cannot be empty".to_string());
    }

    if x.len() != y.len() {
        return Err("X and y must have same number of samples".to_string());
    }

    // Check that y is binary
    for &label in y {
        if label != 0.0 && label != 1.0 {
            return Err("Logistic regression requires binary labels (0 or 1)".to_string());
        }
    }

    let n_samples = x.len();
    let n_features = x[0].len();

    // Initialize weights and bias
    let mut coefficients = vec![0.0; n_features];
    let mut intercept = 0.0;

    // Gradient descent
    for iteration in 0..max_iterations {
        // Compute predictions
        let predictions: Vec<f64> = x
            .iter()
            .map(|sample| {
                let logit = intercept
                    + sample
                        .iter()
                        .zip(&coefficients)
                        .map(|(xi, &coef)| xi * coef)
                        .sum::<f64>();
                sigmoid(logit)
            })
            .collect();

        // Compute gradients
        let mut gradient_coef = vec![0.0; n_features];
        let mut gradient_intercept = 0.0;

        for i in 0..n_samples {
            let error = predictions[i] - y[i];
            gradient_intercept += error;

            for j in 0..n_features {
                gradient_coef[j] += error * x[i][j];
            }
        }

        // Average gradients
        gradient_intercept /= n_samples as f64;
        for g in &mut gradient_coef {
            *g /= n_samples as f64;
        }

        // Update parameters
        let old_intercept = intercept;
        let old_coefficients = coefficients.clone();

        intercept -= learning_rate * gradient_intercept;
        for j in 0..n_features {
            coefficients[j] -= learning_rate * gradient_coef[j];
        }

        // Check convergence
        let param_change = (intercept - old_intercept).abs()
            + coefficients
                .iter()
                .zip(&old_coefficients)
                .map(|(new, old)| (new - old).abs())
                .sum::<f64>();

        if param_change < tolerance {
            break;
        }
    }

    // Compute accuracy
    let predictions: Vec<f64> = x
        .iter()
        .map(|sample| {
            let logit = intercept
                + sample
                    .iter()
                    .zip(&coefficients)
                    .map(|(xi, &coef)| xi * coef)
                    .sum::<f64>();
            if sigmoid(logit) >= 0.5 {
                1.0
            } else {
                0.0
            }
        })
        .collect();

    let accuracy = predictions
        .iter()
        .zip(y)
        .filter(|(pred, actual)| **pred == **actual)
        .count() as f64
        / n_samples as f64;

    Ok(LogisticRegression {
        coefficients,
        intercept,
        accuracy: Some(accuracy),
    })
}

/// Predict probabilities using logistic regression
pub fn logistic_regression_predict_proba(
    model: &LogisticRegression,
    x: &[Vec<f64>],
) -> Vec<f64> {
    x.iter()
        .map(|sample| {
            let logit = model.intercept
                + sample
                    .iter()
                    .zip(&model.coefficients)
                    .map(|(xi, &coef)| xi * coef)
                    .sum::<f64>();
            sigmoid(logit)
        })
        .collect()
}

/// Predict class labels using logistic regression
pub fn logistic_regression_predict(model: &LogisticRegression, x: &[Vec<f64>]) -> Vec<f64> {
    logistic_regression_predict_proba(model, x)
        .into_iter()
        .map(|prob| if prob >= 0.5 { 1.0 } else { 0.0 })
        .collect()
}

/// Sigmoid function: σ(x) = 1 / (1 + e^(-x))
pub(crate) fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// Compute R² (coefficient of determination)
///
/// R² = 1 - (SS_res / SS_tot)
/// where SS_res = Σ(y - ŷ)² and SS_tot = Σ(y - ȳ)²
pub(crate) fn compute_r_squared(y_true: &[f64], y_pred: &[f64]) -> f64 {
    let mean_y: f64 = y_true.iter().sum::<f64>() / y_true.len() as f64;

    let ss_tot: f64 = y_true.iter().map(|&y| (y - mean_y).powi(2)).sum();

    let ss_res: f64 = y_true
        .iter()
        .zip(y_pred)
        .map(|(&y, &y_hat)| (y - y_hat).powi(2))
        .sum();

    1.0 - (ss_res / ss_tot)
}

/// Compute Mean Squared Error
pub fn mean_squared_error(y_true: &[f64], y_pred: &[f64]) -> f64 {
    y_true
        .iter()
        .zip(y_pred)
        .map(|(&y, &y_hat)| (y - y_hat).powi(2))
        .sum::<f64>()
        / y_true.len() as f64
}

/// Compute Mean Absolute Error
pub fn mean_absolute_error(y_true: &[f64], y_pred: &[f64]) -> f64 {
    y_true
        .iter()
        .zip(y_pred)
        .map(|(&y, &y_hat)| (y - y_hat).abs())
        .sum::<f64>()
        / y_true.len() as f64
}

