//! Optimization Module
//!
//! Provides optimization algorithms for:
//! - Unconstrained optimization (gradient descent, Nelder-Mead)
//! - Constrained optimization (linear programming)
//! - Curve fitting and parameter estimation
//! - Sensitivity analysis

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Serialize, Deserialize)]
pub struct OptimizationRequest {
    pub method: String,
    pub objective_function: Option<String>, // For simple expressions
    pub initial_point: Vec<f64>,
    pub gradient: Option<Vec<String>>,
    pub constraints: Option<Vec<Constraint>>,
    pub options: OptimizationOptions,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Constraint {
    pub coefficients: Vec<f64>,
    pub value: f64,
    pub constraint_type: String, // "le", "ge", "eq"
}

#[derive(Debug, Serialize, Deserialize)]
pub struct OptimizationOptions {
    #[serde(default = "default_max_iterations")]
    pub max_iterations: usize,
    #[serde(default = "default_tolerance")]
    pub tolerance: f64,
    #[serde(default = "default_step_size")]
    pub step_size: f64,
    #[serde(default)]
    pub maximize: bool,
}

fn default_max_iterations() -> usize {
    1000
}
fn default_tolerance() -> f64 {
    1e-6
}
fn default_step_size() -> f64 {
    0.01
}

impl Default for OptimizationOptions {
    fn default() -> Self {
        Self {
            max_iterations: 1000,
            tolerance: 1e-6,
            step_size: 0.01,
            maximize: false,
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct OptimizationResult {
    pub success: bool,
    pub variables: HashMap<String, f64>,
    pub objective_value: f64,
    pub iterations: usize,
    pub method: String,
    pub convergence: bool,
    pub error: Option<String>,
}

/// Gradient descent optimization
pub fn gradient_descent(
    objective: fn(&[f64]) -> f64,
    gradient: fn(&[f64]) -> Vec<f64>,
    initial_point: Vec<f64>,
    options: OptimizationOptions,
) -> Result<OptimizationResult, String> {
    let mut x = initial_point.clone();
    let mut iterations = 0;
    let mut prev_objective = objective(&x);

    for iter in 0..options.max_iterations {
        iterations = iter + 1;

        let grad = gradient(&x);
        let grad_norm: f64 = grad.iter().map(|g| g * g).sum::<f64>().sqrt();

        if grad_norm < options.tolerance {
            break;
        }

        // Update variables
        for i in 0..x.len() {
            if options.maximize {
                x[i] += options.step_size * grad[i]; // Maximize: ascend gradient
            } else {
                x[i] -= options.step_size * grad[i]; // Minimize: descend gradient
            }
        }

        let current_objective = objective(&x);

        // Check convergence
        if (current_objective - prev_objective).abs() < options.tolerance {
            break;
        }

        prev_objective = current_objective;
    }

    let mut variables = HashMap::new();
    for (i, &val) in x.iter().enumerate() {
        variables.insert(format!("x{}", i + 1), val);
    }

    Ok(OptimizationResult {
        success: true,
        variables,
        objective_value: objective(&x),
        iterations,
        method: "gradient_descent".to_string(),
        convergence: iterations < options.max_iterations,
        error: None,
    })
}

/// Nelder-Mead simplex method
pub fn nelder_mead(
    objective: fn(&[f64]) -> f64,
    initial_point: Vec<f64>,
    options: OptimizationOptions,
) -> Result<OptimizationResult, String> {
    let n = initial_point.len();
    let alpha = 1.0; // reflection
    let gamma = 2.0; // expansion
    let beta = 0.5; // contraction
    let sigma = 0.5; // shrink

    // Create initial simplex
    let mut simplex: Vec<Vec<f64>> = Vec::new();
    simplex.push(initial_point.clone());

    for i in 0..n {
        let mut point = initial_point.clone();
        point[i] += if point[i] != 0.0 {
            0.05 * point[i]
        } else {
            0.00025
        };
        simplex.push(point);
    }

    let mut iterations = 0;

    for iter in 0..options.max_iterations {
        iterations = iter + 1;

        // Sort simplex by objective value
        if options.maximize {
            // For maximization: best = highest value
            simplex.sort_by(|a, b| objective(b).partial_cmp(&objective(a)).unwrap());
        } else {
            // For minimization: best = lowest value
            simplex.sort_by(|a, b| objective(a).partial_cmp(&objective(b)).unwrap());
        }

        let best = simplex[0].clone();
        let worst = simplex[n].clone();
        let second_worst = simplex[n - 1].clone();

        // Check convergence
        let range = (objective(&worst) - objective(&best)).abs();
        if range < options.tolerance {
            break;
        }

        // Calculate centroid
        let mut centroid = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                centroid[j] += simplex[i][j];
            }
        }
        for j in 0..n {
            centroid[j] /= n as f64;
        }

        // Reflection
        let reflected: Vec<f64> = centroid
            .iter()
            .zip(worst.iter())
            .map(|(c, w)| c + alpha * (c - w))
            .collect();
        let reflected_value = objective(&reflected);

        // Comparisons depend on whether we're maximizing or minimizing
        let is_better = |a: f64, b: f64| {
            if options.maximize {
                a > b
            } else {
                a < b
            }
        };

        let is_between = |val: f64, bound1: f64, bound2: f64| {
            if options.maximize {
                val < bound1 && val >= bound2
            } else {
                val >= bound1 && val < bound2
            }
        };

        if is_between(reflected_value, objective(&best), objective(&second_worst)) {
            simplex[n] = reflected;
            continue;
        }

        // Expansion
        if is_better(reflected_value, objective(&best)) {
            let expanded: Vec<f64> = centroid
                .iter()
                .zip(reflected.iter())
                .map(|(c, r)| c + gamma * (r - c))
                .collect();
            let expanded_value = objective(&expanded);

            if is_better(expanded_value, reflected_value) {
                simplex[n] = expanded;
            } else {
                simplex[n] = reflected;
            }
            continue;
        }

        // Contraction
        let contracted: Vec<f64> = centroid
            .iter()
            .zip(worst.iter())
            .map(|(c, w)| c + beta * (w - c))
            .collect();
        let contracted_value = objective(&contracted);

        if is_better(contracted_value, objective(&worst)) {
            simplex[n] = contracted;
            continue;
        }

        // Shrink
        for i in 1..=n {
            for j in 0..n {
                simplex[i][j] = best[j] + sigma * (simplex[i][j] - best[j]);
            }
        }
    }

    let result = &simplex[0];
    let mut variables = HashMap::new();
    for (i, &val) in result.iter().enumerate() {
        variables.insert(format!("x{}", i + 1), val);
    }

    Ok(OptimizationResult {
        success: true,
        variables,
        objective_value: objective(result),
        iterations,
        method: "nelder_mead".to_string(),
        convergence: iterations < options.max_iterations,
        error: None,
    })
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CurveFitRequest {
    pub x_data: Vec<f64>,
    pub y_data: Vec<f64>,
    pub model: String, // "linear", "quadratic", "exponential", "power"
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CurveFitResult {
    pub coefficients: Vec<f64>,
    pub r_squared: f64,
    pub residuals: Vec<f64>,
    pub model: String,
    pub aic: Option<f64>,
    pub bic: Option<f64>,
    pub aicc: Option<f64>,
}

/// Calculate information criteria for model selection
///
/// # Arguments
/// * `residuals` - Model residuals (observed - predicted)
/// * `n_params` - Number of parameters in the model
/// * `n_samples` - Number of data points
fn calculate_information_criteria(
    residuals: &[f64],
    n_params: usize,
    n_samples: usize,
) -> (f64, f64, f64) {
    let n = n_samples as f64;
    let k = n_params as f64;

    // Calculate RSS (Residual Sum of Squares)
    let rss: f64 = residuals.iter().map(|r| r * r).sum();

    // Estimate variance (MLE)
    let sigma2 = rss / n;

    // Prevent log(0) or negative sigma
    let sigma2 = sigma2.max(1e-10);

    // Calculate log-likelihood (assuming normal errors)
    // L = -n/2 * ln(2π) - n/2 * ln(σ²) - RSS/(2σ²)
    // Simplified: L = -n/2 * ln(σ²) - n/2
    let log_likelihood = -n / 2.0 * sigma2.ln() - n / 2.0;

    // AIC = 2k - 2ln(L)
    let aic = 2.0 * k - 2.0 * log_likelihood;

    // BIC = ln(n)k - 2ln(L)
    let bic = n.ln() * k - 2.0 * log_likelihood;

    // AICc = AIC + 2k(k+1)/(n-k-1)  (corrected for small samples)
    let aicc = if n > k + 1.0 {
        aic + (2.0 * k * (k + 1.0)) / (n - k - 1.0)
    } else {
        aic + 1000.0 // Penalty for too many parameters
    };

    (aic, bic, aicc)
}

/// Curve fitting using least squares
pub fn curve_fitting(request: CurveFitRequest) -> Result<CurveFitResult, String> {
    if request.x_data.len() != request.y_data.len() {
        return Err("X and Y data must have the same length".to_string());
    }

    let n = request.x_data.len();

    match request.model.as_str() {
        "linear" => {
            // y = a + bx
            let sum_x: f64 = request.x_data.iter().sum();
            let sum_y: f64 = request.y_data.iter().sum();
            let sum_xy: f64 = request
                .x_data
                .iter()
                .zip(request.y_data.iter())
                .map(|(x, y)| x * y)
                .sum();
            let sum_x2: f64 = request.x_data.iter().map(|x| x * x).sum();

            let n_f64 = n as f64;
            let b = (n_f64 * sum_xy - sum_x * sum_y) / (n_f64 * sum_x2 - sum_x * sum_x);
            let a = (sum_y - b * sum_x) / n_f64;

            let fitted: Vec<f64> = request.x_data.iter().map(|x| a + b * x).collect();
            let residuals: Vec<f64> = request
                .y_data
                .iter()
                .zip(fitted.iter())
                .map(|(y, f)| y - f)
                .collect();

            let ss_res: f64 = residuals.iter().map(|r| r * r).sum();
            let y_mean = sum_y / n_f64;
            let ss_tot: f64 = request.y_data.iter().map(|y| (y - y_mean).powi(2)).sum();
            let r_squared = 1.0 - (ss_res / ss_tot);

            // Calculate information criteria
            let n_params = 2; // a, b
            let (aic, bic, aicc) = calculate_information_criteria(&residuals, n_params, n);

            Ok(CurveFitResult {
                coefficients: vec![a, b],
                r_squared,
                residuals,
                model: "linear".to_string(),
                aic: Some(aic),
                bic: Some(bic),
                aicc: Some(aicc),
            })
        }
        "quadratic" => {
            // y = a + bx + cx^2
            // Using matrix solution for least squares
            let mut sum_x = 0.0;
            let mut sum_x2 = 0.0;
            let mut sum_x3 = 0.0;
            let mut sum_x4 = 0.0;
            let mut sum_y = 0.0;
            let mut sum_xy = 0.0;
            let mut sum_x2y = 0.0;

            for i in 0..n {
                let x = request.x_data[i];
                let y = request.y_data[i];
                let x2 = x * x;
                sum_x += x;
                sum_x2 += x2;
                sum_x3 += x2 * x;
                sum_x4 += x2 * x2;
                sum_y += y;
                sum_xy += x * y;
                sum_x2y += x2 * y;
            }

            let n_f64 = n as f64;

            // Solve system using Cramer's rule (simplified for 3x3)
            let det = n_f64 * (sum_x2 * sum_x4 - sum_x3 * sum_x3)
                - sum_x * (sum_x * sum_x4 - sum_x3 * sum_x2)
                + sum_x2 * (sum_x * sum_x3 - sum_x2 * sum_x2);

            let a = (sum_y * (sum_x2 * sum_x4 - sum_x3 * sum_x3)
                - sum_x * (sum_xy * sum_x4 - sum_x3 * sum_x2y)
                + sum_x2 * (sum_xy * sum_x3 - sum_x2 * sum_x2y))
                / det;

            let b = (n_f64 * (sum_xy * sum_x4 - sum_x3 * sum_x2y)
                - sum_y * (sum_x * sum_x4 - sum_x3 * sum_x2)
                + sum_x2 * (sum_x * sum_x2y - sum_xy * sum_x2))
                / det;

            let c = (n_f64 * (sum_x2 * sum_x2y - sum_xy * sum_x3)
                - sum_x * (sum_x * sum_x2y - sum_xy * sum_x2)
                + sum_y * (sum_x * sum_x3 - sum_x2 * sum_x2))
                / det;

            let fitted: Vec<f64> = request
                .x_data
                .iter()
                .map(|x| a + b * x + c * x * x)
                .collect();
            let residuals: Vec<f64> = request
                .y_data
                .iter()
                .zip(fitted.iter())
                .map(|(y, f)| y - f)
                .collect();

            let ss_res: f64 = residuals.iter().map(|r| r * r).sum();
            let y_mean = sum_y / n_f64;
            let ss_tot: f64 = request.y_data.iter().map(|y| (y - y_mean).powi(2)).sum();
            let r_squared = 1.0 - (ss_res / ss_tot);

            // Calculate information criteria
            let n_params = 3; // a, b, c
            let (aic, bic, aicc) = calculate_information_criteria(&residuals, n_params, n);

            Ok(CurveFitResult {
                coefficients: vec![a, b, c],
                r_squared,
                residuals,
                model: "quadratic".to_string(),
                aic: Some(aic),
                bic: Some(bic),
                aicc: Some(aicc),
            })
        }
        "trigonometric" | "sinusoidal" => {
            // Trigonometric (linearized): y = a + b*sin(x) + c*cos(x)
            // Note: Equivalent to A*sin(x + φ) where A = √(b²+c²), tan(φ) = c/b
            // This linearized form allows for analytical least-squares solution

            let n_f64 = n as f64;
            let sum_y: f64 = request.y_data.iter().sum();
            let sum_sin: f64 = request.x_data.iter().map(|x| x.sin()).sum();
            let sum_cos: f64 = request.x_data.iter().map(|x| x.cos()).sum();

            let sum_y_sin: f64 = request
                .x_data
                .iter()
                .zip(request.y_data.iter())
                .map(|(x, y)| y * x.sin())
                .sum();
            let sum_y_cos: f64 = request
                .x_data
                .iter()
                .zip(request.y_data.iter())
                .map(|(x, y)| y * x.cos())
                .sum();

            let sum_sin2: f64 = request.x_data.iter().map(|x| x.sin().powi(2)).sum();
            let sum_cos2: f64 = request.x_data.iter().map(|x| x.cos().powi(2)).sum();
            let sum_sin_cos: f64 = request.x_data.iter().map(|x| x.sin() * x.cos()).sum();

            // Solve 3x3 system for [a, b, c]
            let det = n_f64 * (sum_sin2 * sum_cos2 - sum_sin_cos * sum_sin_cos)
                - sum_sin * (sum_sin * sum_cos2 - sum_sin_cos * sum_cos)
                + sum_cos * (sum_sin * sum_sin_cos - sum_sin2 * sum_cos);

            if det.abs() < 1e-10 {
                return Err("Singular matrix in trigonometric fitting".to_string());
            }

            let a = (sum_y * (sum_sin2 * sum_cos2 - sum_sin_cos * sum_sin_cos)
                - sum_sin * (sum_y_sin * sum_cos2 - sum_sin_cos * sum_y_cos)
                + sum_cos * (sum_y_sin * sum_sin_cos - sum_sin2 * sum_y_cos))
                / det;

            let b = (n_f64 * (sum_y_sin * sum_cos2 - sum_sin_cos * sum_y_cos)
                - sum_y * (sum_sin * sum_cos2 - sum_sin_cos * sum_cos)
                + sum_cos * (sum_sin * sum_y_cos - sum_y_sin * sum_cos))
                / det;

            let c = (n_f64 * (sum_sin2 * sum_y_cos - sum_y_sin * sum_sin_cos)
                - sum_sin * (sum_sin * sum_y_cos - sum_y_sin * sum_cos)
                + sum_y * (sum_sin * sum_sin_cos - sum_sin2 * sum_cos))
                / det;

            let fitted: Vec<f64> = request
                .x_data
                .iter()
                .map(|x| a + b * x.sin() + c * x.cos())
                .collect();
            let residuals: Vec<f64> = request
                .y_data
                .iter()
                .zip(fitted.iter())
                .map(|(y, f)| y - f)
                .collect();

            let ss_res: f64 = residuals.iter().map(|r| r * r).sum();
            let y_mean = sum_y / n_f64;
            let ss_tot: f64 = request.y_data.iter().map(|y| (y - y_mean).powi(2)).sum();
            let r_squared = 1.0 - (ss_res / ss_tot);

            // Calculate amplitude and phase for cleaner output
            let amplitude = (b * b + c * c).sqrt();
            let phase = c.atan2(b);

            // Calculate information criteria
            let n_params = 3; // a, b, c (or a, amplitude, phase)
            let (aic, bic, aicc) = calculate_information_criteria(&residuals, n_params, n);

            Ok(CurveFitResult {
                coefficients: vec![a, amplitude, phase],
                r_squared,
                residuals,
                model: "trigonometric".to_string(),
                aic: Some(aic),
                bic: Some(bic),
                aicc: Some(aicc),
            })
        }
        "exponential" => {
            // y = a * e^(b*x)
            // Linearize: ln(y) = ln(a) + b*x

            // Check for non-positive y values
            if request.y_data.iter().any(|&y| y <= 0.0) {
                return Err("Exponential fitting requires all y values to be positive".to_string());
            }

            let ln_y: Vec<f64> = request.y_data.iter().map(|y| y.ln()).collect();

            let sum_x: f64 = request.x_data.iter().sum();
            let sum_ln_y: f64 = ln_y.iter().sum();
            let sum_x_ln_y: f64 = request
                .x_data
                .iter()
                .zip(ln_y.iter())
                .map(|(x, ln_y)| x * ln_y)
                .sum();
            let sum_x2: f64 = request.x_data.iter().map(|x| x * x).sum();

            let n_f64 = n as f64;
            let b = (n_f64 * sum_x_ln_y - sum_x * sum_ln_y) / (n_f64 * sum_x2 - sum_x * sum_x);
            let ln_a = (sum_ln_y - b * sum_x) / n_f64;
            let a = ln_a.exp();

            let fitted: Vec<f64> = request.x_data.iter().map(|x| a * (b * x).exp()).collect();
            let residuals: Vec<f64> = request
                .y_data
                .iter()
                .zip(fitted.iter())
                .map(|(y, f)| y - f)
                .collect();

            let ss_res: f64 = residuals.iter().map(|r| r * r).sum();
            let y_mean = request.y_data.iter().sum::<f64>() / n_f64;
            let ss_tot: f64 = request.y_data.iter().map(|y| (y - y_mean).powi(2)).sum();
            let r_squared = 1.0 - (ss_res / ss_tot);

            let n_params = 2; // a, b
            let (aic, bic, aicc) = calculate_information_criteria(&residuals, n_params, n);

            Ok(CurveFitResult {
                coefficients: vec![a, b],
                r_squared,
                residuals,
                model: "exponential".to_string(),
                aic: Some(aic),
                bic: Some(bic),
                aicc: Some(aicc),
            })
        }
        "logarithmic" => {
            // y = a + b*ln(x)

            // Check for non-positive x values
            if request.x_data.iter().any(|&x| x <= 0.0) {
                return Err("Logarithmic fitting requires all x values to be positive".to_string());
            }

            let ln_x: Vec<f64> = request.x_data.iter().map(|x| x.ln()).collect();

            let sum_ln_x: f64 = ln_x.iter().sum();
            let sum_y: f64 = request.y_data.iter().sum();
            let sum_ln_x_y: f64 = ln_x
                .iter()
                .zip(request.y_data.iter())
                .map(|(ln_x, y)| ln_x * y)
                .sum();
            let sum_ln_x2: f64 = ln_x.iter().map(|ln_x| ln_x * ln_x).sum();

            let n_f64 = n as f64;
            let b =
                (n_f64 * sum_ln_x_y - sum_ln_x * sum_y) / (n_f64 * sum_ln_x2 - sum_ln_x * sum_ln_x);
            let a = (sum_y - b * sum_ln_x) / n_f64;

            let fitted: Vec<f64> = request.x_data.iter().map(|x| a + b * x.ln()).collect();
            let residuals: Vec<f64> = request
                .y_data
                .iter()
                .zip(fitted.iter())
                .map(|(y, f)| y - f)
                .collect();

            let ss_res: f64 = residuals.iter().map(|r| r * r).sum();
            let y_mean = sum_y / n_f64;
            let ss_tot: f64 = request.y_data.iter().map(|y| (y - y_mean).powi(2)).sum();
            let r_squared = 1.0 - (ss_res / ss_tot);

            let n_params = 2; // a, b
            let (aic, bic, aicc) = calculate_information_criteria(&residuals, n_params, n);

            Ok(CurveFitResult {
                coefficients: vec![a, b],
                r_squared,
                residuals,
                model: "logarithmic".to_string(),
                aic: Some(aic),
                bic: Some(bic),
                aicc: Some(aicc),
            })
        }
        "power" => {
            // y = a * x^b
            // Linearize: ln(y) = ln(a) + b*ln(x)

            // Check for non-positive values
            if request.x_data.iter().any(|&x| x <= 0.0) {
                return Err("Power law fitting requires all x values to be positive".to_string());
            }
            if request.y_data.iter().any(|&y| y <= 0.0) {
                return Err("Power law fitting requires all y values to be positive".to_string());
            }

            let ln_x: Vec<f64> = request.x_data.iter().map(|x| x.ln()).collect();
            let ln_y: Vec<f64> = request.y_data.iter().map(|y| y.ln()).collect();

            let sum_ln_x: f64 = ln_x.iter().sum();
            let sum_ln_y: f64 = ln_y.iter().sum();
            let sum_ln_x_ln_y: f64 = ln_x
                .iter()
                .zip(ln_y.iter())
                .map(|(ln_x, ln_y)| ln_x * ln_y)
                .sum();
            let sum_ln_x2: f64 = ln_x.iter().map(|ln_x| ln_x * ln_x).sum();

            let n_f64 = n as f64;
            let b = (n_f64 * sum_ln_x_ln_y - sum_ln_x * sum_ln_y)
                / (n_f64 * sum_ln_x2 - sum_ln_x * sum_ln_x);
            let ln_a = (sum_ln_y - b * sum_ln_x) / n_f64;
            let a = ln_a.exp();

            let fitted: Vec<f64> = request.x_data.iter().map(|x| a * x.powf(b)).collect();
            let residuals: Vec<f64> = request
                .y_data
                .iter()
                .zip(fitted.iter())
                .map(|(y, f)| y - f)
                .collect();

            let ss_res: f64 = residuals.iter().map(|r| r * r).sum();
            let y_mean = request.y_data.iter().sum::<f64>() / n_f64;
            let ss_tot: f64 = request.y_data.iter().map(|y| (y - y_mean).powi(2)).sum();
            let r_squared = 1.0 - (ss_res / ss_tot);

            let n_params = 2; // a, b
            let (aic, bic, aicc) = calculate_information_criteria(&residuals, n_params, n);

            Ok(CurveFitResult {
                coefficients: vec![a, b],
                r_squared,
                residuals,
                model: "power".to_string(),
                aic: Some(aic),
                bic: Some(bic),
                aicc: Some(aicc),
            })
        }
        "rational" => {
            // y = (a + b*x) / (1 + c*x)
            // Simplified rational function fit
            // Rearrange: y + c*x*y = a + b*x
            // This is linear in parameters: y = a - c*x*y + b*x

            // Use iterative approach with initial c=0
            let mut c = 0.0;
            let max_iterations = 10;

            for _ in 0..max_iterations {
                // Adjusted y values
                let y_adjusted: Vec<f64> = request
                    .y_data
                    .iter()
                    .zip(request.x_data.iter())
                    .map(|(y, x)| y * (1.0 + c * x))
                    .collect();

                // Fit linear model: y_adjusted = a + b*x
                let sum_x: f64 = request.x_data.iter().sum();
                let sum_y_adj: f64 = y_adjusted.iter().sum();
                let sum_xy_adj: f64 = request
                    .x_data
                    .iter()
                    .zip(y_adjusted.iter())
                    .map(|(x, y_adj)| x * y_adj)
                    .sum();
                let sum_x2: f64 = request.x_data.iter().map(|x| x * x).sum();

                let n_f64 = n as f64;
                let b = (n_f64 * sum_xy_adj - sum_x * sum_y_adj) / (n_f64 * sum_x2 - sum_x * sum_x);
                let a = (sum_y_adj - b * sum_x) / n_f64;

                // Estimate new c from residuals
                let residuals: Vec<f64> = request
                    .y_data
                    .iter()
                    .zip(request.x_data.iter())
                    .map(|(y, x)| {
                        let predicted = (a + b * x) / (1.0 + c * x);
                        y - predicted
                    })
                    .collect();

                // Simple update to c based on residuals
                let c_update: f64 = residuals
                    .iter()
                    .zip(request.x_data.iter())
                    .zip(request.y_data.iter())
                    .map(|((r, x), y)| r * x * y)
                    .sum::<f64>()
                    / n_f64;

                c += 0.1 * c_update; // Damped update
            }

            // Final coefficients
            let fitted: Vec<f64> = request
                .x_data
                .iter()
                .map(|x| {
                    // Recompute a and b with final c
                    let y_adjusted: Vec<f64> = request
                        .y_data
                        .iter()
                        .zip(request.x_data.iter())
                        .map(|(y, x_i)| y * (1.0 + c * x_i))
                        .collect();

                    let sum_x: f64 = request.x_data.iter().sum();
                    let sum_y_adj: f64 = y_adjusted.iter().sum();
                    let sum_xy_adj: f64 = request
                        .x_data
                        .iter()
                        .zip(y_adjusted.iter())
                        .map(|(x_i, y_adj)| x_i * y_adj)
                        .sum();
                    let sum_x2: f64 = request.x_data.iter().map(|x_i| x_i * x_i).sum();

                    let n_f64 = n as f64;
                    let b =
                        (n_f64 * sum_xy_adj - sum_x * sum_y_adj) / (n_f64 * sum_x2 - sum_x * sum_x);
                    let a = (sum_y_adj - b * sum_x) / n_f64;

                    (a + b * x) / (1.0 + c * x)
                })
                .collect();

            let residuals: Vec<f64> = request
                .y_data
                .iter()
                .zip(fitted.iter())
                .map(|(y, f)| y - f)
                .collect();

            let ss_res: f64 = residuals.iter().map(|r| r * r).sum();
            let y_mean = request.y_data.iter().sum::<f64>() / (n as f64);
            let ss_tot: f64 = request.y_data.iter().map(|y| (y - y_mean).powi(2)).sum();
            let r_squared = 1.0 - (ss_res / ss_tot);

            // Get final a, b coefficients
            let y_adjusted: Vec<f64> = request
                .y_data
                .iter()
                .zip(request.x_data.iter())
                .map(|(y, x)| y * (1.0 + c * x))
                .collect();

            let sum_x: f64 = request.x_data.iter().sum();
            let sum_y_adj: f64 = y_adjusted.iter().sum();
            let sum_xy_adj: f64 = request
                .x_data
                .iter()
                .zip(y_adjusted.iter())
                .map(|(x, y_adj)| x * y_adj)
                .sum();
            let sum_x2: f64 = request.x_data.iter().map(|x| x * x).sum();

            let n_f64 = n as f64;
            let b = (n_f64 * sum_xy_adj - sum_x * sum_y_adj) / (n_f64 * sum_x2 - sum_x * sum_x);
            let a = (sum_y_adj - b * sum_x) / n_f64;

            let n_params = 3; // a, b, c
            let (aic, bic, aicc) = calculate_information_criteria(&residuals, n_params, n);

            Ok(CurveFitResult {
                coefficients: vec![a, b, c],
                r_squared,
                residuals,
                model: "rational".to_string(),
                aic: Some(aic),
                bic: Some(bic),
                aicc: Some(aicc),
            })
        }
        _ => Err(format!("Unsupported model: {}", request.model)),
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct SensitivityRequest {
    pub parameters: HashMap<String, f64>,
    pub perturbation: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct SensitivityResult {
    pub sensitivities: HashMap<String, f64>,
    pub total_sensitivity: f64,
}

/// Simple sensitivity analysis
pub fn sensitivity_analysis(
    objective: fn(&HashMap<String, f64>) -> f64,
    request: SensitivityRequest,
) -> Result<SensitivityResult, String> {
    let base_value = objective(&request.parameters);
    let mut sensitivities = HashMap::new();
    let mut total_sensitivity = 0.0;

    for (param, &value) in &request.parameters {
        let mut perturbed = request.parameters.clone();
        perturbed.insert(param.clone(), value + request.perturbation);

        let perturbed_value = objective(&perturbed);
        let sensitivity = (perturbed_value - base_value) / request.perturbation;

        sensitivities.insert(param.clone(), sensitivity);
        total_sensitivity += sensitivity.abs();
    }

    Ok(SensitivityResult {
        sensitivities,
        total_sensitivity,
    })
}
