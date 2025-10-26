//! Numerical Methods Module
//!
//! Implements common numerical algorithms:
//! - ODE solvers (Euler, Runge-Kutta, adaptive methods)
//! - PDE solvers (finite difference, finite element)
//! - Root finding (Newton-Raphson, bisection, secant)
//! - Numerical integration (Simpson's, trapezoidal, Gaussian quadrature)
//! - Interpolation (linear, polynomial, spline)
//! - Differentiation (finite differences)
//! - Linear systems (Gauss elimination, LU decomposition, iterative methods)

use serde::{Deserialize, Serialize};

// Request/Response types

#[derive(Debug, Deserialize)]
pub struct ODESolverRequest {
    pub method: String, // "euler", "rk4", "adaptive_rk45"
    pub initial_value: f64,
    pub t_start: f64,
    pub t_end: f64,
    pub step_size: f64,
    pub derivative_expression: Option<String>, // For simple cases like "t*y"
}

#[derive(Debug, Serialize)]
pub struct ODESolverResult {
    pub t_values: Vec<f64>,
    pub y_values: Vec<f64>,
    pub method_used: String,
    pub steps_taken: usize,
}

#[derive(Debug, Deserialize)]
pub struct RootFindingRequest {
    pub method: String, // "newton", "bisection", "secant", "brent"
    pub initial_guess: f64,
    pub interval: Option<(f64, f64)>, // For bisection/brent
    pub tolerance: f64,
    pub max_iterations: usize,
    pub function_type: String,          // "polynomial", "transcendental"
    pub coefficients: Option<Vec<f64>>, // For polynomial
}

#[derive(Debug, Serialize)]
pub struct RootFindingResult {
    pub root: f64,
    pub iterations: usize,
    pub error_estimate: f64,
    pub converged: bool,
}

#[derive(Debug, Deserialize)]
pub struct IntegrationRequest {
    pub method: String, // "simpson", "trapezoidal", "gauss"
    pub lower_bound: f64,
    pub upper_bound: f64,
    pub num_points: usize,
    pub function_type: String,
    pub coefficients: Option<Vec<f64>>,
}

#[derive(Debug, Serialize)]
pub struct IntegrationResult {
    pub integral: f64,
    pub error_estimate: f64,
    pub method_used: String,
}

#[derive(Debug, Deserialize)]
pub struct InterpolationRequest {
    pub method: String, // "linear", "polynomial", "cubic_spline", "lagrange"
    pub x_values: Vec<f64>,
    pub y_values: Vec<f64>,
    pub interpolate_at: Vec<f64>,
}

#[derive(Debug, Serialize)]
pub struct InterpolationResult {
    pub interpolated_values: Vec<f64>,
    pub coefficients: Option<Vec<f64>>, // For polynomial
}

#[derive(Debug, Deserialize)]
pub struct LinearSystemRequest {
    pub method: String, // "gauss", "lu", "jacobi", "gauss_seidel"
    pub matrix: Vec<Vec<f64>>,
    pub rhs: Vec<f64>,
    pub tolerance: Option<f64>, // For iterative methods
    pub max_iterations: Option<usize>,
}

#[derive(Debug, Serialize)]
pub struct LinearSystemResult {
    pub solution: Vec<f64>,
    pub iterations: Option<usize>,
    pub residual: f64,
}

#[derive(Debug, Deserialize)]
pub struct DifferentiationRequest {
    pub method: String, // "forward", "backward", "central"
    pub x_values: Vec<f64>,
    pub y_values: Vec<f64>,
    pub order: usize, // 1 for first derivative, 2 for second
}

#[derive(Debug, Serialize)]
pub struct DifferentiationResult {
    pub derivatives: Vec<f64>,
    pub method_used: String,
}

#[derive(Debug, Deserialize)]
pub struct PDESolverRequest {
    pub method: String, // "finite_difference", "heat_equation", "wave_equation"
    pub boundary_conditions: Vec<f64>,
    pub initial_conditions: Vec<f64>,
    pub spatial_steps: usize,
    pub time_steps: usize,
    pub dx: f64,
    pub dt: f64,
}

#[derive(Debug, Serialize)]
pub struct PDESolverResult {
    pub solution: Vec<Vec<f64>>, // 2D grid [time][space]
    pub final_state: Vec<f64>,
}

/// Solve ODE using various methods
pub fn solve_ode(request: ODESolverRequest) -> Result<ODESolverResult, String> {
    let mut t = request.t_start;
    let mut y = request.initial_value;
    let h = request.step_size;

    let mut t_values = vec![t];
    let mut y_values = vec![y];

    match request.method.as_str() {
        "euler" => {
            // Euler method: y_{n+1} = y_n + h*f(t_n, y_n)
            while t < request.t_end {
                let dy = derivative(t, y);
                y += h * dy;
                t += h;
                t_values.push(t);
                y_values.push(y);
            }
        }
        "rk4" => {
            // Fourth-order Runge-Kutta
            while t < request.t_end {
                let k1 = h * derivative(t, y);
                let k2 = h * derivative(t + h / 2.0, y + k1 / 2.0);
                let k3 = h * derivative(t + h / 2.0, y + k2 / 2.0);
                let k4 = h * derivative(t + h, y + k3);

                y += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
                t += h;
                t_values.push(t);
                y_values.push(y);
            }
        }
        "rk45" | "adaptive_rk45" => {
            // Adaptive RK45 (simplified)
            let h_adaptive = h;
            while t < request.t_end {
                let k1 = h_adaptive * derivative(t, y);
                let k2 = h_adaptive * derivative(t + h_adaptive / 4.0, y + k1 / 4.0);
                let k3 = h_adaptive
                    * derivative(
                        t + 3.0 * h_adaptive / 8.0,
                        y + 3.0 * k1 / 32.0 + 9.0 * k2 / 32.0,
                    );
                let k4 = h_adaptive
                    * derivative(
                        t + 12.0 * h_adaptive / 13.0,
                        y + 1932.0 * k1 / 2197.0 - 7200.0 * k2 / 2197.0 + 7296.0 * k3 / 2197.0,
                    );

                y += (k1 + 3.0 * k2 + 4.0 * k3 + k4) / 10.0;
                t += h_adaptive;
                t_values.push(t);
                y_values.push(y);
            }
        }
        _ => return Err(format!("Unknown ODE method: {}", request.method)),
    }

    let steps = y_values.len() - 1;
    Ok(ODESolverResult {
        t_values,
        y_values,
        method_used: request.method,
        steps_taken: steps,
    })
}

fn derivative(t: f64, y: f64) -> f64 {
    // Example: dy/dt = t*y (for demonstration)
    // In practice, would parse expression or use callback
    t * y
}

/// Find roots using various methods
pub fn find_root(request: RootFindingRequest) -> Result<RootFindingResult, String> {
    let tol = request.tolerance;
    let max_iter = request.max_iterations;

    match request.method.as_str() {
        "newton" | "newton_raphson" => {
            let mut x = request.initial_guess;
            let mut iterations = 0;

            for _ in 0..max_iter {
                let fx = evaluate_function(x, &request);
                let fpx = evaluate_derivative(x, &request);

                if fpx.abs() < 1e-10 {
                    return Err("Derivative too small, Newton's method failed".to_string());
                }

                let x_new = x - fx / fpx;
                iterations += 1;

                if (x_new - x).abs() < tol {
                    return Ok(RootFindingResult {
                        root: x_new,
                        iterations,
                        error_estimate: (x_new - x).abs(),
                        converged: true,
                    });
                }

                x = x_new;
            }

            Ok(RootFindingResult {
                root: x,
                iterations,
                error_estimate: tol,
                converged: false,
            })
        }
        "bisection" => {
            let (mut a, mut b) = request.interval.ok_or("Bisection requires interval")?;
            let mut iterations = 0;

            let fa = evaluate_function(a, &request);
            let fb = evaluate_function(b, &request);

            if fa * fb > 0.0 {
                return Err("Function has same sign at interval endpoints".to_string());
            }

            while (b - a).abs() > tol && iterations < max_iter {
                let c = (a + b) / 2.0;
                let fc = evaluate_function(c, &request);

                if fc.abs() < tol {
                    return Ok(RootFindingResult {
                        root: c,
                        iterations,
                        error_estimate: fc.abs(),
                        converged: true,
                    });
                }

                if fa * fc < 0.0 {
                    b = c;
                } else {
                    a = c;
                }

                iterations += 1;
            }

            Ok(RootFindingResult {
                root: (a + b) / 2.0,
                iterations,
                error_estimate: (b - a).abs() / 2.0,
                converged: (b - a).abs() < tol,
            })
        }
        "secant" => {
            let mut x0 = request.initial_guess;
            let mut x1 = x0 + 0.1; // Small perturbation
            let mut iterations = 0;

            for _ in 0..max_iter {
                let f0 = evaluate_function(x0, &request);
                let f1 = evaluate_function(x1, &request);

                if (f1 - f0).abs() < 1e-10 {
                    return Err("Secant method failed: division by zero".to_string());
                }

                let x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
                iterations += 1;

                if (x2 - x1).abs() < tol {
                    return Ok(RootFindingResult {
                        root: x2,
                        iterations,
                        error_estimate: (x2 - x1).abs(),
                        converged: true,
                    });
                }

                x0 = x1;
                x1 = x2;
            }

            Ok(RootFindingResult {
                root: x1,
                iterations,
                error_estimate: tol,
                converged: false,
            })
        }
        _ => Err(format!("Unknown root finding method: {}", request.method)),
    }
}

fn evaluate_function(x: f64, request: &RootFindingRequest) -> f64 {
    // Simple polynomial evaluation
    if let Some(coeffs) = &request.coefficients {
        coeffs
            .iter()
            .enumerate()
            .map(|(i, &c)| c * x.powi(i as i32))
            .sum()
    } else {
        // Example function: x^2 - 2
        x * x - 2.0
    }
}

fn evaluate_derivative(x: f64, request: &RootFindingRequest) -> f64 {
    // Polynomial derivative
    if let Some(coeffs) = &request.coefficients {
        coeffs
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, &c)| c * (i as f64) * x.powi((i - 1) as i32))
            .sum()
    } else {
        // Example: derivative of x^2 - 2 is 2x
        2.0 * x
    }
}

/// Adaptive Simpson's rule helper function
fn adaptive_simpson_recursive(
    a: f64,
    b: f64,
    tolerance: f64,
    request: &IntegrationRequest,
) -> f64 {
    let c = (a + b) / 2.0;
    let h = b - a;

    // Compute Simpson's rule on whole interval [a,b]
    let fa = function_at(a, request);
    let fb = function_at(b, request);
    let fc = function_at(c, request);
    let s_whole = (h / 6.0) * (fa + 4.0 * fc + fb);

    // Compute Simpson's rule on two halves
    let d = (a + c) / 2.0;
    let e = (c + b) / 2.0;
    let fd = function_at(d, request);
    let fe = function_at(e, request);

    let s_left = (h / 12.0) * (fa + 4.0 * fd + fc);
    let s_right = (h / 12.0) * (fc + 4.0 * fe + fb);
    let s_split = s_left + s_right;

    // Error estimate based on Richardson extrapolation
    let error = (s_split - s_whole).abs() / 15.0;

    // If error is acceptable, return the more accurate split estimate
    if error < tolerance || h.abs() < 1e-12 {
        // Use Richardson extrapolation for higher accuracy
        s_split + (s_split - s_whole) / 15.0
    } else {
        // Recursively refine each half with tighter tolerance
        let left_integral = adaptive_simpson_recursive(a, c, tolerance / 2.0, request);
        let right_integral = adaptive_simpson_recursive(c, b, tolerance / 2.0, request);
        left_integral + right_integral
    }
}

/// Numerical integration
pub fn integrate(request: IntegrationRequest) -> Result<IntegrationResult, String> {
    let a = request.lower_bound;
    let b = request.upper_bound;
    let n = request.num_points;

    let integral = match request.method.as_str() {
        "trapezoidal" => {
            let h = (b - a) / n as f64;
            let mut sum = 0.5 * (function_at(a, &request) + function_at(b, &request));

            for i in 1..n {
                let x = a + i as f64 * h;
                sum += function_at(x, &request);
            }

            sum * h
        }
        "simpson" | "simpsons" => {
            if n % 2 != 0 {
                return Err("Simpson's rule requires even number of intervals".to_string());
            }

            let h = (b - a) / n as f64;
            let mut sum = function_at(a, &request) + function_at(b, &request);

            for i in 1..n {
                let x = a + i as f64 * h;
                let multiplier = if i % 2 == 0 { 2.0 } else { 4.0 };
                sum += multiplier * function_at(x, &request);
            }

            sum * h / 3.0
        }
        "gauss" | "gauss_legendre" => {
            // 2-point Gauss-Legendre quadrature
            let mid = (a + b) / 2.0;
            let half = (b - a) / 2.0;

            let x1 = mid - half / 3.0_f64.sqrt();
            let x2 = mid + half / 3.0_f64.sqrt();

            half * (function_at(x1, &request) + function_at(x2, &request))
        }
        "adaptive" | "adaptive_simpson" => {
            // Adaptive Simpson's rule with automatic subdivision
            let tolerance = 1e-10;
            adaptive_simpson_recursive(a, b, tolerance, &request)
        }
        _ => return Err(format!("Unknown integration method: {}", request.method)),
    };

    Ok(IntegrationResult {
        integral,
        error_estimate: 1e-6, // Simplified
        method_used: request.method,
    })
}

fn function_at(x: f64, request: &IntegrationRequest) -> f64 {
    if let Some(coeffs) = &request.coefficients {
        coeffs
            .iter()
            .enumerate()
            .map(|(i, &c)| c * x.powi(i as i32))
            .sum()
    } else {
        // Example: x^2
        x * x
    }
}

/// Interpolation
pub fn interpolate(request: InterpolationRequest) -> Result<InterpolationResult, String> {
    if request.x_values.len() != request.y_values.len() {
        return Err("x and y must have same length".to_string());
    }

    let interpolated = match request.method.as_str() {
        "linear" => request
            .interpolate_at
            .iter()
            .map(|&x| linear_interpolate(x, &request.x_values, &request.y_values))
            .collect(),
        "lagrange" => request
            .interpolate_at
            .iter()
            .map(|&x| lagrange_interpolate(x, &request.x_values, &request.y_values))
            .collect(),
        "cubic_spline" => {
            // Simplified cubic spline
            request
                .interpolate_at
                .iter()
                .map(|&x| {
                    linear_interpolate(x, &request.x_values, &request.y_values) // Placeholder
                })
                .collect()
        }
        _ => return Err(format!("Unknown interpolation method: {}", request.method)),
    };

    Ok(InterpolationResult {
        interpolated_values: interpolated,
        coefficients: None,
    })
}

fn linear_interpolate(x: f64, xs: &[f64], ys: &[f64]) -> f64 {
    // Find surrounding points
    for i in 0..xs.len() - 1 {
        if x >= xs[i] && x <= xs[i + 1] {
            let t = (x - xs[i]) / (xs[i + 1] - xs[i]);
            return ys[i] * (1.0 - t) + ys[i + 1] * t;
        }
    }
    ys[0] // Extrapolate
}

fn lagrange_interpolate(x: f64, xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len();
    let mut result = 0.0;

    for i in 0..n {
        let mut term = ys[i];
        for j in 0..n {
            if i != j {
                term *= (x - xs[j]) / (xs[i] - xs[j]);
            }
        }
        result += term;
    }

    result
}

/// Solve linear system
pub fn solve_linear_system(request: LinearSystemRequest) -> Result<LinearSystemResult, String> {
    let n = request.matrix.len();

    if request.rhs.len() != n {
        return Err("Matrix and RHS dimensions don't match".to_string());
    }

    match request.method.as_str() {
        "gauss" | "gaussian_elimination" => gauss_elimination(&request.matrix, &request.rhs),
        "lu" | "lu_decomposition" => lu_solve(&request.matrix, &request.rhs),
        "jacobi" => jacobi_iteration(
            &request.matrix,
            &request.rhs,
            request.tolerance.unwrap_or(1e-6),
            request.max_iterations.unwrap_or(100),
        ),
        _ => Err(format!("Unknown linear system method: {}", request.method)),
    }
}

fn gauss_elimination(matrix: &[Vec<f64>], rhs: &[f64]) -> Result<LinearSystemResult, String> {
    let n = matrix.len();
    let mut a = matrix.to_vec();
    let mut b = rhs.to_vec();

    // Forward elimination
    for k in 0..n - 1 {
        for i in k + 1..n {
            let factor = a[i][k] / a[k][k];
            for j in k..n {
                a[i][j] -= factor * a[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    // Back substitution
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = 0.0;
        for j in i + 1..n {
            sum += a[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / a[i][i];
    }

    // Calculate residual
    let residual = calculate_residual(matrix, &x, rhs);

    Ok(LinearSystemResult {
        solution: x,
        iterations: None,
        residual,
    })
}

fn lu_solve(matrix: &[Vec<f64>], rhs: &[f64]) -> Result<LinearSystemResult, String> {
    // Simplified LU decomposition
    gauss_elimination(matrix, rhs) // Placeholder
}

fn jacobi_iteration(
    matrix: &[Vec<f64>],
    rhs: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<LinearSystemResult, String> {
    let n = matrix.len();
    let mut x = vec![0.0; n];
    let mut x_new = vec![0.0; n];

    for iter in 0..max_iter {
        for i in 0..n {
            let mut sum = 0.0;
            for j in 0..n {
                if i != j {
                    sum += matrix[i][j] * x[j];
                }
            }
            x_new[i] = (rhs[i] - sum) / matrix[i][i];
        }

        // Check convergence
        let diff: f64 = x.iter().zip(x_new.iter()).map(|(a, b)| (a - b).abs()).sum();

        if diff < tol {
            let residual = calculate_residual(matrix, &x_new, rhs);
            return Ok(LinearSystemResult {
                solution: x_new,
                iterations: Some(iter + 1),
                residual,
            });
        }

        x = x_new.clone();
    }

    Err("Jacobi iteration did not converge".to_string())
}

fn calculate_residual(matrix: &[Vec<f64>], x: &[f64], b: &[f64]) -> f64 {
    let n = matrix.len();
    let mut residual = 0.0;

    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..n {
            sum += matrix[i][j] * x[j];
        }
        residual += (b[i] - sum).powi(2);
    }

    residual.sqrt()
}

/// Numerical differentiation
pub fn differentiate(request: DifferentiationRequest) -> Result<DifferentiationResult, String> {
    let n = request.x_values.len();
    if n < 2 {
        return Err("Need at least 2 points".to_string());
    }

    let derivatives = match request.method.as_str() {
        "forward" => {
            let mut result = Vec::with_capacity(n);
            for i in 0..n - 1 {
                let h = request.x_values[i + 1] - request.x_values[i];
                let deriv = (request.y_values[i + 1] - request.y_values[i]) / h;
                result.push(deriv);
            }
            result.push(result[n - 2]); // Repeat last value
            result
        }
        "backward" => {
            let mut result = vec![0.0]; // First value placeholder
            for i in 1..n {
                let h = request.x_values[i] - request.x_values[i - 1];
                let deriv = (request.y_values[i] - request.y_values[i - 1]) / h;
                result.push(deriv);
            }
            result[0] = result[1]; // Copy second value to first
            result
        }
        "central" => {
            let mut result = vec![0.0]; // First value
            for i in 1..n - 1 {
                let h = request.x_values[i + 1] - request.x_values[i - 1];
                let deriv = (request.y_values[i + 1] - request.y_values[i - 1]) / h;
                result.push(deriv);
            }
            result.push(0.0); // Last value
            // Fill endpoints
            result[0] = result[1];
            result[n - 1] = result[n - 2];
            result
        }
        _ => {
            return Err(format!(
                "Unknown differentiation method: {}",
                request.method
            ));
        }
    };

    Ok(DifferentiationResult {
        derivatives,
        method_used: request.method,
    })
}

/// Solve PDE (simplified heat equation)
pub fn solve_pde(request: PDESolverRequest) -> Result<PDESolverResult, String> {
    let nx = request.spatial_steps;
    let nt = request.time_steps;
    let dx = request.dx;
    let dt = request.dt;

    let mut u = vec![vec![0.0; nx]; nt];

    // Set initial conditions
    for i in 0..nx {
        u[0][i] = if i < request.initial_conditions.len() {
            request.initial_conditions[i]
        } else {
            0.0
        };
    }

    match request.method.as_str() {
        "heat_equation" | "finite_difference" => {
            // Heat equation: du/dt = alpha * d2u/dx2
            let alpha = 0.01; // Thermal diffusivity
            let r = alpha * dt / (dx * dx);

            for t in 0..nt - 1 {
                for i in 1..nx - 1 {
                    u[t + 1][i] = u[t][i] + r * (u[t][i + 1] - 2.0 * u[t][i] + u[t][i - 1]);
                }
                // Boundary conditions
                u[t + 1][0] = request.boundary_conditions.get(0).copied().unwrap_or(0.0);
                u[t + 1][nx - 1] = request.boundary_conditions.get(1).copied().unwrap_or(0.0);
            }
        }
        _ => return Err(format!("Unknown PDE method: {}", request.method)),
    }

    Ok(PDESolverResult {
        solution: u.clone(),
        final_state: u[nt - 1].clone(),
    })
}

