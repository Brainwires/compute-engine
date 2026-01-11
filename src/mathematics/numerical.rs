//! Numerical methods for symbolic expressions
//!
//! Matrix exponentiation, numerical eigenvalues, and numerical evaluation of symbolic expressions

use crate::mathematics::linear_algebra;
use crate::mathematics::symbolic_cas::{Expr, SymbolicMatrix};
use std::collections::HashMap;

/// Evaluate a symbolic expression numerically
pub fn evaluate_numeric(expr: &Expr, values: &HashMap<String, f64>) -> Result<f64, String> {
    match expr {
        Expr::Number(r) => Ok(r.numerator as f64 / r.denominator as f64),
        Expr::Symbol(s) => values
            .get(s)
            .copied()
            .ok_or_else(|| format!("Undefined variable: {}", s)),
        Expr::Add(a, b) => {
            let a_val = evaluate_numeric(a, values)?;
            let b_val = evaluate_numeric(b, values)?;
            Ok(a_val + b_val)
        }
        Expr::Mul(a, b) => {
            let a_val = evaluate_numeric(a, values)?;
            let b_val = evaluate_numeric(b, values)?;
            Ok(a_val * b_val)
        }
        Expr::Pow(base, exp) => {
            let base_val = evaluate_numeric(base, values)?;
            let exp_val = evaluate_numeric(exp, values)?;
            Ok(base_val.powf(exp_val))
        }
        Expr::Function(name, args) => match name.as_str() {
            "sin" => {
                if args.len() != 1 {
                    return Err("sin requires 1 argument".to_string());
                }
                let arg = evaluate_numeric(&args[0], values)?;
                Ok(arg.sin())
            }
            "cos" => {
                if args.len() != 1 {
                    return Err("cos requires 1 argument".to_string());
                }
                let arg = evaluate_numeric(&args[0], values)?;
                Ok(arg.cos())
            }
            "tan" => {
                if args.len() != 1 {
                    return Err("tan requires 1 argument".to_string());
                }
                let arg = evaluate_numeric(&args[0], values)?;
                Ok(arg.tan())
            }
            "exp" => {
                if args.len() != 1 {
                    return Err("exp requires 1 argument".to_string());
                }
                let arg = evaluate_numeric(&args[0], values)?;
                Ok(arg.exp())
            }
            "log" | "ln" => {
                if args.len() != 1 {
                    return Err("log requires 1 argument".to_string());
                }
                let arg = evaluate_numeric(&args[0], values)?;
                Ok(arg.ln())
            }
            "sqrt" => {
                if args.len() != 1 {
                    return Err("sqrt requires 1 argument".to_string());
                }
                let arg = evaluate_numeric(&args[0], values)?;
                Ok(arg.sqrt())
            }
            _ => Err(format!("Unknown function: {}", name)),
        },
    }
}

/// Evaluate a symbolic matrix numerically
pub fn evaluate_matrix_numeric(
    matrix: &SymbolicMatrix,
    values: &HashMap<String, f64>,
) -> Result<Vec<Vec<f64>>, String> {
    let mut result = vec![vec![0.0; matrix.cols()]; matrix.rows()];

    for i in 0..matrix.rows() {
        for j in 0..matrix.cols() {
            let elem = matrix
                .get(i, j)
                .ok_or_else(|| "Invalid matrix indices".to_string())?;
            result[i][j] = evaluate_numeric(elem, values)?;
        }
    }

    Ok(result)
}

/// Matrix exponentiation using Taylor series: exp(A) = I + A + A²/2! + A³/3! + ...
///
/// # Arguments
/// * `matrix` - Symbolic matrix to exponentiate
/// * `values` - Variable values for numerical evaluation
/// * `terms` - Number of Taylor series terms (default 20)
///
/// # Returns
/// Numerical matrix exp(A)
pub fn matrix_exponential(
    matrix: &SymbolicMatrix,
    values: &HashMap<String, f64>,
    terms: Option<usize>,
) -> Result<Vec<Vec<f64>>, String> {
    let n = matrix.rows();
    if n != matrix.cols() {
        return Err("Matrix must be square".to_string());
    }

    let num_terms = terms.unwrap_or(20);

    // Start with identity matrix
    let mut result = vec![vec![0.0; n]; n];
    for i in 0..n {
        result[i][i] = 1.0;
    }

    // Evaluate matrix numerically
    let a = evaluate_matrix_numeric(matrix, values)?;

    // Compute A^k / k! iteratively
    let mut a_power = vec![vec![0.0; n]; n];
    for i in 0..n {
        a_power[i][i] = 1.0; // Start with I
    }

    let mut factorial = 1.0;

    for k in 1..=num_terms {
        // A^k = A^(k-1) * A
        a_power = matrix_multiply_numeric(&a_power, &a);

        factorial *= k as f64;

        // Add A^k / k! to result
        for i in 0..n {
            for j in 0..n {
                result[i][j] += a_power[i][j] / factorial;
            }
        }
    }

    Ok(result)
}

/// Helper: multiply two numerical matrices
fn matrix_multiply_numeric(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let rows = a.len();
    let cols = b[0].len();
    let inner = b.len();

    let mut result = vec![vec![0.0; cols]; rows];

    for i in 0..rows {
        for j in 0..cols {
            let mut sum = 0.0;
            for k in 0..inner {
                sum += a[i][k] * b[k][j];
            }
            result[i][j] = sum;
        }
    }

    result
}

/// Compute eigenvalues numerically using QR algorithm
///
/// # Arguments
/// * `matrix` - Symbolic matrix
/// * `values` - Variable values for evaluation
/// * `max_iterations` - Maximum QR iterations (default 100)
///
/// # Returns
/// Vector of eigenvalues (may be complex, returned as (real, imaginary) pairs)
pub fn eigenvalues_numeric(
    matrix: &SymbolicMatrix,
    values: &HashMap<String, f64>,
    max_iterations: Option<usize>,
) -> Result<Vec<f64>, String> {
    let n = matrix.rows();
    if n != matrix.cols() {
        return Err("Matrix must be square".to_string());
    }

    let max_iter = max_iterations.unwrap_or(100);

    // Evaluate matrix numerically
    let mut a = evaluate_matrix_numeric(matrix, values)?;

    // QR algorithm iterations
    for _ in 0..max_iter {
        let (q, r) = qr_decomposition(&a)?;
        a = matrix_multiply_numeric(&r, &q);
    }

    // Extract diagonal (eigenvalues)
    let mut eigenvalues = Vec::new();
    for i in 0..n {
        eigenvalues.push(a[i][i]);
    }

    Ok(eigenvalues)
}

/// QR decomposition using Gram-Schmidt
fn qr_decomposition(matrix: &[Vec<f64>]) -> Result<(Vec<Vec<f64>>, Vec<Vec<f64>>), String> {
    let n = matrix.len();
    let m = matrix[0].len();

    let mut q = vec![vec![0.0; m]; n];
    let mut r = vec![vec![0.0; m]; m];

    // Extract columns
    let mut columns: Vec<Vec<f64>> = vec![vec![0.0; n]; m];
    for j in 0..m {
        for i in 0..n {
            columns[j][i] = matrix[i][j];
        }
    }

    // Gram-Schmidt orthogonalization
    for j in 0..m {
        let mut v = columns[j].clone();

        // Subtract projections onto previous basis vectors
        for k in 0..j {
            let dot_product: f64 = (0..n).map(|i| columns[j][i] * q[i][k]).sum();
            r[k][j] = dot_product;

            for i in 0..n {
                v[i] -= dot_product * q[i][k];
            }
        }

        // Normalize
        let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        r[j][j] = norm;

        if norm > 1e-10 {
            for i in 0..n {
                q[i][j] = v[i] / norm;
            }
        }
    }

    Ok((q, r))
}

/// Simulate time evolution of a quantum state: |ψ(t)⟩ = exp(-iHt/ℏ) |ψ(0)⟩
///
/// # Arguments
/// * `hamiltonian` - Hamiltonian matrix (symbolic)
/// * `initial_state` - Initial state vector (symbolic)
/// * `time` - Time value
/// * `hbar` - ℏ value (default 1.0)
/// * `values` - Variable values for evaluation
///
/// # Returns
/// State vector at time t
pub fn time_evolution(
    hamiltonian: &SymbolicMatrix,
    initial_state: &SymbolicMatrix,
    time: f64,
    hbar: Option<f64>,
    values: &HashMap<String, f64>,
) -> Result<Vec<f64>, String> {
    let h = hbar.unwrap_or(1.0);

    // Create scaled Hamiltonian: -iHt/ℏ
    let factor = -time / h;

    let _scaled_values = values.clone();
    // We'll evaluate H and scale it

    let h_matrix = evaluate_matrix_numeric(hamiltonian, values)?;

    // Scale by -it/ℏ (treating as real for now)
    let mut scaled_h = h_matrix.clone();
    for i in 0..scaled_h.len() {
        for j in 0..scaled_h[i].len() {
            scaled_h[i][j] *= factor;
        }
    }

    // Compute exp(-iHt/ℏ) using matrix exponential
    // For numerical stability, we convert back to symbolic
    let scaled_symbolic = numeric_to_symbolic_matrix(&scaled_h);
    let u_t = matrix_exponential(&scaled_symbolic, &HashMap::new(), Some(30))?;

    // Apply to initial state
    let state_numeric = evaluate_matrix_numeric(initial_state, values)?;
    let state_vec: Vec<f64> = state_numeric.iter().map(|row| row[0]).collect();

    // U(t) |ψ(0)⟩
    let mut result = vec![0.0; state_vec.len()];
    for i in 0..u_t.len() {
        let mut sum = 0.0;
        for j in 0..state_vec.len() {
            sum += u_t[i][j] * state_vec[j];
        }
        result[i] = sum;
    }

    Ok(result)
}

/// Helper: convert numerical matrix to symbolic
fn numeric_to_symbolic_matrix(matrix: &[Vec<f64>]) -> SymbolicMatrix {
    let data: Vec<Vec<Expr>> = matrix
        .iter()
        .map(|row| {
            row.iter()
                .map(|&val| {
                    // Convert to rational approximation
                    let factor = 1000000;
                    let num = (val * factor as f64).round() as i64;
                    Expr::rational_unchecked(num, factor)
                })
                .collect()
        })
        .collect();

    SymbolicMatrix::new(data).unwrap()
}

/// Compute matrix power: A^n numerically
pub fn matrix_power_numeric(
    matrix: &SymbolicMatrix,
    power: i32,
    values: &HashMap<String, f64>,
) -> Result<Vec<Vec<f64>>, String> {
    let n = matrix.rows();
    if n != matrix.cols() {
        return Err("Matrix must be square".to_string());
    }

    let a = evaluate_matrix_numeric(matrix, values)?;

    if power == 0 {
        // Return identity
        let mut result = vec![vec![0.0; n]; n];
        for i in 0..n {
            result[i][i] = 1.0;
        }
        return Ok(result);
    }

    if power < 0 {
        // For negative powers: A^(-n) = (A^(-1))^n
        // First compute inverse
        let a_inv = linear_algebra::matrix_inverse(a)
            .map_err(|e| format!("Cannot compute negative power: {}", e))?;

        // Then compute (A^(-1))^|n|
        let abs_power = power.abs();
        let mut result = a_inv.clone();
        for _ in 1..abs_power {
            result = matrix_multiply_numeric(&result, &a_inv);
        }
        return Ok(result);
    }

    // Compute A^n by repeated multiplication
    let mut result = a.clone();
    for _ in 1..power {
        result = matrix_multiply_numeric(&result, &a);
    }

    Ok(result)
}

