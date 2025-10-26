//! Linear Algebra Module
//!
//! Provides fundamental linear algebra operations including:
//! - Matrix decompositions (SVD, eigendecomposition, QR, LU, Cholesky)
//! - Dimensionality reduction (PCA)
//! - Matrix operations (rank, pseudoinverse, norms)

use nalgebra::{DMatrix, SVD, Schur};
use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
pub struct MatrixInput {
    pub matrix: Vec<Vec<f64>>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct SVDResult {
    pub u: Vec<Vec<f64>>,
    pub singular_values: Vec<f64>,
    pub v_transpose: Vec<Vec<f64>>,
    pub rank: usize,
    pub condition_number: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct EigenResult {
    pub eigenvalues: Vec<f64>,
    pub eigenvectors: Vec<Vec<f64>>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PCAResult {
    pub principal_components: Vec<Vec<f64>>,
    pub explained_variance: Vec<f64>,
    pub cumulative_variance: Vec<f64>,
    pub transformed_data: Vec<Vec<f64>>,
}

/// Compute Singular Value Decomposition (SVD)
pub fn compute_svd(input: MatrixInput) -> Result<SVDResult, String> {
    let matrix = DMatrix::from_row_slice(
        input.matrix.len(),
        input.matrix[0].len(),
        &input.matrix.concat(),
    );

    let svd = SVD::new(matrix.clone(), true, true);

    let u_mat = svd.u.ok_or("Failed to compute U matrix")?;
    let singular_values = svd.singular_values.as_slice().to_vec();
    let v_t_mat = svd.v_t.ok_or("Failed to compute V^T matrix")?;

    let rank = singular_values.iter().filter(|&&s| s > 1e-10).count();
    let condition_number = if !singular_values.is_empty() {
        singular_values[0] / singular_values.last().unwrap_or(&1e-10)
    } else {
        f64::INFINITY
    };

    let u_rows = u_mat.nrows();
    let u_cols = u_mat.ncols();
    let u_matrix: Vec<Vec<f64>> = (0..u_rows)
        .map(|i| (0..u_cols).map(|j| u_mat[(i, j)]).collect())
        .collect();

    let v_rows = v_t_mat.nrows();
    let v_cols = v_t_mat.ncols();
    let v_t_matrix: Vec<Vec<f64>> = (0..v_rows)
        .map(|i| (0..v_cols).map(|j| v_t_mat[(i, j)]).collect())
        .collect();

    Ok(SVDResult {
        u: u_matrix,
        singular_values,
        v_transpose: v_t_matrix,
        rank,
        condition_number,
    })
}

/// Compute matrix rank
pub fn compute_matrix_rank(input: MatrixInput) -> Result<usize, String> {
    let svd_result = compute_svd(input)?;
    Ok(svd_result.rank)
}

/// Compute Moore-Penrose pseudoinverse
pub fn compute_pseudoinverse(input: MatrixInput) -> Result<Vec<Vec<f64>>, String> {
    let matrix = DMatrix::from_row_slice(
        input.matrix.len(),
        input.matrix[0].len(),
        &input.matrix.concat(),
    );

    let svd = SVD::new(matrix, true, true);
    let pinv = svd
        .pseudo_inverse(1e-10)
        .map_err(|e| format!("Pseudoinverse failed: {}", e))?;

    let rows = pinv.nrows();
    let cols = pinv.ncols();

    // Extract properly using indexing (nalgebra uses column-major storage)
    let result: Vec<Vec<f64>> = (0..rows)
        .map(|i| (0..cols).map(|j| pinv[(i, j)]).collect())
        .collect();

    Ok(result)
}

/// Compute Principal Component Analysis (PCA)
pub fn compute_pca(input: MatrixInput, n_components: Option<usize>) -> Result<PCAResult, String> {
    let n_samples = input.matrix.len();
    let n_features = input.matrix[0].len();

    // Center the data
    let mean: Vec<f64> = (0..n_features)
        .map(|j| input.matrix.iter().map(|row| row[j]).sum::<f64>() / n_samples as f64)
        .collect();

    let centered: Vec<Vec<f64>> = input
        .matrix
        .iter()
        .map(|row| row.iter().zip(&mean).map(|(x, m)| x - m).collect())
        .collect();

    // Compute covariance matrix
    let centered_flat: Vec<f64> = centered.concat();
    let centered_matrix = DMatrix::from_row_slice(n_samples, n_features, &centered_flat);
    let cov = centered_matrix.transpose() * &centered_matrix / (n_samples - 1) as f64;

    // Compute eigendecomposition of covariance matrix
    let eigen = cov.symmetric_eigen();
    let mut eigenvalues: Vec<f64> = eigen.eigenvalues.as_slice().to_vec();
    let mut eigenvectors: Vec<Vec<f64>> = (0..n_features)
        .map(|i| eigen.eigenvectors.column(i).as_slice().to_vec())
        .collect();

    // Sort by eigenvalue (descending)
    let mut pairs: Vec<(f64, Vec<f64>)> = eigenvalues
        .iter()
        .zip(eigenvectors.iter())
        .map(|(val, vec)| (*val, vec.clone()))
        .collect();
    pairs.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

    eigenvalues = pairs.iter().map(|(val, _)| *val).collect();
    eigenvectors = pairs.iter().map(|(_, vec)| vec.clone()).collect();

    // Select top n components (limit to available components)
    let max_comp = n_features.min(n_samples);
    let n_comp = n_components.unwrap_or(max_comp).min(max_comp);
    let principal_components: Vec<Vec<f64>> = eigenvectors.into_iter().take(n_comp).collect();

    // Calculate explained variance
    let total_variance: f64 = eigenvalues.iter().sum();
    let explained_variance: Vec<f64> = eigenvalues
        .iter()
        .take(n_comp)
        .map(|ev| ev / total_variance)
        .collect();

    let cumulative_variance: Vec<f64> = explained_variance
        .iter()
        .scan(0.0, |acc, &x| {
            *acc += x;
            Some(*acc)
        })
        .collect();

    // Transform data
    let transformed_data: Vec<Vec<f64>> = centered
        .iter()
        .map(|row| {
            principal_components
                .iter()
                .map(|pc| row.iter().zip(pc).map(|(a, b)| a * b).sum())
                .collect()
        })
        .collect();

    Ok(PCAResult {
        principal_components,
        explained_variance,
        cumulative_variance,
        transformed_data,
    })
}

/// Compute eigendecomposition (for symmetric matrices)
pub fn compute_eigendecomposition(input: MatrixInput) -> Result<EigenResult, String> {
    let matrix = DMatrix::from_row_slice(
        input.matrix.len(),
        input.matrix[0].len(),
        &input.matrix.concat(),
    );

    // Check if matrix is symmetric
    if !matrix.is_square() {
        return Err("Matrix must be square for eigendecomposition".to_string());
    }

    let eigen = matrix.symmetric_eigen();
    let eigenvalues = eigen.eigenvalues.as_slice().to_vec();

    let n = eigen.eigenvectors.nrows();
    let eigenvectors: Vec<Vec<f64>> = (0..n)
        .map(|i| eigen.eigenvectors.column(i).as_slice().to_vec())
        .collect();

    Ok(EigenResult {
        eigenvalues,
        eigenvectors,
    })
}

/// Matrix operations (determinant, trace, norm)
#[derive(Debug, Serialize, Deserialize)]
pub struct MatrixOpsResult {
    pub determinant: Option<f64>,
    pub trace: Option<f64>,
    pub frobenius_norm: f64,
    pub max_norm: f64,
}

pub fn matrix_operations(input: MatrixInput) -> Result<MatrixOpsResult, String> {
    let matrix = DMatrix::from_row_slice(
        input.matrix.len(),
        input.matrix[0].len(),
        &input.matrix.concat(),
    );

    let determinant = if matrix.is_square() {
        matrix.determinant().into()
    } else {
        None
    };

    let trace = if matrix.is_square() {
        Some(matrix.trace())
    } else {
        None
    };

    let frobenius_norm = matrix.norm();
    let max_norm = matrix.amax();

    Ok(MatrixOpsResult {
        determinant,
        trace,
        frobenius_norm,
        max_norm,
    })
}

// Additional matrix operations

#[derive(Debug, Deserialize)]
pub struct MatrixNormRequest {
    pub matrix: Vec<Vec<f64>>,
    pub norm_type: String, // "frobenius", "1", "inf", "2"
}

#[derive(Debug, Serialize)]
pub struct MatrixNormResult {
    pub norm: f64,
    pub norm_type: String,
}

#[derive(Debug, Deserialize)]
pub struct MatrixPowerRequest {
    pub matrix: Vec<Vec<f64>>,
    pub power: i32,
}

#[derive(Debug, Serialize)]
pub struct MatrixPowerResult {
    pub result: Vec<Vec<f64>>,
}

#[derive(Debug, Deserialize)]
pub struct MatrixExpRequest {
    pub matrix: Vec<Vec<f64>>,
    pub terms: Option<usize>, // Number of Taylor series terms
}

#[derive(Debug, Serialize)]
pub struct MatrixExpResult {
    pub result: Vec<Vec<f64>>,
}

#[derive(Debug, Deserialize)]
pub struct CholeskyRequest {
    pub matrix: Vec<Vec<f64>>,
}

#[derive(Debug, Serialize)]
pub struct CholeskyResult {
    pub lower: Vec<Vec<f64>>,
}

#[derive(Debug, Deserialize)]
pub struct SchurRequest {
    pub matrix: Vec<Vec<f64>>,
}

#[derive(Debug, Serialize)]
pub struct SchurResult {
    pub t: Vec<Vec<f64>>, // Upper triangular matrix
    pub u: Vec<Vec<f64>>, // Unitary matrix
}

#[derive(Debug, Deserialize)]
pub struct LURequest {
    pub matrix: Vec<Vec<f64>>,
}

#[derive(Debug, Serialize)]
pub struct LUResult {
    pub l: Vec<Vec<f64>>, // Lower triangular matrix
    pub u: Vec<Vec<f64>>, // Upper triangular matrix
    pub p: Vec<Vec<f64>>, // Permutation matrix
}

pub fn matrix_norm(request: MatrixNormRequest) -> Result<MatrixNormResult, String> {
    let mat = &request.matrix;
    let m = mat.len();
    let n = mat[0].len();

    let norm = match request.norm_type.as_str() {
        "frobenius" | "fro" => {
            // Frobenius norm: sqrt of sum of squares of all elements
            mat.iter()
                .flat_map(|row| row.iter())
                .map(|&x| x * x)
                .sum::<f64>()
                .sqrt()
        }
        "1" => {
            // Maximum absolute column sum
            (0..n)
                .map(|j| (0..m).map(|i| mat[i][j].abs()).sum::<f64>())
                .fold(0.0, f64::max)
        }
        "inf" | "infinity" => {
            // Maximum absolute row sum
            mat.iter()
                .map(|row| row.iter().map(|x| x.abs()).sum::<f64>())
                .fold(0.0, f64::max)
        }
        _ => return Err(format!("Unknown norm type: {}", request.norm_type)),
    };

    Ok(MatrixNormResult {
        norm,
        norm_type: request.norm_type,
    })
}

pub fn matrix_power(request: MatrixPowerRequest) -> Result<MatrixPowerResult, String> {
    let n = request.matrix.len();
    if request.matrix[0].len() != n {
        return Err("Matrix must be square".to_string());
    }

    let mut result = identity_matrix(n);

    if request.power == 0 {
        return Ok(MatrixPowerResult { result });
    }

    let mut base = request.matrix.clone();
    let mut exp = request.power.abs();

    while exp > 0 {
        if exp % 2 == 1 {
            result = matrix_multiply(&result, &base);
        }
        base = matrix_multiply(&base, &base);
        exp /= 2;
    }

    Ok(MatrixPowerResult { result })
}

pub fn matrix_exp(request: MatrixExpRequest) -> Result<MatrixExpResult, String> {
    let n = request.matrix.len();
    let terms = request.terms.unwrap_or(20);

    let mut result = identity_matrix(n);
    let mut term = identity_matrix(n);
    let mut fact = 1.0_f64;

    for k in 1..terms {
        term = matrix_multiply(&term, &request.matrix);
        fact *= k as f64;
        let scalar = 1.0 / fact;
        let scaled_term = matrix_scale(&term, scalar);
        result = matrix_add(&result, &scaled_term);
    }

    Ok(MatrixExpResult { result })
}

pub fn cholesky_decomposition(request: CholeskyRequest) -> Result<CholeskyResult, String> {
    let a = &request.matrix;
    let n = a.len();

    let mut l = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in 0..=i {
            let mut sum = 0.0;
            for k in 0..j {
                sum += l[i][k] * l[j][k];
            }

            if i == j {
                let val = a[i][i] - sum;
                if val < 0.0 {
                    return Err("Matrix is not positive definite".to_string());
                }
                l[i][j] = val.sqrt();
            } else {
                l[i][j] = (a[i][j] - sum) / l[j][j];
            }
        }
    }

    Ok(CholeskyResult { lower: l })
}

/// LU decomposition with partial pivoting: PA = LU
pub fn lu_decomposition(request: LURequest) -> Result<LUResult, String> {
    let a = &request.matrix;
    let n = a.len();

    if a.is_empty() || a[0].len() != n {
        return Err("Matrix must be square".to_string());
    }

    // Initialize L, U, and permutation matrix P
    let mut l = identity_matrix(n);
    let mut u = a.clone();
    let mut p = identity_matrix(n);

    // Perform LU decomposition with partial pivoting
    for k in 0..n {
        // Find pivot (row with largest absolute value in column k)
        let mut max_val = u[k][k].abs();
        let mut pivot_row = k;

        for i in (k + 1)..n {
            let val = u[i][k].abs();
            if val > max_val {
                max_val = val;
                pivot_row = i;
            }
        }

        // Swap rows if needed
        if pivot_row != k {
            u.swap(k, pivot_row);
            p.swap(k, pivot_row);

            // Swap already computed parts of L
            if k > 0 {
                for j in 0..k {
                    let temp = l[k][j];
                    l[k][j] = l[pivot_row][j];
                    l[pivot_row][j] = temp;
                }
            }
        }

        // Check for singular matrix
        if u[k][k].abs() < 1e-14 {
            return Err("Matrix is singular or nearly singular".to_string());
        }

        // Elimination
        for i in (k + 1)..n {
            let factor = u[i][k] / u[k][k];
            l[i][k] = factor;

            for j in k..n {
                u[i][j] -= factor * u[k][j];
            }
        }
    }

    Ok(LUResult { l, u, p })
}

/// Schur decomposition: A = U T U^H
/// where T is upper triangular (Schur form) and U is unitary
pub fn schur_decomposition(request: SchurRequest) -> Result<SchurResult, String> {
    let a = &request.matrix;
    let n = a.len();

    if a.is_empty() || a[0].len() != n {
        return Err("Matrix must be square".to_string());
    }

    // Convert to nalgebra matrix
    let matrix = DMatrix::from_row_slice(n, n, &a.concat());

    // Compute Schur decomposition
    let schur = Schur::new(matrix);

    // Unpack returns (T, U) where T is upper triangular and U is unitary
    let (t_mat, u_mat) = schur.unpack();

    // Extract T (upper triangular Schur form)
    let t_rows = t_mat.nrows();
    let t_cols = t_mat.ncols();
    let t: Vec<Vec<f64>> = (0..t_rows)
        .map(|i| (0..t_cols).map(|j| t_mat[(i, j)]).collect())
        .collect();

    // Extract U (unitary matrix)
    let u_rows = u_mat.nrows();
    let u_cols = u_mat.ncols();
    let u: Vec<Vec<f64>> = (0..u_rows)
        .map(|i| (0..u_cols).map(|j| u_mat[(i, j)]).collect())
        .collect();

    Ok(SchurResult { t, u })
}

/// Compute matrix inverse using nalgebra
pub fn matrix_inverse(matrix: Vec<Vec<f64>>) -> Result<Vec<Vec<f64>>, String> {
    let n = matrix.len();

    if matrix.is_empty() || matrix[0].len() != n {
        return Err("Matrix must be square".to_string());
    }

    // Convert to nalgebra matrix
    let mat = DMatrix::from_row_slice(n, n, &matrix.concat());

    // Compute inverse
    let inv = mat
        .try_inverse()
        .ok_or("Matrix is singular and cannot be inverted".to_string())?;

    // Convert back to Vec<Vec<f64>>
    // nalgebra uses column-major storage, so we need to extract properly
    let result: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| inv[(i, j)]).collect())
        .collect();

    Ok(result)
}

// Helper functions
fn identity_matrix(n: usize) -> Vec<Vec<f64>> {
    let mut mat = vec![vec![0.0; n]; n];
    for i in 0..n {
        mat[i][i] = 1.0;
    }
    mat
}

fn matrix_multiply(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let m = a.len();
    let n = b[0].len();
    let p = b.len();

    let mut result = vec![vec![0.0; n]; m];

    for i in 0..m {
        for j in 0..n {
            for k in 0..p {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    result
}

fn matrix_add(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let m = a.len();
    let n = a[0].len();

    let mut result = vec![vec![0.0; n]; m];

    for i in 0..m {
        for j in 0..n {
            result[i][j] = a[i][j] + b[i][j];
        }
    }

    result
}

fn matrix_scale(a: &[Vec<f64>], scalar: f64) -> Vec<Vec<f64>> {
    a.iter()
        .map(|row| row.iter().map(|&x| x * scalar).collect())
        .collect()
}

fn factorial(n: usize) -> usize {
    (1..=n).product()
}

// Test module - part of this module, can access private functions
#[cfg(test)]
#[path = "../../../tests/unit/mathematics/linear_algebra_tests.rs"]
mod tests;
