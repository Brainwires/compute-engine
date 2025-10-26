//! Unit tests for mathematics::linear_algebra

use crate::mathematics::linear_algebra::*;

// Test SVD computation
#[test]
fn test_svd_basic() {
    let input = MatrixInput {
        matrix: vec![vec![1.0, 2.0], vec![3.0, 4.0]],
    };

    let result = compute_svd(input).unwrap();

    assert_eq!(result.u.len(), 2);
    assert_eq!(result.u[0].len(), 2);
    assert_eq!(result.singular_values.len(), 2);
    assert_eq!(result.v_transpose.len(), 2);
    assert!(result.singular_values[0] > result.singular_values[1]);
}

#[test]
fn test_svd_rank_deficient() {
    let input = MatrixInput {
        matrix: vec![vec![1.0, 2.0], vec![2.0, 4.0]], // Rank 1
    };

    let result = compute_svd(input).unwrap();

    assert_eq!(result.rank, 1);
    assert!(result.singular_values[1] < 1e-10);
}

#[test]
fn test_svd_identity() {
    let input = MatrixInput {
        matrix: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
    };

    let result = compute_svd(input).unwrap();

    assert_eq!(result.rank, 2);
    assert!((result.singular_values[0] - 1.0).abs() < 1e-10);
    assert!((result.singular_values[1] - 1.0).abs() < 1e-10);
}

// Test matrix rank
#[test]
fn test_matrix_rank_full() {
    let input = MatrixInput {
        matrix: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
    };

    let rank = compute_matrix_rank(input).unwrap();
    assert_eq!(rank, 2);
}

#[test]
fn test_matrix_rank_deficient() {
    let input = MatrixInput {
        matrix: vec![vec![1.0, 2.0, 3.0], vec![2.0, 4.0, 6.0], vec![3.0, 6.0, 9.0]],
    };

    let rank = compute_matrix_rank(input).unwrap();
    assert_eq!(rank, 1);
}

// Test pseudoinverse
#[test]
fn test_pseudoinverse_full_rank() {
    let input = MatrixInput {
        matrix: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
    };

    let pinv = compute_pseudoinverse(input).unwrap();

    assert_eq!(pinv.len(), 2);
    assert_eq!(pinv[0].len(), 2);
    assert!((pinv[0][0] - 1.0).abs() < 1e-10);
    assert!((pinv[1][1] - 1.0).abs() < 1e-10);
}

#[test]
fn test_pseudoinverse_rank_deficient() {
    let input = MatrixInput {
        matrix: vec![vec![1.0, 2.0], vec![2.0, 4.0]],
    };

    let pinv = compute_pseudoinverse(input).unwrap();

    assert_eq!(pinv.len(), 2);
    assert_eq!(pinv[0].len(), 2);
}

// Test PCA
#[test]
fn test_pca_basic() {
    let input = MatrixInput {
        matrix: vec![
            vec![2.5, 2.4],
            vec![0.5, 0.7],
            vec![2.2, 2.9],
            vec![1.9, 2.2],
            vec![3.1, 3.0],
        ],
    };

    let result = compute_pca(input, Some(2)).unwrap();

    assert_eq!(result.principal_components.len(), 2);
    assert_eq!(result.explained_variance.len(), 2);
    assert_eq!(result.cumulative_variance.len(), 2);
    assert_eq!(result.transformed_data.len(), 5);
    assert_eq!(result.transformed_data[0].len(), 2);
}

#[test]
fn test_pca_variance_sums_to_one() {
    let input = MatrixInput {
        matrix: vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ],
    };

    let result = compute_pca(input, None).unwrap();

    let total: f64 = result.explained_variance.iter().sum();
    assert!((total - 1.0).abs() < 1e-10);
}

#[test]
fn test_pca_cumulative_variance() {
    let input = MatrixInput {
        matrix: vec![
            vec![1.0, 2.0],
            vec![3.0, 4.0],
            vec![5.0, 6.0],
        ],
    };

    let result = compute_pca(input, None).unwrap();

    assert!((result.cumulative_variance.last().unwrap() - 1.0).abs() < 1e-10);
    assert!(result.cumulative_variance[0] <= result.cumulative_variance[1]);
}

// Test eigendecomposition
#[test]
fn test_eigendecomposition_symmetric() {
    let input = MatrixInput {
        matrix: vec![vec![4.0, 1.0], vec![1.0, 3.0]],
    };

    let result = compute_eigendecomposition(input).unwrap();

    assert_eq!(result.eigenvalues.len(), 2);
    assert_eq!(result.eigenvectors.len(), 2);
}

#[test]
fn test_eigendecomposition_identity() {
    let input = MatrixInput {
        matrix: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
    };

    let result = compute_eigendecomposition(input).unwrap();

    assert!((result.eigenvalues[0] - 1.0).abs() < 1e-10);
    assert!((result.eigenvalues[1] - 1.0).abs() < 1e-10);
}

#[test]
fn test_eigendecomposition_non_square() {
    let input = MatrixInput {
        matrix: vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]],
    };

    let result = compute_eigendecomposition(input);
    assert!(result.is_err());
}

// Test matrix operations
#[test]
fn test_matrix_operations_square() {
    let input = MatrixInput {
        matrix: vec![vec![1.0, 2.0], vec![3.0, 4.0]],
    };

    let result = matrix_operations(input).unwrap();

    assert!(result.determinant.is_some());
    assert!((result.determinant.unwrap() - (-2.0)).abs() < 1e-10);
    assert!(result.trace.is_some());
    assert!((result.trace.unwrap() - 5.0).abs() < 1e-10);
    assert!(result.frobenius_norm > 0.0);
    assert!(result.max_norm > 0.0);
}

#[test]
fn test_matrix_operations_non_square() {
    let input = MatrixInput {
        matrix: vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]],
    };

    let result = matrix_operations(input).unwrap();

    assert!(result.determinant.is_none());
    assert!(result.trace.is_none());
    assert!(result.frobenius_norm > 0.0);
}

// Test matrix norm
#[test]
fn test_matrix_norm_frobenius() {
    let request = MatrixNormRequest {
        matrix: vec![vec![1.0, 2.0], vec![3.0, 4.0]],
        norm_type: "frobenius".to_string(),
    };

    let result = matrix_norm(request).unwrap();

    let expected = (1.0 + 4.0 + 9.0 + 16.0_f64).sqrt();
    assert!((result.norm - expected).abs() < 1e-10);
}

#[test]
fn test_matrix_norm_1() {
    let request = MatrixNormRequest {
        matrix: vec![vec![1.0, -2.0], vec![3.0, -4.0]],
        norm_type: "1".to_string(),
    };

    let result = matrix_norm(request).unwrap();

    assert!((result.norm - 6.0).abs() < 1e-10); // max column sum: max(|1|+|3|, |-2|+|-4|) = 6
}

#[test]
fn test_matrix_norm_inf() {
    let request = MatrixNormRequest {
        matrix: vec![vec![1.0, -2.0], vec![3.0, -4.0]],
        norm_type: "inf".to_string(),
    };

    let result = matrix_norm(request).unwrap();

    assert!((result.norm - 7.0).abs() < 1e-10); // max row sum: max(|1|+|-2|, |3|+|-4|) = 7
}

// Test matrix power
#[test]
fn test_matrix_power_zero() {
    let request = MatrixPowerRequest {
        matrix: vec![vec![1.0, 2.0], vec![3.0, 4.0]],
        power: 0,
    };

    let result = matrix_power(request).unwrap();

    assert_eq!(result.result, vec![vec![1.0, 0.0], vec![0.0, 1.0]]);
}

#[test]
fn test_matrix_power_one() {
    let mat = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
    let request = MatrixPowerRequest {
        matrix: mat.clone(),
        power: 1,
    };

    let result = matrix_power(request).unwrap();

    assert_eq!(result.result, mat);
}

#[test]
fn test_matrix_power_two() {
    let request = MatrixPowerRequest {
        matrix: vec![vec![2.0, 0.0], vec![0.0, 3.0]],
        power: 2,
    };

    let result = matrix_power(request).unwrap();

    assert!((result.result[0][0] - 4.0).abs() < 1e-10);
    assert!((result.result[1][1] - 9.0).abs() < 1e-10);
}

#[test]
fn test_matrix_power_non_square() {
    let request = MatrixPowerRequest {
        matrix: vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]],
        power: 2,
    };

    let result = matrix_power(request);
    assert!(result.is_err());
}

// Test matrix exponential
#[test]
fn test_matrix_exp_zero() {
    let request = MatrixExpRequest {
        matrix: vec![vec![0.0, 0.0], vec![0.0, 0.0]],
        terms: Some(10),
    };

    let result = matrix_exp(request).unwrap();

    // exp(0) = I
    assert!((result.result[0][0] - 1.0).abs() < 1e-10);
    assert!((result.result[1][1] - 1.0).abs() < 1e-10);
    assert!(result.result[0][1].abs() < 1e-10);
    assert!(result.result[1][0].abs() < 1e-10);
}

#[test]
fn test_matrix_exp_diagonal() {
    let request = MatrixExpRequest {
        matrix: vec![vec![1.0, 0.0], vec![0.0, 2.0]],
        terms: Some(20),
    };

    let result = matrix_exp(request).unwrap();

    // exp([[1,0],[0,2]]) ≈ [[e,0],[0,e^2]]
    assert!((result.result[0][0] - 1.0_f64.exp()).abs() < 1e-6);
    assert!((result.result[1][1] - 2.0_f64.exp()).abs() < 1e-6);
}

// Test Cholesky decomposition
#[test]
fn test_cholesky_basic() {
    let request = CholeskyRequest {
        matrix: vec![vec![4.0, 2.0], vec![2.0, 3.0]],
    };

    let result = cholesky_decomposition(request).unwrap();

    assert_eq!(result.lower.len(), 2);
    assert!((result.lower[0][0] - 2.0).abs() < 1e-10);
}

#[test]
fn test_cholesky_identity() {
    let request = CholeskyRequest {
        matrix: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
    };

    let result = cholesky_decomposition(request).unwrap();

    assert!((result.lower[0][0] - 1.0).abs() < 1e-10);
    assert!((result.lower[1][1] - 1.0).abs() < 1e-10);
}

#[test]
fn test_cholesky_not_positive_definite() {
    let request = CholeskyRequest {
        matrix: vec![vec![-1.0, 0.0], vec![0.0, 1.0]],
    };

    let result = cholesky_decomposition(request);
    assert!(result.is_err());
}

// Test LU decomposition
#[test]
fn test_lu_basic() {
    let request = LURequest {
        matrix: vec![vec![2.0, 1.0], vec![4.0, 3.0]],
    };

    let result = lu_decomposition(request).unwrap();

    assert_eq!(result.l.len(), 2);
    assert_eq!(result.u.len(), 2);
    assert_eq!(result.p.len(), 2);
}

#[test]
fn test_lu_identity() {
    let request = LURequest {
        matrix: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
    };

    let result = lu_decomposition(request).unwrap();

    // L should be identity
    assert!((result.l[0][0] - 1.0).abs() < 1e-10);
    assert!((result.l[1][1] - 1.0).abs() < 1e-10);
    // U should be identity
    assert!((result.u[0][0] - 1.0).abs() < 1e-10);
    assert!((result.u[1][1] - 1.0).abs() < 1e-10);
}

#[test]
fn test_lu_non_square() {
    let request = LURequest {
        matrix: vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]],
    };

    let result = lu_decomposition(request);
    assert!(result.is_err());
}

// Test Schur decomposition
#[test]
fn test_schur_basic() {
    let request = SchurRequest {
        matrix: vec![vec![1.0, 2.0], vec![3.0, 4.0]],
    };

    let result = schur_decomposition(request).unwrap();

    assert_eq!(result.t.len(), 2);
    assert_eq!(result.u.len(), 2);
}

#[test]
fn test_schur_identity() {
    let request = SchurRequest {
        matrix: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
    };

    let result = schur_decomposition(request).unwrap();

    // T should be identity (already in Schur form)
    assert!((result.t[0][0] - 1.0).abs() < 1e-10);
    assert!((result.t[1][1] - 1.0).abs() < 1e-10);
}

#[test]
fn test_schur_non_square() {
    let request = SchurRequest {
        matrix: vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]],
    };

    let result = schur_decomposition(request);
    assert!(result.is_err());
}

// Test matrix inverse
#[test]
fn test_matrix_inverse_identity() {
    let matrix = vec![vec![1.0, 0.0], vec![0.0, 1.0]];

    let result = matrix_inverse(matrix.clone()).unwrap();

    assert_eq!(result, matrix);
}

#[test]
fn test_matrix_inverse_basic() {
    let matrix = vec![vec![1.0, 2.0], vec![3.0, 4.0]];

    let result = matrix_inverse(matrix).unwrap();

    // Inverse of [[1,2],[3,4]] is [[-2,1],[1.5,-0.5]]
    assert!((result[0][0] - (-2.0)).abs() < 1e-10);
    assert!((result[0][1] - 1.0).abs() < 1e-10);
    assert!((result[1][0] - 1.5).abs() < 1e-10);
    assert!((result[1][1] - (-0.5)).abs() < 1e-10);
}

#[test]
fn test_matrix_inverse_singular() {
    let matrix = vec![vec![1.0, 2.0], vec![2.0, 4.0]]; // Singular matrix

    let result = matrix_inverse(matrix);
    assert!(result.is_err());
}

#[test]
fn test_matrix_inverse_non_square() {
    let matrix = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];

    let result = matrix_inverse(matrix);
    assert!(result.is_err());
}
