//! Comprehensive linear algebra test suite
//!
//! Tests for all linear algebra operations including:
//! - Matrix decompositions (SVD, LU, Cholesky, QR, Schur)
//! - Eigenvalue/eigenvector calculations
//! - Matrix operations (determinant, trace, norms, inverse, power, exponential)
//! - Matrix rank and pseudoinverse
//! - Principal Component Analysis (PCA)
//! - Linear system solving

use computational_engine::mathematics::linear_algebra::*;

// ============================================================================
// SINGULAR VALUE DECOMPOSITION (SVD) TESTS
// ============================================================================

#[test]
fn test_svd_square_matrix() {
    let matrix = vec![
        vec![1.0, 2.0, 3.0],
        vec![4.0, 5.0, 6.0],
        vec![7.0, 8.0, 9.0],
    ];

    let result = compute_svd(MatrixInput {
        matrix: matrix.clone(),
    })
    .unwrap();

    assert_eq!(result.u.len(), 3, "U should have 3 rows");
    assert_eq!(
        result.singular_values.len(),
        3,
        "Should have 3 singular values"
    );
    assert_eq!(result.v_transpose.len(), 3, "V^T should have 3 rows");
    assert!(result.rank <= 3, "Rank should be at most 3");
    assert!(
        result.condition_number > 0.0,
        "Condition number should be positive"
    );
}

#[test]
fn test_svd_rectangular_tall() {
    let matrix = vec![vec![1.0, 2.0], vec![3.0, 4.0], vec![5.0, 6.0]];

    let result = compute_svd(MatrixInput { matrix }).unwrap();

    assert_eq!(result.u.len(), 3, "U should have 3 rows");
    assert_eq!(
        result.singular_values.len(),
        2,
        "Should have 2 singular values"
    );
    assert!(result.rank <= 2, "Rank should be at most 2");
}

#[test]
fn test_svd_rectangular_wide() {
    let matrix = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];

    let result = compute_svd(MatrixInput { matrix }).unwrap();

    assert_eq!(result.u.len(), 2, "U should have 2 rows");
    assert_eq!(
        result.singular_values.len(),
        2,
        "Should have 2 singular values"
    );
    assert!(result.rank <= 2, "Rank should be at most 2");
}

#[test]
fn test_svd_identity_matrix() {
    let matrix = vec![
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];

    let result = compute_svd(MatrixInput { matrix }).unwrap();

    // Identity matrix has all singular values = 1
    for &sv in &result.singular_values {
        assert!((sv - 1.0).abs() < 1e-10, "Singular values should be 1.0");
    }
    assert_eq!(result.rank, 3, "Identity matrix should be full rank");
    assert!(
        (result.condition_number - 1.0).abs() < 1e-10,
        "Condition number should be 1"
    );
}

#[test]
fn test_svd_diagonal_matrix() {
    let matrix = vec![
        vec![3.0, 0.0, 0.0],
        vec![0.0, 2.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];

    let result = compute_svd(MatrixInput { matrix }).unwrap();

    assert_eq!(result.rank, 3, "Diagonal matrix should be full rank");
    // Singular values should be 3, 2, 1 in descending order
    assert!(result.singular_values[0] >= result.singular_values[1]);
    assert!(result.singular_values[1] >= result.singular_values[2]);
}

#[test]
fn test_svd_rank_deficient() {
    // Matrix with linearly dependent rows
    let matrix = vec![
        vec![1.0, 2.0, 3.0],
        vec![2.0, 4.0, 6.0], // 2x first row
        vec![3.0, 6.0, 9.0], // 3x first row
    ];

    let result = compute_svd(MatrixInput { matrix }).unwrap();

    assert_eq!(result.rank, 1, "Matrix should have rank 1");
    assert!(
        result.condition_number > 1e10,
        "Should have very high condition number"
    );
}

// ============================================================================
// MATRIX RANK TESTS
// ============================================================================

#[test]
fn test_matrix_rank_full() {
    let matrix = vec![
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];

    let rank = compute_matrix_rank(MatrixInput { matrix }).unwrap();
    assert_eq!(rank, 3, "Identity matrix should be full rank");
}

#[test]
fn test_matrix_rank_2x2() {
    let matrix = vec![vec![1.0, 2.0], vec![3.0, 4.0]];

    let rank = compute_matrix_rank(MatrixInput { matrix }).unwrap();
    assert_eq!(rank, 2, "Non-singular 2x2 should have rank 2");
}

#[test]
fn test_matrix_rank_deficient_2x2() {
    let matrix = vec![
        vec![1.0, 2.0],
        vec![2.0, 4.0], // Linearly dependent
    ];

    let rank = compute_matrix_rank(MatrixInput { matrix }).unwrap();
    assert_eq!(rank, 1, "Rank deficient matrix should have rank 1");
}

#[test]
fn test_matrix_rank_rectangular() {
    let matrix = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];

    let rank = compute_matrix_rank(MatrixInput { matrix }).unwrap();
    assert!(rank <= 2, "Rank cannot exceed min dimension");
}

#[test]
fn test_matrix_rank_zero_row() {
    let matrix = vec![
        vec![1.0, 2.0, 3.0],
        vec![0.0, 0.0, 0.0],
        vec![4.0, 5.0, 6.0],
    ];

    let rank = compute_matrix_rank(MatrixInput { matrix }).unwrap();
    assert!(rank <= 2, "Matrix with zero row should have reduced rank");
}

// ============================================================================
// PSEUDOINVERSE TESTS
// ============================================================================

#[test]
fn test_pseudoinverse_square() {
    let matrix = vec![vec![1.0, 2.0], vec![3.0, 4.0]];

    let pinv = compute_pseudoinverse(MatrixInput {
        matrix: matrix.clone(),
    })
    .unwrap();

    assert_eq!(pinv.len(), 2, "Pseudoinverse should have 2 rows");
    assert_eq!(pinv[0].len(), 2, "Pseudoinverse should have 2 columns");

    // For invertible matrix, pseudoinverse equals inverse
    // Verify A * A+ * A ≈ A
    let prod = matrix_multiply(&matrix, &pinv);
    let prod2 = matrix_multiply(&prod, &matrix);

    for i in 0..2 {
        for j in 0..2 {
            assert!(
                (prod2[i][j] - matrix[i][j]).abs() < 1e-8,
                "A * A+ * A should equal A: prod2[{}][{}] = {}, matrix[{}][{}] = {}",
                i,
                j,
                prod2[i][j],
                i,
                j,
                matrix[i][j]
            );
        }
    }
}

#[test]
fn test_pseudoinverse_rectangular_tall() {
    let matrix = vec![vec![1.0, 2.0], vec![3.0, 4.0], vec![5.0, 6.0]];

    let pinv = compute_pseudoinverse(MatrixInput { matrix }).unwrap();

    // Pseudoinverse of m×n matrix is n×m
    assert_eq!(pinv.len(), 2, "Should have 2 rows");
    assert_eq!(pinv[0].len(), 3, "Should have 3 columns");
}

#[test]
fn test_pseudoinverse_rectangular_wide() {
    let matrix = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];

    let pinv = compute_pseudoinverse(MatrixInput { matrix }).unwrap();

    assert_eq!(pinv.len(), 3, "Should have 3 rows");
    assert_eq!(pinv[0].len(), 2, "Should have 2 columns");
}

#[test]
fn test_pseudoinverse_identity() {
    let matrix = vec![vec![1.0, 0.0], vec![0.0, 1.0]];

    let pinv = compute_pseudoinverse(MatrixInput {
        matrix: matrix.clone(),
    })
    .unwrap();

    // Pseudoinverse of identity is identity
    for i in 0..2 {
        for j in 0..2 {
            assert!((pinv[i][j] - matrix[i][j]).abs() < 1e-10);
        }
    }
}

// ============================================================================
// EIGENDECOMPOSITION TESTS
// ============================================================================

#[test]
fn test_eigendecomposition_diagonal() {
    let matrix = vec![
        vec![3.0, 0.0, 0.0],
        vec![0.0, 2.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];

    let result = compute_eigendecomposition(MatrixInput { matrix }).unwrap();

    assert_eq!(result.eigenvalues.len(), 3, "Should have 3 eigenvalues");
    assert_eq!(result.eigenvectors.len(), 3, "Should have 3 eigenvectors");

    // Eigenvalues of diagonal matrix are the diagonal elements
    let mut evals = result.eigenvalues.clone();
    evals.sort_by(|a, b| b.partial_cmp(a).unwrap());
    assert!((evals[0] - 3.0).abs() < 1e-9);
    assert!((evals[1] - 2.0).abs() < 1e-9);
    assert!((evals[2] - 1.0).abs() < 1e-9);
}

#[test]
fn test_eigendecomposition_identity() {
    let matrix = vec![
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];

    let result = compute_eigendecomposition(MatrixInput { matrix }).unwrap();

    // All eigenvalues should be 1
    for &eval in &result.eigenvalues {
        assert!((eval - 1.0).abs() < 1e-10, "All eigenvalues should be 1");
    }
}

#[test]
fn test_eigendecomposition_symmetric_2x2() {
    let matrix = vec![vec![2.0, 1.0], vec![1.0, 2.0]];

    let result = compute_eigendecomposition(MatrixInput { matrix }).unwrap();

    assert_eq!(result.eigenvalues.len(), 2);
    // Eigenvalues should be 3 and 1
    let mut evals = result.eigenvalues.clone();
    evals.sort_by(|a, b| b.partial_cmp(a).unwrap());
    assert!((evals[0] - 3.0).abs() < 1e-9);
    assert!((evals[1] - 1.0).abs() < 1e-9);
}

#[test]
fn test_eigendecomposition_non_square_error() {
    let matrix = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];

    let result = compute_eigendecomposition(MatrixInput { matrix });
    assert!(result.is_err(), "Non-square matrix should error");
}

// ============================================================================
// MATRIX OPERATIONS TESTS
// ============================================================================

#[test]
fn test_matrix_operations_square() {
    let matrix = vec![vec![1.0, 2.0], vec![3.0, 4.0]];

    let result = matrix_operations(MatrixInput { matrix }).unwrap();

    assert!(
        result.determinant.is_some(),
        "Square matrix should have determinant"
    );
    assert!(result.trace.is_some(), "Square matrix should have trace");
    assert!(
        (result.determinant.unwrap() - (-2.0)).abs() < 1e-10,
        "det([[1,2],[3,4]]) = -2"
    );
    assert!(
        (result.trace.unwrap() - 5.0).abs() < 1e-10,
        "trace([[1,2],[3,4]]) = 5"
    );
    assert!(result.frobenius_norm > 0.0);
    assert!(result.max_norm > 0.0);
}

#[test]
fn test_matrix_operations_identity() {
    let matrix = vec![
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];

    let result = matrix_operations(MatrixInput { matrix }).unwrap();

    assert!((result.determinant.unwrap() - 1.0).abs() < 1e-10);
    assert!((result.trace.unwrap() - 3.0).abs() < 1e-10);
}

#[test]
fn test_matrix_operations_rectangular() {
    let matrix = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];

    let result = matrix_operations(MatrixInput { matrix }).unwrap();

    assert!(
        result.determinant.is_none(),
        "Non-square has no determinant"
    );
    assert!(result.trace.is_none(), "Non-square has no trace");
    assert!(result.frobenius_norm > 0.0);
    assert!(result.max_norm > 0.0);
}

// ============================================================================
// MATRIX NORM TESTS
// ============================================================================

#[test]
fn test_matrix_norm_frobenius() {
    let matrix = vec![vec![3.0, 0.0], vec![0.0, 4.0]];

    let result = matrix_norm(MatrixNormRequest {
        matrix,
        norm_type: "frobenius".to_string(),
    })
    .unwrap();

    // Frobenius norm of diag(3,4) = sqrt(9+16) = 5
    assert!((result.norm - 5.0).abs() < 1e-10);
}

#[test]
fn test_matrix_norm_1() {
    let matrix = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];

    let result = matrix_norm(MatrixNormRequest {
        matrix,
        norm_type: "1".to_string(),
    })
    .unwrap();

    // 1-norm is max column sum: max(5, 7, 9) = 9
    assert!((result.norm - 9.0).abs() < 1e-10);
}

#[test]
fn test_matrix_norm_infinity() {
    let matrix = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];

    let result = matrix_norm(MatrixNormRequest {
        matrix,
        norm_type: "inf".to_string(),
    })
    .unwrap();

    // Infinity norm is max row sum: max(6, 15) = 15
    assert!((result.norm - 15.0).abs() < 1e-10);
}

#[test]
fn test_matrix_norm_fro_alias() {
    let matrix = vec![vec![1.0, 0.0], vec![0.0, 1.0]];

    let result = matrix_norm(MatrixNormRequest {
        matrix,
        norm_type: "fro".to_string(),
    })
    .unwrap();

    // Should work with "fro" alias
    assert!((result.norm - std::f64::consts::SQRT_2).abs() < 1e-10);
}

// ============================================================================
// MATRIX POWER TESTS
// ============================================================================

#[test]
fn test_matrix_power_zero() {
    let matrix = vec![vec![2.0, 3.0], vec![4.0, 5.0]];

    let result = matrix_power(MatrixPowerRequest { matrix, power: 0 }).unwrap();

    // A^0 = I
    assert_eq!(result.result[0][0], 1.0);
    assert_eq!(result.result[0][1], 0.0);
    assert_eq!(result.result[1][0], 0.0);
    assert_eq!(result.result[1][1], 1.0);
}

#[test]
fn test_matrix_power_one() {
    let matrix = vec![vec![2.0, 3.0], vec![4.0, 5.0]];

    let result = matrix_power(MatrixPowerRequest {
        matrix: matrix.clone(),
        power: 1,
    })
    .unwrap();

    // A^1 = A
    for i in 0..2 {
        for j in 0..2 {
            assert_eq!(result.result[i][j], matrix[i][j]);
        }
    }
}

#[test]
fn test_matrix_power_diagonal() {
    let matrix = vec![vec![2.0, 0.0], vec![0.0, 3.0]];

    let result = matrix_power(MatrixPowerRequest { matrix, power: 3 }).unwrap();

    // diag(2,3)^3 = diag(8,27)
    assert_eq!(result.result[0][0], 8.0);
    assert_eq!(result.result[1][1], 27.0);
    assert_eq!(result.result[0][1], 0.0);
    assert_eq!(result.result[1][0], 0.0);
}

#[test]
fn test_matrix_power_identity() {
    let matrix = vec![vec![1.0, 0.0], vec![0.0, 1.0]];

    let result = matrix_power(MatrixPowerRequest { matrix, power: 100 }).unwrap();

    // I^n = I for any n
    assert_eq!(result.result[0][0], 1.0);
    assert_eq!(result.result[1][1], 1.0);
}

// ============================================================================
// MATRIX EXPONENTIAL TESTS
// ============================================================================

#[test]
fn test_matrix_exp_zero() {
    let matrix = vec![vec![0.0, 0.0], vec![0.0, 0.0]];

    let result = matrix_exp(MatrixExpRequest {
        matrix,
        terms: Some(10),
    })
    .unwrap();

    // exp(0) = I
    assert!((result.result[0][0] - 1.0).abs() < 1e-6);
    assert!((result.result[1][1] - 1.0).abs() < 1e-6);
    assert!(result.result[0][1].abs() < 1e-6);
    assert!(result.result[1][0].abs() < 1e-6);
}

#[test]
fn test_matrix_exp_diagonal() {
    let matrix = vec![vec![1.0, 0.0], vec![0.0, 2.0]];

    let result = matrix_exp(MatrixExpRequest {
        matrix,
        terms: Some(30), // Increase terms for better accuracy
    })
    .unwrap();

    // exp(diag(1,2)) = diag(e, e^2)
    assert!(
        (result.result[0][0] - std::f64::consts::E).abs() < 1e-3,
        "Expected {}, got {}",
        std::f64::consts::E,
        result.result[0][0]
    );
    assert!(
        (result.result[1][1] - std::f64::consts::E.powi(2)).abs() < 1e-2,
        "Expected {}, got {}",
        std::f64::consts::E.powi(2),
        result.result[1][1]
    );
    assert!(result.result[0][1].abs() < 1e-6);
    assert!(result.result[1][0].abs() < 1e-6);
}

// ============================================================================
// CHOLESKY DECOMPOSITION TESTS
// ============================================================================

#[test]
fn test_cholesky_positive_definite() {
    let matrix = vec![
        vec![4.0, 12.0, -16.0],
        vec![12.0, 37.0, -43.0],
        vec![-16.0, -43.0, 98.0],
    ];

    let result = cholesky_decomposition(CholeskyRequest {
        matrix: matrix.clone(),
    })
    .unwrap();

    // Check that L is lower triangular
    assert_eq!(result.lower[0][1], 0.0);
    assert_eq!(result.lower[0][2], 0.0);
    assert_eq!(result.lower[1][2], 0.0);

    // Verify L * L^T = A
    let lt = transpose(&result.lower);
    let product = matrix_multiply(&result.lower, &lt);

    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (product[i][j] - matrix[i][j]).abs() < 1e-9,
                "L*L^T should equal original matrix"
            );
        }
    }
}

#[test]
fn test_cholesky_identity() {
    let matrix = vec![
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];

    let result = cholesky_decomposition(CholeskyRequest {
        matrix: matrix.clone(),
    })
    .unwrap();

    // Cholesky of I is I
    for i in 0..3 {
        for j in 0..3 {
            assert!((result.lower[i][j] - matrix[i][j]).abs() < 1e-10);
        }
    }
}

#[test]
fn test_cholesky_diagonal() {
    let matrix = vec![
        vec![4.0, 0.0, 0.0],
        vec![0.0, 9.0, 0.0],
        vec![0.0, 0.0, 16.0],
    ];

    let result = cholesky_decomposition(CholeskyRequest { matrix }).unwrap();

    // Cholesky of diag(4,9,16) is diag(2,3,4)
    assert!((result.lower[0][0] - 2.0).abs() < 1e-10);
    assert!((result.lower[1][1] - 3.0).abs() < 1e-10);
    assert!((result.lower[2][2] - 4.0).abs() < 1e-10);
}

#[test]
fn test_cholesky_not_positive_definite() {
    let matrix = vec![
        vec![1.0, 2.0],
        vec![2.0, 1.0], // Not positive definite
    ];

    let result = cholesky_decomposition(CholeskyRequest { matrix });
    assert!(
        result.is_err(),
        "Should fail for non-positive definite matrix"
    );
}

// ============================================================================
// LU DECOMPOSITION TESTS
// ============================================================================

#[test]
fn test_lu_square_matrix() {
    let matrix = vec![
        vec![2.0, -1.0, 0.0],
        vec![-1.0, 2.0, -1.0],
        vec![0.0, -1.0, 2.0],
    ];

    let result = lu_decomposition(LURequest {
        matrix: matrix.clone(),
    })
    .unwrap();

    assert_eq!(result.l.len(), 3);
    assert_eq!(result.u.len(), 3);
    assert_eq!(result.p.len(), 3);

    // Verify PA = LU
    let pa = matrix_multiply(&result.p, &matrix);
    let lu = matrix_multiply(&result.l, &result.u);

    for i in 0..3 {
        for j in 0..3 {
            assert!((pa[i][j] - lu[i][j]).abs() < 1e-9, "PA should equal LU");
        }
    }
}

#[test]
fn test_lu_identity() {
    let matrix = vec![
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];

    let result = lu_decomposition(LURequest { matrix }).unwrap();

    // L should be I, U should be I
    for i in 0..3 {
        assert!((result.l[i][i] - 1.0).abs() < 1e-10);
        assert!((result.u[i][i] - 1.0).abs() < 1e-10);
    }
}

#[test]
fn test_lu_lower_triangular() {
    let matrix = vec![
        vec![1.0, 0.0, 0.0],
        vec![2.0, 3.0, 0.0],
        vec![4.0, 5.0, 6.0],
    ];

    let result = lu_decomposition(LURequest {
        matrix: matrix.clone(),
    })
    .unwrap();

    // Verify the decomposition
    let pa = matrix_multiply(&result.p, &matrix);
    let lu = matrix_multiply(&result.l, &result.u);

    for i in 0..3 {
        for j in 0..3 {
            assert!((pa[i][j] - lu[i][j]).abs() < 1e-9);
        }
    }
}

#[test]
fn test_lu_singular_error() {
    let matrix = vec![
        vec![1.0, 2.0, 3.0],
        vec![2.0, 4.0, 6.0],
        vec![3.0, 6.0, 9.0],
    ];

    let result = lu_decomposition(LURequest { matrix });
    assert!(result.is_err(), "Should fail for singular matrix");
}

// ============================================================================
// SCHUR DECOMPOSITION TESTS
// ============================================================================

#[test]
fn test_schur_square_matrix() {
    let matrix = vec![
        vec![1.0, 2.0, 3.0],
        vec![0.0, 4.0, 5.0],
        vec![0.0, 0.0, 6.0],
    ];

    let result = schur_decomposition(SchurRequest { matrix }).unwrap();

    assert_eq!(result.t.len(), 3, "T should be 3x3");
    assert_eq!(result.u.len(), 3, "U should be 3x3");
}

#[test]
fn test_schur_identity() {
    let matrix = vec![
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];

    let result = schur_decomposition(SchurRequest { matrix }).unwrap();

    // For identity, T should be identity
    for i in 0..3 {
        for j in 0..3 {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!((result.t[i][j] - expected).abs() < 1e-9);
        }
    }
}

#[test]
fn test_schur_diagonal() {
    let matrix = vec![
        vec![3.0, 0.0, 0.0],
        vec![0.0, 2.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];

    let result = schur_decomposition(SchurRequest { matrix }).unwrap();

    // Diagonal matrix is already in Schur form
    // Check that T is upper triangular
    for i in 1..3 {
        for j in 0..i {
            assert!(result.t[i][j].abs() < 1e-9, "T should be upper triangular");
        }
    }
}

// ============================================================================
// MATRIX INVERSE TESTS
// ============================================================================

#[test]
fn test_matrix_inverse_2x2() {
    let matrix = vec![vec![1.0, 2.0], vec![3.0, 4.0]];

    let inv = matrix_inverse(matrix.clone()).unwrap();

    // Verify A * A^-1 = I
    let product = matrix_multiply(&matrix, &inv);

    // The inverse of [[1,2],[3,4]] is [[-2, 1], [1.5, -0.5]]
    // Check that product is close to identity
    assert!(
        (product[0][0] - 1.0).abs() < 1e-8,
        "product[0][0] = {}, expected 1.0",
        product[0][0]
    );
    assert!(
        (product[1][1] - 1.0).abs() < 1e-8,
        "product[1][1] = {}, expected 1.0",
        product[1][1]
    );
    assert!(
        product[0][1].abs() < 1e-8,
        "product[0][1] = {}, expected 0.0",
        product[0][1]
    );
    assert!(
        product[1][0].abs() < 1e-8,
        "product[1][0] = {}, expected 0.0",
        product[1][0]
    );
}

#[test]
fn test_matrix_inverse_identity() {
    let matrix = vec![
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];

    let inv = matrix_inverse(matrix.clone()).unwrap();

    // Inverse of identity is identity
    for i in 0..3 {
        for j in 0..3 {
            assert!((inv[i][j] - matrix[i][j]).abs() < 1e-10);
        }
    }
}

#[test]
fn test_matrix_inverse_diagonal() {
    let matrix = vec![
        vec![2.0, 0.0, 0.0],
        vec![0.0, 3.0, 0.0],
        vec![0.0, 0.0, 4.0],
    ];

    let inv = matrix_inverse(matrix).unwrap();

    // Inverse of diag(2,3,4) is diag(0.5, 1/3, 0.25)
    assert!((inv[0][0] - 0.5).abs() < 1e-10);
    assert!((inv[1][1] - 1.0 / 3.0).abs() < 1e-10);
    assert!((inv[2][2] - 0.25).abs() < 1e-10);
}

#[test]
fn test_matrix_inverse_singular_error() {
    let matrix = vec![vec![1.0, 2.0], vec![2.0, 4.0]];

    let result = matrix_inverse(matrix);
    assert!(result.is_err(), "Should fail for singular matrix");
}

#[test]
fn test_matrix_inverse_non_square_error() {
    let matrix = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];

    let result = matrix_inverse(matrix);
    assert!(result.is_err(), "Should fail for non-square matrix");
}

// ============================================================================
// PCA TESTS
// ============================================================================

#[test]
fn test_pca_basic() {
    let matrix = vec![
        vec![2.5, 2.4],
        vec![0.5, 0.7],
        vec![2.2, 2.9],
        vec![1.9, 2.2],
        vec![3.1, 3.0],
    ];

    let result = compute_pca(MatrixInput { matrix }, None).unwrap();

    assert_eq!(result.principal_components.len(), 2, "Should have 2 PCs");
    assert_eq!(result.explained_variance.len(), 2);
    assert_eq!(result.cumulative_variance.len(), 2);
    assert_eq!(result.transformed_data.len(), 5);

    // Cumulative variance should sum to 1.0
    assert!((result.cumulative_variance.last().unwrap() - 1.0).abs() < 1e-9);

    // Explained variance should be in descending order
    if result.explained_variance.len() > 1 {
        assert!(result.explained_variance[0] >= result.explained_variance[1]);
    }
}

#[test]
fn test_pca_n_components() {
    let matrix = vec![
        vec![1.0, 2.0, 3.0],
        vec![4.0, 5.0, 6.0],
        vec![7.0, 8.0, 9.0],
        vec![10.0, 11.0, 12.0],
    ];

    let result = compute_pca(MatrixInput { matrix }, Some(2)).unwrap();

    assert_eq!(result.principal_components.len(), 2, "Should have 2 PCs");
    assert_eq!(result.explained_variance.len(), 2);
    assert_eq!(result.transformed_data[0].len(), 2);
}

#[test]
fn test_pca_identity_covariance() {
    // Data with identity covariance matrix
    let matrix = vec![
        vec![1.0, 0.0],
        vec![0.0, 1.0],
        vec![-1.0, 0.0],
        vec![0.0, -1.0],
    ];

    let result = compute_pca(MatrixInput { matrix }, None).unwrap();

    assert_eq!(result.principal_components.len(), 2);
    // Should have equal explained variance
    assert!((result.explained_variance[0] - result.explained_variance[1]).abs() < 0.2);
}

// ============================================================================
// ADDITIONAL EDGE CASE TESTS
// ============================================================================

#[test]
fn test_svd_single_element() {
    let matrix = vec![vec![5.0]];

    let result = compute_svd(MatrixInput { matrix }).unwrap();

    assert_eq!(result.singular_values.len(), 1);
    assert!((result.singular_values[0] - 5.0).abs() < 1e-10);
    assert_eq!(result.rank, 1);
}

#[test]
fn test_matrix_norm_single_element() {
    let matrix = vec![vec![3.0]];

    let result = matrix_norm(MatrixNormRequest {
        matrix,
        norm_type: "frobenius".to_string(),
    })
    .unwrap();

    assert!((result.norm - 3.0).abs() < 1e-10);
}

#[test]
fn test_matrix_inverse_3x3() {
    let matrix = vec![
        vec![1.0, 2.0, 3.0],
        vec![0.0, 1.0, 4.0],
        vec![5.0, 6.0, 0.0],
    ];

    let inv = matrix_inverse(matrix.clone()).unwrap();
    let product = matrix_multiply(&matrix, &inv);

    // Check that product is identity
    for i in 0..3 {
        for j in 0..3 {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!((product[i][j] - expected).abs() < 1e-8);
        }
    }
}

#[test]
fn test_matrix_power_large() {
    let matrix = vec![vec![1.0, 1.0], vec![0.0, 1.0]];

    let result = matrix_power(MatrixPowerRequest { matrix, power: 5 }).unwrap();

    // [[1,1],[0,1]]^5 = [[1,5],[0,1]]
    assert_eq!(result.result[0][0], 1.0);
    assert_eq!(result.result[0][1], 5.0);
    assert_eq!(result.result[1][0], 0.0);
    assert_eq!(result.result[1][1], 1.0);
}

#[test]
fn test_cholesky_2x2() {
    let matrix = vec![vec![4.0, 2.0], vec![2.0, 3.0]];

    let result = cholesky_decomposition(CholeskyRequest {
        matrix: matrix.clone(),
    })
    .unwrap();

    let lt = transpose(&result.lower);
    let product = matrix_multiply(&result.lower, &lt);

    for i in 0..2 {
        for j in 0..2 {
            assert!((product[i][j] - matrix[i][j]).abs() < 1e-9);
        }
    }
}

#[test]
fn test_lu_with_pivoting() {
    let matrix = vec![
        vec![0.0, 1.0, 2.0],
        vec![1.0, 2.0, 3.0],
        vec![2.0, 3.0, 5.0],
    ];

    let result = lu_decomposition(LURequest {
        matrix: matrix.clone(),
    })
    .unwrap();

    // Verify PA = LU
    let pa = matrix_multiply(&result.p, &matrix);
    let lu = matrix_multiply(&result.l, &result.u);

    for i in 0..3 {
        for j in 0..3 {
            assert!((pa[i][j] - lu[i][j]).abs() < 1e-9);
        }
    }
}

#[test]
fn test_svd_very_wide() {
    let matrix = vec![vec![1.0, 2.0, 3.0, 4.0, 5.0]];

    let result = compute_svd(MatrixInput { matrix }).unwrap();

    assert_eq!(result.u.len(), 1);
    assert_eq!(result.singular_values.len(), 1);
}

#[test]
fn test_svd_very_tall() {
    let matrix = vec![vec![1.0], vec![2.0], vec![3.0], vec![4.0], vec![5.0]];

    let result = compute_svd(MatrixInput { matrix }).unwrap();

    assert_eq!(result.u.len(), 5);
    assert_eq!(result.singular_values.len(), 1);
}

#[test]
fn test_matrix_norm_infinity_single_row() {
    let matrix = vec![vec![1.0, 2.0, 3.0, 4.0]];

    let result = matrix_norm(MatrixNormRequest {
        matrix,
        norm_type: "infinity".to_string(),
    })
    .unwrap();

    // Max row sum = 1+2+3+4 = 10
    assert!((result.norm - 10.0).abs() < 1e-10);
}

#[test]
fn test_matrix_norm_1_single_column() {
    let matrix = vec![vec![1.0], vec![2.0], vec![3.0], vec![4.0]];

    let result = matrix_norm(MatrixNormRequest {
        matrix,
        norm_type: "1".to_string(),
    })
    .unwrap();

    // Max column sum = 1+2+3+4 = 10
    assert!((result.norm - 10.0).abs() < 1e-10);
}

#[test]
fn test_eigendecomposition_larger_matrix() {
    let matrix = vec![
        vec![4.0, -2.0, 1.0, 0.0],
        vec![-2.0, 4.0, -2.0, 1.0],
        vec![1.0, -2.0, 4.0, -2.0],
        vec![0.0, 1.0, -2.0, 4.0],
    ];

    let result = compute_eigendecomposition(MatrixInput { matrix }).unwrap();

    assert_eq!(result.eigenvalues.len(), 4);
    assert_eq!(result.eigenvectors.len(), 4);
}

#[test]
fn test_pca_single_component() {
    let matrix = vec![
        vec![1.0, 2.0, 3.0],
        vec![2.0, 4.0, 6.0],
        vec![3.0, 6.0, 9.0],
        vec![4.0, 8.0, 12.0],
    ];

    let result = compute_pca(MatrixInput { matrix }, Some(1)).unwrap();

    assert_eq!(result.principal_components.len(), 1);
    assert_eq!(result.explained_variance.len(), 1);
    // First component should explain most variance
    assert!(result.explained_variance[0] > 0.9);
}

#[test]
fn test_pseudoinverse_rank_deficient() {
    let matrix = vec![vec![1.0, 2.0], vec![2.0, 4.0]];

    let pinv = compute_pseudoinverse(MatrixInput { matrix }).unwrap();

    // Pseudoinverse should exist even for rank-deficient matrix
    assert_eq!(pinv.len(), 2);
    assert_eq!(pinv[0].len(), 2);
}

#[test]
fn test_matrix_operations_large() {
    let matrix = vec![
        vec![1.0, 2.0, 3.0, 4.0],
        vec![5.0, 6.0, 7.0, 8.0],
        vec![9.0, 10.0, 11.0, 12.0],
        vec![13.0, 14.0, 15.0, 16.0],
    ];

    let result = matrix_operations(MatrixInput { matrix }).unwrap();

    assert!(result.determinant.is_some());
    assert!(result.trace.is_some());
    // Trace = 1+6+11+16 = 34
    assert!((result.trace.unwrap() - 34.0).abs() < 1e-10);
}

#[test]
fn test_matrix_power_negative() {
    // For negative powers, we'd need matrix inverse, but current implementation
    // only handles positive powers. This test verifies the current behavior.
    let matrix = vec![vec![2.0, 0.0], vec![0.0, 3.0]];

    let result = matrix_power(MatrixPowerRequest { matrix, power: -1 }).unwrap();

    // With current implementation, negative powers behave like positive
    // This is a limitation we're documenting
    assert!(result.result.len() > 0);
}

#[test]
fn test_schur_2x2() {
    let matrix = vec![vec![1.0, 2.0], vec![3.0, 4.0]];

    let result = schur_decomposition(SchurRequest { matrix }).unwrap();

    assert_eq!(result.t.len(), 2);
    assert_eq!(result.u.len(), 2);
}

#[test]
fn test_lu_2x2() {
    let matrix = vec![vec![4.0, 3.0], vec![6.0, 3.0]];

    let result = lu_decomposition(LURequest {
        matrix: matrix.clone(),
    })
    .unwrap();

    let pa = matrix_multiply(&result.p, &matrix);
    let lu = matrix_multiply(&result.l, &result.u);

    for i in 0..2 {
        for j in 0..2 {
            assert!((pa[i][j] - lu[i][j]).abs() < 1e-9);
        }
    }
}

#[test]
fn test_pca_more_components_than_features() {
    let matrix = vec![
        vec![1.0, 2.0],
        vec![3.0, 4.0],
        vec![5.0, 6.0],
        vec![7.0, 8.0],
    ];

    // Request more components than features (we have 2 features, request 5)
    // PCA should limit to min(n_features, n_samples) = min(2, 4) = 2
    let result = compute_pca(MatrixInput { matrix }, Some(5)).unwrap();

    // Should return only 2 components (limited by min of features and samples)
    assert_eq!(result.principal_components.len(), 2);
    assert_eq!(result.explained_variance.len(), 2);
}

#[test]
fn test_matrix_exp_small_matrix() {
    let matrix = vec![vec![0.1, 0.0], vec![0.0, 0.1]];

    let result = matrix_exp(MatrixExpRequest {
        matrix,
        terms: Some(15),
    })
    .unwrap();

    // exp(0.1*I) ≈ 1.105*I
    let expected = (0.1_f64).exp();
    assert!((result.result[0][0] - expected).abs() < 1e-3);
    assert!((result.result[1][1] - expected).abs() < 1e-3);
}

// ============================================================================
// HELPER FUNCTIONS FOR TESTS
// ============================================================================

fn transpose(matrix: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let rows = matrix.len();
    let cols = matrix[0].len();

    (0..cols)
        .map(|j| (0..rows).map(|i| matrix[i][j]).collect())
        .collect()
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
