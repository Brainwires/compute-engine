use computational_engine::linear_algebra::*;

// compute_qr function not implemented yet
// #[test]
// fn test_qr_decomposition() {
//     let matrix = vec![
//         vec![12.0, -51.0, 4.0],
//         vec![6.0, 167.0, -68.0],
//         vec![-4.0, 24.0, -41.0],
//     ];
//
//     let result = compute_qr(MatrixInput { matrix }).unwrap();
//
//     assert_eq!(result.q.len(), 3);
//     assert_eq!(result.r.len(), 3);
//
//     // Q should be orthogonal (Q^T Q = I)
//     // Just check dimensions for now
//     assert_eq!(result.q[0].len(), 3);
// }

#[test]
fn test_matrix_rank() {
    let matrix = vec![
        vec![1.0, 2.0, 3.0],
        vec![4.0, 5.0, 6.0],
        vec![7.0, 8.0, 9.0],
    ];

    let rank = compute_matrix_rank(MatrixInput { matrix }).unwrap();

    // This matrix has rank 2
    assert!(rank >= 2);
}

#[test]
fn test_matrix_norm() {
    let matrix = vec![
        vec![3.0, 0.0],
        vec![0.0, 4.0],
    ];

    let result = matrix_norm(MatrixNormRequest {
        matrix,
        norm_type: "frobenius".to_string(),
    }).unwrap();

    // Frobenius norm of diag(3,4) = sqrt(9+16) = 5
    assert!((result.norm - 5.0).abs() < 1e-10);
}

#[test]
fn test_matrix_power() {
    let matrix = vec![
        vec![2.0, 0.0],
        vec![0.0, 3.0],
    ];

    let result = matrix_power(MatrixPowerRequest {
        matrix,
        power: 2,
    }).unwrap();

    // [2,0; 0,3]^2 = [4,0; 0,9]
    assert_eq!(result.result[0][0], 4.0);
    assert_eq!(result.result[1][1], 9.0);
}

#[test]
fn test_cholesky() {
    let matrix = vec![
        vec![4.0, 12.0, -16.0],
        vec![12.0, 37.0, -43.0],
        vec![-16.0, -43.0, 98.0],
    ];

    let result = cholesky_decomposition(CholeskyRequest { matrix }).unwrap();

    // Check that L is lower triangular
    assert_eq!(result.lower[0][1], 0.0);
    assert_eq!(result.lower[0][2], 0.0);
    assert_eq!(result.lower[1][2], 0.0);
}

#[test]
fn test_pca() {
    let matrix = vec![
        vec![2.5, 2.4],
        vec![0.5, 0.7],
        vec![2.2, 2.9],
        vec![1.9, 2.2],
    ];

    let result = compute_pca(MatrixInput { matrix }, Some(1)).unwrap();

    assert!(result.explained_variance.len() > 0);
    assert!(result.principal_components.len() > 0);
}
