//! Unit tests for lqr module

use super::lqr;

#[test]
fn test_lqr_basic() {
    let a = vec![vec![0.0, 1.0], vec![-1.0, -0.5]];
    let b = vec![vec![0.0], vec![1.0]];
    let q = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    let r = vec![vec![1.0]];

    let result = lqr(&a, &b, &q, &r);
    assert!(result.is_ok());
}

#[test]
fn test_lqr_result_dimensions() {
    let a = vec![vec![0.0, 1.0], vec![-1.0, -0.5]];
    let b = vec![vec![0.0], vec![1.0]];
    let q = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    let r = vec![vec![1.0]];

    let result = lqr(&a, &b, &q, &r).unwrap();
    assert_eq!(result.k.len(), 1); // m rows (num inputs)
    assert_eq!(result.k[0].len(), 2); // n columns (num states)
    assert_eq!(result.p.len(), 2); // n x n
    assert_eq!(result.p[0].len(), 2);
}

#[test]
fn test_lqr_siso_system() {
    let a = vec![vec![-1.0]];
    let b = vec![vec![1.0]];
    let q = vec![vec![1.0]];
    let r = vec![vec![1.0]];

    let result = lqr(&a, &b, &q, &r);
    assert!(result.is_ok());
}
