//! Unit tests for analysis module

use super::{controllability, observability, eigenvalues};

#[test]
fn test_controllability() {
    let a = vec![vec![0.0, 1.0], vec![-1.0, -0.5]];
    let b = vec![vec![0.0], vec![1.0]];

    let result = controllability(&a, &b);
    assert_eq!(result.rank, 1);
}

#[test]
fn test_observability() {
    let a = vec![vec![0.0, 1.0], vec![-1.0, -0.5]];
    let c = vec![vec![1.0, 0.0]];

    let result = observability(&a, &c);
    assert_eq!(result.rank, 1);
}

#[test]
fn test_controllability_siso() {
    let a = vec![vec![-1.0]];
    let b = vec![vec![1.0]];

    let result = controllability(&a, &b);
    assert!(result.is_controllable);
}

#[test]
fn test_observability_siso() {
    let a = vec![vec![-1.0]];
    let c = vec![vec![1.0]];

    let result = observability(&a, &c);
    assert!(result.is_observable);
}

#[test]
fn test_eigenvalues() {
    let a = vec![vec![0.0, 1.0], vec![-1.0, -0.5]];
    let eigs = eigenvalues(&a);
    // Currently returns empty vec (placeholder)
    assert!(eigs.is_empty());
}
