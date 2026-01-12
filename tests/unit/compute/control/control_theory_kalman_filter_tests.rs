//! Unit tests for kalman_filter module

use super::KalmanFilter;

#[test]
fn test_kalman_filter_creation() {
    let kf = KalmanFilter::new(
        vec![0.0, 0.0],
        vec![vec![1.0, 0.0], vec![0.0, 1.0]],
        vec![vec![0.1, 0.0], vec![0.0, 0.1]],
        vec![vec![1.0]],
    );
    assert_eq!(kf.x.len(), 2);
}

#[test]
fn test_kalman_predict() {
    let mut kf = KalmanFilter::new(
        vec![1.0, 0.0],
        vec![vec![1.0, 0.0], vec![0.0, 1.0]],
        vec![vec![0.1, 0.0], vec![0.0, 0.1]],
        vec![vec![1.0]],
    );

    let a = vec![vec![1.0, 0.1], vec![0.0, 1.0]];
    kf.predict(&a);
    assert_eq!(kf.x.len(), 2);
}

#[test]
fn test_kalman_update() {
    let mut kf = KalmanFilter::new(
        vec![1.0, 0.0],
        vec![vec![1.0, 0.0], vec![0.0, 1.0]],
        vec![vec![0.1, 0.0], vec![0.0, 0.1]],
        vec![vec![1.0]],
    );

    let z = vec![1.5];
    let h = vec![vec![1.0, 0.0]];
    kf.update(&z, &h);
    assert_eq!(kf.x.len(), 2);
}

#[test]
fn test_kalman_predict_update_cycle() {
    let mut kf = KalmanFilter::new(
        vec![0.0, 0.0],
        vec![vec![1.0, 0.0], vec![0.0, 1.0]],
        vec![vec![0.01, 0.0], vec![0.0, 0.01]],
        vec![vec![0.1]],
    );

    let a = vec![vec![1.0, 0.1], vec![0.0, 1.0]];
    let h = vec![vec![1.0, 0.0]];

    // Predict
    kf.predict(&a);
    // Update
    kf.update(&[1.0], &h);

    assert!(kf.x[0].is_finite());
}
