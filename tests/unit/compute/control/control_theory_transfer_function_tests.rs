//! Unit tests for transfer_function module

use super::TransferFunction;
use num_complex::Complex64;

#[test]
fn test_transfer_function_creation() {
    // G(s) = 1 / (s + 1)
    let tf = TransferFunction::new(vec![1.0], vec![1.0, 1.0]);
    assert!(tf.is_ok());
}

#[test]
fn test_transfer_function_evaluate() {
    let tf = TransferFunction::new(vec![1.0], vec![1.0, 1.0]).unwrap();
    let s = Complex64::new(0.0, 1.0);
    let result = tf.evaluate(s);
    assert!(result.norm() > 0.0);
}

#[test]
fn test_frequency_response() {
    let tf = TransferFunction::new(vec![1.0], vec![1.0, 1.0]).unwrap();
    let h = tf.frequency_response(1.0);
    assert!(h.norm() > 0.0);
}

#[test]
fn test_bode_plot() {
    let tf = TransferFunction::new(vec![1.0], vec![1.0, 1.0]).unwrap();
    let bode = tf.bode_plot(0.1, 100.0, 50);
    assert_eq!(bode.frequencies.len(), 50);
    assert_eq!(bode.magnitudes.len(), 50);
    assert_eq!(bode.phases.len(), 50);
}

#[test]
fn test_dc_gain() {
    // G(s) = 2 / (s + 1) has DC gain = 2
    let tf = TransferFunction::new(vec![2.0], vec![1.0, 1.0]).unwrap();
    let gain = tf.dc_gain();
    assert!((gain - 2.0).abs() < 1e-10);
}
