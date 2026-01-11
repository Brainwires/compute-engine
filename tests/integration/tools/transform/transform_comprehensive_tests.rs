//! Comprehensive TRANSFORM tool test suite
//!
//! Tests for COMPUTE tool transform operations including:
//! - FFT (Forward, Inverse)
//! - Fourier (Continuous Forward, Inverse)
//! - Laplace (Forward, Inverse)
//! - Wavelets (Haar, Daubechies, Morlet, Mexican Hat)
//! - Filters (LowPass, HighPass, BandPass, BandStop)

use computational_engine::create_default_dispatcher;
use computational_engine::engine::equations::*;
use computational_engine::engine::*;
use serde_json::json;
use std::collections::HashMap;

// ============================================================================
// FFT TESTS
// ============================================================================

#[test]
fn test_fft_forward() {
    let dispatcher = create_default_dispatcher();

    // Simple sine wave
    let data: Vec<f64> = (0..64)
        .map(|i| (2.0 * std::f64::consts::PI * i as f64 / 64.0).sin())
        .collect();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Transform(TransformType::FFT(FFTType::Forward)),
        data: json!({
            "data": data,
            "sampling_rate": 1000.0
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform forward FFT: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        assert!(output.result.get("real").is_some(), "Should have real component");
        assert!(output.result.get("imag").is_some(), "Should have imaginary component");
    }
}

#[test]
fn test_fft_inverse() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 0.5, 0.25, 0.1];

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Transform(TransformType::FFT(FFTType::Inverse)),
        data: json!({
            "data": data
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform inverse FFT: {:?}", result);
}

// ============================================================================
// FOURIER TRANSFORM TESTS
// ============================================================================

#[test]
fn test_fourier_forward() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Transform(TransformType::Fourier(FourierTransform::Forward)),
        data: json!({
            "expression": "exp(-x^2)",
            "variable": "x"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform forward Fourier transform: {:?}", result);
}

#[test]
fn test_fourier_inverse() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Transform(TransformType::Fourier(FourierTransform::Inverse)),
        data: json!({
            "expression": "exp(-k^2)",
            "variable": "k"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform inverse Fourier transform: {:?}", result);
}

// ============================================================================
// LAPLACE TRANSFORM TESTS
// ============================================================================

#[test]
fn test_laplace_forward() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Transform(TransformType::Laplace(LaplaceTransform::Forward)),
        data: json!({
            "expression": "exp(-a*t)",
            "variable": "t"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform forward Laplace transform: {:?}", result);
}

#[test]
fn test_laplace_inverse() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Transform(TransformType::Laplace(LaplaceTransform::Inverse)),
        data: json!({
            "expression": "1/(s+a)",
            "variable": "s"
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform inverse Laplace transform: {:?}", result);
}

// ============================================================================
// WAVELET TRANSFORM TESTS
// ============================================================================

#[test]
fn test_wavelet_haar() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..64).map(|i| (i as f64 * 0.1).sin()).collect();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Transform(TransformType::Wavelet(WaveletType::Haar)),
        data: json!({
            "data": data
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform Haar wavelet transform: {:?}", result);
}

#[test]
fn test_wavelet_daubechies() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..64).map(|i| (i as f64 * 0.1).sin()).collect();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Transform(TransformType::Wavelet(WaveletType::Daubechies)),
        data: json!({
            "data": data,
            "order": 4
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform Daubechies wavelet transform: {:?}", result);
}

// ============================================================================
// FILTER TESTS
// ============================================================================

#[test]
fn test_filter_lowpass() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..128).map(|i| (i as f64 * 0.1).sin() + (i as f64 * 0.5).sin()).collect();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Transform(TransformType::Filter(FilterType::LowPass)),
        data: json!({
            "data": data,
            "cutoff": 0.2
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply lowpass filter: {:?}", result);
}

#[test]
fn test_filter_highpass() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..128).map(|i| (i as f64 * 0.1).sin() + (i as f64 * 0.5).sin()).collect();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Transform(TransformType::Filter(FilterType::HighPass)),
        data: json!({
            "data": data,
            "cutoff": 0.3
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply highpass filter: {:?}", result);
}

// ============================================================================
// WINDOW FUNCTION TESTS
// ============================================================================

#[test]
fn test_window_hamming() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Transform(TransformType::Window(WindowType::Hamming)),
        data: json!({
            "length": 64
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should generate Hamming window: {:?}", result);
}

#[test]
fn test_window_blackman() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Transform(TransformType::Window(WindowType::Blackman)),
        data: json!({
            "length": 64
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should generate Blackman window: {:?}", result);
}

// ============================================================================
// CONFORMAL MAPPING TESTS
// ============================================================================

#[test]
fn test_conformal_mapping() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Transform(TransformType::Conformal),
        data: json!({
            "mapping": "joukowski",
            "point": {"re": 1.0, "im": 0.5}
        }),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform conformal mapping: {:?}", result);
}
