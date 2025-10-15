//! Comprehensive TRANSFORM tool test suite
//!
//! Tests for all TRANSFORM operations including:
//! - FFT (Forward, Inverse)
//! - Fourier (Continuous Forward, Inverse)
//! - Laplace (Forward, Inverse)
//! - Wavelets (Haar, Daubechies, Morlet, Mexican Hat)
//! - Filters (LowPass, HighPass, BandPass, BandStop)
//! - Window Functions (Hamming, Hanning, Blackman, Kaiser)
//! - Conformal Mappings

use computational_engine::engine::*;
use computational_engine::create_default_dispatcher;
use std::collections::HashMap;

// ============================================================================
// FFT TESTS
// ============================================================================

#[test]
fn test_fft_forward() {
    let dispatcher = create_default_dispatcher();

    // Simple sine wave
    let data: Vec<f64> = (0..64).map(|i| (2.0 * std::f64::consts::PI * i as f64 / 64.0).sin()).collect();

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::FFT(FFTType::Forward),
        data,
        sampling_rate: Some(1000.0),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform forward FFT: {:?}", result);

    if let Ok(ToolResponse::Transform(output)) = result {
        assert!(output.frequencies.is_some(), "Should have frequencies");
        assert!(output.magnitude.is_some(), "Should have magnitude spectrum");
    }
}

#[test]
fn test_fft_inverse() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 0.5, 0.25, 0.1];
    let mut params = HashMap::new();
    params.insert("frequency_data".to_string(), serde_json::json!(data.clone()));

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::FFT(FFTType::Inverse),
        data,
        sampling_rate: Some(1000.0),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform inverse FFT: {:?}", result);
}

// ============================================================================
// FOURIER TRANSFORM TESTS (Continuous)
// ============================================================================

#[test]
fn test_fourier_forward() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 2.0, 3.0, 4.0];

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Fourier(FourierTransform::Forward),
        data,
        sampling_rate: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should handle Fourier forward transform: {:?}", result);
}

#[test]
fn test_fourier_inverse() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 0.5, 0.25, 0.125];

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Fourier(FourierTransform::Inverse),
        data,
        sampling_rate: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should handle Fourier inverse transform: {:?}", result);
}

// ============================================================================
// LAPLACE TRANSFORM TESTS
// ============================================================================

#[test]
fn test_laplace_forward() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 0.9, 0.8, 0.7]; // Decaying signal

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Laplace(LaplaceTransform::Forward),
        data,
        sampling_rate: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should handle Laplace forward transform: {:?}", result);
}

#[test]
fn test_laplace_inverse() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 0.5, 0.25];

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Laplace(LaplaceTransform::Inverse),
        data,
        sampling_rate: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should handle Laplace inverse transform: {:?}", result);
}

// ============================================================================
// WAVELET TRANSFORM TESTS
// ============================================================================

#[test]
fn test_wavelet_haar() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..32).map(|i| (i as f64 / 32.0).sin()).collect();

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Wavelet(WaveletType::Haar),
        data,
        sampling_rate: Some(1000.0),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute Haar wavelet transform: {:?}", result);
}

#[test]
fn test_wavelet_daubechies() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..32).map(|i| (i as f64 / 32.0).cos()).collect();

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Wavelet(WaveletType::Daubechies),
        data,
        sampling_rate: Some(1000.0),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute Daubechies wavelet transform: {:?}", result);
}

#[test]
fn test_wavelet_morlet() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..64).map(|i| {
        let t = i as f64 / 64.0;
        (2.0 * std::f64::consts::PI * 5.0 * t).sin()
    }).collect();

    let mut params = HashMap::new();
    params.insert("scales".to_string(), serde_json::json!([1.0, 2.0, 4.0, 8.0]));

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Wavelet(WaveletType::Morlet),
        data,
        sampling_rate: Some(1000.0),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute Morlet wavelet transform: {:?}", result);
}

#[test]
fn test_wavelet_mexican_hat() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..64).map(|i| {
        let t = i as f64 / 64.0;
        (-t * t).exp()
    }).collect();

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Wavelet(WaveletType::Mexican),
        data,
        sampling_rate: Some(1000.0),
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute Mexican Hat wavelet transform: {:?}", result);
}

// ============================================================================
// FILTER TESTS
// ============================================================================

#[test]
fn test_filter_lowpass() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..100).map(|i| {
        (2.0 * std::f64::consts::PI * 5.0 * i as f64 / 100.0).sin() +
        0.5 * (2.0 * std::f64::consts::PI * 50.0 * i as f64 / 100.0).sin()
    }).collect();

    let mut params = HashMap::new();
    params.insert("cutoff".to_string(), serde_json::json!(10.0));

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Filter(FilterType::LowPass),
        data,
        sampling_rate: Some(1000.0),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply lowpass filter: {:?}", result);
}

#[test]
fn test_filter_highpass() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..100).map(|i| {
        (2.0 * std::f64::consts::PI * 5.0 * i as f64 / 100.0).sin()
    }).collect();

    let mut params = HashMap::new();
    params.insert("cutoff".to_string(), serde_json::json!(50.0));

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Filter(FilterType::HighPass),
        data,
        sampling_rate: Some(1000.0),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply highpass filter: {:?}", result);
}

#[test]
fn test_filter_bandpass() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..100).map(|i| (i as f64).sin()).collect();

    let mut params = HashMap::new();
    params.insert("cutoff".to_string(), serde_json::json!(20.0));
    params.insert("bandwidth".to_string(), serde_json::json!(10.0));

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Filter(FilterType::BandPass),
        data,
        sampling_rate: Some(1000.0),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply bandpass filter: {:?}", result);
}

#[test]
fn test_filter_bandstop() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..100).map(|i| (i as f64 / 10.0).cos()).collect();

    let mut params = HashMap::new();
    params.insert("cutoff".to_string(), serde_json::json!(30.0));

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Filter(FilterType::BandStop),
        data,
        sampling_rate: Some(1000.0),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply bandstop filter: {:?}", result);
}

// ============================================================================
// WINDOW FUNCTION TESTS
// ============================================================================

#[test]
fn test_window_hamming() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..64).map(|i| 1.0).collect();

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Window(WindowType::Hamming),
        data,
        sampling_rate: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply Hamming window: {:?}", result);

    if let Ok(ToolResponse::Transform(output)) = result {
        assert_eq!(output.result.len(), 64, "Should have same length as input");
        // Window should taper at edges
        assert!(output.result[0] < 1.0, "Window should taper at start");
        assert!(output.result[63] < 1.0, "Window should taper at end");
    }
}

#[test]
fn test_window_hanning() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..64).map(|_| 1.0).collect();

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Window(WindowType::Hanning),
        data,
        sampling_rate: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply Hanning window: {:?}", result);
}

#[test]
fn test_window_blackman() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..64).map(|_| 1.0).collect();

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Window(WindowType::Blackman),
        data,
        sampling_rate: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply Blackman window: {:?}", result);
}

#[test]
fn test_window_kaiser() {
    let dispatcher = create_default_dispatcher();

    let data: Vec<f64> = (0..64).map(|_| 1.0).collect();

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Window(WindowType::Kaiser),
        data,
        sampling_rate: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should apply Kaiser window: {:?}", result);
}

// ============================================================================
// CONFORMAL MAPPING TESTS
// ============================================================================

#[test]
fn test_conformal_mobius() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 0.5, 0.0, -0.5, -1.0];
    let mut params = HashMap::new();
    params.insert("mapping".to_string(), serde_json::json!("mobius"));

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Conformal,
        data,
        sampling_rate: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should handle Möbius conformal mapping: {:?}", result);
}

#[test]
fn test_conformal_joukowski() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.5, 1.0, 0.5];
    let mut params = HashMap::new();
    params.insert("mapping".to_string(), serde_json::json!("joukowski"));

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Conformal,
        data,
        sampling_rate: None,
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should handle Joukowski conformal mapping: {:?}", result);
}

// ============================================================================
// ERROR HANDLING TESTS
// ============================================================================

#[test]
fn test_fft_requires_sampling_rate() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 2.0, 3.0, 4.0];

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::FFT(FFTType::Forward),
        data,
        sampling_rate: None, // Missing required parameter
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "FFT should require sampling_rate");
}

#[test]
fn test_filter_requires_cutoff() {
    let dispatcher = create_default_dispatcher();

    let data = vec![1.0, 2.0, 3.0, 4.0];

    let request = ToolRequest::Transform(TransformInput {
        transform_type: TransformType::Filter(FilterType::LowPass),
        data,
        sampling_rate: Some(1000.0),
        parameters: HashMap::new(), // Missing cutoff parameter
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Filter should require cutoff frequency");
}
