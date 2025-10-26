//! Unit tests for signal_processing API handler

use crate::api::handlers::signal_processing::handle;
use crate::api::types::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};
use std::collections::HashMap;

/// Helper to create a request
fn create_request(operation: &str, parameters: HashMap<String, Value>) -> ComputationRequest {
    ComputationRequest {
        module: "signal_processing".to_string(),
        operation: operation.to_string(),
        parameters,
    }
}

/// Helper to extract result value
fn extract_result(response: ComputationResponse) -> Value {
    assert!(
        response.success,
        "Expected success but got error: {:?}",
        response.error
    );
    response.result.unwrap()
}

/// Helper to create a test signal (sine wave)
fn create_sine_wave(frequency: f64, sample_rate: f64, duration: f64) -> Vec<f64> {
    let num_samples = (duration * sample_rate) as usize;
    (0..num_samples)
        .map(|i| {
            let t = i as f64 / sample_rate;
            (2.0 * std::f64::consts::PI * frequency * t).sin()
        })
        .collect()
}

// ============================================================================
// FFT Tests
// ============================================================================

#[test]
fn test_fft_with_rectangular_window() {
    let signal = create_sine_wave(10.0, 100.0, 1.0);

    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!(signal));
    params.insert("sample_rate".to_string(), json!(100.0));
    params.insert("window_type".to_string(), json!("rectangular"));

    let request = create_request("fft", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert!(result["frequencies"].is_array());
    assert!(result["magnitudes"].is_array());
    assert!(result["phases"].is_array());
    assert!(result["computation_time"].is_number());

    let frequencies = result["frequencies"].as_array().unwrap();
    let magnitudes = result["magnitudes"].as_array().unwrap();
    assert_eq!(frequencies.len(), magnitudes.len());
}

#[test]
fn test_fft_with_hanning_window() {
    let signal = create_sine_wave(5.0, 50.0, 2.0);

    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!(signal));
    params.insert("sample_rate".to_string(), json!(50.0));
    params.insert("window_type".to_string(), json!("hanning"));

    let request = create_request("fft", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert!(result["frequencies"].is_array());
    assert!(result["magnitudes"].is_array());
    assert!(result["phases"].is_array());
}

#[test]
fn test_fft_with_hamming_window() {
    let signal = create_sine_wave(20.0, 200.0, 0.5);

    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!(signal));
    params.insert("sample_rate".to_string(), json!(200.0));
    params.insert("window_type".to_string(), json!("hamming"));

    let request = create_request("fft", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert!(result["magnitudes"].is_array());
}

#[test]
fn test_fft_with_blackman_window() {
    let signal = create_sine_wave(15.0, 100.0, 1.0);

    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!(signal));
    params.insert("sample_rate".to_string(), json!(100.0));
    params.insert("window_type".to_string(), json!("blackman"));

    let request = create_request("fft", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert!(result["frequencies"].is_array());
}

#[test]
fn test_fft_invalid_parameters() {
    let mut params = HashMap::new();
    params.insert("invalid_param".to_string(), json!("test"));

    let request = create_request("fft", params);
    let response = handle(&request);

    assert!(!response.success);
    assert!(response.error.is_some());
}

// ============================================================================
// Filter Tests
// ============================================================================

#[test]
fn test_lowpass_filter() {
    let signal = create_sine_wave(10.0, 100.0, 1.0);

    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!(signal.clone()));
    params.insert("filter_type".to_string(), json!("lowpass"));
    params.insert("cutoff_frequency".to_string(), json!(15.0));
    params.insert("sample_rate".to_string(), json!(100.0));
    params.insert("order".to_string(), json!(2));

    let request = create_request("filter", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert!(result["filtered_signal"].is_array());
    assert!(result["computation_time"].is_number());

    let filtered = result["filtered_signal"].as_array().unwrap();
    assert_eq!(filtered.len(), signal.len());
}

#[test]
fn test_highpass_filter() {
    let signal = create_sine_wave(10.0, 100.0, 1.0);

    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!(signal));
    params.insert("filter_type".to_string(), json!("highpass"));
    params.insert("cutoff_frequency".to_string(), json!(5.0));
    params.insert("sample_rate".to_string(), json!(100.0));
    params.insert("order".to_string(), json!(3));

    let request = create_request("filter", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert!(result["filtered_signal"].is_array());
}

#[test]
fn test_bandpass_filter() {
    let signal = create_sine_wave(10.0, 100.0, 1.0);

    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!(signal));
    params.insert("filter_type".to_string(), json!("bandpass"));
    params.insert("cutoff_frequency".to_string(), json!(10.0));
    params.insert("sample_rate".to_string(), json!(100.0));
    params.insert("order".to_string(), json!(4));

    let request = create_request("filter", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert!(result["filtered_signal"].is_array());
}

#[test]
fn test_bandstop_filter() {
    let signal = create_sine_wave(10.0, 100.0, 1.0);

    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!(signal));
    params.insert("filter_type".to_string(), json!("bandstop"));
    params.insert("cutoff_frequency".to_string(), json!(10.0));
    params.insert("sample_rate".to_string(), json!(100.0));
    params.insert("order".to_string(), json!(4));

    let request = create_request("filter", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert!(result["filtered_signal"].is_array());
}

#[test]
fn test_filter_unsupported_type() {
    let signal = create_sine_wave(10.0, 100.0, 1.0);

    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!(signal));
    params.insert("filter_type".to_string(), json!("notch")); // Unsupported
    params.insert("cutoff_frequency".to_string(), json!(10.0));
    params.insert("sample_rate".to_string(), json!(100.0));
    params.insert("order".to_string(), json!(2));

    let request = create_request("filter", params);
    let response = handle(&request);

    assert!(!response.success);
    assert!(response.error.unwrap().contains("Unsupported filter type"));
}

#[test]
fn test_filter_invalid_parameters() {
    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!([1.0, 2.0, 3.0]));
    params.insert("filter_type".to_string(), json!("lowpass"));
    // Missing cutoff_frequency and sample_rate

    let request = create_request("filter", params);
    let response = handle(&request);

    assert!(!response.success);
}

// ============================================================================
// Spectrogram Tests
// ============================================================================

#[test]
fn test_spectrogram_basic() {
    let signal = create_sine_wave(10.0, 100.0, 2.0);

    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!(signal));
    params.insert("sample_rate".to_string(), json!(100.0));
    params.insert("window_size".to_string(), json!(32));
    params.insert("hop_size".to_string(), json!(16));
    params.insert("window_type".to_string(), json!("hanning"));

    let request = create_request("spectrogram", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert!(result["frequencies"].is_array());
    assert!(result["times"].is_array());
    assert!(result["magnitudes"].is_array());
    assert!(result["computation_time"].is_number());
}

#[test]
fn test_spectrogram_with_rectangular_window() {
    let signal = create_sine_wave(20.0, 200.0, 1.0);

    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!(signal));
    params.insert("sample_rate".to_string(), json!(200.0));
    params.insert("window_size".to_string(), json!(64));
    params.insert("hop_size".to_string(), json!(32));
    params.insert("window_type".to_string(), json!("rectangular"));

    let request = create_request("spectrogram", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert!(result["magnitudes"].is_array());
    let magnitudes = result["magnitudes"].as_array().unwrap();
    assert!(!magnitudes.is_empty());
}

#[test]
fn test_spectrogram_window_size_too_large() {
    let signal = create_sine_wave(10.0, 100.0, 0.5);

    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!(signal.clone()));
    params.insert("sample_rate".to_string(), json!(100.0));
    params.insert("window_size".to_string(), json!(signal.len() + 10)); // Too large
    params.insert("hop_size".to_string(), json!(16));
    params.insert("window_type".to_string(), json!("hanning"));

    let request = create_request("spectrogram", params);
    let response = handle(&request);

    assert!(!response.success);
    assert!(response.error.unwrap().contains("Window size cannot be larger"));
}

#[test]
fn test_spectrogram_different_hop_size() {
    let signal = create_sine_wave(15.0, 150.0, 1.5);

    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!(signal));
    params.insert("sample_rate".to_string(), json!(150.0));
    params.insert("window_size".to_string(), json!(64));
    params.insert("hop_size".to_string(), json!(8)); // Smaller hop for more time resolution
    params.insert("window_type".to_string(), json!("hamming"));

    let request = create_request("spectrogram", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert!(result["times"].is_array());
    let times = result["times"].as_array().unwrap();
    assert!(times.len() > 1);
}

// ============================================================================
// Mock Operations Tests (IFFT, Convolution)
// ============================================================================

#[test]
fn test_ifft_mock_operation() {
    let mut params = HashMap::new();
    params.insert("dummy".to_string(), json!("data"));

    let request = create_request("ifft", params);
    let response = handle(&request);
    let result = extract_result(response);

    // Mock operation returns generic success
    assert_eq!(result["result"], "computed");
    assert!(result["output"].is_array());
}

#[test]
fn test_convolution_mock_operation() {
    let mut params = HashMap::new();
    params.insert("dummy".to_string(), json!("data"));

    let request = create_request("convolution", params);
    let response = handle(&request);
    let result = extract_result(response);

    // Mock operation returns generic success
    assert_eq!(result["result"], "computed");
    assert!(result["output"].is_array());
}

// ============================================================================
// Unknown Operation Tests
// ============================================================================

#[test]
fn test_unknown_operation() {
    let mut params = HashMap::new();
    params.insert("dummy".to_string(), json!("data"));

    let request = create_request("unknown_operation", params);
    let response = handle(&request);

    assert!(!response.success);
    assert!(response.error.unwrap().contains("Unknown operation"));
}

#[test]
fn test_invalid_operation_name() {
    let mut params = HashMap::new();
    params.insert("test".to_string(), json!(123));

    let request = create_request("not_a_real_operation", params);
    let response = handle(&request);

    assert!(!response.success);
}

// ============================================================================
// Edge Cases and Validation Tests
// ============================================================================

#[test]
fn test_fft_with_small_signal() {
    let signal = vec![1.0, 2.0, 3.0, 4.0];

    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!(signal));
    params.insert("sample_rate".to_string(), json!(100.0));
    params.insert("window_type".to_string(), json!("rectangular"));

    let request = create_request("fft", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert!(result["frequencies"].is_array());
    assert!(result["magnitudes"].is_array());
}

#[test]
fn test_filter_with_high_order() {
    let signal = create_sine_wave(10.0, 100.0, 1.0);

    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!(signal));
    params.insert("filter_type".to_string(), json!("lowpass"));
    params.insert("cutoff_frequency".to_string(), json!(20.0));
    params.insert("sample_rate".to_string(), json!(100.0));
    params.insert("order".to_string(), json!(10)); // High order

    let request = create_request("filter", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert!(result["filtered_signal"].is_array());
}

#[test]
fn test_spectrogram_with_blackman_window() {
    let signal = create_sine_wave(10.0, 100.0, 1.0);

    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!(signal));
    params.insert("sample_rate".to_string(), json!(100.0));
    params.insert("window_size".to_string(), json!(32));
    params.insert("hop_size".to_string(), json!(16));
    params.insert("window_type".to_string(), json!("blackman"));

    let request = create_request("spectrogram", params);
    let response = handle(&request);
    let result = extract_result(response);

    assert!(result["frequencies"].is_array());
}

#[test]
fn test_fft_preserves_signal_length_relationship() {
    let signal = create_sine_wave(5.0, 64.0, 1.0);

    let mut params = HashMap::new();
    params.insert("signal".to_string(), json!(signal.clone()));
    params.insert("sample_rate".to_string(), json!(64.0));
    params.insert("window_type".to_string(), json!("hanning"));

    let request = create_request("fft", params);
    let response = handle(&request);
    let result = extract_result(response);

    let frequencies = result["frequencies"].as_array().unwrap();
    let magnitudes = result["magnitudes"].as_array().unwrap();
    let phases = result["phases"].as_array().unwrap();

    // All output arrays should have the same length
    assert_eq!(frequencies.len(), magnitudes.len());
    assert_eq!(frequencies.len(), phases.len());
}
