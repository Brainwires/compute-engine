use computational_engine::signal_processing::*;

#[test]
fn test_fft() {
    let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];

    let result = compute_fft(FFTRequest {
        signal,
        sample_rate: 100.0,
        window_type: "hann".to_string(),
    }).unwrap();

    assert!(result.frequencies.len() > 0);
    assert_eq!(result.frequencies.len(), result.magnitudes.len());
    assert_eq!(result.frequencies.len(), result.phases.len());
}

#[test]
fn test_filter() {
    let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];

    let result = apply_filter(FilterRequest {
        signal,
        filter_type: "lowpass".to_string(),
        cutoff_frequency: 0.3,
        sample_rate: 100.0,
        order: 2,
    }).unwrap();

    assert!(result.filtered_signal.len() > 0);
}

// spectral_analysis function not implemented yet
// #[test]
// fn test_spectral_analysis() {
//     let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
//
//     let result = spectral_analysis(SpectralRequest {
//         signal,
//         sampling_rate: 100.0,
//         window: Some("hann".to_string()),
//     }).unwrap();
//
//     assert!(result.frequencies.len() > 0);
//     assert!(result.power_spectrum.len() > 0);
//     assert!(result.peak_frequency >= 0.0);
// }

#[test]
fn test_peak_detection() {
    let signal = vec![0.0, 1.0, 0.0, 2.0, 0.0, 3.0, 0.0];

    let result = detect_peaks(PeakDetectionRequest {
        signal,
        height: Some(0.5),
        distance: Some(1),
        prominence: None,
    }).unwrap();

    // Should detect 3 peaks at positions 1, 3, 5
    assert_eq!(result.peak_indices.len(), 3);
    assert_eq!(result.peak_heights.len(), 3);
}

// cross_correlation function not implemented yet
// #[test]
// fn test_cross_correlation() {
//     let signal1 = vec![1.0, 2.0, 3.0, 4.0];
//     let signal2 = vec![1.0, 2.0, 3.0, 4.0];
//
//     let result = cross_correlation(CrossCorrelationRequest {
//         signal1,
//         signal2,
//     }).unwrap();
//
//     assert!(result.correlation.len() > 0);
//     // Max correlation should be at zero lag for identical signals
//     assert!(result.max_correlation > 0.0);
// }

#[test]
fn test_wavelet_transform() {
    let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];

    let result = wavelet_transform(WaveletTransformRequest {
        signal,
        wavelet_type: "morlet".to_string(),
        scales: vec![1.0, 2.0, 4.0],
        sample_rate: 100.0,
    }).unwrap();

    assert!(result.coefficients.len() > 0);
    assert!(result.scales.len() > 0);
}
