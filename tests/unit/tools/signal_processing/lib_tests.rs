// Unit tests for compute::transforms (signal processing)
use computational_engine::signal_processing::*;
use std::f64::consts::PI;

    #[test]
    fn test_fft_sine_wave() {
        // Generate 10 Hz sine wave at 100 Hz sample rate
        let sample_rate = 100.0;
        let freq = 10.0;
        let duration = 1.0;
        let samples = (sample_rate * duration) as usize;

        let signal: Vec<f64> = (0..samples)
            .map(|i| (2.0 * PI * freq * i as f64 / sample_rate).sin())
            .collect();

        let result = compute_fft(FFTRequest {
            signal,
            sample_rate,
            window_type: "rectangular".to_string(),
        })
        .unwrap();

        assert_eq!(result.frequencies.len(), result.magnitudes.len());
        assert!(result.computation_time >= 0.0);

        // Peak should be near 10 Hz
        let peak_idx = result
            .magnitudes
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .unwrap();

        assert!((result.frequencies[peak_idx] - 10.0).abs() < 2.0);
    }

    #[test]
    fn test_fft_hanning_window() {
        let signal = vec![1.0; 64];
        let result = compute_fft(FFTRequest {
            signal,
            sample_rate: 100.0,
            window_type: "hanning".to_string(),
        })
        .unwrap();

        assert!(result.magnitudes.len() > 0);
    }

    #[test]
    fn test_fft_hamming_window() {
        let signal = vec![1.0; 64];
        let result = compute_fft(FFTRequest {
            signal,
            sample_rate: 100.0,
            window_type: "hamming".to_string(),
        })
        .unwrap();

        assert!(result.magnitudes.len() > 0);
    }

    #[test]
    fn test_fft_blackman_window() {
        let signal = vec![1.0; 64];
        let result = compute_fft(FFTRequest {
            signal,
            sample_rate: 100.0,
            window_type: "blackman".to_string(),
        })
        .unwrap();

        assert!(result.magnitudes.len() > 0);
    }

    #[test]
    fn test_lowpass_filter() {
        let signal: Vec<f64> = (0..100).map(|i| (i as f64).sin()).collect();

        let result = apply_filter(FilterRequest {
            signal: signal.clone(),
            filter_type: "lowpass".to_string(),
            cutoff_frequency: 5.0,
            sample_rate: 100.0,
            order: 2,
        })
        .unwrap();

        assert_eq!(result.filtered_signal.len(), signal.len());
    }

    #[test]
    fn test_highpass_filter() {
        let signal: Vec<f64> = (0..100).map(|i| (i as f64).sin()).collect();

        let result = apply_filter(FilterRequest {
            signal: signal.clone(),
            filter_type: "highpass".to_string(),
            cutoff_frequency: 5.0,
            sample_rate: 100.0,
            order: 2,
        })
        .unwrap();

        assert_eq!(result.filtered_signal.len(), signal.len());
    }

    #[test]
    fn test_bandpass_filter() {
        let signal: Vec<f64> = (0..100).map(|i| (i as f64).sin()).collect();

        let result = apply_filter(FilterRequest {
            signal: signal.clone(),
            filter_type: "bandpass".to_string(),
            cutoff_frequency: 10.0,
            sample_rate: 100.0,
            order: 4,
        })
        .unwrap();

        assert_eq!(result.filtered_signal.len(), signal.len());
    }

    #[test]
    fn test_bandstop_filter() {
        let signal: Vec<f64> = (0..100).map(|i| (i as f64).sin()).collect();

        let result = apply_filter(FilterRequest {
            signal: signal.clone(),
            filter_type: "bandstop".to_string(),
            cutoff_frequency: 10.0,
            sample_rate: 100.0,
            order: 4,
        })
        .unwrap();

        assert_eq!(result.filtered_signal.len(), signal.len());
    }

    #[test]
    fn test_invalid_filter_type() {
        let signal = vec![1.0; 10];

        let result = apply_filter(FilterRequest {
            signal,
            filter_type: "invalid".to_string(),
            cutoff_frequency: 5.0,
            sample_rate: 100.0,
            order: 2,
        });

        assert!(result.is_err());
    }

    #[test]
    fn test_spectrogram() {
        let signal: Vec<f64> = (0..1000).map(|i| (i as f64 * 0.1).sin()).collect();

        let result = compute_spectrogram(SpectrogramRequest {
            signal,
            sample_rate: 100.0,
            window_size: 128,
            hop_size: 64,
            window_type: "hanning".to_string(),
        })
        .unwrap();

        assert!(result.frequencies.len() > 0);
        assert!(result.times.len() > 0);
        assert_eq!(result.magnitudes.len(), result.times.len());
    }

    #[test]
    fn test_spectrogram_window_too_large() {
        let signal = vec![1.0; 10];

        let result = compute_spectrogram(SpectrogramRequest {
            signal,
            sample_rate: 100.0,
            window_size: 20,
            hop_size: 10,
            window_type: "hanning".to_string(),
        });

        assert!(result.is_err());
    }

    #[test]
    fn test_psd_welch() {
        let signal: Vec<f64> = (0..1000).map(|i| (i as f64 * 0.1).sin()).collect();

        let result = compute_psd(PSDRequest {
            signal,
            sample_rate: 100.0,
            method: "welch".to_string(),
            window_size: 256,
        })
        .unwrap();

        assert_eq!(result.frequencies.len(), result.psd.len());
        assert!(result.total_power > 0.0);
        assert!(result.frequency_resolution > 0.0);
    }

    #[test]
    fn test_psd_periodogram() {
        let signal: Vec<f64> = (0..256).map(|i| (i as f64 * 0.1).sin()).collect();

        let result = compute_psd(PSDRequest {
            signal,
            sample_rate: 100.0,
            method: "periodogram".to_string(),
            window_size: 256,
        })
        .unwrap();

        assert_eq!(result.frequencies.len(), result.psd.len());
        assert!(result.total_power > 0.0);
    }

    #[test]
    fn test_psd_invalid_method() {
        let signal = vec![1.0; 100];

        let result = compute_psd(PSDRequest {
            signal,
            sample_rate: 100.0,
            method: "invalid".to_string(),
            window_size: 64,
        });

        assert!(result.is_err());
    }

    // Note: apply_window tests removed because apply_window is a private function.
    // Window functionality is tested indirectly through compute_fft with different window types.

    #[test]
    fn test_filter_order_1() {
        let signal = vec![1.0; 20];
        let result = apply_filter(FilterRequest {
            signal,
            filter_type: "lowpass".to_string(),
            cutoff_frequency: 10.0,
            sample_rate: 100.0,
            order: 1,
        })
        .unwrap();

        assert!(result.filtered_signal.len() == 20);
    }

    #[test]
    fn test_filter_high_order() {
        let signal = vec![1.0; 20];
        let result = apply_filter(FilterRequest {
            signal,
            filter_type: "lowpass".to_string(),
            cutoff_frequency: 10.0,
            sample_rate: 100.0,
            order: 8,
        })
        .unwrap();

        assert!(result.filtered_signal.len() == 20);
    }

    #[test]
    fn test_fft_power_of_two_padding() {
        let signal = vec![1.0; 100]; // Not a power of 2
        let result = compute_fft(FFTRequest {
            signal,
            sample_rate: 100.0,
            window_type: "rectangular".to_string(),
        })
        .unwrap();

        // Should be padded to 128 (next power of 2)
        assert!(result.magnitudes.len() == 64); // Half of 128
    }

    #[test]
    fn test_spectrogram_different_windows() {
        let signal = vec![1.0; 512];

        for window_type in &["hanning", "hamming", "blackman", "rectangular"] {
            let result = compute_spectrogram(SpectrogramRequest {
                signal: signal.clone(),
                sample_rate: 100.0,
                window_size: 64,
                hop_size: 32,
                window_type: window_type.to_string(),
            })
            .unwrap();

            assert!(result.magnitudes.len() > 0);
        }
    }

    #[test]
    fn test_psd_welch_small_signal() {
        let signal = vec![1.0; 100];
        let result = compute_psd(PSDRequest {
            signal,
            sample_rate: 100.0,
            method: "welch".to_string(),
            window_size: 32,
        })
        .unwrap();

        assert!(result.frequencies.len() > 0);
    }
