use rustfft::{FftPlanner, num_complex::Complex};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use std::time::Instant;

#[derive(Deserialize, Serialize)]
pub struct FFTRequest {
    pub signal: Vec<f64>,
    pub sample_rate: f64,
    pub window_type: String,
}

#[derive(Deserialize, Serialize)]
pub struct FilterRequest {
    pub signal: Vec<f64>,
    pub filter_type: String,
    pub cutoff_frequency: f64,
    pub sample_rate: f64,
    pub order: usize,
}

#[derive(Deserialize, Serialize)]
pub struct SpectrogramRequest {
    pub signal: Vec<f64>,
    pub sample_rate: f64,
    pub window_size: usize,
    pub hop_size: usize,
    pub window_type: String,
}

#[derive(Deserialize, Serialize)]
pub struct PSDRequest {
    pub signal: Vec<f64>,
    pub sample_rate: f64,
    pub method: String,
    pub window_size: usize,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct FFTResult {
    pub frequencies: Vec<f64>,
    pub magnitudes: Vec<f64>,
    pub phases: Vec<f64>,
    pub computation_time: f64,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct FilterResult {
    pub filtered_signal: Vec<f64>,
    pub computation_time: f64,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct SpectrogramResult {
    pub frequencies: Vec<f64>,
    pub times: Vec<f64>,
    pub magnitudes: Vec<Vec<f64>>,
    pub computation_time: f64,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct PSDResult {
    pub frequencies: Vec<f64>,
    pub psd: Vec<f64>,
    pub total_power: f64,
    pub frequency_resolution: f64,
    pub computation_time: f64,
}

fn apply_window(signal: &[f64], window_type: &str) -> Vec<f64> {
    let n = signal.len();
    let windowed: Vec<f64> = match window_type {
        "hanning" => signal
            .iter()
            .enumerate()
            .map(|(i, &x)| {
                let window_val = 0.5 * (1.0 - (2.0 * PI * i as f64 / (n - 1) as f64).cos());
                x * window_val
            })
            .collect(),
        "hamming" => signal
            .iter()
            .enumerate()
            .map(|(i, &x)| {
                let window_val = 0.54 - 0.46 * (2.0 * PI * i as f64 / (n - 1) as f64).cos();
                x * window_val
            })
            .collect(),
        "blackman" => signal
            .iter()
            .enumerate()
            .map(|(i, &x)| {
                let window_val = 0.42 
                    - 0.5 * (2.0 * PI * i as f64 / (n - 1) as f64).cos()
                    + 0.08 * (4.0 * PI * i as f64 / (n - 1) as f64).cos();
                x * window_val
            })
            .collect(),
        _ => signal.to_vec(), // rectangular window (no modification)
    };
    windowed
}

pub fn compute_fft(request: FFTRequest) -> Result<FFTResult, String> {
    let start_time = Instant::now();
    
    // Apply window function
    let windowed_signal = apply_window(&request.signal, &request.window_type);
    
    // Pad to next power of 2 for efficiency
    let original_len = windowed_signal.len();
    let fft_len = original_len.next_power_of_two();
    
    // Convert to complex numbers and pad with zeros
    let mut complex_signal: Vec<Complex<f64>> = windowed_signal
        .into_iter()
        .map(|x| Complex::new(x, 0.0))
        .collect();
    complex_signal.resize(fft_len, Complex::new(0.0, 0.0));
    
    // Perform FFT
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(fft_len);
    fft.process(&mut complex_signal);
    
    // Extract results (only positive frequencies)
    let half_len = fft_len / 2;
    let frequencies: Vec<f64> = (0..half_len)
        .map(|i| i as f64 * request.sample_rate / fft_len as f64)
        .collect();
    
    let magnitudes: Vec<f64> = complex_signal[0..half_len]
        .iter()
        .map(|c| c.norm() / original_len as f64)
        .collect();
    
    let phases: Vec<f64> = complex_signal[0..half_len]
        .iter()
        .map(|c| c.arg())
        .collect();
    
    let computation_time = start_time.elapsed().as_secs_f64();
    
    Ok(FFTResult {
        frequencies,
        magnitudes,
        phases,
        computation_time,
    })
}

fn apply_lowpass_filter(signal: &[f64], cutoff: f64, sample_rate: f64, order: usize) -> Vec<f64> {
    // Simple RC low-pass filter cascade (approximation of Butterworth)
    let rc = 1.0 / (2.0 * PI * cutoff);
    let dt = 1.0 / sample_rate;
    let alpha = dt / (rc + dt);
    
    let mut filtered = signal.to_vec();
    
    // Apply multiple stages for higher order
    for _ in 0..order {
        let mut prev_output = filtered[0];
        for i in 1..filtered.len() {
            prev_output = prev_output + alpha * (filtered[i] - prev_output);
            filtered[i] = prev_output;
        }
    }
    
    filtered
}

fn apply_highpass_filter(signal: &[f64], cutoff: f64, sample_rate: f64, order: usize) -> Vec<f64> {
    // Simple RC high-pass filter cascade
    let rc = 1.0 / (2.0 * PI * cutoff);
    let dt = 1.0 / sample_rate;
    let alpha = rc / (rc + dt);
    
    let mut filtered = signal.to_vec();
    
    // Apply multiple stages for higher order
    for _ in 0..order {
        let mut prev_input = filtered[0];
        let mut prev_output = filtered[0];
        
        for i in 1..filtered.len() {
            let current_input = filtered[i];
            let current_output = alpha * (prev_output + current_input - prev_input);
            filtered[i] = current_output;
            prev_input = current_input;
            prev_output = current_output;
        }
    }
    
    filtered
}

pub fn apply_filter(request: FilterRequest) -> Result<FilterResult, String> {
    let start_time = Instant::now();
    
    let filtered_signal = match request.filter_type.as_str() {
        "lowpass" => apply_lowpass_filter(&request.signal, request.cutoff_frequency, request.sample_rate, request.order),
        "highpass" => apply_highpass_filter(&request.signal, request.cutoff_frequency, request.sample_rate, request.order),
        "bandpass" => {
            // Simple bandpass: high-pass followed by low-pass
            // For a proper bandpass, cutoff_frequency should be center frequency
            let bandwidth = request.cutoff_frequency * 0.2; // 20% bandwidth
            let low_cutoff = request.cutoff_frequency - bandwidth / 2.0;
            let high_cutoff = request.cutoff_frequency + bandwidth / 2.0;
            
            let high_passed = apply_highpass_filter(&request.signal, low_cutoff, request.sample_rate, request.order / 2);
            apply_lowpass_filter(&high_passed, high_cutoff, request.sample_rate, request.order / 2)
        },
        "bandstop" => {
            // Simple bandstop: parallel high-pass and low-pass, then add
            let bandwidth = request.cutoff_frequency * 0.2;
            let low_cutoff = request.cutoff_frequency - bandwidth / 2.0;
            let high_cutoff = request.cutoff_frequency + bandwidth / 2.0;
            
            let low_passed = apply_lowpass_filter(&request.signal, low_cutoff, request.sample_rate, request.order / 2);
            let high_passed = apply_highpass_filter(&request.signal, high_cutoff, request.sample_rate, request.order / 2);
            
            low_passed.iter().zip(high_passed.iter()).map(|(a, b)| a + b).collect()
        },
        _ => return Err(format!("Unsupported filter type: {}", request.filter_type)),
    };
    
    let computation_time = start_time.elapsed().as_secs_f64();
    
    Ok(FilterResult {
        filtered_signal,
        computation_time,
    })
}

pub fn compute_spectrogram(request: SpectrogramRequest) -> Result<SpectrogramResult, String> {
    let start_time = Instant::now();
    
    if request.window_size > request.signal.len() {
        return Err("Window size cannot be larger than signal length".to_string());
    }
    
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(request.window_size);
    
    let num_windows = (request.signal.len() - request.window_size) / request.hop_size + 1;
    let mut magnitudes = Vec::with_capacity(num_windows);
    let mut times = Vec::with_capacity(num_windows);
    
    let frequencies: Vec<f64> = (0..request.window_size / 2)
        .map(|i| i as f64 * request.sample_rate / request.window_size as f64)
        .collect();
    
    for window_idx in 0..num_windows {
        let start = window_idx * request.hop_size;
        let end = start + request.window_size;
        
        if end > request.signal.len() {
            break;
        }
        
        // Extract and window the signal segment
        let window_signal = &request.signal[start..end];
        let windowed = apply_window(window_signal, &request.window_type);
        
        // Convert to complex and perform FFT
        let mut complex_signal: Vec<Complex<f64>> = windowed
            .into_iter()
            .map(|x| Complex::new(x, 0.0))
            .collect();
        
        fft.process(&mut complex_signal);
        
        // Extract magnitudes (only positive frequencies)
        let window_magnitudes: Vec<f64> = complex_signal[0..request.window_size / 2]
            .iter()
            .map(|c| c.norm() / request.window_size as f64)
            .collect();
        
        magnitudes.push(window_magnitudes);
        times.push(start as f64 / request.sample_rate);
    }
    
    let computation_time = start_time.elapsed().as_secs_f64();
    
    Ok(SpectrogramResult {
        frequencies,
        times,
        magnitudes,
        computation_time,
    })
}

fn compute_psd_welch(signal: &[f64], sample_rate: f64, window_size: usize) -> (Vec<f64>, Vec<f64>, f64) {
    let overlap = window_size / 2;
    let hop_size = window_size - overlap;
    let num_windows = (signal.len() - overlap) / hop_size;
    
    if num_windows == 0 {
        return (vec![], vec![], 0.0);
    }
    
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(window_size);
    
    let frequencies: Vec<f64> = (0..window_size / 2)
        .map(|i| i as f64 * sample_rate / window_size as f64)
        .collect();
    
    let mut psd_sum = vec![0.0; window_size / 2];
    
    for window_idx in 0..num_windows {
        let start = window_idx * hop_size;
        let end = start + window_size;
        
        if end > signal.len() {
            break;
        }
        
        // Extract and window the signal segment
        let window_signal = &signal[start..end];
        let windowed = apply_window(window_signal, "hanning");
        
        // Convert to complex and perform FFT
        let mut complex_signal: Vec<Complex<f64>> = windowed
            .into_iter()
            .map(|x| Complex::new(x, 0.0))
            .collect();
        
        fft.process(&mut complex_signal);
        
        // Accumulate power spectral density
        for (i, c) in complex_signal[0..window_size / 2].iter().enumerate() {
            let magnitude_squared = c.norm_sqr();
            psd_sum[i] += magnitude_squared;
        }
    }
    
    // Average and normalize
    let psd: Vec<f64> = psd_sum
        .into_iter()
        .map(|sum| sum / (num_windows as f64 * sample_rate))
        .collect();
    
    let total_power = psd.iter().sum::<f64>() * (sample_rate / window_size as f64);
    
    (frequencies, psd, total_power)
}

pub fn compute_psd(request: PSDRequest) -> Result<PSDResult, String> {
    let start_time = Instant::now();

    let (frequencies, psd, total_power) = match request.method.as_str() {
        "welch" => compute_psd_welch(&request.signal, request.sample_rate, request.window_size),
        "periodogram" => {
            // Simple periodogram (single FFT)
            let windowed = apply_window(&request.signal, "hanning");
            let mut complex_signal: Vec<Complex<f64>> = windowed
                .into_iter()
                .map(|x| Complex::new(x, 0.0))
                .collect();

            let fft_len = complex_signal.len();
            let mut planner = FftPlanner::new();
            let fft = planner.plan_fft_forward(fft_len);
            fft.process(&mut complex_signal);

            let frequencies: Vec<f64> = (0..fft_len / 2)
                .map(|i| i as f64 * request.sample_rate / fft_len as f64)
                .collect();

            let psd: Vec<f64> = complex_signal[0..fft_len / 2]
                .iter()
                .map(|c| c.norm_sqr() / (fft_len as f64 * request.sample_rate))
                .collect();

            let total_power = psd.iter().sum::<f64>() * (request.sample_rate / fft_len as f64);

            (frequencies, psd, total_power)
        },
        _ => return Err(format!("Unsupported PSD method: {}", request.method)),
    };

    let frequency_resolution = request.sample_rate / request.window_size as f64;
    let computation_time = start_time.elapsed().as_secs_f64();

    Ok(PSDResult {
        frequencies,
        psd,
        total_power,
        frequency_resolution,
        computation_time,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

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
        }).unwrap();

        assert_eq!(result.frequencies.len(), result.magnitudes.len());
        assert!(result.computation_time >= 0.0);

        // Peak should be near 10 Hz
        let peak_idx = result.magnitudes.iter()
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
        }).unwrap();

        assert!(result.magnitudes.len() > 0);
    }

    #[test]
    fn test_fft_hamming_window() {
        let signal = vec![1.0; 64];
        let result = compute_fft(FFTRequest {
            signal,
            sample_rate: 100.0,
            window_type: "hamming".to_string(),
        }).unwrap();

        assert!(result.magnitudes.len() > 0);
    }

    #[test]
    fn test_fft_blackman_window() {
        let signal = vec![1.0; 64];
        let result = compute_fft(FFTRequest {
            signal,
            sample_rate: 100.0,
            window_type: "blackman".to_string(),
        }).unwrap();

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
        }).unwrap();

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
        }).unwrap();

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
        }).unwrap();

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
        }).unwrap();

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
        }).unwrap();

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
        }).unwrap();

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
        }).unwrap();

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

    #[test]
    fn test_apply_window_rectangular() {
        let signal = vec![1.0; 10];
        let windowed = apply_window(&signal, "rectangular");
        assert_eq!(windowed.len(), signal.len());
        assert_eq!(windowed, signal);
    }

    #[test]
    fn test_apply_window_hanning() {
        let signal = vec![1.0; 10];
        let windowed = apply_window(&signal, "hanning");
        assert_eq!(windowed.len(), signal.len());
        // Edges should be attenuated
        assert!(windowed[0] < 0.5);
        assert!(windowed[9] < 0.5);
    }

    #[test]
    fn test_filter_order_1() {
        let signal = vec![1.0; 20];
        let result = apply_filter(FilterRequest {
            signal,
            filter_type: "lowpass".to_string(),
            cutoff_frequency: 10.0,
            sample_rate: 100.0,
            order: 1,
        }).unwrap();

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
        }).unwrap();

        assert!(result.filtered_signal.len() == 20);
    }

    #[test]
    fn test_fft_power_of_two_padding() {
        let signal = vec![1.0; 100]; // Not a power of 2
        let result = compute_fft(FFTRequest {
            signal,
            sample_rate: 100.0,
            window_type: "rectangular".to_string(),
        }).unwrap();

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
            }).unwrap();

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
        }).unwrap();

        assert!(result.frequencies.len() > 0);
    }
}

