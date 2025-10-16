use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Deserialize, Serialize)]
pub struct PeakDetectionRequest {
    pub signal: Vec<f64>,
    pub height: Option<f64>,
    pub distance: Option<usize>,
    pub prominence: Option<f64>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct PeakDetectionResult {
    pub peak_indices: Vec<usize>,
    pub peak_heights: Vec<f64>,
    pub peak_prominences: Vec<f64>,
    pub total_peaks_found: usize,
}

pub fn detect_peaks(request: PeakDetectionRequest) -> Result<PeakDetectionResult, String> {
    if request.signal.is_empty() {
        return Err("Signal cannot be empty".to_string());
    }

    let mut peaks = Vec::new();
    let mut peak_heights = Vec::new();

    // Find local maxima
    for i in 1..request.signal.len() - 1 {
        let current = request.signal[i];
        let prev = request.signal[i - 1];
        let next = request.signal[i + 1];

        if current > prev && current > next {
            // Apply height filter if specified
            if let Some(min_height) = request.height {
                if current < min_height {
                    continue;
                }
            }

            // Apply distance filter
            if let Some(min_distance) = request.distance {
                if min_distance > 1 && !peaks.is_empty() {
                    let last_peak = peaks[peaks.len() - 1];
                    if i - last_peak < min_distance {
                        continue;
                    }
                }
            }

            peaks.push(i);
            peak_heights.push(current);
        }
    }

    // Calculate prominence
    let prominences: Vec<f64> = peaks
        .iter()
        .map(|&peak_idx| {
            let peak_value = request.signal[peak_idx];
            let mut left_min = peak_value;
            let mut right_min = peak_value;

            // Find minimum to the left
            for i in (0..peak_idx).rev() {
                left_min = left_min.min(request.signal[i]);
            }

            // Find minimum to the right
            for i in (peak_idx + 1)..request.signal.len() {
                right_min = right_min.min(request.signal[i]);
            }

            peak_value - left_min.max(right_min)
        })
        .collect();

    // Filter by prominence if specified
    if let Some(min_prominence) = request.prominence {
        let filtered: Vec<_> = peaks
            .iter()
            .zip(peak_heights.iter())
            .zip(prominences.iter())
            .filter(|((_idx, _height), prom)| **prom >= min_prominence)
            .map(|((idx, height), prom)| (*idx, *height, *prom))
            .collect();

        let (filtered_peaks, (filtered_heights, filtered_prominences)): (Vec<_>, (Vec<_>, Vec<_>)) =
            filtered
                .into_iter()
                .map(|(idx, height, prom)| (idx, (height, prom)))
                .unzip();

        Ok(PeakDetectionResult {
            total_peaks_found: filtered_peaks.len(),
            peak_indices: filtered_peaks,
            peak_heights: filtered_heights,
            peak_prominences: filtered_prominences,
        })
    } else {
        Ok(PeakDetectionResult {
            total_peaks_found: peaks.len(),
            peak_indices: peaks,
            peak_heights,
            peak_prominences: prominences,
        })
    }
}

#[derive(Deserialize, Serialize)]
pub struct FourierSeriesRequest {
    pub period: f64,
    pub coefficients: usize,
    pub samples: usize,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct FourierSeriesResult {
    pub time: Vec<f64>,
    pub signal: Vec<f64>,
    pub a0: f64,
    pub an: Vec<f64>,
    pub bn: Vec<f64>,
}

pub fn compute_fourier_series(
    request: FourierSeriesRequest,
) -> Result<FourierSeriesResult, String> {
    // Generate time points
    let time: Vec<f64> = (0..request.samples)
        .map(|i| i as f64 * request.period / request.samples as f64)
        .collect();

    // For demonstration, create a sample signal (square wave)
    // In practice, you'd compute coefficients from an actual signal
    let a0 = 0.0;
    let mut an = vec![0.0; request.coefficients];
    let mut bn = Vec::new();

    // Square wave Fourier series coefficients
    for n in 1..=request.coefficients {
        an[n - 1] = 0.0; // Square wave has no cosine terms
        let bn_val = if n % 2 == 1 {
            4.0 / (PI * n as f64)
        } else {
            0.0
        };
        bn.push(bn_val);
    }

    // Reconstruct signal from coefficients
    let signal: Vec<f64> = time
        .iter()
        .map(|&t| {
            let mut value = a0;
            for n in 1..=request.coefficients {
                let omega = 2.0 * PI * n as f64 / request.period;
                value += an[n - 1] * (omega * t).cos();
                value += bn[n - 1] * (omega * t).sin();
            }
            value
        })
        .collect();

    Ok(FourierSeriesResult {
        time,
        signal,
        a0,
        an,
        bn,
    })
}

#[derive(Deserialize, Serialize)]
pub struct WaveletTransformRequest {
    pub signal: Vec<f64>,
    pub wavelet_type: String,
    pub scales: Vec<f64>,
    pub sample_rate: f64,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct WaveletTransformResult {
    pub coefficients: Vec<Vec<f64>>,
    pub scales: Vec<f64>,
    pub time: Vec<f64>,
}

fn morlet_wavelet(t: f64, scale: f64) -> f64 {
    let sigma = 1.0;
    let scaled_t = t / scale;
    let envelope = (-scaled_t * scaled_t / (2.0 * sigma * sigma)).exp();
    let oscillation = (5.0 * scaled_t).cos();
    envelope * oscillation / scale.sqrt()
}

fn haar_wavelet(t: f64, scale: f64) -> f64 {
    let scaled_t = t / scale;
    if scaled_t >= 0.0 && scaled_t < 0.5 {
        1.0 / scale.sqrt()
    } else if scaled_t >= 0.5 && scaled_t < 1.0 {
        -1.0 / scale.sqrt()
    } else {
        0.0
    }
}

fn daubechies_wavelet(t: f64, scale: f64) -> f64 {
    // Daubechies 4 wavelet (db4) approximation
    let scaled_t = t / scale;
    let h0 = 0.6830127;
    let h1 = 1.1830127;
    let h2 = 0.3169873;
    let h3 = -0.1830127;

    if scaled_t >= 0.0 && scaled_t < 1.0 {
        h0 / scale.sqrt()
    } else if scaled_t >= 1.0 && scaled_t < 2.0 {
        h1 / scale.sqrt()
    } else if scaled_t >= 2.0 && scaled_t < 3.0 {
        h2 / scale.sqrt()
    } else if scaled_t >= 3.0 && scaled_t < 4.0 {
        h3 / scale.sqrt()
    } else {
        0.0
    }
}

fn mexican_hat_wavelet(t: f64, scale: f64) -> f64 {
    // Mexican Hat (Ricker) wavelet: second derivative of Gaussian
    let scaled_t = t / scale;
    let sigma2: f64 = 1.0;
    let norm = 2.0 / ((3.0_f64 * sigma2).sqrt() * PI.powf(0.25));
    let gaussian = (-scaled_t * scaled_t / (2.0 * sigma2)).exp();
    norm * (1.0 - scaled_t * scaled_t / sigma2) * gaussian / scale.sqrt()
}

pub fn wavelet_transform(
    request: WaveletTransformRequest,
) -> Result<WaveletTransformResult, String> {
    if request.signal.is_empty() {
        return Err("Signal cannot be empty".to_string());
    }

    let n = request.signal.len();
    let time: Vec<f64> = (0..n).map(|i| i as f64 / request.sample_rate).collect();

    let mut coefficients = Vec::new();

    for &scale in &request.scales {
        let mut scale_coeffs = Vec::new();

        for shift in 0..n {
            let mut coeff = 0.0;

            // Convolution with wavelet
            for i in 0..n {
                let t = (i as f64 - shift as f64) / request.sample_rate;
                let wavelet_val = match request.wavelet_type.as_str() {
                    "morlet" => morlet_wavelet(t, scale),
                    "haar" => haar_wavelet(t, scale),
                    "daubechies" => daubechies_wavelet(t, scale),
                    "mexican_hat" => mexican_hat_wavelet(t, scale),
                    _ => {
                        return Err(format!(
                            "Unsupported wavelet type: {}",
                            request.wavelet_type
                        ));
                    }
                };

                coeff += request.signal[i] * wavelet_val;
            }

            scale_coeffs.push(coeff);
        }

        coefficients.push(scale_coeffs);
    }

    Ok(WaveletTransformResult {
        coefficients,
        scales: request.scales.clone(),
        time,
    })
}

#[derive(Deserialize, Serialize)]
pub struct WindowFunctionRequest {
    pub size: usize,
    pub window_type: String,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct WindowFunctionResult {
    pub window: Vec<f64>,
}

// Modified Bessel function of the first kind of order zero (I0)
// Used for Kaiser window computation
fn bessel_i0(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 3.75 {
        let y = (x / 3.75) * (x / 3.75);
        1.0 + y
            * (3.5156229
                + y * (3.0899424
                    + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))))
    } else {
        let y = 3.75 / ax;
        (ax.exp() / ax.sqrt())
            * (0.39894228
                + y * (0.1328592e-1
                    + y * (0.225319e-2
                        + y * (-0.157565e-2
                            + y * (0.916281e-2
                                + y * (-0.2057706e-1
                                    + y * (0.2635537e-1
                                        + y * (-0.1647633e-1 + y * 0.392377e-2))))))))
    }
}

pub fn windowing_functions(request: WindowFunctionRequest) -> Result<WindowFunctionResult, String> {
    let n = request.size;

    let window: Vec<f64> = match request.window_type.as_str() {
        "hanning" => (0..n)
            .map(|i| 0.5 * (1.0 - (2.0 * PI * i as f64 / (n - 1) as f64).cos()))
            .collect(),
        "hamming" => (0..n)
            .map(|i| 0.54 - 0.46 * (2.0 * PI * i as f64 / (n - 1) as f64).cos())
            .collect(),
        "blackman" => (0..n)
            .map(|i| {
                0.42 - 0.5 * (2.0 * PI * i as f64 / (n - 1) as f64).cos()
                    + 0.08 * (4.0 * PI * i as f64 / (n - 1) as f64).cos()
            })
            .collect(),
        "kaiser" => {
            // Kaiser window with beta = 5.0 (adjustable parameter)
            let beta = 5.0;
            let alpha = (n - 1) as f64 / 2.0;
            (0..n)
                .map(|i| {
                    let x = (i as f64 - alpha) / alpha;
                    let arg = beta * (1.0 - x * x).max(0.0).sqrt();
                    bessel_i0(arg) / bessel_i0(beta)
                })
                .collect()
        }
        "bartlett" => (0..n)
            .map(|i| {
                let mid = (n - 1) as f64 / 2.0;
                1.0 - ((i as f64 - mid).abs() / mid)
            })
            .collect(),
        "rectangular" => vec![1.0; n],
        _ => return Err(format!("Unsupported window type: {}", request.window_type)),
    };

    Ok(WindowFunctionResult { window })
}
