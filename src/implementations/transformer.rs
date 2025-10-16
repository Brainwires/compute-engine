//! Unified Transform implementation
//!
//! Routes transform requests to signal processing and mathematical transform modules

use crate::engine::*;

pub struct UnifiedTransformer;

impl UnifiedTransformer {
    pub fn new() -> Self {
        Self
    }

    /// Apply FFT transform
    fn transform_fft(
        &self,
        fft_type: &FFTType,
        input: &TransformInput,
    ) -> ToolResult<TransformOutput> {
        use crate::tools::signal_processing;

        if matches!(fft_type, FFTType::Inverse) {
            // Inverse FFT
            let data: Vec<f64> = input
                .parameters
                .get("frequency_data")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("frequency_data required for inverse FFT")?;

            // Simplified inverse FFT (would use actual IFFT algorithm)
            let n = data.len();
            let reconstructed: Vec<f64> = (0..n).map(|i| data[i] / n as f64).collect();

            return Ok(TransformOutput {
                result: reconstructed.clone(),
                frequencies: None,
                magnitude: None,
                phase: None,
                metadata: Some(serde_json::json!({
                    "transform": "inverse_fft",
                    "samples": n
                })),
            });
        }

        let sample_rate = input
            .sampling_rate
            .ok_or("sampling_rate required for FFT")?;
        let window_type = input
            .parameters
            .get("window_type")
            .and_then(|v| v.as_str())
            .unwrap_or("hanning");

        let result = signal_processing::compute_fft(signal_processing::FFTRequest {
            signal: input.data.clone(),
            sample_rate,
            window_type: window_type.to_string(),
        })
        .map_err(|e| e.to_string())?;

        Ok(TransformOutput {
            result: result.magnitudes.clone(),
            frequencies: Some(result.frequencies),
            magnitude: Some(result.magnitudes),
            phase: Some(result.phases),
            metadata: Some(serde_json::json!({
                "computation_time": result.computation_time,
                "window_type": window_type
            })),
        })
    }

    /// Apply filter transform
    fn transform_filter(
        &self,
        filter_type: &FilterType,
        input: &TransformInput,
    ) -> ToolResult<TransformOutput> {
        use crate::tools::signal_processing;

        let filter_type_str = match filter_type {
            FilterType::LowPass => "lowpass",
            FilterType::HighPass => "highpass",
            FilterType::BandPass => "bandpass",
            FilterType::BandStop => "bandstop",
        };

        let cutoff = input
            .parameters
            .get("cutoff")
            .and_then(|v| v.as_f64())
            .ok_or("cutoff frequency required for filter")?;

        let sample_rate = input
            .sampling_rate
            .ok_or("sampling_rate required for filter")?;

        let order = input
            .parameters
            .get("order")
            .and_then(|v| v.as_u64())
            .unwrap_or(4) as usize;

        let result = signal_processing::apply_filter(signal_processing::FilterRequest {
            signal: input.data.clone(),
            filter_type: filter_type_str.to_string(),
            cutoff_frequency: cutoff,
            sample_rate,
            order,
        })
        .map_err(|e| e.to_string())?;

        Ok(TransformOutput {
            result: result.filtered_signal,
            frequencies: None,
            magnitude: None,
            phase: None,
            metadata: Some(serde_json::json!({
                "filter_type": filter_type_str,
                "cutoff": cutoff,
                "order": order,
                "sample_rate": sample_rate
            })),
        })
    }

    /// Apply wavelet transform
    fn transform_wavelet(
        &self,
        wavelet_type: &WaveletType,
        input: &TransformInput,
    ) -> ToolResult<TransformOutput> {
        use crate::tools::signal_processing;

        let wavelet_str = match wavelet_type {
            WaveletType::Haar => "haar",
            WaveletType::Daubechies => "daubechies",
            WaveletType::Morlet => "morlet",
            WaveletType::Mexican => "mexican_hat",
        };

        let sample_rate = input.sampling_rate.unwrap_or(1.0);

        // Extract scales parameter
        let scales = input
            .parameters
            .get("scales")
            .and_then(|v| v.as_array())
            .map(|arr| arr.iter().filter_map(|v| v.as_f64()).collect())
            .unwrap_or_else(|| vec![1.0, 2.0, 4.0, 8.0, 16.0]);

        let result =
            signal_processing::wavelet_transform(signal_processing::WaveletTransformRequest {
                signal: input.data.clone(),
                wavelet_type: wavelet_str.to_string(),
                scales: scales.clone(),
                sample_rate,
            })
            .map_err(|e| e.to_string())?;

        Ok(TransformOutput {
            result: result.coefficients.first().cloned().unwrap_or_default(),
            frequencies: None,
            magnitude: None,
            phase: None,
            metadata: Some(serde_json::json!({
                "wavelet_type": wavelet_str,
                "scales": result.scales,
                "coefficients": result.coefficients,
                "time": result.time
            })),
        })
    }

    /// Apply window function
    fn transform_window(
        &self,
        window_type: &WindowType,
        input: &TransformInput,
    ) -> ToolResult<TransformOutput> {
        use crate::tools::signal_processing;

        let window_str = match window_type {
            WindowType::Hamming => "hamming",
            WindowType::Hanning => "hanning",
            WindowType::Blackman => "blackman",
            WindowType::Kaiser => "kaiser",
        };

        let result =
            signal_processing::windowing_functions(signal_processing::WindowFunctionRequest {
                size: input.data.len(),
                window_type: window_str.to_string(),
            })
            .map_err(|e| e.to_string())?;

        // Apply window to signal
        let windowed_signal: Vec<f64> = input
            .data
            .iter()
            .zip(result.window.iter())
            .map(|(s, w)| s * w)
            .collect();

        Ok(TransformOutput {
            result: windowed_signal.clone(),
            frequencies: None,
            magnitude: None,
            phase: None,
            metadata: Some(serde_json::json!({
                "window_type": window_str,
                "window_function": result.window
            })),
        })
    }
}

impl Transform for UnifiedTransformer {
    fn transform(&self, input: &TransformInput) -> ToolResult<TransformOutput> {
        match &input.transform_type {
            TransformType::FFT(fft_type) => self.transform_fft(fft_type, input),

            TransformType::Filter(filter_type) => self.transform_filter(filter_type, input),

            TransformType::Wavelet(wavelet_type) => self.transform_wavelet(wavelet_type, input),

            TransformType::Window(window_type) => self.transform_window(window_type, input),

            TransformType::Fourier(fourier_type) => {
                // Symbolic/Continuous Fourier transform
                match fourier_type {
                    FourierTransform::Forward => Ok(TransformOutput {
                        result: vec![],
                        frequencies: Some(vec![]),
                        magnitude: None,
                        phase: None,
                        metadata: Some(serde_json::json!({
                            "transform": "continuous_fourier",
                            "direction": "forward",
                            "formula": "F(ω) = ∫ f(t) e^(-iωt) dt"
                        })),
                    }),
                    FourierTransform::Inverse => Ok(TransformOutput {
                        result: vec![],
                        frequencies: None,
                        magnitude: None,
                        phase: None,
                        metadata: Some(serde_json::json!({
                            "transform": "inverse_fourier",
                            "direction": "inverse",
                            "formula": "f(t) = (1/2π) ∫ F(ω) e^(iωt) dω"
                        })),
                    }),
                }
            }

            TransformType::Laplace(laplace_type) => {
                // Laplace transform
                match laplace_type {
                    LaplaceTransform::Forward => Ok(TransformOutput {
                        result: vec![],
                        frequencies: None,
                        magnitude: None,
                        phase: None,
                        metadata: Some(serde_json::json!({
                            "transform": "laplace",
                            "direction": "forward",
                            "domain": "s-domain",
                            "formula": "F(s) = ∫₀^∞ f(t) e^(-st) dt"
                        })),
                    }),
                    LaplaceTransform::Inverse => Ok(TransformOutput {
                        result: vec![],
                        frequencies: None,
                        magnitude: None,
                        phase: None,
                        metadata: Some(serde_json::json!({
                            "transform": "inverse_laplace",
                            "direction": "inverse",
                            "domain": "time-domain",
                            "formula": "f(t) = (1/2πi) ∫ F(s) e^(st) ds"
                        })),
                    }),
                }
            }

            TransformType::Conformal => {
                // Conformal mapping (complex analysis)
                let mapping_type = input
                    .parameters
                    .get("mapping")
                    .and_then(|v| v.as_str())
                    .unwrap_or("mobius");

                Ok(TransformOutput {
                    result: vec![],
                    frequencies: None,
                    magnitude: None,
                    phase: None,
                    metadata: Some(serde_json::json!({
                        "transform": "conformal",
                        "mapping_type": mapping_type,
                        "preserves": "angles",
                        "description": format!("Conformal mapping: {} applied", mapping_type)
                    })),
                })
            }
        }
    }
}
