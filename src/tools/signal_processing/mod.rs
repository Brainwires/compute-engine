//! Signal processing module
//!
//! Provides digital signal processing capabilities:
//! - Fast Fourier Transform (FFT) and inverse FFT
//! - Digital filters (low-pass, high-pass, band-pass)
//! - Spectrogram generation
//! - Power spectral density estimation
//! - Window functions
//! - Peak detection
//! - Fourier series
//! - Wavelet transform

mod additional;
mod lib;

pub use lib::{
    FFTRequest, FFTResult, FilterRequest, FilterResult, PSDRequest, PSDResult, SpectrogramRequest,
    SpectrogramResult, apply_filter, compute_fft, compute_psd, compute_spectrogram,
};

pub use additional::{
    FourierSeriesRequest, FourierSeriesResult, PeakDetectionRequest, PeakDetectionResult,
    WaveletTransformRequest, WaveletTransformResult, WindowFunctionRequest, WindowFunctionResult,
    compute_fourier_series, detect_peaks, wavelet_transform, windowing_functions,
};
