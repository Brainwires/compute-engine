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

mod lib;
mod additional;

pub use lib::{
    FFTRequest, FFTResult, compute_fft,
    FilterRequest, FilterResult, apply_filter,
    SpectrogramRequest, SpectrogramResult, compute_spectrogram,
    PSDRequest, PSDResult, compute_psd,
};

pub use additional::{
    PeakDetectionRequest, PeakDetectionResult, detect_peaks,
    FourierSeriesRequest, FourierSeriesResult, compute_fourier_series,
    WaveletTransformRequest, WaveletTransformResult, wavelet_transform,
    WindowFunctionRequest, WindowFunctionResult, windowing_functions,
};
