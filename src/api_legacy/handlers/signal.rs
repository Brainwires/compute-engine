use super::super::{ComputationRequest, ComputationResponse};
use serde_json::json;

pub fn process_signal_request(request: &ComputationRequest) -> ComputationResponse {
    use crate::signal_processing::*;

    let result = match request.operation.as_str() {
        "fft" => {
            let fft_request: FFTRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse {
                    success: false,
                    module: request.module.clone(),
                    operation: request.operation.clone(),
                    result: None,
                    error: Some(format!("Invalid parameters: {}", e)),
                },
            };

            match compute_fft(fft_request) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "filter" => {
            let filter_request: FilterRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse {
                    success: false,
                    module: request.module.clone(),
                    operation: request.operation.clone(),
                    result: None,
                    error: Some(format!("Invalid parameters: {}", e)),
                },
            };

            match apply_filter(filter_request) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "spectrogram" => {
            let spec_request: SpectrogramRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse {
                    success: false,
                    module: request.module.clone(),
                    operation: request.operation.clone(),
                    result: None,
                    error: Some(format!("Invalid parameters: {}", e)),
                },
            };

            match compute_spectrogram(spec_request) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "psd" | "analyze_power_spectrum" => {
            let psd_request: PSDRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse {
                    success: false,
                    module: request.module.clone(),
                    operation: request.operation.clone(),
                    result: None,
                    error: Some(format!("Invalid parameters: {}", e)),
                },
            };

            match compute_psd(psd_request) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "detect_peaks" => {
            let req: PeakDetectionRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse {
                    success: false,
                    module: request.module.clone(),
                    operation: request.operation.clone(),
                    result: None,
                    error: Some(format!("Invalid parameters: {}", e)),
                },
            };

            match detect_peaks(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "compute_fourier_series" => {
            let req: FourierSeriesRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse {
                    success: false,
                    module: request.module.clone(),
                    operation: request.operation.clone(),
                    result: None,
                    error: Some(format!("Invalid parameters: {}", e)),
                },
            };

            match compute_fourier_series(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "wavelet_transform" => {
            let req: WaveletTransformRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse {
                    success: false,
                    module: request.module.clone(),
                    operation: request.operation.clone(),
                    result: None,
                    error: Some(format!("Invalid parameters: {}", e)),
                },
            };

            match wavelet_transform(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "windowing_functions" => {
            let req: WindowFunctionRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse {
                    success: false,
                    module: request.module.clone(),
                    operation: request.operation.clone(),
                    result: None,
                    error: Some(format!("Invalid parameters: {}", e)),
                },
            };

            match windowing_functions(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        _ => Err(format!("Unknown signal operation: {}", request.operation)),
    };

    match result {
        Ok(value) => ComputationResponse {
            success: true,
            module: request.module.clone(),
            operation: request.operation.clone(),
            result: Some(value),
            error: None,
        },
        Err(e) => ComputationResponse {
            success: false,
            module: request.module.clone(),
            operation: request.operation.clone(),
            result: None,
            error: Some(e),
        },
    }
}
