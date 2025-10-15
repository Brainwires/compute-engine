//! Electromagnetism request handler

use crate::api::types::{ComputationRequest, ComputationResponse};
use crate::electromagnetism::*;
use serde_json::json;

pub fn handle(request: &ComputationRequest) -> ComputationResponse {
    let result = match request.operation.as_str() {
        "maxwell" | "maxwell_equations" => {
            let req: MaxwellRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match maxwell_equations(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "em_wave" | "wave" => {
            let req: WaveRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match em_wave(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "antenna" => {
            let req: AntennaRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match antenna_analysis(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "transmission_line" => {
            let req: TransmissionLineRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match transmission_line(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "waveguide" => {
            let req: WaveguideRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match waveguide(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "scattering" => {
            let req: ScatteringRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match scattering(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "poynting" | "poynting_vector" => {
            let req: PoyntingRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match poynting_vector(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        "skin_effect" => {
            let req: SkinEffectRequest = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(r) => r,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid parameters: {}", e),
                ),
            };
            match skin_effect(req) {
                Ok(result) => Ok(json!(result)),
                Err(e) => Err(e),
            }
        },
        _ => Err(format!("Unknown electromagnetism operation: {}", request.operation)),
    };

    match result {
        Ok(value) => ComputationResponse::success(
            request.module.clone(),
            request.operation.clone(),
            value,
        ),
        Err(e) => ComputationResponse::error(
            request.module.clone(),
            request.operation.clone(),
            e,
        ),
    }
}
