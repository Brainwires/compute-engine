use crate::api::types::ComputationResponse;
use crate::specialized::chemistry as chem;
use serde_json::{Value, json};

pub fn handle(request: &crate::api::types::ComputationRequest) -> ComputationResponse {
    let params_map: serde_json::Map<String, Value> =
        request.parameters.clone().into_iter().collect();

    match request.operation.as_str() {
        "molar_mass" => {
            let input: chem::MolarMassRequest = match serde_json::from_value(Value::Object(params_map)) {
                Ok(i) => i,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid molar_mass request: {}", e),
                ),
            };
            match chem::molar_mass(input) {
                Ok(result) => ComputationResponse::success(
                    request.module.clone(),
                    request.operation.clone(),
                    json!(result),
                ),
                Err(e) => ComputationResponse::error(request.module.clone(), request.operation.clone(), e),
            }
        }
        "gas_law" => {
            let input: chem::GasLawRequest = match serde_json::from_value(Value::Object(params_map)) {
                Ok(i) => i,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid gas_law request: {}", e),
                ),
            };
            match chem::gas_law(input) {
                Ok(result) => ComputationResponse::success(
                    request.module.clone(),
                    request.operation.clone(),
                    json!(result),
                ),
                Err(e) => ComputationResponse::error(request.module.clone(), request.operation.clone(), e),
            }
        }
        "thermodynamics" => {
            let input: chem::ThermodynamicsRequest = match serde_json::from_value(Value::Object(params_map)) {
                Ok(i) => i,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid thermodynamics request: {}", e),
                ),
            };
            match chem::thermodynamics(input) {
                Ok(result) => ComputationResponse::success(
                    request.module.clone(),
                    request.operation.clone(),
                    json!(result),
                ),
                Err(e) => ComputationResponse::error(request.module.clone(), request.operation.clone(), e),
            }
        }
        "acid_base" => {
            let input: chem::AcidBaseRequest = match serde_json::from_value(Value::Object(params_map)) {
                Ok(i) => i,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid acid_base request: {}", e),
                ),
            };
            match chem::acid_base(input) {
                Ok(result) => ComputationResponse::success(
                    request.module.clone(),
                    request.operation.clone(),
                    json!(result),
                ),
                Err(e) => ComputationResponse::error(request.module.clone(), request.operation.clone(), e),
            }
        }
        "kinetics" => {
            let input: chem::KineticsRequest = match serde_json::from_value(Value::Object(params_map)) {
                Ok(i) => i,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid kinetics request: {}", e),
                ),
            };
            match chem::kinetics(input) {
                Ok(result) => ComputationResponse::success(
                    request.module.clone(),
                    request.operation.clone(),
                    json!(result),
                ),
                Err(e) => ComputationResponse::error(request.module.clone(), request.operation.clone(), e),
            }
        }
        "electrochemistry" => {
            let input: chem::ElectrochemistryRequest = match serde_json::from_value(Value::Object(params_map)) {
                Ok(i) => i,
                Err(e) => return ComputationResponse::error(
                    request.module.clone(),
                    request.operation.clone(),
                    format!("Invalid electrochemistry request: {}", e),
                ),
            };
            match chem::electrochemistry(input) {
                Ok(result) => ComputationResponse::success(
                    request.module.clone(),
                    request.operation.clone(),
                    json!(result),
                ),
                Err(e) => ComputationResponse::error(request.module.clone(), request.operation.clone(), e),
            }
        }
        _ => ComputationResponse::error(
            request.module.clone(),
            request.operation.clone(),
            format!("Unknown chemistry operation: {}", request.operation),
        ),
    }
}

