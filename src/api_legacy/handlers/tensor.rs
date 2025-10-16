use super::super::{ComputationRequest, ComputationResponse};
use serde_json::json;

pub fn process_tensor_request(request: &ComputationRequest) -> ComputationResponse {
    use crate::tensor_calculus::*;

    let result = match request.operation.as_str() {
        "solve_vacuum_einstein" => {
            let coords = request.parameters.get("coordinates")
                .and_then(|v| serde_json::from_value::<Vec<String>>(v.clone()).ok())
                .unwrap_or_else(|| vec!["t".to_string(), "r".to_string(), "theta".to_string(), "phi".to_string()]);
            let system = request.parameters.get("system")
                .and_then(|v| v.as_str())
                .unwrap_or("spherical");

            match solve_vacuum_einstein_equations(&coords, system, &[]) {
                Ok(solutions) => Ok(json!({
                    "solutions": solutions,
                    "count": solutions.len()
                })),
                Err(e) => Err(format!("Einstein equation error: {}", e)),
            }
        },
        // TODO: Add remaining 9 tensor operations with proper type conversion
        _ => Err(format!("Unknown tensor operation: {}", request.operation)),
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
