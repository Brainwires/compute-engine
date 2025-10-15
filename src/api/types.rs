//! Core API types shared across all handlers

use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;

#[derive(Debug, Serialize, Deserialize)]
pub struct ComputationRequest {
    pub module: String,
    pub operation: String,
    pub parameters: HashMap<String, Value>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ComputationResponse {
    pub success: bool,
    pub module: String,
    pub operation: String,
    pub result: Option<Value>,
    pub error: Option<String>,
}

impl ComputationResponse {
    pub fn success(module: String, operation: String, result: Value) -> Self {
        Self {
            success: true,
            module,
            operation,
            result: Some(result),
            error: None,
        }
    }

    pub fn error(module: String, operation: String, error: String) -> Self {
        Self {
            success: false,
            module,
            operation,
            result: None,
            error: Some(error),
        }
    }
}
