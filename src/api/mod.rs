//! Unified API for all computational engine operations
//!
//! This module provides a clean, modular JSON-based interface to all computational functions.
//! Handlers are organized by domain for better maintainability.

pub mod handlers;
pub mod operations;
pub mod types;

pub use operations::list_all_operations;
pub use types::{ComputationRequest, ComputationResponse};

/// Process a computation request and route it to the appropriate module handler
pub fn process_request(request: &ComputationRequest) -> ComputationResponse {
    handlers::route_request(request)
}

/// Process JSON string requests
pub fn process_json_request(json_str: &str) -> String {
    let request: ComputationRequest = match serde_json::from_str(json_str) {
        Ok(req) => req,
        Err(e) => {
            return serde_json::to_string(&ComputationResponse {
                success: false,
                module: "unknown".to_string(),
                operation: "parse".to_string(),
                result: None,
                error: Some(format!("Invalid JSON: {}", e)),
            })
            .unwrap();
        }
    };

    let response = process_request(&request);
    serde_json::to_string(&response).unwrap()
}

