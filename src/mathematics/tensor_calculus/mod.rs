//! Tensor calculus module for symbolic tensor operations and Einstein field equations
//!
//! This module provides comprehensive tools for:
//! - Symbolic tensor calculus operations
//! - Einstein field equation solving
//! - Spacetime metric computations
//! - General relativity calculations
//! - Quantum tensor computations

pub mod einstein;
pub mod quantum_tensors;
pub mod symbolic;
pub mod tensor;

// Re-export commonly used types and functions
pub use einstein::*;
pub use quantum_tensors::*;
pub use symbolic::*;
pub use tensor::*;

use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum TensorError {
    #[error("JSON parsing error: {0}")]
    JsonError(#[from] serde_json::Error),
    #[error("Invalid metric tensor: {0}")]
    InvalidMetric(String),
    #[error("Computation error: {0}")]
    ComputationError(String),
}

#[derive(Serialize, Deserialize, Debug)]
pub struct TensorResult {
    pub result_type: String,
    pub data: serde_json::Value,
    pub coordinates: Vec<String>,
    pub success: bool,
    pub error: Option<String>,
}
