//! Mathematics Domain
//!
//! Core mathematical computation modules

pub mod tensor_calculus;
pub mod calculus;
pub mod linear_algebra;
pub mod symbolic_regression;
pub mod special_functions;
pub mod symbolic_cas;
pub mod numerical;

// Backwards compatibility
pub use calculus as advanced_calculus;
