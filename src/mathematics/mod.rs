//! Mathematics Domain
//!
//! Core mathematical computation modules

pub mod calculus;
pub mod linear_algebra;
pub mod numerical;
pub mod special_functions;
pub mod symbolic_cas;
pub mod symbolic_regression;
pub mod tensor_calculus;

// Backwards compatibility
pub use calculus as advanced_calculus;
