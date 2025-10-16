//! Legacy API handlers organized by domain

pub mod tensor;
pub mod calculus;
pub mod fluid;
pub mod signal;
pub mod stochastic;
pub mod crypto;
pub mod regression;
pub mod dimensional;
pub mod validator;
pub mod geometry;
pub mod quantum;
pub mod linear_algebra;
pub mod statistics;
pub mod optimization;
pub mod graph_theory;
pub mod information_theory;

// Re-export all handler functions
pub use tensor::process_tensor_request;
pub use calculus::process_calculus_request;
pub use fluid::process_fluid_request;
pub use signal::process_signal_request;
pub use stochastic::process_stochastic_request;
pub use crypto::process_crypto_request;
pub use regression::process_regression_request;
pub use dimensional::process_dimensional_request;
pub use validator::process_validator_request;
pub use geometry::process_geometry_request;
pub use quantum::process_quantum_request;
pub use linear_algebra::process_linear_algebra_request;
pub use statistics::process_statistics_request;
pub use optimization::process_optimization_request;
pub use graph_theory::process_graph_theory_request;
pub use information_theory::process_information_theory_request;
