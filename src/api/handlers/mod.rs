//! Request handlers for all computational modules
//!
//! Each module has its own handler function that processes requests
//! and routes them to the appropriate operations.

use super::types::{ComputationRequest, ComputationResponse};

// Mathematics domain handlers
pub mod tensor_calculus;
pub mod advanced_calculus;
pub mod linear_algebra;
pub mod special_functions;

// Physics domain handlers
pub mod fluid_dynamics;
pub mod quantum_physics;
pub mod electromagnetism;

// Tools domain handlers
pub mod signal_processing;
pub mod dimensional_analysis;
pub mod equation_validation;
pub mod computational_geometry;
pub mod numerical_methods;
pub mod function_approximator;

// Specialized domain handlers
pub mod stochastic_processes;
pub mod cryptographic_mathematics;
pub mod symbolic_regression;
pub mod statistics;
pub mod optimization;
pub mod graph_theory;
pub mod information_theory;
pub mod chemistry;

// Scientific formula modules (2025 expansion)
pub mod biology;
pub mod thermodynamics;
pub mod optics;
pub mod geophysics;
pub mod engineering;
pub mod datetime;

/// Route a request to the appropriate module handler
pub fn route_request(request: &ComputationRequest) -> ComputationResponse {
    match request.module.as_str() {
        // Mathematics
        "tensor_calculus" | "tensor" => tensor_calculus::handle(request),
        "advanced_calculus" | "calculus" => advanced_calculus::handle(request),
        "linear_algebra" | "linalg" => linear_algebra::handle(request),
        "special_functions" | "special" => special_functions::handle(request),

        // Physics
        "fluid_dynamics" | "fluid" => fluid_dynamics::handle(request),
        "quantum_physics" | "quantum" => quantum_physics::handle(request),
        "electromagnetism" | "em" => electromagnetism::handle(request),

        // Tools
        "signal_processing" | "signal" => signal_processing::handle(request),
        "dimensional_analysis" | "dimensional" => dimensional_analysis::handle(request),
        "equation_validation" | "validator" => equation_validation::handle(request),
        "computational_geometry" | "geometry" => computational_geometry::handle(request),
        "numerical_methods" | "numerical" => numerical_methods::handle(request),
        "function_approximator" | "approximate" => function_approximator::handle(request),

        // Specialized
        "stochastic_processes" | "stochastic" => stochastic_processes::handle(request),
        "cryptographic_mathematics" | "crypto" => cryptographic_mathematics::handle(request),
        "symbolic_regression" | "regression" => symbolic_regression::handle(request),
        "statistics" | "stats" => statistics::handle(request),
        "optimization" | "optimize" => optimization::handle(request),
        "graph_theory" | "graph" => graph_theory::handle(request),
        "information_theory" | "info" => information_theory::handle(request),
        "chemistry" | "chem" => chemistry::handle(request),

        // Scientific formulas (2025 expansion)
        "biology" | "bio" => biology::handle(request),
        "thermodynamics" | "thermo" | "heat" => thermodynamics::handle(request),
        "optics" | "optical" => optics::handle(request),
        "geophysics" | "geo" => geophysics::handle(request),
        "engineering" | "eng" => engineering::handle(request),
        "datetime" | "date" | "time" => datetime::handle(request),

        _ => ComputationResponse::error(
            request.module.clone(),
            request.operation.clone(),
            format!("Unknown module: {}", request.module),
        ),
    }
}
