//! # Computational Engine
//!
//! A unified computational engine built around 5 core tools with a clean, modular architecture.
//!
//! ## Architecture
//!
//! ### 5 Core Tools (NEW - Recommended)
//! The engine is organized around 5 powerful, flexible tools:
//!
//! 1. **Solve** - Equations, systems, optimization, root finding
//! 2. **Differentiate** - Derivatives, gradients, Jacobians, Hessians
//! 3. **Integrate** - Definite/indefinite, single/multivariable integrals
//! 4. **Analyze** - Simplify, expand, factor, transform, evaluate
//! 5. **Simulate** - ODEs, PDEs, physics, mechanics, optimization
//!
//! ### Quick Start (New API)
//!
//! ```rust
//! use computational_engine::core::{ToolRequest, SolveInput};
//! use computational_engine::implementations::create_default_dispatcher;
//!
//! // Create dispatcher
//! let dispatcher = create_default_dispatcher();
//!
//! // Solve a quadratic equation
//! let request = ToolRequest::Solve(SolveInput {
//!     equations: vec!["x^2 + 2*x - 8 = 0".to_string()],
//!     variables: None,
//!     initial_guess: None,
//!     domain: None,
//!     method: None,
//! });
//!
//! let response = dispatcher.dispatch(request).unwrap();
//! ```
//!
//! ### JSON API
//!
//! ```json
//! {
//!   "tool": "solve",
//!   "input": {
//!     "equations": ["x^2 - 4 = 0"]
//!   }
//! }
//! ```
//!
//! ## Legacy Domain-Specific Modules (Still Available)
//!
//! ### Mathematics Domain
//! - **Tensor Calculus**: Symbolic tensor operations and Einstein field equations
//! - **Advanced Calculus**: Fractional, variational, and stochastic calculus
//! - **Linear Algebra**: Matrix operations, decompositions, PCA
//! - **Symbolic Regression**: Automated function discovery
//!
//! ### Physics Domain
//! - **Fluid Dynamics**: Navier-Stokes solvers and flow analysis
//! - **Quantum Physics**: Quantum simulations and particle physics
//!
//! ### Engineering Tools
//! - **Signal Processing**: FFT, filters, and signal analysis
//! - **Dimensional Analysis**: Unit checking and dimensional validation
//! - **Equation Validation**: Mathematical equation verification
//! - **Computational Geometry**: Geometric algorithms and spatial computations
//!
//! ### Specialized Domains
//! - **Stochastic Processes**: Random walks, Monte Carlo, and stochastic simulations
//! - **Cryptographic Mathematics**: Number theory and cryptographic primitives
//! - **Statistics**: Statistical analysis and probability

// Core 10-tool engine (NEW - Recommended)
pub mod engine;
pub mod implementations;

// MCP Server (Official MCP SDK integration - primary interface)
pub mod mcp_server;

// Scientific formula modules (NEW - 2025 Expansion)
pub mod chemistry;
pub mod biology;
pub mod thermodynamics;
pub mod optics;
pub mod geophysics;
pub mod engineering;
pub mod datetime;

// WebAssembly bindings (when compiled with wasm feature)
#[cfg(feature = "wasm")]
pub mod wasm;

// Domain modules (Legacy - for implementation)
pub mod mathematics;
pub mod physics;
pub mod tools;
pub mod specialized;

// Backwards compatibility - re-export modules at top level
pub use mathematics::tensor_calculus;
pub use mathematics::advanced_calculus;
pub use mathematics::linear_algebra;
pub use mathematics::symbolic_regression;
pub use mathematics::special_functions;
pub use physics::fluid_dynamics;
pub use physics::quantum_physics;
pub use physics::electromagnetism;
pub use tools::signal_processing;
pub use tools::dimensional_analysis;
pub use tools::equation_validation;
pub use tools::computational_geometry;
pub use tools::numerical_methods;
pub use specialized::stochastic_processes;
pub use specialized::cryptographic_mathematics;
pub use specialized::statistics;
pub use specialized::optimization;
pub use specialized::graph_theory;
pub use specialized::information_theory;
// Note: specialized::chemistry has been replaced by the new chemistry module

// Unified API for all operations
pub mod api;

// Re-export new engine types (Recommended)
pub use engine::{
    // Traits
    Solve, Differentiate, Integrate, Analyze, Simulate,
    Compute, Transform, FieldTheory, Sample, Optimize,
    // Types
    SolveInput, SolveOutput,
    DifferentiateInput, DifferentiateOutput,
    IntegrateInput, IntegrateOutput,
    AnalyzeInput, AnalyzeOutput,
    SimulateInput, SimulateOutput,
    ComputeInput, ComputeOutput,
    TransformInput, TransformOutput,
    FieldTheoryInput, FieldTheoryOutput,
    SampleInput, SampleOutput,
    OptimizeInput, OptimizeOutput,
    // Dispatcher
    ToolDispatcher, ToolRequest, ToolResponse, ToolResult,
    // Equation types
    EquationType, EinsteinEquation, FluidEquation, DifferentialEquation,
    EMEquation, ChemicalEquation, DifferentiationOp, IntegrationType,
    AnalysisOp, SimulationModel, ComputeOp, TransformType, FieldType,
    SamplingMethod, OptimizationMethod,
};

// Re-export implementation helper
pub use implementations::create_default_dispatcher;

// Re-export legacy API (Still supported for backwards compatibility)
pub use tensor_calculus::TensorError;
pub use api::{ComputationRequest, ComputationResponse, process_request, process_json_request, list_all_operations};

/// Common result type for computational operations
pub type ComputationalResult<T> = Result<T, Box<dyn std::error::Error>>;

/// Version information
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub const NAME: &str = env!("CARGO_PKG_NAME");

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_version() {
        assert!(!VERSION.is_empty());
        assert_eq!(NAME, "computational-engine");
    }
}
