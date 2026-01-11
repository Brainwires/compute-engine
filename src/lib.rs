//! # Computational Engine
//!
//! A unified computational engine built around 8 tools with a clean, modular architecture.
//!
//! ## 8 Tools
//!
//! 1. **Solve** - Equations, systems, optimization, root finding
//! 2. **Compute** - Matrix ops, calculus, transforms, field theory, sampling
//! 3. **Analyze** - Series, limits, stability analysis, simplification
//! 4. **Simulate** - Time evolution, stochastic processes, fluid dynamics
//! 5. **ML** - Machine learning (clustering, neural nets, regression)
//! 6. **Chaos** - Chaos theory (fractals, attractors, Lyapunov exponents)
//! 7. **Units** - Dimensional analysis and unit conversion
//! 8. **Validate** - Equation and physics validation
//!
//! ## Quick Start
//!
//! ```rust,ignore
//! use computational_engine::{ToolRequest, SolveInput, create_default_dispatcher};
//!
//! // Create dispatcher
//! let dispatcher = create_default_dispatcher();
//!
//! // Solve a quadratic equation
//! let request = ToolRequest::Solve(SolveInput {
//!     equation_type: Default::default(),
//!     equations: vec!["x^2 + 2*x - 8 = 0".to_string()],
//!     ..Default::default()
//! });
//!
//! let response = dispatcher.dispatch(request).unwrap();
//! ```
//!
//! ## JSON API
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
//! ## Domain-Specific Modules
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

// Core engine
pub mod engine;
pub mod help;
pub mod help_auto;
pub mod implementations;

// MCP Server (Official MCP SDK integration)
#[cfg(feature = "rmcp")]
pub mod mcp_server;

// Scientific formula modules
pub mod biology;
pub mod chemistry;
pub mod datetime;
pub mod electrical;
pub mod engineering;
pub mod geophysics;
pub mod materials_science;
pub mod optics;
pub mod thermodynamics;

// WebAssembly bindings (when compiled with wasm feature)
#[cfg(feature = "wasm")]
pub mod wasm;

// Domain modules (implementation details)
pub mod mathematics;
pub mod physics;
pub mod specialized;
pub mod tools;

// Re-export modules at top level for convenience
pub use mathematics::advanced_calculus;
pub use mathematics::linear_algebra;
pub use mathematics::special_functions;
pub use mathematics::symbolic_regression;
pub use mathematics::tensor_calculus;
pub use physics::electromagnetism;
pub use physics::fluid_dynamics;
pub use physics::quantum_physics;
pub use specialized::cryptographic_mathematics;
pub use specialized::graph_theory;
pub use specialized::information_theory;
pub use specialized::optimization;
pub use specialized::statistics;
pub use specialized::stochastic_processes;
pub use tools::computational_geometry;
pub use tools::dimensional_analysis;
pub use tools::equation_validation;
pub use tools::numerical_methods;
pub use tools::signal_processing;

// Re-export engine types
pub use engine::{
    // 8 Tools - Traits
    Analyze, Chaos, Compute, MachineLearning, Simulate, Solve, Units, Validate,
    // 8 Tools - Input/Output types
    AnalyzeInput, AnalyzeOutput, ChaosInput, ChaosOutput, ComputeInput, ComputeOutput, MLInput,
    MLOutput, SimulateInput, SimulateOutput, SolveInput, SolveOutput, UnitsInput, UnitsOutput,
    ValidateInput, ValidateOutput, ValidationResult,
    // Operation enums
    AnalysisOp, ChaosOp, ComputeOp, MLOp, SimulationModel, UnitsOp, ValidateOp,
    // ML sub-enums
    ClassificationMethod, ClusteringMethod, DimReductionMethod, NeuralNetworkOp, RegressionMethod,
    // Chaos sub-enums
    AttractorType, BifurcationType, DimensionMethod, FractalType, LyapunovMethod,
    // Equation types
    ChemicalEquation, DifferentialEquation, EinsteinEquation, EMEquation, EquationType,
    FluidEquation,
    // Dispatcher
    ToolDispatcher, ToolRequest, ToolResponse, ToolResult,
    // Common types
    Domain, Method,
};

// Re-export implementation helper
pub use implementations::create_default_dispatcher;

// Re-export operation listing from engine
pub use engine::{list_all_operations, list_tools};

pub use tensor_calculus::TensorError;

/// Common result type for computational operations
pub type ComputationalResult<T> = Result<T, Box<dyn std::error::Error>>;

/// Version information
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub const NAME: &str = env!("CARGO_PKG_NAME");

