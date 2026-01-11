//! # Computational Engine
//!
//! A unified computational engine built around 8 primary tools with a clean, modular architecture.
//!
//! ## Architecture
//!
//! ### 8 Primary Tools (NEW - Recommended)
//! The engine is organized around 8 powerful, flexible tools:
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
//! ### Legacy Tools (Backward Compatible)
//! - **Differentiate** - Routes to Compute with differentiate operation
//! - **Integrate** - Routes to Compute with integrate operation
//! - **Transform** - Routes to Compute with transform operation
//! - **FieldTheory** - Routes to Compute with field operation
//! - **Sample** - Routes to Compute with sample operation
//! - **Optimize** - Routes to Solve with optimize equation type
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
pub mod help;
pub mod help_auto;

// MCP Server (Official MCP SDK integration - primary interface)
#[cfg(feature = "rmcp")]
pub mod mcp_server;

// Scientific formula modules (NEW - 2025 Expansion)
pub mod biology;
pub mod chemistry;
pub mod datetime;
pub mod electrical; // Circuit analysis + NEC calculations
pub mod engineering;
pub mod geophysics;
pub mod materials_science; // Crystal structures, band theory, mechanical properties
pub mod optics;
pub mod thermodynamics;

// WebAssembly bindings (when compiled with wasm feature)
#[cfg(feature = "wasm")]
pub mod wasm;

// Domain modules (Legacy - for implementation)
pub mod mathematics;
pub mod physics;
pub mod specialized;
pub mod tools;

// Backwards compatibility - re-export modules at top level
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
// Note: specialized::chemistry has been replaced by the new chemistry module

// Re-export new engine types (Recommended)
pub use engine::{
    // Primary 8 Tools - Traits
    Analyze,
    Chaos,
    Compute,
    MachineLearning,
    Simulate,
    Solve,
    Units,
    Validate,
    // Primary 8 Tools - Input/Output types
    AnalyzeInput,
    AnalyzeOutput,
    ChaosInput,
    ChaosOutput,
    ComputeInput,
    ComputeOutput,
    MLInput,
    MLOutput,
    SimulateInput,
    SimulateOutput,
    SolveInput,
    SolveOutput,
    UnitsInput,
    UnitsOutput,
    ValidateInput,
    ValidateOutput,
    // Primary 8 Tools - Operation enums
    AnalysisOp,
    ChaosOp,
    ComputeOp,
    MLOp,
    SimulationModel,
    UnitsOp,
    ValidateOp,
    // ML sub-enums
    ClusteringMethod,
    NeuralNetworkOp,
    RegressionMethod,
    DimReductionMethod,
    ClassificationMethod,
    // Chaos sub-enums
    FractalType,
    AttractorType,
    LyapunovMethod,
    BifurcationType,
    DimensionMethod,
    // Legacy Tools - Traits
    Differentiate,
    FieldTheory,
    Integrate,
    Optimize,
    Sample,
    Transform,
    // Legacy Tools - Input/Output types
    DifferentiateInput,
    DifferentiateOutput,
    FieldTheoryInput,
    FieldTheoryOutput,
    IntegrateInput,
    IntegrateOutput,
    OptimizeInput,
    OptimizeOutput,
    SampleInput,
    SampleOutput,
    TransformInput,
    TransformOutput,
    // Legacy Tools - Operation enums
    DifferentiationOp,
    FieldType,
    IntegrationType,
    OptimizationMethod,
    SamplingMethod,
    TransformType,
    // Equation types
    ChemicalEquation,
    DifferentialEquation,
    EinsteinEquation,
    EMEquation,
    EquationType,
    FluidEquation,
    // Dispatcher
    ToolDispatcher,
    ToolRequest,
    ToolResponse,
    ToolResult,
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

