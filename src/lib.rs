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
//! All domain functionality is organized under the 8 tool modules:
//!
//! ### Solve Tool (`solve/`)
//! - Equation solving, optimization, root finding
//! - Specialized: Game theory, linear programming
//!
//! ### Compute Tool (`compute/`)
//! - Matrix operations, tensor calculus, linear algebra
//! - Calculus: Differentiation, integration
//! - Transforms: FFT, Laplace, wavelets
//! - Physics: EM, quantum, relativity, nuclear
//! - Scientific: Chemistry, biology, thermodynamics, optics
//! - Geometry, graph algorithms, information theory
//!
//! ### Analyze Tool (`analyze/`)
//! - Symbolic simplification, expansion, factorization
//! - Series: Taylor, Laurent
//! - Stability analysis
//!
//! ### Simulate Tool (`simulate/`)
//! - ODE solvers: Euler, Runge-Kutta
//! - Stochastic: Brownian motion, Monte Carlo
//! - Fluid dynamics: Navier-Stokes
//! - Finance: Black-Scholes, Heston
//!
//! ### ML Tool (`ml/`)
//! - Clustering, regression, classification
//! - Neural networks, dimensionality reduction
//!
//! ### Chaos Tool (`chaos/`)
//! - Fractals: Mandelbrot, Julia
//! - Attractors: Lorenz, Rossler
//! - Lyapunov exponents, bifurcation analysis
//!
//! ### Units Tool (`units/`)
//! - Unit conversion, dimensional analysis
//!
//! ### Validate Tool (`validate/`)
//! - Equation validation, physics consistency checks

// Core engine
pub mod engine;
pub mod help;
pub mod help_auto;

// MCP Server (Official MCP SDK integration)
#[cfg(feature = "rmcp")]
pub mod mcp_server;

// ============================================================================
// 8 TOOL MODULES - All computational functionality organized by tool
// ============================================================================
pub mod analyze;
pub mod chaos;
pub mod compute;
pub mod ml;
pub mod simulate;
pub mod solve;
pub mod units;
pub mod validate;

// WebAssembly bindings (when compiled with wasm feature)
#[cfg(feature = "wasm")]
pub mod wasm;

// ============================================================================
// BACKWARD COMPATIBILITY RE-EXPORTS
// ============================================================================
// These provide the old module paths for backward compatibility.
// New code should use the 8 tool modules directly.

// Mathematics domain -> compute/
pub use compute::calculus as advanced_calculus;
pub use compute::matrix::linear_algebra;
pub use compute::special_functions;
pub use compute::symbolic_regression;
pub use compute::tensor as tensor_calculus;

// Physics domain -> compute/physics/ and simulate/fluids/
pub use compute::physics::electromagnetism;
pub use simulate::fluids as fluid_dynamics;
pub use compute::physics::quantum as quantum_physics;

// Specialized domain -> various tools
pub use compute::cryptographic_mathematics;
pub use compute::graph as graph_theory;
pub use compute::information as information_theory;
pub use solve::optimization;
pub use compute::sampling as statistics;
pub use simulate::stochastic as stochastic_processes;

// Tools domain -> compute/
pub use compute::geometry as computational_geometry;
pub use compute::numerical_methods;
pub use compute::transforms as signal_processing;

// Units and validation
pub use units::dimensional_analysis;
pub use validate as equation_validation;

// Scientific formula modules at top level
pub use compute::biology;
pub use compute::chemistry;
pub use compute::datetime;
pub use compute::electrical;
pub use compute::engineering;
pub use compute::geophysics;
pub use compute::materials_science;
pub use compute::optics;
pub use compute::thermodynamics;

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

// Re-export operation listing from engine
pub use engine::{list_all_operations, list_tools};

pub use tensor_calculus::TensorError;

/// Common result type for computational operations
pub type ComputationalResult<T> = Result<T, Box<dyn std::error::Error>>;

/// Version information
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub const NAME: &str = env!("CARGO_PKG_NAME");

// ============================================================================
// DISPATCHER FACTORY
// ============================================================================

use analyze::UnifiedAnalyzer;
use chaos::UnifiedChaos;
use compute::UnifiedComputer;
use ml::UnifiedML;
use simulate::UnifiedSimulator;
use solve::UnifiedSolver;
use units::UnifiedUnits;
use validate::UnifiedValidator;

/// Create a fully configured ToolDispatcher with all 8 tool implementations
pub fn create_default_dispatcher() -> ToolDispatcher {
    ToolDispatcher::new(
        Box::new(UnifiedSolver::new()),
        Box::new(UnifiedComputer::new()),
        Box::new(UnifiedAnalyzer::new()),
        Box::new(UnifiedSimulator::new()),
        Box::new(UnifiedML::new()),
        Box::new(UnifiedChaos::new()),
        Box::new(UnifiedUnits::new()),
        Box::new(UnifiedValidator::new()),
    )
}
