//! 8-Tool Computational Engine Architecture
//!
//! This module provides the unified interface for all 350+ mathematical operations,
//! organized around 8 primary tools that cover all domains of mathematics and physics.
//!
//! ## 8 Tools
//! - Solve: Equations, systems, optimization, root finding
//! - Compute: Matrix ops, calculus, transforms, field theory, sampling
//! - Analyze: Series, limits, stability analysis, simplification
//! - Simulate: Time evolution, stochastic processes, fluid dynamics
//! - ML: Machine learning (clustering, neural nets, regression)
//! - Chaos: Chaos theory (fractals, attractors, Lyapunov)
//! - Units: Dimensional analysis and unit conversion
//! - Validate: Equation and physics validation

pub mod dispatcher;
pub mod equations;
pub mod operations;
pub mod traits;
pub mod types;

// 8 tool traits
pub use traits::{Analyze, Chaos, Compute, MachineLearning, Simulate, Solve, Units, Validate};

// 8 tools - Input/Output types
pub use types::{
    AnalyzeInput, AnalyzeOutput, ChaosInput, ChaosOutput, ComputeInput, ComputeOutput, MLInput,
    MLOutput, SimulateInput, SimulateOutput, SolveInput, SolveOutput, UnitsInput, UnitsOutput,
    ValidateInput, ValidateOutput,
    // Common types
    Domain, Method, ToolResult, ValidationResult,
};

pub use dispatcher::{TimedResponse, ToolDispatcher, ToolRequest, ToolResponse};
pub use equations::*;
pub use operations::{list_all_operations, list_tools};
