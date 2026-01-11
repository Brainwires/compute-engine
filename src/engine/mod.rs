//! 8-Tool Computational Engine Architecture
//!
//! This module provides the unified interface for all 350+ mathematical operations,
//! organized around 8 primary tools that cover all domains of mathematics and physics.
//!
//! ## Primary 8 Tools
//! - Solve: Equations, systems, optimization, root finding
//! - Compute: Matrix ops, calculus, transforms, field theory, sampling
//! - Analyze: Series, limits, stability analysis, simplification
//! - Simulate: Time evolution, stochastic processes, fluid dynamics
//! - ML: Machine learning (clustering, neural nets, regression)
//! - Chaos: Chaos theory (fractals, attractors, Lyapunov)
//! - Units: Dimensional analysis and unit conversion
//! - Validate: Equation and physics validation
//!
//! ## Legacy Tools (Backward Compatible)
//! Differentiate, Integrate, Transform, FieldTheory, Sample, Optimize

pub mod dispatcher;
pub mod equations;
pub mod operations;
pub mod traits;
pub mod types;

// Primary 8 tool traits
pub use traits::{
    Analyze, Chaos, Compute, MachineLearning, Simulate, Solve, Units, Validate,
};

// Legacy tool traits (backward compatibility)
pub use traits::{Differentiate, FieldTheory, Integrate, Optimize, Sample, Transform};

pub use types::{
    // Primary 8 tools - Input/Output
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
    // Legacy tools - Input/Output
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
    // Common types
    Domain,
    Method,
    ToolResult,
};

pub use dispatcher::{TimedResponse, ToolDispatcher, ToolRequest, ToolResponse};
pub use equations::*;
pub use operations::{list_all_operations, list_tools};
