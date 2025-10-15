//! 10-Tool Computational Engine Architecture
//!
//! This module provides the unified interface for all 180+ mathematical operations,
//! organized around 10 core tools that cover all domains of mathematics and physics.

pub mod traits;
pub mod types;
pub mod equations;
pub mod dispatcher;

pub use traits::{
    Solve, Differentiate, Integrate, Analyze, Simulate,
    Compute, Transform, FieldTheory, Sample, Optimize
};

pub use types::{
    // Solve
    SolveInput, SolveOutput,
    // Differentiate
    DifferentiateInput, DifferentiateOutput,
    // Integrate
    IntegrateInput, IntegrateOutput,
    // Analyze
    AnalyzeInput, AnalyzeOutput,
    // Simulate
    SimulateInput, SimulateOutput,
    // Compute
    ComputeInput, ComputeOutput,
    // Transform
    TransformInput, TransformOutput,
    // FieldTheory
    FieldTheoryInput, FieldTheoryOutput,
    // Sample
    SampleInput, SampleOutput,
    // Optimize
    OptimizeInput, OptimizeOutput,
    // Common
    ToolResult, Domain, Method,
};

pub use equations::*;
pub use dispatcher::{ToolDispatcher, ToolRequest, ToolResponse};
