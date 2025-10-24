//! 10-Tool Computational Engine Architecture
//!
//! This module provides the unified interface for all 180+ mathematical operations,
//! organized around 10 core tools that cover all domains of mathematics and physics.

pub mod dispatcher;
pub mod equations;
pub mod traits;
pub mod types;

pub use traits::{
    Analyze, Compute, Differentiate, FieldTheory, Integrate, Optimize, Sample, Simulate, Solve,
    Transform,
};

pub use types::{
    // Analyze
    AnalyzeInput,
    AnalyzeOutput,
    // Compute
    ComputeInput,
    ComputeOutput,
    // Differentiate
    DifferentiateInput,
    DifferentiateOutput,
    Domain,
    // FieldTheory
    FieldTheoryInput,
    FieldTheoryOutput,
    // Integrate
    IntegrateInput,
    IntegrateOutput,
    Method,
    // Optimize
    OptimizeInput,
    OptimizeOutput,
    // Sample
    SampleInput,
    SampleOutput,
    // Simulate
    SimulateInput,
    SimulateOutput,
    // Solve
    SolveInput,
    SolveOutput,
    // Common
    ToolResult,
    // Transform
    TransformInput,
    TransformOutput,
};

pub use dispatcher::{TimedResponse, ToolDispatcher, ToolRequest, ToolResponse};
pub use equations::*;
