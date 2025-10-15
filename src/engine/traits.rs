//! Core traits defining the 10-tool interface

use super::types::*;

/// Tool 1: Solve equations, systems, and mathematical problems
pub trait Solve {
    fn solve(&self, input: &SolveInput) -> ToolResult<SolveOutput>;
}

/// Tool 2: Compute derivatives, gradients, and differential operators
pub trait Differentiate {
    fn differentiate(&self, input: &DifferentiateInput) -> ToolResult<DifferentiateOutput>;
}

/// Tool 3: Compute integrals (definite, indefinite, line, surface, volume)
pub trait Integrate {
    fn integrate(&self, input: &IntegrateInput) -> ToolResult<IntegrateOutput>;
}

/// Tool 4: Analyze expressions (validate, simplify, parse, check)
pub trait Analyze {
    fn analyze(&self, input: &AnalyzeInput) -> ToolResult<AnalyzeOutput>;
}

/// Tool 5: Simulate time evolution and stochastic processes
pub trait Simulate {
    fn simulate(&self, input: &SimulateInput) -> ToolResult<SimulateOutput>;
}

/// Tool 6: Compute tensor operations, matrix decompositions, special functions
pub trait Compute {
    fn compute(&self, input: &ComputeInput) -> ToolResult<ComputeOutput>;
}

/// Tool 7: Apply transforms (Fourier, Laplace, wavelet, filters)
pub trait Transform {
    fn transform(&self, input: &TransformInput) -> ToolResult<TransformOutput>;
}

/// Tool 8: Work with physics fields (EM, gravity, quantum)
pub trait FieldTheory {
    fn field_theory(&self, input: &FieldTheoryInput) -> ToolResult<FieldTheoryOutput>;
}

/// Tool 9: Sample and analyze statistical data
pub trait Sample {
    fn sample(&self, input: &SampleInput) -> ToolResult<SampleOutput>;
}

/// Tool 10: Optimize and fit data to models
pub trait Optimize {
    fn optimize(&self, input: &OptimizeInput) -> ToolResult<OptimizeOutput>;
}
