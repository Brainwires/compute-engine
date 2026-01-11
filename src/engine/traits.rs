//! Core traits defining the 8-tool interface
//!
//! These 8 tools provide a complete interface for mathematical and scientific computing:
//! - Solve: Equations, systems, optimization
//! - Compute: Calculus, transforms, fields, matrices, sampling
//! - Analyze: Simplify, series, limits, stability
//! - Simulate: ODEs, stochastic, fluid dynamics
//! - MachineLearning: Clustering, neural nets, regression
//! - Chaos: Fractals, attractors, Lyapunov
//! - Units: Dimensional analysis, conversion
//! - Validate: Equation and physics validation

use super::types::*;

/// Tool 1: Solve equations, systems, and optimization problems
pub trait Solve {
    fn solve(&self, input: &SolveInput) -> ToolResult<SolveOutput>;
}

/// Tool 2: Compute mathematical operations (calculus, transforms, fields, matrices, sampling)
pub trait Compute {
    fn compute(&self, input: &ComputeInput) -> ToolResult<ComputeOutput>;
}

/// Tool 3: Analyze expressions (simplify, series, limits, stability)
pub trait Analyze {
    fn analyze(&self, input: &AnalyzeInput) -> ToolResult<AnalyzeOutput>;
}

/// Tool 4: Simulate time evolution and stochastic processes
pub trait Simulate {
    fn simulate(&self, input: &SimulateInput) -> ToolResult<SimulateOutput>;
}

/// Tool 5: Machine Learning operations (clustering, neural networks, regression)
pub trait MachineLearning {
    fn ml(&self, input: &MLInput) -> ToolResult<MLOutput>;
}

/// Tool 6: Chaos theory (fractals, attractors, Lyapunov exponents)
pub trait Chaos {
    fn chaos(&self, input: &ChaosInput) -> ToolResult<ChaosOutput>;
}

/// Tool 7: Units and dimensional analysis
pub trait Units {
    fn units(&self, input: &UnitsInput) -> ToolResult<UnitsOutput>;
}

/// Tool 8: Equation and physics validation
pub trait Validate {
    fn validate(&self, input: &ValidateInput) -> ToolResult<ValidateOutput>;
}
