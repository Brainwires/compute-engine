//! Concrete implementations of the 8 tools
//!
//! 8 Tools:
//! - Solve: Equations, optimization, root finding
//! - Compute: Matrix, calculus, transforms, field theory, sampling
//! - Analyze: Series, limits, stability analysis
//! - Simulate: Time evolution, stochastic, fluid dynamics
//! - ML: Machine learning (clustering, neural nets, regression)
//! - Chaos: Chaos theory (fractals, attractors, Lyapunov)
//! - Units: Dimensional analysis and unit conversion
//! - Validate: Equation and physics validation

pub mod analyzer;
pub mod chaos;
pub mod compute;
pub mod computer;
pub mod ml;
pub mod simulator;
pub mod solver;
pub mod units;
pub mod validator;

pub use analyzer::UnifiedAnalyzer;
pub use chaos::UnifiedChaos;
pub use computer::UnifiedComputer;
pub use ml::UnifiedML;
pub use simulator::UnifiedSimulator;
pub use solver::UnifiedSolver;
pub use units::UnifiedUnits;
pub use validator::UnifiedValidator;

use crate::engine::*;

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
