//! Concrete implementations of the 8 primary tools + 6 legacy tools
//!
//! Primary 8 Tools:
//! - Solve: Equations, optimization, root finding
//! - Compute: Matrix, calculus, transforms, field theory, sampling
//! - Analyze: Series, limits, stability analysis
//! - Simulate: Time evolution, stochastic, fluid dynamics
//! - ML: Machine learning (clustering, neural nets, regression)
//! - Chaos: Chaos theory (fractals, attractors, Lyapunov)
//! - Units: Dimensional analysis and unit conversion
//! - Validate: Equation and physics validation
//!
//! Legacy Tools (backward compatibility, route to primary):
//! - Differentiate, Integrate, Transform, FieldTheory, Sample, Optimize

pub mod analyzer;
pub mod chaos;
pub mod compute;
pub mod computer;
pub mod differentiator;
pub mod field_solver;
pub mod integrator;
pub mod ml;
pub mod optimizer;
pub mod sampler;
pub mod simulator;
pub mod solver;
pub mod transformer;
pub mod units;
pub mod validator;

// Primary 8 tools
pub use analyzer::UnifiedAnalyzer;
pub use chaos::UnifiedChaos;
pub use computer::UnifiedComputer;
pub use ml::UnifiedML;
pub use simulator::UnifiedSimulator;
pub use solver::UnifiedSolver;
pub use units::UnifiedUnits;
pub use validator::UnifiedValidator;

// Legacy tools (backward compatibility)
pub use differentiator::UnifiedDifferentiator;
pub use field_solver::UnifiedFieldSolver;
pub use integrator::UnifiedIntegrator;
pub use optimizer::UnifiedOptimizer;
pub use sampler::UnifiedSampler;
pub use transformer::UnifiedTransformer;

use crate::engine::*;

/// Create a fully configured ToolDispatcher with all 8 primary + 6 legacy tool implementations
pub fn create_default_dispatcher() -> ToolDispatcher {
    ToolDispatcher::new_consolidated(
        // Primary 8 tools
        Box::new(UnifiedSolver::new()),     // Solve
        Box::new(UnifiedComputer::new()),   // Compute
        Box::new(UnifiedAnalyzer::new()),   // Analyze
        Box::new(UnifiedSimulator::new()),  // Simulate
        Box::new(UnifiedML::new()),         // ML
        Box::new(UnifiedChaos::new()),      // Chaos
        Box::new(UnifiedUnits::new()),      // Units
        Box::new(UnifiedValidator::new()),  // Validate
        // Legacy tools (backward compatibility)
        Box::new(UnifiedDifferentiator::new()), // Differentiate
        Box::new(UnifiedIntegrator::new()),     // Integrate
        Box::new(UnifiedTransformer::new()),    // Transform
        Box::new(UnifiedFieldSolver::new()),    // FieldTheory
        Box::new(UnifiedSampler::new()),        // Sample
        Box::new(UnifiedOptimizer::new()),      // Optimize
    )
}
