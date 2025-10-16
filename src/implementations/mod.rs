//! Concrete implementations of the 10 core tools
//!
//! These implementations wire the new 10-tool API to the existing
//! 180+ domain-specific functions across all modules.

pub mod analyzer;
pub mod compute;
pub mod computer;
pub mod differentiator;
pub mod field_solver;
pub mod integrator;
pub mod optimizer;
pub mod sampler;
pub mod simulator;
pub mod solver;
pub mod transformer;

pub use analyzer::UnifiedAnalyzer;
pub use computer::UnifiedComputer;
pub use differentiator::UnifiedDifferentiator;
pub use field_solver::UnifiedFieldSolver;
pub use integrator::UnifiedIntegrator;
pub use optimizer::UnifiedOptimizer;
pub use sampler::UnifiedSampler;
pub use simulator::UnifiedSimulator;
pub use solver::UnifiedSolver;
pub use transformer::UnifiedTransformer;

use crate::engine::*;

/// Create a fully configured ToolDispatcher with all 10 tool implementations
pub fn create_default_dispatcher() -> ToolDispatcher {
    ToolDispatcher::new(
        Box::new(UnifiedSolver::new()),         // Solve
        Box::new(UnifiedDifferentiator::new()), // Differentiate
        Box::new(UnifiedIntegrator::new()),     // Integrate
        Box::new(UnifiedAnalyzer::new()),       // Analyze
        Box::new(UnifiedSimulator::new()),      // Simulate
        Box::new(UnifiedComputer::new()),       // Compute
        Box::new(UnifiedTransformer::new()),    // Transform
        Box::new(UnifiedFieldSolver::new()),    // FieldTheory
        Box::new(UnifiedSampler::new()),        // Sample
        Box::new(UnifiedOptimizer::new()),      // Optimize
    )
}
