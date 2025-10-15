//! Stub implementations for deprecated/legacy tools
//!
//! These stubs remain for backwards compatibility but are no longer used
//! in the default dispatcher. All functionality has been moved to the
//! main unified tool implementations.

use crate::engine::*;

// Legacy stubs - kept for backwards compatibility only
// Not used in create_default_dispatcher()

pub struct StubSimulator;
impl Simulate for StubSimulator {
    fn simulate(&self, _input: &SimulateInput) -> ToolResult<SimulateOutput> {
        Err("Use UnifiedSimulator instead".to_string())
    }
}

pub struct StubTransformer;
impl Transform for StubTransformer {
    fn transform(&self, _input: &TransformInput) -> ToolResult<TransformOutput> {
        Err("Use UnifiedTransformer instead".to_string())
    }
}

pub struct StubSampler;
impl Sample for StubSampler {
    fn sample(&self, _input: &SampleInput) -> ToolResult<SampleOutput> {
        Err("Use UnifiedSampler instead".to_string())
    }
}

pub struct StubOptimizer;
impl Optimize for StubOptimizer {
    fn optimize(&self, _input: &OptimizeInput) -> ToolResult<OptimizeOutput> {
        Err("Use UnifiedOptimizer instead".to_string())
    }
}
