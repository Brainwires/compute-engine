//! Central dispatcher for routing tool requests

use super::traits::*;
use super::types::*;
use serde::{Deserialize, Serialize};

/// Unified tool request envelope
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "tool", content = "input")]
#[serde(rename_all = "lowercase")]
pub enum ToolRequest {
    Solve(SolveInput),
    Differentiate(DifferentiateInput),
    Integrate(IntegrateInput),
    Analyze(AnalyzeInput),
    Simulate(SimulateInput),
    Compute(ComputeInput),
    Transform(TransformInput),
    FieldTheory(FieldTheoryInput),
    Sample(SampleInput),
    Optimize(OptimizeInput),
}

/// Unified tool response envelope
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "tool", content = "output")]
#[serde(rename_all = "lowercase")]
pub enum ToolResponse {
    Solve(SolveOutput),
    Differentiate(DifferentiateOutput),
    Integrate(IntegrateOutput),
    Analyze(AnalyzeOutput),
    Simulate(SimulateOutput),
    Compute(ComputeOutput),
    Transform(TransformOutput),
    FieldTheory(FieldTheoryOutput),
    Sample(SampleOutput),
    Optimize(OptimizeOutput),
}

/// Main dispatcher that routes requests to appropriate tool implementations
pub struct ToolDispatcher {
    solver: Box<dyn Solve + Send + Sync>,
    differentiator: Box<dyn Differentiate + Send + Sync>,
    integrator: Box<dyn Integrate + Send + Sync>,
    analyzer: Box<dyn Analyze + Send + Sync>,
    simulator: Box<dyn Simulate + Send + Sync>,
    computer: Box<dyn Compute + Send + Sync>,
    transformer: Box<dyn Transform + Send + Sync>,
    field_theory: Box<dyn FieldTheory + Send + Sync>,
    sampler: Box<dyn Sample + Send + Sync>,
    optimizer: Box<dyn Optimize + Send + Sync>,
}

impl ToolDispatcher {
    /// Create a new dispatcher with the given tool implementations
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        solver: Box<dyn Solve + Send + Sync>,
        differentiator: Box<dyn Differentiate + Send + Sync>,
        integrator: Box<dyn Integrate + Send + Sync>,
        analyzer: Box<dyn Analyze + Send + Sync>,
        simulator: Box<dyn Simulate + Send + Sync>,
        computer: Box<dyn Compute + Send + Sync>,
        transformer: Box<dyn Transform + Send + Sync>,
        field_theory: Box<dyn FieldTheory + Send + Sync>,
        sampler: Box<dyn Sample + Send + Sync>,
        optimizer: Box<dyn Optimize + Send + Sync>,
    ) -> Self {
        Self {
            solver,
            differentiator,
            integrator,
            analyzer,
            simulator,
            computer,
            transformer,
            field_theory,
            sampler,
            optimizer,
        }
    }

    /// Process a tool request and return the appropriate response
    pub fn dispatch(&self, request: ToolRequest) -> ToolResult<ToolResponse> {
        match request {
            ToolRequest::Solve(input) => self.solver.solve(&input).map(ToolResponse::Solve),
            ToolRequest::Differentiate(input) => self
                .differentiator
                .differentiate(&input)
                .map(ToolResponse::Differentiate),
            ToolRequest::Integrate(input) => self
                .integrator
                .integrate(&input)
                .map(ToolResponse::Integrate),
            ToolRequest::Analyze(input) => self.analyzer.analyze(&input).map(ToolResponse::Analyze),
            ToolRequest::Simulate(input) => {
                self.simulator.simulate(&input).map(ToolResponse::Simulate)
            }
            ToolRequest::Compute(input) => self.computer.compute(&input).map(ToolResponse::Compute),
            ToolRequest::Transform(input) => self
                .transformer
                .transform(&input)
                .map(ToolResponse::Transform),
            ToolRequest::FieldTheory(input) => self
                .field_theory
                .field_theory(&input)
                .map(ToolResponse::FieldTheory),
            ToolRequest::Sample(input) => self.sampler.sample(&input).map(ToolResponse::Sample),
            ToolRequest::Optimize(input) => {
                self.optimizer.optimize(&input).map(ToolResponse::Optimize)
            }
        }
    }

    /// Process a JSON request string
    pub fn dispatch_json(&self, json_str: &str) -> String {
        let request: ToolRequest = match serde_json::from_str(json_str) {
            Ok(req) => req,
            Err(e) => {
                return serde_json::json!({
                    "success": false,
                    "error": format!("Invalid JSON request: {}", e)
                })
                .to_string();
            }
        };

        match self.dispatch(request) {
            Ok(response) => serde_json::json!({
                "success": true,
                "response": response
            })
            .to_string(),
            Err(e) => serde_json::json!({
                "success": false,
                "error": e
            })
            .to_string(),
        }
    }
}
