//! Central dispatcher for routing tool requests

use super::traits::*;
use super::types::*;
#[cfg(feature = "schemars")]
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::time::Instant;

/// Unified tool request envelope
///
/// All tool names are **lowercase** in JSON: solve, differentiate, integrate,
/// analyze, simulate, compute, transform, fieldtheory, sample, optimize.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
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
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
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

/// Response wrapper with timing and metadata for billing
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimedResponse {
    /// Success status
    pub success: bool,

    /// The actual response (if successful)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub response: Option<ToolResponse>,

    /// Error message (if failed)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,

    /// Compute time in milliseconds
    pub compute_time_ms: f64,

    /// Compute time in microseconds (for high-precision billing)
    pub compute_time_us: u128,

    /// Compute time in nanoseconds (most precise)
    pub compute_time_ns: u128,
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

    /// Process a JSON request string with timing information
    pub fn dispatch_json(&self, json_str: &str) -> String {
        let request: ToolRequest = match serde_json::from_str(json_str) {
            Ok(req) => req,
            Err(e) => {
                // Parsing error - no compute time
                let response = TimedResponse {
                    success: false,
                    response: None,
                    error: Some(format!("Invalid JSON request: {}", e)),
                    compute_time_ms: 0.0,
                    compute_time_us: 0,
                    compute_time_ns: 0,
                };
                return serde_json::to_string(&response).unwrap_or_else(|_| {
                    r#"{"success":false,"error":"Failed to serialize error response"}"#.to_string()
                });
            }
        };

        // Start timing
        let start = Instant::now();

        // Execute the request
        let result = self.dispatch(request);

        // Measure elapsed time
        let elapsed = start.elapsed();
        let compute_time_ns = elapsed.as_nanos();
        let compute_time_us = elapsed.as_micros();
        let compute_time_ms = elapsed.as_secs_f64() * 1000.0;

        // Build timed response
        let timed_response = match result {
            Ok(response) => TimedResponse {
                success: true,
                response: Some(response),
                error: None,
                compute_time_ms,
                compute_time_us,
                compute_time_ns,
            },
            Err(e) => TimedResponse {
                success: false,
                response: None,
                error: Some(e),
                compute_time_ms,
                compute_time_us,
                compute_time_ns,
            },
        };

        serde_json::to_string(&timed_response).unwrap_or_else(|e| {
            serde_json::json!({
                "success": false,
                "error": format!("Failed to serialize response: {}", e)
            })
            .to_string()
        })
    }

    /// Process a JSON request string (legacy method without timing - for backward compatibility)
    #[deprecated(
        since = "0.1.0",
        note = "Use dispatch_json for timing information. This method is kept for backward compatibility."
    )]
    pub fn dispatch_json_legacy(&self, json_str: &str) -> String {
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
