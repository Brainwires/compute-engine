//! Central dispatcher for routing tool requests
//!
//! 8-tool architecture:
//! - Solve: Equations, optimization, root finding, curve fitting
//! - Compute: All numerical calculations (calculus, transforms, physics, statistics)
//! - Analyze: Series, limits, stability analysis, simplification
//! - Simulate: Time evolution, stochastic processes, fluid dynamics
//! - ML: Machine learning (clustering, neural networks, regression, dimensionality reduction)
//! - Chaos: Chaos theory (fractals, attractors, Lyapunov exponents, bifurcation)
//! - Units: Dimensional analysis and unit conversion
//! - Validate: Equation and physics validation

use super::traits::*;
use super::types::*;
#[cfg(feature = "schemars")]
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::time::Instant;

/// Unified tool request envelope
///
/// Tools (lowercase in JSON): solve, compute, analyze, simulate, ml, chaos, units, validate
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(tag = "tool", content = "input")]
#[serde(rename_all = "lowercase")]
pub enum ToolRequest {
    Solve(SolveInput),
    Compute(ComputeInput),
    Analyze(AnalyzeInput),
    Simulate(SimulateInput),
    /// Machine Learning: clustering, neural networks, regression, dimensionality reduction
    #[serde(alias = "ML")]
    Ml(MLInput),
    /// Chaos Theory: fractals, attractors, Lyapunov exponents, bifurcation
    Chaos(ChaosInput),
    /// Dimensional Analysis: unit conversion and validation
    Units(UnitsInput),
    /// Validation: equation/physics validation
    Validate(ValidateInput),
}

/// Unified tool response envelope
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(tag = "tool", content = "output")]
#[serde(rename_all = "lowercase")]
pub enum ToolResponse {
    Solve(SolveOutput),
    Compute(ComputeOutput),
    Analyze(AnalyzeOutput),
    Simulate(SimulateOutput),
    /// Machine Learning response
    #[serde(alias = "ML")]
    Ml(MLOutput),
    /// Chaos Theory response
    Chaos(ChaosOutput),
    /// Dimensional Analysis response
    Units(UnitsOutput),
    /// Validation response
    Validate(ValidateOutput),
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
///
/// 8-tool architecture.
pub struct ToolDispatcher {
    solver: Box<dyn Solve + Send + Sync>,
    computer: Box<dyn Compute + Send + Sync>,
    analyzer: Box<dyn Analyze + Send + Sync>,
    simulator: Box<dyn Simulate + Send + Sync>,
    ml: Box<dyn MachineLearning + Send + Sync>,
    chaos: Box<dyn Chaos + Send + Sync>,
    units: Box<dyn Units + Send + Sync>,
    validator: Box<dyn Validate + Send + Sync>,
}

impl ToolDispatcher {
    /// Create a new dispatcher with all 8 tool implementations
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        solver: Box<dyn Solve + Send + Sync>,
        computer: Box<dyn Compute + Send + Sync>,
        analyzer: Box<dyn Analyze + Send + Sync>,
        simulator: Box<dyn Simulate + Send + Sync>,
        ml: Box<dyn MachineLearning + Send + Sync>,
        chaos: Box<dyn Chaos + Send + Sync>,
        units: Box<dyn Units + Send + Sync>,
        validator: Box<dyn Validate + Send + Sync>,
    ) -> Self {
        Self {
            solver,
            computer,
            analyzer,
            simulator,
            ml,
            chaos,
            units,
            validator,
        }
    }

    /// Process a tool request and return the appropriate response
    pub fn dispatch(&self, request: ToolRequest) -> ToolResult<ToolResponse> {
        match request {
            ToolRequest::Solve(input) => self.solver.solve(&input).map(ToolResponse::Solve),
            ToolRequest::Compute(input) => self.computer.compute(&input).map(ToolResponse::Compute),
            ToolRequest::Analyze(input) => self.analyzer.analyze(&input).map(ToolResponse::Analyze),
            ToolRequest::Simulate(input) => {
                self.simulator.simulate(&input).map(ToolResponse::Simulate)
            }
            ToolRequest::Ml(input) => self.ml.ml(&input).map(ToolResponse::Ml),
            ToolRequest::Chaos(input) => self.chaos.chaos(&input).map(ToolResponse::Chaos),
            ToolRequest::Units(input) => self.units.units(&input).map(ToolResponse::Units),
            ToolRequest::Validate(input) => {
                self.validator.validate(&input).map(ToolResponse::Validate)
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
}

/// Get the list of available tools (8-tool architecture)
pub fn available_tools() -> Vec<&'static str> {
    vec![
        "solve", "compute", "analyze", "simulate", "ml", "chaos", "units", "validate",
    ]
}
