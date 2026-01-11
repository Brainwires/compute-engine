//! Central dispatcher for routing tool requests
//!
//! Expanded 8-tool architecture:
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
/// Primary tools (lowercase in JSON): solve, compute, analyze, simulate, ml, chaos, units, validate
///
/// Legacy tools (for backward compatibility, routed to primary tools):
/// - differentiate -> compute
/// - integrate -> compute
/// - transform -> compute
/// - fieldtheory -> compute
/// - sample -> compute
/// - optimize -> solve
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(tag = "tool", content = "input")]
#[serde(rename_all = "lowercase")]
pub enum ToolRequest {
    // ===== Primary 8 Tools =====
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

    // ===== Legacy Tools (backward compatibility) =====
    /// Legacy: Use compute with operation: {differentiate: ...}
    Differentiate(DifferentiateInput),
    /// Legacy: Use compute with operation: {integrate: ...}
    Integrate(IntegrateInput),
    /// Legacy: Use compute with operation: {transform: ...}
    Transform(TransformInput),
    /// Legacy: Use compute with operation: {field: ...}
    FieldTheory(FieldTheoryInput),
    /// Legacy: Use compute with operation: {sample: ...}
    Sample(SampleInput),
    /// Legacy: Use solve with equation_type: {optimize: ...}
    Optimize(OptimizeInput),
}

/// Unified tool response envelope
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(tag = "tool", content = "output")]
#[serde(rename_all = "lowercase")]
pub enum ToolResponse {
    // Primary 8 tool responses
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

    // Legacy responses (for backward compatibility)
    Differentiate(DifferentiateOutput),
    Integrate(IntegrateOutput),
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
///
/// Expanded 8-tool architecture. Legacy tools are routed internally.
pub struct ToolDispatcher {
    // ===== Primary 8 Tools =====
    solver: Box<dyn Solve + Send + Sync>,
    computer: Box<dyn Compute + Send + Sync>,
    analyzer: Box<dyn Analyze + Send + Sync>,
    simulator: Box<dyn Simulate + Send + Sync>,
    ml: Box<dyn MachineLearning + Send + Sync>,
    chaos: Box<dyn Chaos + Send + Sync>,
    units: Box<dyn Units + Send + Sync>,
    validator: Box<dyn Validate + Send + Sync>,

    // Legacy tool implementations (for backward compatibility routing)
    differentiator: Box<dyn Differentiate + Send + Sync>,
    integrator: Box<dyn Integrate + Send + Sync>,
    transformer: Box<dyn Transform + Send + Sync>,
    field_theory: Box<dyn FieldTheory + Send + Sync>,
    sampler: Box<dyn Sample + Send + Sync>,
    optimizer: Box<dyn Optimize + Send + Sync>,
}

impl ToolDispatcher {
    /// Create a new dispatcher with the 8 primary tool implementations
    #[allow(clippy::too_many_arguments)]
    pub fn new_consolidated(
        solver: Box<dyn Solve + Send + Sync>,
        computer: Box<dyn Compute + Send + Sync>,
        analyzer: Box<dyn Analyze + Send + Sync>,
        simulator: Box<dyn Simulate + Send + Sync>,
        ml: Box<dyn MachineLearning + Send + Sync>,
        chaos: Box<dyn Chaos + Send + Sync>,
        units: Box<dyn Units + Send + Sync>,
        validator: Box<dyn Validate + Send + Sync>,
        // Legacy implementations for backward compatibility
        differentiator: Box<dyn Differentiate + Send + Sync>,
        integrator: Box<dyn Integrate + Send + Sync>,
        transformer: Box<dyn Transform + Send + Sync>,
        field_theory: Box<dyn FieldTheory + Send + Sync>,
        sampler: Box<dyn Sample + Send + Sync>,
        optimizer: Box<dyn Optimize + Send + Sync>,
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
            differentiator,
            integrator,
            transformer,
            field_theory,
            sampler,
            optimizer,
        }
    }

    /// Legacy constructor - kept for backward compatibility
    ///
    /// **Deprecated**: Use `new_consolidated` instead.
    /// This constructor now requires all 14 tool implementations (8 primary + 6 legacy).
    #[deprecated(
        since = "0.2.0",
        note = "Use new_consolidated() instead. This constructor signature has changed to include ML, Chaos, Units, and Validate tools."
    )]
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
        // New required tools
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
            differentiator,
            integrator,
            transformer,
            field_theory,
            sampler,
            optimizer,
        }
    }

    /// Process a tool request and return the appropriate response
    pub fn dispatch(&self, request: ToolRequest) -> ToolResult<ToolResponse> {
        match request {
            // ===== Primary 8 Tools =====
            ToolRequest::Solve(input) => self.solver.solve(&input).map(ToolResponse::Solve),
            ToolRequest::Compute(input) => self.computer.compute(&input).map(ToolResponse::Compute),
            ToolRequest::Analyze(input) => self.analyzer.analyze(&input).map(ToolResponse::Analyze),
            ToolRequest::Simulate(input) => self.simulator.simulate(&input).map(ToolResponse::Simulate),
            ToolRequest::Ml(input) => self.ml.ml(&input).map(ToolResponse::Ml),
            ToolRequest::Chaos(input) => self.chaos.chaos(&input).map(ToolResponse::Chaos),
            ToolRequest::Units(input) => self.units.units(&input).map(ToolResponse::Units),
            ToolRequest::Validate(input) => self.validator.validate(&input).map(ToolResponse::Validate),

            // ===== Legacy Tools (backward compatibility) =====
            ToolRequest::Differentiate(input) => self
                .differentiator
                .differentiate(&input)
                .map(ToolResponse::Differentiate),
            ToolRequest::Integrate(input) => self
                .integrator
                .integrate(&input)
                .map(ToolResponse::Integrate),
            ToolRequest::Transform(input) => self
                .transformer
                .transform(&input)
                .map(ToolResponse::Transform),
            ToolRequest::FieldTheory(input) => self
                .field_theory
                .field_theory(&input)
                .map(ToolResponse::FieldTheory),
            ToolRequest::Sample(input) => self.sampler.sample(&input).map(ToolResponse::Sample),
            ToolRequest::Optimize(input) => self
                .optimizer
                .optimize(&input)
                .map(ToolResponse::Optimize),
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

/// Get the list of primary tools (8-tool architecture)
pub fn primary_tools() -> Vec<&'static str> {
    vec!["solve", "compute", "analyze", "simulate", "ml", "chaos", "units", "validate"]
}

/// Get the list of legacy tools (for backward compatibility)
pub fn legacy_tools() -> Vec<&'static str> {
    vec!["differentiate", "integrate", "transform", "fieldtheory", "sample", "optimize"]
}

/// Get all available tools (primary + legacy)
pub fn all_tools() -> Vec<&'static str> {
    let mut tools = primary_tools();
    tools.extend(legacy_tools());
    tools
}
