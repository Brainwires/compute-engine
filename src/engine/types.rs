//! Type definitions for the 10-tool architecture

use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;
use super::equations::*;

/// Common result type for all tool operations
pub type ToolResult<T> = Result<T, String>;

/// Mathematical domain (real or complex)
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Domain {
    Real,
    Complex,
}

/// Computation method preference
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Method {
    Symbolic,
    Numeric,
    Auto,
}

// ============================================================================
// TOOL 1: SOLVE
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolveInput {
    /// Type of equation to solve
    pub equation_type: EquationType,

    /// Equations (as strings or structured data)
    #[serde(default)]
    pub equations: Vec<String>,

    /// Variables to solve for
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub variables: Option<Vec<String>>,

    /// Initial guess for numeric solving
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub initial_guess: Option<HashMap<String, f64>>,

    /// Boundary conditions (for BVPs, PDEs)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub boundary_conditions: Option<HashMap<String, Value>>,

    /// Domain (real or complex)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub domain: Option<Domain>,

    /// Method preference
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub method: Option<Method>,

    /// Additional parameters
    #[serde(default)]
    pub parameters: HashMap<String, Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolveOutput {
    /// Solution sets
    pub solutions: Vec<HashMap<String, Value>>,

    /// Symbolic representation
    #[serde(skip_serializing_if = "Option::is_none")]
    pub symbolic: Option<String>,

    /// Numeric solutions
    #[serde(skip_serializing_if = "Option::is_none")]
    pub numeric: Option<Vec<HashMap<String, f64>>>,

    /// Step-by-step solution
    #[serde(skip_serializing_if = "Option::is_none")]
    pub steps: Option<Vec<String>>,

    /// Additional metadata
    #[serde(skip_serializing_if = "Option::is_none")]
    pub metadata: Option<Value>,
}

// ============================================================================
// TOOL 2: DIFFERENTIATE
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DifferentiateInput {
    /// Type of differentiation
    pub operation: DifferentiationOp,

    /// Expression to differentiate
    pub expression: String,

    /// Variables to differentiate with respect to
    pub variables: Vec<String>,

    /// Order of derivative (per variable)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub order: Option<Vec<usize>>,

    /// Point at which to evaluate
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub evaluate_at: Option<HashMap<String, f64>>,

    /// Additional parameters
    #[serde(default)]
    pub parameters: HashMap<String, Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DifferentiateOutput {
    /// Derivatives (symbolic or evaluated)
    pub derivatives: HashMap<String, Value>,

    /// LaTeX representations
    #[serde(skip_serializing_if = "Option::is_none")]
    pub latex: Option<HashMap<String, String>>,

    /// Additional results (Jacobian, Hessian, etc.)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub metadata: Option<Value>,
}

// ============================================================================
// TOOL 3: INTEGRATE
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntegrateInput {
    /// Type of integration
    pub integration_type: IntegrationType,

    /// Expression to integrate
    pub expression: String,

    /// Integration variables (order matters)
    pub variables: Vec<String>,

    /// Integration limits
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub limits: Option<Vec<[f64; 2]>>,

    /// Path (for line integrals)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path: Option<Value>,

    /// Method preference
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub method: Option<Method>,

    /// Additional parameters
    #[serde(default)]
    pub parameters: HashMap<String, Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntegrateOutput {
    /// Result (symbolic expression or numeric value)
    pub result: Value,

    /// Symbolic form
    #[serde(skip_serializing_if = "Option::is_none")]
    pub symbolic: Option<String>,

    /// LaTeX representation
    #[serde(skip_serializing_if = "Option::is_none")]
    pub latex: Option<String>,

    /// Error estimate (for numeric integration)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error_estimate: Option<f64>,

    /// Additional metadata
    #[serde(skip_serializing_if = "Option::is_none")]
    pub metadata: Option<Value>,
}

// ============================================================================
// TOOL 4: ANALYZE
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalyzeInput {
    /// Operation to perform
    pub operation: AnalysisOp,

    /// Expression or data to analyze
    pub expression: String,

    /// Context-specific options
    #[serde(default)]
    pub options: HashMap<String, Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalyzeOutput {
    /// Result (type depends on operation)
    pub result: Value,

    /// LaTeX representation
    #[serde(skip_serializing_if = "Option::is_none")]
    pub latex: Option<String>,

    /// Validation results
    #[serde(skip_serializing_if = "Option::is_none")]
    pub validation: Option<ValidationResult>,

    /// Additional details
    #[serde(skip_serializing_if = "Option::is_none")]
    pub details: Option<Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationResult {
    pub is_valid: bool,
    pub errors: Vec<String>,
    pub warnings: Vec<String>,
}

// ============================================================================
// TOOL 5: SIMULATE
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulateInput {
    /// Type of simulation
    pub model: SimulationModel,

    /// Equations defining the system
    pub equations: Vec<String>,

    /// System variables
    pub variables: Vec<String>,

    /// System parameters
    #[serde(default)]
    pub parameters: HashMap<String, f64>,

    /// Initial conditions
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub initial_conditions: Option<HashMap<String, f64>>,

    /// Time/parameter range [start, end]
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub range: Option<[f64; 2]>,

    /// Number of steps
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub steps: Option<usize>,

    /// Solver method
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub method: Option<String>,

    /// Number of paths (for stochastic)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub num_paths: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulateOutput {
    /// Simulation results (variable name -> values over time)
    pub results: HashMap<String, Vec<f64>>,

    /// Time points
    #[serde(skip_serializing_if = "Option::is_none")]
    pub time: Option<Vec<f64>>,

    /// Statistical moments (for stochastic)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub moments: Option<HashMap<String, f64>>,

    /// Plot data or URLs
    #[serde(skip_serializing_if = "Option::is_none")]
    pub plots: Option<Vec<String>>,

    /// Additional metadata
    #[serde(skip_serializing_if = "Option::is_none")]
    pub metadata: Option<Value>,
}

// ============================================================================
// TOOL 6: COMPUTE
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComputeInput {
    /// Operation to perform
    pub operation: ComputeOp,

    /// Input data (matrices, tensors, values, etc.)
    pub data: Value,

    /// Additional parameters
    #[serde(default)]
    pub parameters: HashMap<String, Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComputeOutput {
    /// Computed result
    pub result: Value,

    /// Additional outputs (e.g., eigenvalues AND eigenvectors)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub additional: Option<HashMap<String, Value>>,

    /// Metadata (convergence info, condition number, etc.)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub metadata: Option<Value>,
}

// ============================================================================
// TOOL 7: TRANSFORM
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransformInput {
    /// Type of transform
    pub transform_type: TransformType,

    /// Input signal/function
    pub data: Vec<f64>,

    /// Sampling rate (if applicable)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sampling_rate: Option<f64>,

    /// Additional parameters
    #[serde(default)]
    pub parameters: HashMap<String, Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransformOutput {
    /// Transformed data
    pub result: Vec<f64>,

    /// Frequencies (for FFT, etc.)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub frequencies: Option<Vec<f64>>,

    /// Magnitude spectrum
    #[serde(skip_serializing_if = "Option::is_none")]
    pub magnitude: Option<Vec<f64>>,

    /// Phase spectrum
    #[serde(skip_serializing_if = "Option::is_none")]
    pub phase: Option<Vec<f64>>,

    /// Additional metadata
    #[serde(skip_serializing_if = "Option::is_none")]
    pub metadata: Option<Value>,
}

// ============================================================================
// TOOL 8: FIELDTHEORY
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FieldTheoryInput {
    /// Type of field
    pub field_type: FieldType,

    /// Field configuration
    pub configuration: HashMap<String, Value>,

    /// Spatial points to evaluate
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub points: Option<Vec<Vec<f64>>>,

    /// Additional parameters
    #[serde(default)]
    pub parameters: HashMap<String, Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FieldTheoryOutput {
    /// Field values at requested points
    pub field_values: Vec<Value>,

    /// Vector/tensor field components
    #[serde(skip_serializing_if = "Option::is_none")]
    pub components: Option<HashMap<String, Vec<f64>>>,

    /// Additional metadata
    #[serde(skip_serializing_if = "Option::is_none")]
    pub metadata: Option<Value>,
}

// ============================================================================
// TOOL 9: SAMPLE
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleInput {
    /// Sampling method
    pub method: SamplingMethod,

    /// Input data (for statistical analysis)
    #[serde(default)]
    pub data: Vec<f64>,

    /// Number of samples to generate
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub num_samples: Option<usize>,

    /// Additional parameters
    #[serde(default)]
    pub parameters: HashMap<String, Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleOutput {
    /// Generated samples or computed statistics
    pub result: Value,

    /// Statistical moments
    #[serde(skip_serializing_if = "Option::is_none")]
    pub moments: Option<HashMap<String, f64>>,

    /// Confidence intervals
    #[serde(skip_serializing_if = "Option::is_none")]
    pub confidence_intervals: Option<HashMap<String, [f64; 2]>>,

    /// Additional metadata
    #[serde(skip_serializing_if = "Option::is_none")]
    pub metadata: Option<Value>,
}

// ============================================================================
// TOOL 10: OPTIMIZE
// ============================================================================

/// Default optimization method - Auto selection with AICc criteria
fn default_optimization_method() -> OptimizationMethod {
    OptimizationMethod::Auto {
        criteria: SelectionCriteria::Aicc,
        candidates: vec![],
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptimizeInput {
    /// Optimization method (defaults to Auto if not specified)
    #[serde(default = "default_optimization_method")]
    pub method: OptimizationMethod,

    /// Data to fit (x, y pairs)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub data: Option<(Vec<f64>, Vec<f64>)>,

    /// Objective function (as string)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub objective: Option<String>,

    /// Initial guess
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub initial_guess: Option<Vec<f64>>,

    /// Constraints
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub constraints: Option<Vec<String>>,

    /// Additional parameters
    #[serde(default)]
    pub parameters: HashMap<String, Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptimizeOutput {
    /// Optimized parameters
    pub parameters: Vec<f64>,

    /// Fitted function (as string)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub function: Option<String>,

    /// Error/residuals
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<f64>,

    /// R-squared (for fitting)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub r_squared: Option<f64>,

    /// Akaike Information Criterion
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aic: Option<f64>,

    /// Bayesian Information Criterion
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bic: Option<f64>,

    /// Corrected AIC (for small samples)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aicc: Option<f64>,

    /// Convergence info
    #[serde(skip_serializing_if = "Option::is_none")]
    pub convergence: Option<Value>,

    /// Additional metadata
    #[serde(skip_serializing_if = "Option::is_none")]
    pub metadata: Option<Value>,
}
