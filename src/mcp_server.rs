//! MCP Server Implementation
//!
//! Official Model Context Protocol server using rmcp SDK v0.8.1
//! Exposes 4 primary tools: solve, compute, analyze, simulate

pub mod server {
    use crate::engine::*;
    use crate::implementations::create_default_dispatcher;
    use rmcp::{
        handler::server::{ServerHandler, tool::ToolRouter, wrapper::Parameters},
        model::{Implementation, ProtocolVersion, ServerCapabilities, ServerInfo},
        service::ServiceExt,
        tool, tool_handler, tool_router,
    };
    use schemars::JsonSchema;
    use serde::{Deserialize, Serialize};
    use serde_json::Value;
    use std::collections::HashMap;
    use std::sync::Arc;

    // =========================================================================
    // SOLVE TOOL - Equations, systems, optimization
    // =========================================================================

    #[derive(Debug, Serialize, Deserialize, JsonSchema)]
    pub struct SolveRequest {
        /// Type of equation to solve
        /// Options: "root_finding", "linear_system", or nested object like {"einstein": "schwarzschild"}, {"fluid": "navier_stokes"}, {"differential": "ode"}, {"optimize": {"fit": "polynomial"}}
        pub equation_type: Value,

        /// Equations to solve (array of strings)
        /// Example: ["x^2 - 4 = 0"] or ["2x + 3y = 7", "x - y = 1"]
        pub equations: Vec<String>,

        /// Variables to solve for (optional)
        #[serde(default, skip_serializing_if = "Option::is_none")]
        pub variables: Option<Vec<String>>,

        /// Initial guess for iterative methods (optional)
        /// Map of variable name to initial value
        #[serde(default, skip_serializing_if = "Option::is_none")]
        pub initial_guess: Option<HashMap<String, f64>>,
    }

    // =========================================================================
    // COMPUTE TOOL - Calculus, transforms, fields, sampling, matrix ops
    // =========================================================================

    #[derive(Debug, Serialize, Deserialize, JsonSchema)]
    pub struct ComputeRequest {
        /// Operation to perform. Common operations:
        /// - Matrix: {"matrix": "determinant"}, {"matrix": "inverse"}, {"matrix_decomp": "svd"}
        /// - Calculus: {"differentiate": "symbolic"}, {"integrate": "definite"}
        /// - Transforms: {"transform": {"fft": "forward"}}, {"transform": {"laplace": "forward"}}
        /// - Fields: {"field": "decoherence_scale"}, {"field": "green_function"}, {"field": "bohm_potential"}
        /// - Sampling: {"sample": "monte_carlo"}, {"sample": {"distribution": "normal"}}
        /// - Physics: {"physics": {"quantum": "wave_function"}}, {"physics": {"relativity": "lorentz_factor"}}
        pub operation: Value,

        /// Input data for the operation
        /// For matrix: {"matrix": [[1,2],[3,4]]}
        /// For calculus: {"expression": "x^2", "variable": "x"}
        /// For fields: {"mass": 3e-26, "temperature": 300}
        pub data: Value,

        /// Additional parameters (optional)
        #[serde(default)]
        pub parameters: HashMap<String, Value>,
    }

    // =========================================================================
    // ANALYZE TOOL - Series, limits, stability, simplification
    // =========================================================================

    #[derive(Debug, Serialize, Deserialize, JsonSchema)]
    pub struct AnalyzeRequest {
        /// Analysis operation to perform
        /// Options: "simplify", "expand", "factor", "series", "limit", "asymptotic", "stability", "roots", "extrema", "validate"
        pub operation: String,

        /// Expression to analyze
        pub expression: String,

        /// Additional options (optional)
        /// For series: {"variable": "x", "point": 0, "order": 5}
        /// For limit: {"variable": "x", "point": "infinity"}
        #[serde(default)]
        pub options: HashMap<String, Value>,
    }

    // =========================================================================
    // SIMULATE TOOL - Time evolution, stochastic, fluid dynamics
    // =========================================================================

    #[derive(Debug, Serialize, Deserialize, JsonSchema)]
    pub struct SimulateRequest {
        /// Simulation model type
        /// Options:
        /// - Time evolution: {"time_evolution": "euler"}, {"time_evolution": "runge_kutta4"}
        /// - Stochastic: {"stochastic": "brownian_motion"}, {"stochastic": "geometric_brownian"}
        /// - Finance: {"finance": "black_scholes"}, {"finance": "heston"}
        /// - Fluid dynamics: {"fluid_dynamics": "navier_stokes_2d"}, {"fluid_dynamics": "quantum_navier_stokes_1d"}
        pub model: Value,

        /// System equations (array of strings)
        /// Example: ["dx/dt = -k*x", "dy/dt = k*x - m*y"]
        pub equations: Vec<String>,

        /// Variable names
        pub variables: Vec<String>,

        /// System parameters (name -> value)
        #[serde(default)]
        pub parameters: HashMap<String, f64>,

        /// Initial conditions (variable -> value)
        #[serde(default, skip_serializing_if = "Option::is_none")]
        pub initial_conditions: Option<HashMap<String, f64>>,

        /// Time range [start, end]
        #[serde(default, skip_serializing_if = "Option::is_none")]
        pub range: Option<[f64; 2]>,

        /// Number of time steps
        #[serde(default, skip_serializing_if = "Option::is_none")]
        pub steps: Option<usize>,
    }

    // =========================================================================
    // MCP SERVER
    // =========================================================================

    /// Computational Engine MCP Server - 4 Tool Architecture
    #[derive(Clone)]
    pub struct ComputationalEngine {
        dispatcher: Arc<ToolDispatcher>,
        tool_router: ToolRouter<Self>,
    }

    #[tool_router(router = tool_router)]
    impl ComputationalEngine {
        pub fn new() -> Self {
            Self {
                dispatcher: Arc::new(create_default_dispatcher()),
                tool_router: Self::tool_router(),
            }
        }

        /// Solve equations, systems, and optimization problems
        #[tool(
            description = r#"Solve equations and mathematical problems.

EQUATION TYPES:
- "root_finding" - Find roots of equations: {"equations": ["x^2 - 4 = 0"]}
- "linear_system" - Solve Ax=b systems: {"equations": ["2x + y = 5", "x - y = 1"]}
- {"differential": "ode"} - Solve ODEs
- {"einstein": "schwarzschild"} - Einstein field equations
- {"fluid": "navier_stokes"} - Fluid dynamics
- {"optimize": {"fit": "polynomial"}} - Curve fitting (pass data in initial_guess)

EXAMPLES:
1. Quadratic: equation_type="root_finding", equations=["x^2 - 4 = 0"]
2. System: equation_type="linear_system", equations=["2x + y = 5", "x - y = 1"], variables=["x", "y"]"#
        )]
        async fn solve(&self, Parameters(req): Parameters<SolveRequest>) -> Result<String, String> {
            // Convert to internal format
            let equation_type: EquationType = serde_json::from_value(req.equation_type.clone())
                .map_err(|e| format!("Invalid equation_type: {}. Use 'root_finding', 'linear_system', or object like {{\"optimize\": ...}}", e))?;

            let input = SolveInput {
                equation_type,
                equations: req.equations,
                variables: req.variables,
                initial_guess: req.initial_guess,
                boundary_conditions: None,
                domain: None,
                method: None,
                parameters: HashMap::new(),
            };

            let request = ToolRequest::Solve(input);
            match self.dispatcher.dispatch(request) {
                Ok(response) => serde_json::to_string_pretty(&response)
                    .map_err(|e| format!("Serialization error: {}", e)),
                Err(e) => Err(e),
            }
        }

        /// Compute mathematical operations (calculus, transforms, fields, matrices)
        #[tool(
            description = r#"Compute mathematical operations including calculus, transforms, field theory, and matrix operations.

OPERATIONS:
- Matrix: {"matrix": "determinant"}, {"matrix": "inverse"}, {"matrix_decomp": "svd"}
- Calculus: {"differentiate": "symbolic"}, {"integrate": "definite"}
- Transforms: {"transform": {"fft": "forward"}}, {"transform": {"laplace": "forward"}}
- Fields: {"field": "decoherence_scale"}, {"field": "bohm_potential"}, {"field": "green_function"}
- Sampling: {"sample": "monte_carlo"}, {"sample": {"distribution": "normal"}}
- Physics: {"physics": {"quantum": "energy_levels"}}, {"chemistry": "molar_mass"}

DATA FORMAT (depends on operation):
- Matrix ops: {"matrix": [[1,2],[3,4]]}
- Calculus: {"expression": "x^2", "variable": "x", "lower": 0, "upper": 1}
- Fields: {"mass": 3e-26, "temperature": 300}

EXAMPLES:
1. Determinant: operation={"matrix": "determinant"}, data={"matrix": [[1,2],[3,4]]}
2. Decoherence: operation={"field": "decoherence_scale"}, data={"mass": 3e-26, "temperature": 300}
3. Derivative: operation={"differentiate": "symbolic"}, data={"expression": "x^3", "variable": "x"}"#
        )]
        async fn compute(
            &self,
            Parameters(req): Parameters<ComputeRequest>,
        ) -> Result<String, String> {
            // Convert to internal format
            let operation: ComputeOp = serde_json::from_value(req.operation.clone())
                .map_err(|e| format!("Invalid operation: {}. See tool description for valid operations.", e))?;

            let input = ComputeInput {
                operation,
                data: req.data,
                parameters: req.parameters,
            };

            let request = ToolRequest::Compute(input);
            match self.dispatcher.dispatch(request) {
                Ok(response) => serde_json::to_string_pretty(&response)
                    .map_err(|e| format!("Serialization error: {}", e)),
                Err(e) => Err(e),
            }
        }

        /// Analyze mathematical expressions (simplify, series, limits, stability)
        #[tool(
            description = r#"Analyze mathematical expressions - simplify, expand, series expansion, limits, stability analysis.

OPERATIONS:
- "simplify" - Simplify algebraic expressions
- "expand" - Expand products and powers
- "factor" - Factor polynomials
- "series" - Taylor/Laurent series (use options: {"variable": "x", "point": 0, "order": 5})
- "limit" - Compute limits (use options: {"variable": "x", "point": "infinity"})
- "stability" - Analyze system stability
- "roots" - Find roots of polynomials
- "extrema" - Find critical points

EXAMPLES:
1. Simplify: operation="simplify", expression="(x+1)^2 - x^2 - 2*x - 1"
2. Series: operation="series", expression="sin(x)", options={"variable": "x", "point": 0, "order": 5}
3. Limit: operation="limit", expression="sin(x)/x", options={"variable": "x", "point": 0}"#
        )]
        async fn analyze(
            &self,
            Parameters(req): Parameters<AnalyzeRequest>,
        ) -> Result<String, String> {
            // Convert operation string to AnalysisOp
            let operation: AnalysisOp = serde_json::from_value(serde_json::json!(req.operation))
                .map_err(|e| format!("Invalid operation '{}': {}. Use: simplify, expand, factor, series, limit, stability, roots, extrema", req.operation, e))?;

            let input = AnalyzeInput {
                operation,
                expression: req.expression,
                options: req.options,
            };

            let request = ToolRequest::Analyze(input);
            match self.dispatcher.dispatch(request) {
                Ok(response) => serde_json::to_string_pretty(&response)
                    .map_err(|e| format!("Serialization error: {}", e)),
                Err(e) => Err(e),
            }
        }

        /// Simulate dynamic systems (ODEs, stochastic processes, fluid dynamics)
        #[tool(
            description = r#"Simulate time evolution of dynamic systems including ODEs, stochastic processes, and fluid dynamics.

MODELS:
- Time evolution: {"time_evolution": "euler"}, {"time_evolution": "runge_kutta4"}, {"time_evolution": "adaptive_step"}
- Stochastic: {"stochastic": "brownian_motion"}, {"stochastic": "geometric_brownian"}, {"stochastic": "ornstein_uhlenbeck"}
- Finance: {"finance": "black_scholes"}, {"finance": "heston"}
- Fluids: {"fluid_dynamics": "navier_stokes_2d"}, {"fluid_dynamics": "quantum_navier_stokes_1d"}

EXAMPLES:
1. Harmonic oscillator:
   model={"time_evolution": "runge_kutta4"}
   equations=["dx/dt = v", "dv/dt = -omega^2 * x"]
   variables=["x", "v"]
   parameters={"omega": 1.0}
   initial_conditions={"x": 1.0, "v": 0.0}
   range=[0.0, 10.0]
   steps=100

2. Brownian motion:
   model={"stochastic": "brownian_motion"}
   equations=["dX = mu*dt + sigma*dW"]
   variables=["X"]
   parameters={"mu": 0.0, "sigma": 1.0}
   initial_conditions={"X": 0.0}
   range=[0.0, 1.0]
   steps=1000"#
        )]
        async fn simulate(
            &self,
            Parameters(req): Parameters<SimulateRequest>,
        ) -> Result<String, String> {
            // Convert model to SimulationModel
            let model: SimulationModel = serde_json::from_value(req.model.clone())
                .map_err(|e| format!("Invalid model: {}. Use {{\"time_evolution\": \"euler\"}}, {{\"stochastic\": \"brownian_motion\"}}, etc.", e))?;

            let input = SimulateInput {
                model,
                equations: req.equations,
                variables: req.variables,
                parameters: req.parameters,
                initial_conditions: req.initial_conditions,
                range: req.range,
                steps: req.steps,
                method: None,
                num_paths: None,
            };

            let request = ToolRequest::Simulate(input);
            match self.dispatcher.dispatch(request) {
                Ok(response) => serde_json::to_string_pretty(&response)
                    .map_err(|e| format!("Serialization error: {}", e)),
                Err(e) => Err(e),
            }
        }
    }

    #[tool_handler(router = self.tool_router)]
    impl ServerHandler for ComputationalEngine {
        fn get_info(&self) -> ServerInfo {
            ServerInfo {
                protocol_version: ProtocolVersion::default(),
                capabilities: ServerCapabilities::builder().enable_tools().build(),
                server_info: Implementation {
                    name: "brainwires-compute-engine".into(),
                    title: Some("Brainwires Compute Engine".into()),
                    version: env!("CARGO_PKG_VERSION").into(),
                    icons: None,
                    website_url: None,
                },
                instructions: Some(
                    "Computational engine with 229+ operations across 4 tools:\n\
                    - solve: Equations, systems, optimization\n\
                    - compute: Matrix ops, calculus, transforms, field theory, sampling\n\
                    - analyze: Simplify, series, limits, stability\n\
                    - simulate: ODEs, stochastic processes, fluid dynamics"
                        .into(),
                ),
            }
        }
    }

    /// Start the MCP server using stdio transport
    pub async fn serve_stdio() -> Result<(), Box<dyn std::error::Error>> {
        let engine = ComputationalEngine::new();
        let transport = rmcp::transport::io::stdio();

        eprintln!(
            "Starting Computational Engine MCP Server v{}",
            env!("CARGO_PKG_VERSION")
        );
        eprintln!("Protocol: MCP 2024-11-05");
        eprintln!("Tools: solve, compute, analyze, simulate");
        eprintln!("Operations: 229+ mathematical and physics computations");
        eprintln!("");

        engine.serve(transport).await?.waiting().await?;

        Ok(())
    }
}
