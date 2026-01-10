//! MCP Server Implementation
//!
//! Official Model Context Protocol server using rmcp SDK v0.8.1

pub mod server {
    use crate::engine::*;
    use crate::implementations::create_default_dispatcher;
    use rmcp::{
        handler::server::{ServerHandler, tool::ToolRouter, wrapper::Parameters},
        model::{Implementation, ProtocolVersion, ServerCapabilities, ServerInfo},
        service::ServiceExt,
        tool, tool_handler, tool_router,
    };
    use schemars::{schema_for, JsonSchema};
    use serde::{Deserialize, Serialize};
    use std::sync::Arc;

    /// Request parameters for compute_json tool
    #[derive(Debug, Serialize, Deserialize, JsonSchema)]
    pub struct ComputeJsonRequest {
        /// JSON string containing ToolRequest data
        pub request_json: String,
    }

    /// Generate the full JSON schema for ToolRequest
    pub fn get_tool_request_schema() -> String {
        let schema = schema_for!(ToolRequest);
        serde_json::to_string_pretty(&schema).unwrap_or_else(|_| "{}".to_string())
    }

    /// Get examples for each tool type
    pub fn get_tool_examples() -> String {
        r#"{
  "examples": {
    "solve": {
      "description": "Solve equations (algebraic, differential, linear systems)",
      "example": {"tool":"solve","input":{"equation_type":"root_finding","equations":["x^2 - 4 = 0"]}}
    },
    "differentiate": {
      "description": "Compute derivatives (symbolic, numeric, partial, gradient)",
      "example": {"tool":"differentiate","input":{"operation":"symbolic","expression":"x^3 + 2*x","variables":["x"]}}
    },
    "integrate": {
      "description": "Compute integrals (definite, indefinite, multiple, contour)",
      "example": {"tool":"integrate","input":{"integration_type":"symbolic","expression":"x^2","variables":["x"]}}
    },
    "analyze": {
      "description": "Analyze expressions (simplify, expand, series, limits)",
      "example": {"tool":"analyze","input":{"operation":"simplify","expression":"(x+1)^2 - x^2 - 2*x"}}
    },
    "simulate": {
      "description": "Run simulations (ODE, stochastic, Monte Carlo)",
      "example": {"tool":"simulate","input":{"model":{"time_evolution":"runge_kutta4"},"equations":["dx/dt = -x"],"variables":["x"],"initial_conditions":{"x":1.0},"range":[0.0,10.0],"steps":100}}
    },
    "compute": {
      "description": "Compute operations (matrix, tensor, special functions)",
      "example": {"tool":"compute","input":{"operation":{"matrix":"determinant"},"data":[[1,2],[3,4]]}}
    },
    "transform": {
      "description": "Apply transforms (FFT, Fourier, Laplace, wavelets)",
      "example": {"tool":"transform","input":{"transform_type":{"fft":"forward"},"data":[1.0,0.0,-1.0,0.0]}}
    },
    "fieldtheory": {
      "description": "Field theory calculations (EM fields, Green's functions)",
      "example": {"tool":"fieldtheory","input":{"field_type":"green_function","configuration":{"source":[0,0,0],"field_point":[1,0,0]}}}
    },
    "sample": {
      "description": "Statistical sampling (Monte Carlo, MCMC, distributions)",
      "example": {"tool":"sample","input":{"method":"path_generation","num_samples":1000,"parameters":{"process":"brownian_motion"}}}
    },
    "optimize": {
      "description": "Optimization (curve fitting, minimization, regression)",
      "example": {"tool":"optimize","input":{"method":{"fit":"polynomial"},"data":[[0,1,2,3,4],[0,1,4,9,16]],"parameters":{"degree":2}}}
    }
  },
  "important_notes": [
    "Tool names MUST be lowercase in JSON (solve, analyze, NOT Solve, Analyze)",
    "Each tool has a required operation/type field - check the schema",
    "Use get_schema tool to see the complete JSON schema for all types"
  ]
}"#.to_string()
    }

    /// Computational Engine MCP Server
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

        /// Execute a computational request using JSON
        #[tool(
            description = r#"Execute computational operations via JSON.

IMPORTANT: Tool names MUST be lowercase: solve, differentiate, integrate, analyze, simulate, compute, transform, fieldtheory, sample, optimize

Quick examples:
- solve: {"tool":"solve","input":{"equation_type":"root_finding","equations":["x^2-4=0"]}}
- analyze: {"tool":"analyze","input":{"operation":"simplify","expression":"(x+1)^2"}}
- differentiate: {"tool":"differentiate","input":{"operation":"symbolic","expression":"x^2","variables":["x"]}}
- compute: {"tool":"compute","input":{"operation":{"matrix":"determinant"},"data":[[1,2],[3,4]]}}

Use get_schema tool for complete JSON schema, or get_examples for all tool examples."#
        )]
        async fn compute_json(
            &self,
            Parameters(req): Parameters<ComputeJsonRequest>,
        ) -> Result<String, String> {
            // Parse the JSON request
            let tool_request: ToolRequest = serde_json::from_str(&req.request_json)
                .map_err(|e| format!("JSON parse error: {}", e))?;

            // Dispatch the request
            match self.dispatcher.dispatch(tool_request) {
                Ok(response) => {
                    let result_json = serde_json::to_string_pretty(&response)
                        .map_err(|e| format!("Serialization failed: {}", e))?;

                    Ok(result_json)
                }
                Err(e) => Err(format!("Computation error: {}", e)),
            }
        }

        /// Get the complete JSON schema for ToolRequest
        #[tool(
            description = "Get the complete JSON schema for all computational tools. Returns the full type definitions for ToolRequest including all 10 tools (solve, differentiate, integrate, analyze, simulate, compute, transform, fieldtheory, sample, optimize) and their input types."
        )]
        async fn get_schema(&self) -> Result<String, String> {
            Ok(get_tool_request_schema())
        }

        /// Get examples for all tools
        #[tool(
            description = "Get working examples for all 10 computational tools. Each example shows the correct JSON format including required fields. Use this to understand how to construct requests for each tool type."
        )]
        async fn get_examples(&self) -> Result<String, String> {
            Ok(get_tool_examples())
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
                    title: Some("Brainwires Compute Engine - Mathematics & Physics".into()),
                    version: env!("CARGO_PKG_VERSION").into(),
                    icons: None,
                    website_url: None,
                },
                instructions: Some(
                    "Computational engine with 229+ operations. Use compute_json with lowercase tool names \
                    (solve, analyze, etc). Call get_schema for full JSON schema or get_examples for working examples. \
                    IMPORTANT: Tool names in JSON must be lowercase (e.g., 'analyze' not 'Analyze')."
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
        eprintln!("Tool: compute_json");
        eprintln!(
            "Operations: Solve, Differentiate, Integrate, Analyze, Simulate, Compute, Transform, FieldTheory, Sample, Optimize"
        );
        eprintln!("Total: 229+ mathematical and physics operations");
        eprintln!("");

        engine.serve(transport).await?.waiting().await?;

        Ok(())
    }
}
