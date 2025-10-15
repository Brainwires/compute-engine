//! MCP Server Implementation
//!
//! Official Model Context Protocol server using rmcp SDK v0.8.1

pub mod server {
    use crate::engine::*;
    use crate::implementations::create_default_dispatcher;
    use rmcp::{
        handler::server::{ServerHandler, tool::ToolRouter, wrapper::Parameters},
        model::{ServerInfo, Implementation, ProtocolVersion, ServerCapabilities},
        tool, tool_router, tool_handler,
        service::ServiceExt,
    };
    use schemars::JsonSchema;
    use serde::{Deserialize, Serialize};
    use std::sync::Arc;

    /// Request parameters for compute_json tool
    #[derive(Debug, Serialize, Deserialize, JsonSchema)]
    pub struct ComputeJsonRequest {
        /// JSON string containing ToolRequest data
        pub request_json: String,
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
        #[tool(description = "Execute computational operations via JSON. Accepts ToolRequest JSON format with 10 tools: Solve, Differentiate, Integrate, Analyze, Simulate, Compute, Transform, FieldTheory, Sample, Optimize")]
        async fn compute_json(&self, Parameters(req): Parameters<ComputeJsonRequest>) -> Result<String, String> {
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
                Err(e) => {
                    Err(format!("Computation error: {}", e))
                }
            }
        }
    }

    #[tool_handler(router = self.tool_router)]
    impl ServerHandler for ComputationalEngine {
        fn get_info(&self) -> ServerInfo {
            ServerInfo {
                protocol_version: ProtocolVersion::default(),
                capabilities: ServerCapabilities::builder()
                    .enable_tools()
                    .build(),
                server_info: Implementation {
                    name: "computational-engine".into(),
                    title: Some("Computational Engine - Mathematics & Physics".into()),
                    version: env!("CARGO_PKG_VERSION").into(),
                    icons: None,
                    website_url: None,
                },
                instructions: Some(
                    "Use compute_json tool with ToolRequest JSON. Example: \
                    {\"tool\":\"solve\",\"input\":{\"equations\":[\"x^2-4=0\"]}}. \
                    Supports 229+ operations across 10 tools.".into()
                ),
            }
        }
    }

    /// Start the MCP server using stdio transport
    pub async fn serve_stdio() -> Result<(), Box<dyn std::error::Error>> {
        let engine = ComputationalEngine::new();
        let transport = rmcp::transport::io::stdio();

        eprintln!("Starting Computational Engine MCP Server v{}", env!("CARGO_PKG_VERSION"));
        eprintln!("Protocol: MCP 2024-11-05");
        eprintln!("Tool: compute_json");
        eprintln!("Operations: Solve, Differentiate, Integrate, Analyze, Simulate, Compute, Transform, FieldTheory, Sample, Optimize");
        eprintln!("Total: 229+ mathematical and physics operations");
        eprintln!("");

        engine.serve(transport).await?.waiting().await?;

        Ok(())
    }
}
