//! MCP Server Integration Tests
//!
//! Tests for the 4-tool MCP architecture: solve, compute, analyze, simulate

use computational_engine::mcp_server::server::ComputationalEngine;
use rmcp::{ClientHandler, model::CallToolRequestParam, model::ClientInfo, service::ServiceExt};

#[derive(Debug, Clone, Default)]
struct DummyClientHandler {}

impl ClientHandler for DummyClientHandler {
    fn get_info(&self) -> ClientInfo {
        ClientInfo::default()
    }
}

#[tokio::test]
async fn test_mcp_tool_list() -> anyhow::Result<()> {
    let (server_transport, client_transport) = tokio::io::duplex(4096);

    // Server setup
    let server = ComputationalEngine::new();
    let server_handle = tokio::spawn(async move {
        server.serve(server_transport).await?.waiting().await?;
        anyhow::Ok(())
    });

    // Client setup
    let client_handler = DummyClientHandler::default();
    let client = client_handler.serve(client_transport).await?;

    // List tools
    let tools = client.list_tools(None).await?;
    println!("Tools: {:?}", tools);

    // Should have exactly 8 tools: solve, compute, analyze, simulate, ml, chaos, units, validate
    assert_eq!(tools.tools.len(), 8, "Should have exactly 8 tools");

    let tool_names: Vec<String> = tools.tools.iter().map(|t| t.name.to_string()).collect();
    assert!(tool_names.iter().any(|n| n == "solve"), "Should have solve tool");
    assert!(tool_names.iter().any(|n| n == "compute"), "Should have compute tool");
    assert!(tool_names.iter().any(|n| n == "analyze"), "Should have analyze tool");
    assert!(tool_names.iter().any(|n| n == "simulate"), "Should have simulate tool");
    assert!(tool_names.iter().any(|n| n == "ml"), "Should have ml tool");
    assert!(tool_names.iter().any(|n| n == "chaos"), "Should have chaos tool");
    assert!(tool_names.iter().any(|n| n == "units"), "Should have units tool");
    assert!(tool_names.iter().any(|n| n == "validate"), "Should have validate tool");

    client.cancel().await?;
    server_handle.await??;
    Ok(())
}

#[tokio::test]
async fn test_mcp_solve_tool() -> anyhow::Result<()> {
    let (server_transport, client_transport) = tokio::io::duplex(4096);

    // Server setup
    let server = ComputationalEngine::new();
    let server_handle = tokio::spawn(async move {
        server.serve(server_transport).await?.waiting().await?;
        anyhow::Ok(())
    });

    // Client setup
    let client_handler = DummyClientHandler::default();
    let client = client_handler.serve(client_transport).await?;

    // Call solve tool - find roots of x^2 - 4 = 0
    let result = client
        .call_tool(CallToolRequestParam {
            name: "solve".into(),
            arguments: Some(
                serde_json::json!({
                    "equation_type": "root_finding",
                    "equations": ["x^2 - 4 = 0"]
                })
                .as_object()
                .unwrap()
                .clone(),
            ),
        })
        .await?;

    println!("Solve Result: {:?}", result);

    let result_text = result
        .content
        .first()
        .and_then(|content| content.raw.as_text())
        .map(|text| text.text.as_str())
        .expect("Expected text content");

    println!("Result text: {}", result_text);

    // Parse the result JSON
    let response: serde_json::Value = serde_json::from_str(result_text)?;
    assert!(response.is_object(), "Should return a JSON object");

    client.cancel().await?;
    server_handle.await??;
    Ok(())
}

#[tokio::test]
async fn test_mcp_compute_tool() -> anyhow::Result<()> {
    let (server_transport, client_transport) = tokio::io::duplex(4096);

    // Server setup
    let server = ComputationalEngine::new();
    let server_handle = tokio::spawn(async move {
        server.serve(server_transport).await?.waiting().await?;
        anyhow::Ok(())
    });

    // Client setup
    let client_handler = DummyClientHandler::default();
    let client = client_handler.serve(client_transport).await?;

    // Call compute tool - matrix determinant
    let result = client
        .call_tool(CallToolRequestParam {
            name: "compute".into(),
            arguments: Some(
                serde_json::json!({
                    "operation": {"matrix": "determinant"},
                    "data": {"matrix": [[1.0, 2.0], [3.0, 4.0]]}
                })
                .as_object()
                .unwrap()
                .clone(),
            ),
        })
        .await?;

    println!("Compute Result: {:?}", result);

    let result_text = result
        .content
        .first()
        .and_then(|content| content.raw.as_text())
        .map(|text| text.text.as_str())
        .expect("Expected text content");

    println!("Result text: {}", result_text);

    // Parse the result JSON
    let response: serde_json::Value = serde_json::from_str(result_text)?;
    assert!(response.is_object(), "Should return a JSON object");

    client.cancel().await?;
    server_handle.await??;
    Ok(())
}

#[tokio::test]
async fn test_mcp_analyze_tool() -> anyhow::Result<()> {
    let (server_transport, client_transport) = tokio::io::duplex(4096);

    // Server setup
    let server = ComputationalEngine::new();
    let server_handle = tokio::spawn(async move {
        server.serve(server_transport).await?.waiting().await?;
        anyhow::Ok(())
    });

    // Client setup
    let client_handler = DummyClientHandler::default();
    let client = client_handler.serve(client_transport).await?;

    // Call analyze tool - simplify expression
    let result = client
        .call_tool(CallToolRequestParam {
            name: "analyze".into(),
            arguments: Some(
                serde_json::json!({
                    "operation": "simplify",
                    "expression": "x + x"
                })
                .as_object()
                .unwrap()
                .clone(),
            ),
        })
        .await?;

    println!("Analyze Result: {:?}", result);

    let result_text = result
        .content
        .first()
        .and_then(|content| content.raw.as_text())
        .map(|text| text.text.as_str())
        .expect("Expected text content");

    println!("Result text: {}", result_text);

    // Parse the result JSON
    let response: serde_json::Value = serde_json::from_str(result_text)?;
    assert!(response.is_object(), "Should return a JSON object");

    client.cancel().await?;
    server_handle.await??;
    Ok(())
}

#[tokio::test]
async fn test_mcp_simulate_tool() -> anyhow::Result<()> {
    let (server_transport, client_transport) = tokio::io::duplex(4096);

    // Server setup
    let server = ComputationalEngine::new();
    let server_handle = tokio::spawn(async move {
        server.serve(server_transport).await?.waiting().await?;
        anyhow::Ok(())
    });

    // Client setup
    let client_handler = DummyClientHandler::default();
    let client = client_handler.serve(client_transport).await?;

    // Call simulate tool - simple ODE (exponential decay)
    let result = client
        .call_tool(CallToolRequestParam {
            name: "simulate".into(),
            arguments: Some(
                serde_json::json!({
                    "model": {"time_evolution": "euler"},
                    "equations": ["dx/dt = -x"],
                    "variables": ["x"],
                    "parameters": {},
                    "initial_conditions": {"x": 1.0},
                    "range": [0.0, 1.0],
                    "steps": 10
                })
                .as_object()
                .unwrap()
                .clone(),
            ),
        })
        .await?;

    println!("Simulate Result: {:?}", result);

    let result_text = result
        .content
        .first()
        .and_then(|content| content.raw.as_text())
        .map(|text| text.text.as_str())
        .expect("Expected text content");

    println!("Result text: {}", result_text);

    // Parse the result JSON
    let response: serde_json::Value = serde_json::from_str(result_text)?;
    assert!(response.is_object(), "Should return a JSON object");

    client.cancel().await?;
    server_handle.await??;
    Ok(())
}

#[tokio::test]
async fn test_mcp_ml_tool() -> anyhow::Result<()> {
    let (server_transport, client_transport) = tokio::io::duplex(4096);

    // Server setup
    let server = ComputationalEngine::new();
    let server_handle = tokio::spawn(async move {
        server.serve(server_transport).await?.waiting().await?;
        anyhow::Ok(())
    });

    // Client setup
    let client_handler = DummyClientHandler::default();
    let client = client_handler.serve(client_transport).await?;

    // Call ml tool - k-means clustering
    let result = client
        .call_tool(CallToolRequestParam {
            name: "ml".into(),
            arguments: Some(
                serde_json::json!({
                    "operation": {"clustering": "k_means"},
                    "data": [[1.0, 2.0], [1.5, 1.8], [5.0, 8.0], [8.0, 8.0]],
                    "parameters": {"k": 2}
                })
                .as_object()
                .unwrap()
                .clone(),
            ),
        })
        .await?;

    println!("ML Result: {:?}", result);

    let result_text = result
        .content
        .first()
        .and_then(|content| content.raw.as_text())
        .map(|text| text.text.as_str())
        .expect("Expected text content");

    println!("Result text: {}", result_text);

    // Parse the result JSON
    let response: serde_json::Value = serde_json::from_str(result_text)?;
    assert!(response.is_object(), "Should return a JSON object");

    client.cancel().await?;
    server_handle.await??;
    Ok(())
}

#[tokio::test]
async fn test_mcp_chaos_tool() -> anyhow::Result<()> {
    let (server_transport, client_transport) = tokio::io::duplex(4096);

    // Server setup
    let server = ComputationalEngine::new();
    let server_handle = tokio::spawn(async move {
        server.serve(server_transport).await?.waiting().await?;
        anyhow::Ok(())
    });

    // Client setup
    let client_handler = DummyClientHandler::default();
    let client = client_handler.serve(client_transport).await?;

    // Call chaos tool - Lorenz attractor
    let result = client
        .call_tool(CallToolRequestParam {
            name: "chaos".into(),
            arguments: Some(
                serde_json::json!({
                    "operation": {"attractor": "lorenz"},
                    "parameters": {"sigma": 10.0, "rho": 28.0, "beta": 2.666}
                })
                .as_object()
                .unwrap()
                .clone(),
            ),
        })
        .await?;

    println!("Chaos Result: {:?}", result);

    let result_text = result
        .content
        .first()
        .and_then(|content| content.raw.as_text())
        .map(|text| text.text.as_str())
        .expect("Expected text content");

    println!("Result text: {}", result_text);

    // Parse the result JSON
    let response: serde_json::Value = serde_json::from_str(result_text)?;
    assert!(response.is_object(), "Should return a JSON object");

    client.cancel().await?;
    server_handle.await??;
    Ok(())
}

#[tokio::test]
async fn test_mcp_units_tool() -> anyhow::Result<()> {
    let (server_transport, client_transport) = tokio::io::duplex(4096);

    // Server setup
    let server = ComputationalEngine::new();
    let server_handle = tokio::spawn(async move {
        server.serve(server_transport).await?.waiting().await?;
        anyhow::Ok(())
    });

    // Client setup
    let client_handler = DummyClientHandler::default();
    let client = client_handler.serve(client_transport).await?;

    // Call units tool - check unit compatibility (N and kg*m/s^2 are both force)
    let result = client
        .call_tool(CallToolRequestParam {
            name: "units".into(),
            arguments: Some(
                serde_json::json!({
                    "operation": "check_compatibility",
                    "from_unit": "N",
                    "to_unit": "kg*m/s^2"
                })
                .as_object()
                .unwrap()
                .clone(),
            ),
        })
        .await?;

    println!("Units Result: {:?}", result);

    let result_text = result
        .content
        .first()
        .and_then(|content| content.raw.as_text())
        .map(|text| text.text.as_str())
        .expect("Expected text content");

    println!("Result text: {}", result_text);

    // Parse the result JSON
    let response: serde_json::Value = serde_json::from_str(result_text)?;
    assert!(response.is_object(), "Should return a JSON object");

    client.cancel().await?;
    server_handle.await??;
    Ok(())
}

#[tokio::test]
async fn test_mcp_validate_tool() -> anyhow::Result<()> {
    let (server_transport, client_transport) = tokio::io::duplex(4096);

    // Server setup
    let server = ComputationalEngine::new();
    let server_handle = tokio::spawn(async move {
        server.serve(server_transport).await?.waiting().await?;
        anyhow::Ok(())
    });

    // Client setup
    let client_handler = DummyClientHandler::default();
    let client = client_handler.serve(client_transport).await?;

    // Call validate tool - check equation validity
    let result = client
        .call_tool(CallToolRequestParam {
            name: "validate".into(),
            arguments: Some(
                serde_json::json!({
                    "operation": "equation",
                    "expression": "E = m*c^2",
                    "parameters": {"domain": "relativity"}
                })
                .as_object()
                .unwrap()
                .clone(),
            ),
        })
        .await?;

    println!("Validate Result: {:?}", result);

    let result_text = result
        .content
        .first()
        .and_then(|content| content.raw.as_text())
        .map(|text| text.text.as_str())
        .expect("Expected text content");

    println!("Result text: {}", result_text);

    // Parse the result JSON
    let response: serde_json::Value = serde_json::from_str(result_text)?;
    assert!(response.is_object(), "Should return a JSON object");

    client.cancel().await?;
    server_handle.await??;
    Ok(())
}
