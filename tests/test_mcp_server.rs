use computational_engine::mcp_server::server::ComputationalEngine;
use rmcp::{
    ClientHandler, model::ClientInfo,
    model::CallToolRequestParam,
    service::ServiceExt,
};

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

    assert!(!tools.tools.is_empty(), "Should have at least one tool");

    client.cancel().await?;
    server_handle.await??;
    Ok(())
}

#[tokio::test]
async fn test_mcp_tool_call() -> anyhow::Result<()> {
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

    // Call compute_json tool
    let result = client
        .call_tool(CallToolRequestParam {
            name: "compute_json".into(),
            arguments: Some(
                serde_json::json!({
                    "request_json": "{\"tool\":\"solve\",\"input\":{\"equation_type\":\"root_finding\",\"equations\":[\"x^2-4=0\"]}}"
                })
                .as_object()
                .unwrap()
                .clone(),
            ),
        })
        .await?;

    println!("Result: {:?}", result);

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
