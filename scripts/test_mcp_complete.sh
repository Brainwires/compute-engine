#!/bin/bash

# Test complete MCP interaction
(
  echo '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{},"clientInfo":{"name":"test-client","version":"1.0"}}}'
  echo '{"jsonrpc":"2.0","method":"notifications/initialized"}'
  echo '{"jsonrpc":"2.0","id":2,"method":"tools/list"}'
  echo '{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"compute_json","arguments":{"request_json":"{\"tool\":\"solve\",\"input\":{\"equations\":[\"x^2-4=0\"]}}"}}}'
) | ./target/release/brainwires-compute-engine mcp-server 2>&1
