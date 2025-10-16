#!/bin/bash

# Test MCP server with proper JSON-RPC messages

# Start the server in background
cargo run -- mcp-server 2>/dev/null &
SERVER_PID=$!

sleep 1

# Send initialize request
echo '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{},"clientInfo":{"name":"test","version":"1.0"}}}' | cargo run -- mcp-server 2>/dev/null &
sleep 2

# Kill the server
kill $SERVER_PID 2>/dev/null
