#!/usr/bin/env node

// Simple MCP client to test the brainwires-compute-engine MCP server (release binary)
const { spawn } = require('child_process');
const path = require('path');

async function testMCPServer() {
    console.log('Starting MCP server (release)...');

    const serverPath = path.join(__dirname, 'target/release/brainwires-compute-engine');
    const server = spawn(serverPath, ['mcp-server'], {
        stdio: ['pipe', 'pipe', 'pipe']
    });

    let responseData = '';
    let stderrData = '';

    server.stdout.on('data', (data) => {
        responseData += data.toString();
        console.log('=== Server Response ===');
        console.log(data.toString());
    });

    server.stderr.on('data', (data) => {
        stderrData += data.toString();
    });

    server.on('error', (err) => {
        console.error('Failed to start server:', err);
        process.exit(1);
    });

    // Wait for server to start
    await new Promise(resolve => setTimeout(resolve, 500));

    // Send initialize request
    console.log('\n>>> Sending initialize request...');
    const initRequest = {
        jsonrpc: '2.0',
        id: 1,
        method: 'initialize',
        params: {
            protocolVersion: '2024-11-05',
            capabilities: {},
            clientInfo: { name: 'test-client', version: '1.0' }
        }
    };
    server.stdin.write(JSON.stringify(initRequest) + '\n');

    await new Promise(resolve => setTimeout(resolve, 500));

    // Send initialized notification (required by MCP protocol)
    console.log('\n>>> Sending initialized notification...');
    const initializedNotif = {
        jsonrpc: '2.0',
        method: 'notifications/initialized',
        params: {}
    };
    server.stdin.write(JSON.stringify(initializedNotif) + '\n');

    await new Promise(resolve => setTimeout(resolve, 200));

    // Send tools/list request
    console.log('\n>>> Sending tools/list request...');
    const toolsRequest = {
        jsonrpc: '2.0',
        id: 2,
        method: 'tools/list',
        params: {}
    };
    server.stdin.write(JSON.stringify(toolsRequest) + '\n');

    await new Promise(resolve => setTimeout(resolve, 1000));

    // Send a compute request
    console.log('\n>>> Sending tools/call request (solve x^2-4=0)...');
    const computeRequest = {
        jsonrpc: '2.0',
        id: 3,
        method: 'tools/call',
        params: {
            name: 'compute_json',
            arguments: {
                request_json: JSON.stringify({
                    tool: 'solve',
                    input: {
                        equations: ['x^2 - 4 = 0']
                    }
                })
            }
        }
    };
    server.stdin.write(JSON.stringify(computeRequest) + '\n');

    await new Promise(resolve => setTimeout(resolve, 2000));

    // Clean up
    server.kill();

    console.log('\n=== Server stderr ===');
    console.log(stderrData);

    console.log('\n=== Test Complete ===');
}

testMCPServer().catch(console.error);
