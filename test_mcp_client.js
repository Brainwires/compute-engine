#!/usr/bin/env node

// Simple MCP client to test the computational-engine MCP server
const { spawn } = require('child_process');

async function testMCPServer() {
    console.log('Starting MCP server...');

    const server = spawn('cargo', ['run', '--', 'mcp-server'], {
        cwd: __dirname,
        stdio: ['pipe', 'pipe', 'pipe']
    });

    let responseData = '';

    server.stdout.on('data', (data) => {
        responseData += data.toString();
        console.log('Server stdout:', data.toString());
    });

    server.stderr.on('data', (data) => {
        console.error('Server stderr:', data.toString());
    });

    // Wait for server to start
    await new Promise(resolve => setTimeout(resolve, 2000));

    // Send initialize request
    console.log('\nSending initialize request...');
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

    await new Promise(resolve => setTimeout(resolve, 1000));

    // Send tools/list request
    console.log('\nSending tools/list request...');
    const toolsRequest = {
        jsonrpc: '2.0',
        id: 2,
        method: 'tools/list',
        params: {}
    };
    server.stdin.write(JSON.stringify(toolsRequest) + '\n');

    await new Promise(resolve => setTimeout(resolve, 1000));

    // Send a compute request
    console.log('\nSending tools/call request (solve x^2-4=0)...');
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

    console.log('\n=== Test Complete ===');
}

testMCPServer().catch(console.error);
