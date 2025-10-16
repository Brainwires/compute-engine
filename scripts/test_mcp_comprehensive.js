#!/usr/bin/env node

const { spawn } = require('child_process');
const path = require('path');

async function testMCP() {
    const serverPath = path.join(__dirname, 'target/debug/brainwires-compute-engine');
    console.log('=== Starting MCP Server Test ===\n');

    const server = spawn(serverPath, ['mcp-server'], {
        stdio: ['pipe', 'pipe', 'pipe']
    });

    let allStdout = '';
    let allStderr = '';
    let responseCount = 0;

    server.stdout.on('data', (data) => {
        const text = data.toString();
        allStdout += text;

        // Split by lines to handle multiple JSON responses
        const lines = text.split('\n').filter(l => l.trim());
        lines.forEach(line => {
            if (line.trim()) {
                responseCount++;
                console.log(`\n[RESPONSE ${responseCount}]`);
                try {
                    const json = JSON.parse(line);
                    console.log(JSON.stringify(json, null, 2));
                } catch {
                    console.log(line);
                }
            }
        });
    });

    server.stderr.on('data', (data) => {
        allStderr += data.toString();
        process.stderr.write(`[STDERR] ${data}`);
    });

    server.on('error', (err) => {
        console.error('\n[ERROR] Failed to start server:', err);
        process.exit(1);
    });

    const send = (msg, label) => {
        console.log(`\n[REQUEST ${label}]`);
        console.log(JSON.stringify(msg, null, 2));
        server.stdin.write(JSON.stringify(msg) + '\n');
    };

    await new Promise(r => setTimeout(r, 500));

    // 1. Initialize
    send({
        jsonrpc: '2.0',
        id: 1,
        method: 'initialize',
        params: {
            protocolVersion: '2024-11-05',
            capabilities: {},
            clientInfo: { name: 'test', version: '1.0' }
        }
    }, '1-INITIALIZE');

    await new Promise(r => setTimeout(r, 500));

    // 2. Initialized notification
    send({
        jsonrpc: '2.0',
        method: 'notifications/initialized',
        params: {}
    }, '2-INITIALIZED');

    await new Promise(r => setTimeout(r, 500));

    // 3. List tools
    send({
        jsonrpc: '2.0',
        id: 2,
        method: 'tools/list',
        params: {}
    }, '3-TOOLS/LIST');

    await new Promise(r => setTimeout(r, 1000));

    // 4. Call compute_json tool
    send({
        jsonrpc: '2.0',
        id: 3,
        method: 'tools/call',
        params: {
            name: 'compute_json',
            arguments: {
                request_json: '{"tool":"solve","input":{"equation_type":"root_finding","equations":["x^2-4=0"]}}'
            }
        }
    }, '4-TOOLS/CALL');

    await new Promise(r => setTimeout(r, 2000));

    server.kill();

    console.log('\n\n=== Test Summary ===');
    console.log(`Total responses received: ${responseCount}`);
    console.log(`\nFull stdout length: ${allStdout.length} bytes`);
    console.log(`Full stderr length: ${allStderr.length} bytes`);
}

testMCP().catch(console.error);
