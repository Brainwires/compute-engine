#!/usr/bin/env node

// MCP test using dev build
const { spawn } = require('child_process');
const path = require('path');
const readline = require('readline');

async function testMCP() {
    const serverPath = path.join(__dirname, 'target/debug/brainwires-compute-engine');
    console.log('Starting server:', serverPath);

    const server = spawn(serverPath, ['mcp-server'], {
        stdio: ['pipe', 'pipe', 'pipe']
    });

    // Read line-by-line from stdout
    const rl = readline.createInterface({
        input: server.stdout,
        crlfDelay: Infinity
    });

    rl.on('line', (line) => {
        console.log('<<< Response:', line);
    });

    server.stderr.on('data', (data) => {
        process.stderr.write(data);
    });

    await new Promise(resolve => setTimeout(resolve, 500));

    // Initialize
    const init = {jsonrpc:'2.0',id:1,method:'initialize',params:{protocolVersion:'2024-11-05',capabilities:{},clientInfo:{name:'test',version:'1.0'}}};
    console.log('>>> Request:', JSON.stringify(init));
    server.stdin.write(JSON.stringify(init) + '\n');

    await new Promise(resolve => setTimeout(resolve, 500));

    // Initialized notification
    const initialized = {jsonrpc:'2.0',method:'notifications/initialized',params:{}};
    console.log('>>> Notification:', JSON.stringify(initialized));
    server.stdin.write(JSON.stringify(initialized) + '\n');

    await new Promise(resolve => setTimeout(resolve, 500));

    // List tools
    const list = {jsonrpc:'2.0',id:2,method:'tools/list',params:{}};
    console.log('>>> Request:', JSON.stringify(list));
    server.stdin.write(JSON.stringify(list) + '\n');

    await new Promise(resolve => setTimeout(resolve, 1000));

    // Call tool
    const call = {
        jsonrpc:'2.0',
        id:3,
        method:'tools/call',
        params:{
            name:'compute_json',
            arguments:{request_json:'{"tool":"solve","input":{"equations":["x^2-4=0"]}}'}
        }
    };
    console.log('>>> Request:', JSON.stringify(call));
    server.stdin.write(JSON.stringify(call) + '\n');

    await new Promise(resolve => setTimeout(resolve, 2000));

    server.kill();
    console.log('\nTest complete');
}

testMCP().catch(console.error);
