// Node.js example for using the computational engine WASM module

const { ComputationalEngine } = require('../../pkg/nodejs/computational_engine.js');

async function main() {
    console.log('🧮 Computational Engine - Node.js Example\n');

    // Create engine instance
    const engine = new ComputationalEngine();
    console.log('Version:', engine.version());
    console.log('');

    // Example 1: Solve an equation
    console.log('1. Solving x^2 - 4 = 0:');
    try {
        const solveResult = engine.solve({
            equations: ['x^2 - 4 = 0']
        });
        console.log(JSON.stringify(solveResult, null, 2));
    } catch (e) {
        console.error('Error:', e.toString());
    }
    console.log('');

    // Example 2: Differentiate
    console.log('2. Differentiating x^3 + 2*x^2 - 5*x + 1:');
    try {
        const diffResult = engine.differentiate({
            expression: 'x^3 + 2*x^2 - 5*x + 1',
            variables: ['x']
        });
        console.log(JSON.stringify(diffResult, null, 2));
    } catch (e) {
        console.error('Error:', e.toString());
    }
    console.log('');

    // Example 3: Integrate
    console.log('3. Integrating x^2:');
    try {
        const intResult = engine.integrate({
            expression: 'x^2',
            variable: 'x'
        });
        console.log(JSON.stringify(intResult, null, 2));
    } catch (e) {
        console.error('Error:', e.toString());
    }
    console.log('');

    // Example 4: Analyze (expand)
    console.log('4. Expanding (x + 1)^2:');
    try {
        const analyzeResult = engine.analyze({
            expression: '(x + 1)^2',
            operation: 'Expand'
        });
        console.log(JSON.stringify(analyzeResult, null, 2));
    } catch (e) {
        console.error('Error:', e.toString());
    }
    console.log('');

    // Example 5: Raw JSON request
    console.log('5. Raw JSON request:');
    const jsonRequest = JSON.stringify({
        tool: 'solve',
        input: {
            equations: ['2*x + 3 = 7']
        }
    });

    try {
        const jsonResult = engine.processJson(jsonRequest);
        console.log(JSON.parse(jsonResult));
    } catch (e) {
        console.error('Error:', e.toString());
    }
    console.log('');

    // Example 6: List operations
    console.log('6. Listing available operations:');
    try {
        const ops = engine.listOperations();
        console.log(JSON.stringify(ops, null, 2));
    } catch (e) {
        console.error('Error:', e.toString());
    }
}

main().catch(console.error);
