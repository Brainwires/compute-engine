# WebAssembly Integration

This document describes how to build and use the computational engine as a WebAssembly module for browser and Node.js environments.

## Overview

The computational engine provides two interfaces:

1. **MCP Server Interface** (binary: `computational-engine`)
   - CLI tool for tool-calling/MCP integration
   - Runs as a standalone server process
   - JSON-based request/response over stdin/stdout

2. **WASM Library** (feature: `wasm`)
   - Browser and Node.js compatible
   - Direct JavaScript/TypeScript API
   - On-device computation (no server required)

Both interfaces share the same core engine and provide access to all computational capabilities.

## Building for WASM

### Prerequisites

```bash
# Install wasm-pack
cargo install wasm-pack

# Or use the build script (installs wasm-pack if needed)
./build-wasm.sh
```

### Build Commands

```bash
# Build all targets
./build-wasm.sh

# Or build individual targets
npm run build:bundler    # For webpack/rollup/vite
npm run build:nodejs     # For Node.js
npm run build:web        # For ES modules in browser
npm run build:no-modules # For classic <script> tags
```

### Output

Build artifacts are generated in `pkg/`:

```
pkg/
├── bundler/           # For bundlers (webpack, rollup, vite)
├── nodejs/            # For Node.js
├── web/               # For ES modules
└── no-modules/        # For <script> tags
```

## Usage

### Browser (ES Modules)

```html
<!DOCTYPE html>
<html>
<head>
    <script type="module">
        import init, { ComputationalEngine } from './pkg/web/computational_engine.js';

        async function run() {
            await init();
            const engine = new ComputationalEngine();

            // Solve an equation
            const result = engine.solve({
                equations: ["x^2 - 4 = 0"]
            });
            console.log(result);
        }

        run();
    </script>
</head>
<body>
    <h1>Computational Engine Demo</h1>
</body>
</html>
```

### Browser (Classic Script Tag)

```html
<script src="./pkg/no-modules/computational_engine.js"></script>
<script>
    wasm_bindgen('./pkg/no-modules/computational_engine_bg.wasm')
        .then(() => {
            const engine = new wasm_bindgen.ComputationalEngine();
            const result = engine.solve({
                equations: ["x^2 - 4 = 0"]
            });
            console.log(result);
        });
</script>
```

### Node.js

```javascript
const { ComputationalEngine } = require('./pkg/nodejs/computational_engine.js');

const engine = new ComputationalEngine();

// Solve equations
const result = engine.solve({
    equations: ["x^2 - 4 = 0"]
});
console.log(result);
```

### Bundlers (Webpack, Vite, etc.)

```javascript
import init, { ComputationalEngine } from '@brainwires/computational-engine';

await init();
const engine = new ComputationalEngine();

const result = engine.solve({
    equations: ["x^2 - 4 = 0"]
});
```

## API Reference

### ComputationalEngine Class

#### Constructor

```javascript
const engine = new ComputationalEngine();
```

#### Methods

##### `solve(input: SolveInput): SolveOutput`

Solve equations, optimization problems, and systems.

```javascript
const result = engine.solve({
    equations: ["x^2 - 4 = 0"],
    variables: ["x"],        // optional
    initial_guess: null,     // optional
    domain: null,           // optional
    method: null            // optional
});
```

##### `differentiate(input: DifferentiateInput): DifferentiateOutput`

Compute derivatives, gradients, Jacobians, and Hessians.

```javascript
const result = engine.differentiate({
    expression: "x^3 + 2*x^2 - 5*x + 1",
    variables: ["x"],
    order: 1,               // optional
    point: null            // optional
});
```

##### `integrate(input: IntegrateInput): IntegrateOutput`

Compute definite and indefinite integrals.

```javascript
const result = engine.integrate({
    expression: "x^2",
    variable: "x",
    bounds: null,          // optional [lower, upper]
    method: null          // optional
});
```

##### `analyze(input: AnalyzeInput): AnalyzeOutput`

Simplify, expand, factor, and transform expressions.

```javascript
const result = engine.analyze({
    expression: "(x + 1)^2",
    operation: "Expand",   // "Simplify", "Expand", "Factor"
    variables: null       // optional
});
```

##### `simulate(input: SimulateInput): SimulateOutput`

Simulate ODEs, PDEs, physics systems, and optimization.

```javascript
const result = engine.simulate({
    model: "ODE",
    equations: ["dy/dt = -k*y"],
    initial_conditions: { y: 1.0 },
    parameters: { k: 0.1 },
    time_span: [0.0, 10.0],
    options: null
});
```

##### `processJson(jsonRequest: string): string`

Process a raw JSON request (compatible with MCP interface).

```javascript
const jsonRequest = JSON.stringify({
    tool: "solve",
    input: { equations: ["x^2 - 4 = 0"] }
});

const jsonResponse = engine.processJson(jsonRequest);
const result = JSON.parse(jsonResponse);
```

##### `version(): string`

Get the engine version.

```javascript
console.log(engine.version()); // "0.1.0"
```

##### `listOperations(): object`

List all available operations across all modules.

```javascript
const ops = engine.listOperations();
console.log(ops);
// {
//   "Advanced Calculus": [...],
//   "Tensor Calculus": [...],
//   ...
// }
```

## TypeScript Support

The package includes automatically generated TypeScript definitions.

```typescript
import init, { ComputationalEngine } from '@brainwires/computational-engine';

await init();

const engine: ComputationalEngine = new ComputationalEngine();

interface SolveInput {
    equations: string[];
    variables?: string[] | null;
    initial_guess?: number[] | null;
    domain?: [number, number][] | null;
    method?: string | null;
}

const result = engine.solve({
    equations: ["x^2 - 4 = 0"]
});
```

## Examples

See the `examples/wasm/` directory for complete examples:

- `index.html` - Interactive browser demo
- `node-example.js` - Node.js usage
- `README.md` - Detailed instructions

To run the examples:

```bash
# Build WASM
./build-wasm.sh

# Run browser example
cd examples/wasm
npm run dev
# Open http://localhost:8000

# Run Node.js example
npm run node
```

## Performance Considerations

### WASM vs MCP Server

- **WASM**: Runs in-browser/in-process, instant startup, no network overhead
- **MCP Server**: Separate process, better for long-running computations, easier to update

### Bundle Size

The WASM module is approximately:
- Gzipped: ~200-300 KB
- Uncompressed: ~1-2 MB

Use code splitting to load only when needed:

```javascript
// Lazy load the engine
async function loadEngine() {
    const { default: init, ComputationalEngine } = await import(
        '@brainwires/computational-engine'
    );
    await init();
    return new ComputationalEngine();
}

// Use when needed
button.addEventListener('click', async () => {
    const engine = await loadEngine();
    // ...
});
```

### Optimization

For production builds, ensure:

1. `wasm-opt` is installed (part of binaryen)
2. Use release builds: `wasm-pack build --release`
3. Enable LTO in Cargo.toml (already configured)

## Troubleshooting

### "RuntimeError: unreachable executed"

This usually indicates a panic in Rust code. Set up the panic hook:

```javascript
import init, { ComputationalEngine } from '@brainwires/computational-engine';

await init();
// Panic messages will now appear in the browser console
```

### Module not found

Ensure you've built the correct target:

- Browser ES modules: `npm run build:web`
- Node.js: `npm run build:nodejs`
- Bundlers: `npm run build:bundler`

### CORS errors

When using the `web` target, serve files from a local server:

```bash
python3 -m http.server 8000
```

## Publishing

To publish the WASM package to npm:

```bash
# Build all targets
npm run build

# Publish (requires npm login)
npm publish --access public
```

## Comparison: WASM vs MCP Server

| Feature | WASM | MCP Server |
|---------|------|------------|
| Environment | Browser, Node.js | CLI, Server |
| Startup | Instant (~10ms) | Process spawn (~100ms) |
| Overhead | None | IPC/stdio |
| Updates | Rebuild required | Hot reload |
| Integration | Direct API calls | JSON over stdio |
| Use Cases | Web apps, interactive tools | LLM integration, CLI tools |

Both use the same core engine, so choose based on your deployment needs.

## License

MIT OR Apache-2.0
