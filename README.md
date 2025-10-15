# Computational Engine

A unified computational engine for mathematical and scientific computing in Rust, built around **5 powerful, flexible tools**.

## 🚀 Quick Start

```rust
use computational_engine::core::{ToolRequest, SolveInput};
use computational_engine::implementations::create_default_dispatcher;

// Create dispatcher
let dispatcher = create_default_dispatcher();

// Solve a quadratic equation
let request = ToolRequest::Solve(SolveInput {
    equations: vec!["x^2 + 2*x - 8 = 0".to_string()],
    variables: None,
    initial_guess: None,
    domain: None,
    method: None,
});

let response = dispatcher.dispatch(request).unwrap();
println!("{:#?}", response);
```

## 🧮 The 5 Core Tools

| Tool | Purpose | Example |
|------|---------|---------|
| **Solve** | Equations, systems, optimization | Solve `x^2 - 4 = 0` |
| **Differentiate** | Derivatives, gradients, Jacobians | `∂/∂x (sin(x)*e^x)` |
| **Integrate** | Definite/indefinite integrals | `∫ sin(x) dx` from 0 to π |
| **Analyze** | Simplify, expand, factor, transform | Simplify `(x^2-1)/(x-1)` |
| **Simulate** | ODEs, PDEs, physics simulations | Solve `dy/dt = -k*y` |

## 📚 Examples

### 1. Solve an Equation

```rust
let request = ToolRequest::Solve(SolveInput {
    equations: vec!["x^2 - 4 = 0".to_string()],
    variables: None,
    initial_guess: None,
    domain: None,
    method: None,
});
```

### 2. Differentiate

```rust
let request = ToolRequest::Differentiate(DifferentiateInput {
    expression: "x^2*y".to_string(),
    variables: vec!["x".to_string(), "y".to_string()],
    order: None,
    evaluate_at: None,
});
```

### 3. Integrate

```rust
let request = ToolRequest::Integrate(IntegrateInput {
    expression: "sin(x)".to_string(),
    variables: vec!["x".to_string()],
    limits: Some(vec![[0.0, std::f64::consts::PI]]),
    method: None,
});
```

### 4. Simplify

```rust
let request = ToolRequest::Analyze(AnalyzeInput {
    operation: AnalyzeOperation::Simplify,
    expression: "(x^2 - 1)/(x - 1)".to_string(),
    options: HashMap::new(),
});
```

### 5. Simulate ODE

```rust
let mut initial_conditions = HashMap::new();
initial_conditions.insert("y".to_string(), 10.0);

let request = ToolRequest::Simulate(SimulateInput {
    model: SimulationModel::ODE,
    equations: vec!["dy/dt = -k*y".to_string()],
    variables: vec!["y".to_string()],
    parameters: HashMap::from([("k".to_string(), 0.1)]),
    initial_conditions: Some(initial_conditions),
    range: Some([0.0, 10.0]),
    steps: Some(100),
    method: None,
});
```

## 🔧 Multiple Interfaces

### 1. Rust API (Native)

Direct Rust integration with full type safety:

```rust
let dispatcher = create_default_dispatcher();
let response = dispatcher.dispatch(request).unwrap();
```

### 2. JSON API (MCP Server)

JSON interface for AI agents and CLI tools:

```json
{
  "tool": "solve",
  "input": {
    "equations": ["x^2 - 4 = 0"]
  }
}
```

Run as MCP server:
```bash
echo '{"tool":"solve","input":{"equations":["x^2 - 4 = 0"]}}' | computational-engine stdin
```

### 3. WebAssembly (Browser/Node.js)

Use in JavaScript/TypeScript for on-device computation:

```javascript
import init, { ComputationalEngine } from '@brainwires/computational-engine';

await init();
const engine = new ComputationalEngine();

const result = engine.solve({
    equations: ["x^2 - 4 = 0"]
});
```

See [WASM.md](WASM.md) for complete WebAssembly documentation.

## 🏗️ Architecture

The engine follows a clean layered architecture:

```
JSON API ─→ ToolDispatcher ─→ Core Traits ─→ Implementations ─→ Domain Modules
```

- **Core Traits**: Define the 5-tool interface
- **Implementations**: Default implementations of each tool
- **Domain Modules**: Specialized math/physics modules (legacy, still supported)

See [ARCHITECTURE.md](ARCHITECTURE.md) for detailed documentation.

## 📦 Features

- ✅ **Equation Solving**: Polynomial, linear systems, nonlinear
- ✅ **Calculus**: Symbolic differentiation and integration
- ✅ **Algebra**: Simplification, expansion, factoring
- ✅ **Simulation**: ODE/PDE solvers, physics engines
- ✅ **Linear Algebra**: Matrix operations, decompositions
- ✅ **Signal Processing**: FFT, filters, transforms
- ✅ **Tensor Calculus**: Einstein summation, symbolic tensors
- ✅ **Quantum Physics**: Quantum simulations
- ✅ **Fluid Dynamics**: Navier-Stokes solvers
- ✅ **Optimization**: Gradient descent, numerical optimization
- ✅ **Statistics**: Probability distributions, statistical analysis
- ✅ **Cryptography**: Number theory, cryptographic primitives

## 🚧 Domain Modules (Legacy API)

The original domain-specific API is still supported for backwards compatibility:

### Mathematics
- Tensor Calculus
- Advanced Calculus (fractional, variational, stochastic)
- Linear Algebra
- Symbolic Regression
- Special Functions

### Physics
- Fluid Dynamics
- Quantum Physics
- Electromagnetism

### Tools
- Signal Processing
- Dimensional Analysis
- Equation Validation
- Computational Geometry
- Numerical Methods

### Specialized
- Stochastic Processes
- Cryptographic Mathematics
- Statistics
- Optimization
- Graph Theory
- Information Theory
- Chemistry

## 🔬 Building & Running

### Native Build

```bash
# Build the library and CLI
cargo build --release

# Run the MCP server
./target/release/computational-engine stdin

# Run tests
cargo test

# Build documentation
cargo doc --open
```

### WebAssembly Build

```bash
# Build for all WASM targets
./build-wasm.sh

# Or build specific targets
npm run build:web        # Browser ES modules
npm run build:nodejs     # Node.js
npm run build:bundler    # Webpack/Vite/Rollup

# Run browser example
cd examples/wasm && npm run dev

# Run Node.js example
cd examples/wasm && npm run node
```

See [WASM.md](WASM.md) for detailed WASM build instructions.

## 🎯 Use Cases

- Mathematical modeling and simulation
- Scientific computing
- AI/ML mathematical operations
- Educational tools
- Research and prototyping
- MCP servers for AI agents
- Web services requiring math capabilities

## 📖 Documentation

- [ARCHITECTURE.md](ARCHITECTURE.md) - Detailed architecture documentation
- [WASM.md](WASM.md) - WebAssembly build and usage guide
- [examples/](examples/) - Working code examples
- [examples/wasm/](examples/wasm/) - Browser and Node.js WASM examples
- API documentation: `cargo doc --open`

## 🛠️ Development

```bash
# Check code
cargo check

# Run tests
cargo test

# Format code
cargo fmt

# Lint
cargo clippy
```

## 📝 License

MIT OR Apache-2.0

## 🤝 Contributing

Contributions welcome! The 5-tool architecture makes it easy to:

1. Add new solving methods
2. Implement new analysis operations
3. Create new simulation models
4. Integrate specialized libraries

See the implementation files in `src/implementations/` for examples.

## 🌟 Design Philosophy

- **Simplicity**: 5 intuitive tools instead of dozens of functions
- **Power**: Each tool handles a broad category of operations
- **Flexibility**: Rich input schemas with sensible defaults
- **Extensibility**: Easy to add new capabilities
- **Type Safety**: Full Rust type checking
- **JSON Ready**: Seamless JSON serialization
