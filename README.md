# Brainwires Compute Engine

A comprehensive computational engine for mathematical and scientific computing in Rust, providing **194+ operations** through a clean **10-tool API**.

[![Tests](https://img.shields.io/badge/tests-523%20passing-brightgreen)](https://github.com/nightness/brainwires-compute-engine)
[![Coverage](https://img.shields.io/badge/coverage-1500%2B%20tests-blue)](https://github.com/nightness/brainwires-compute-engine)
[![Rust](https://img.shields.io/badge/rust-2024%20edition-orange)](https://www.rust-lang.org/)

## 🚀 Quick Start

```rust
use computational_engine::engine::{ToolRequest, SolveInput};
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

## 🧮 The 10 Core Tools

| Tool | Purpose | Example Operations |
|------|---------|-------------------|
| **Solve** | Equations & systems | Polynomial, linear systems, differential equations, Einstein field equations |
| **Differentiate** | Derivatives | Symbolic, numeric, partial, implicit, parametric, directional |
| **Integrate** | Integration | Definite, indefinite, numeric, improper, multiple, contour |
| **Analyze** | Analysis & transforms | Series, limits, roots, extrema, Fourier, wavelets, stability |
| **Simulate** | Dynamic systems | ODEs, PDEs, stochastic processes, Monte Carlo, cellular automata |
| **Compute** | Matrix & tensor ops | Linear algebra, tensor calculus, special functions, scientific formulas |
| **Transform** | Signal transforms | FFT, Fourier, Laplace, wavelets, filters, window functions |
| **FieldTheory** | Physics fields | Electromagnetic, gravitational, quantum fields, Green's functions |
| **Sample** | Statistical sampling | Monte Carlo, MCMC, distributions, bootstrap, signal analysis |
| **Optimize** | Optimization | Curve fitting, interpolation, minimization, symbolic regression |

## ✨ Key Features

### Core Capabilities
- ✅ **523 unit tests** with 100% pass rate
- ✅ **1500+ total tests** (comprehensive integration + unit tests)
- ✅ **Custom Computer Algebra System (CAS)** - no Python dependencies
- ✅ **Multiple interfaces**: Native Rust, JSON/MCP, WebAssembly
- ✅ **Rust 2024 edition** with full type safety

### Mathematics
- **Calculus**: Symbolic differentiation, integration, series analysis
- **Linear Algebra**: Matrix operations, decompositions (SVD, LU, QR, Cholesky), eigenvalues, PCA
- **Tensor Calculus**: Einstein summation, Christoffel symbols, Riemann curvature, metric tensors
- **Symbolic CAS**: Expression parsing, simplification, symbolic differentiation/integration
- **Special Functions**: Bessel, gamma, beta, error functions, elliptic integrals, orthogonal polynomials
- **Numerical Methods**: Root finding, interpolation, ODE/PDE solvers, numerical integration

### Physics
- **Quantum Mechanics**: Wave functions, operators, perturbation theory
- **Relativity**: Special and general relativity calculations
- **Electromagnetism**: Maxwell's equations, EM waves, antennas, waveguides
- **Nuclear Physics**: Radioactive decay, binding energy, fission/fusion
- **Fluid Dynamics**: Navier-Stokes solvers, boundary conditions, flow analysis
- **Control Systems**: Transfer functions, stability analysis, PID tuning

### Science Formulas
- **Chemistry**: Gas laws, pH calculations, molar mass, thermochemistry (23 tests)
- **Biology**: Michaelis-Menten kinetics, Hardy-Weinberg equilibrium, pharmacokinetics (19 tests)
- **Thermodynamics**: Heat transfer, entropy, thermodynamic cycles (16 tests)
- **Optics**: Thin lens, diffraction, interference, polarization (14 tests)
- **Engineering**: Acoustics, materials science, fluid mechanics, control theory (20 tests)
- **Geophysics**: Seismology, atmospheric physics, radiometric dating, planetary science (40 tests)
- **DateTime**: Date arithmetic, business days, leap years, time zones (29 tests)

### Specialized Modules
- **Signal Processing**: FFT, filters, spectrograms, wavelets, window functions
- **Statistics**: Distributions, hypothesis testing, regression, MCMC
- **Optimization**: Gradient descent, Nelder-Mead, curve fitting, symbolic regression
- **Graph Theory**: Shortest paths, MST, connected components, topological sort
- **Information Theory**: Entropy, mutual information, channel capacity, Huffman coding
- **Cryptography**: RSA, prime generation, modular arithmetic
- **Computational Geometry**: Convex hull, Delaunay triangulation, Voronoi diagrams

## 📦 Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
brainwires-compute-engine = { git = "https://github.com/nightness/brainwires-compute-engine" }
```

## 🔧 Multiple Interfaces

### 1. Rust API (Native)

Direct Rust integration with full type safety:

```rust
use computational_engine::engine::{ToolRequest, DifferentiateInput};
use computational_engine::implementations::create_default_dispatcher;

let dispatcher = create_default_dispatcher();

let request = ToolRequest::Differentiate(DifferentiateInput {
    expression: "x^2*sin(x)".to_string(),
    variables: vec!["x".to_string()],
    order: None,
    evaluate_at: None,
});

let response = dispatcher.dispatch(request).unwrap();
```

### 2. JSON API (MCP Server)

JSON interface for AI agents and CLI tools:

```bash
# Run as MCP server (reads JSON from stdin)
echo '{"tool":"solve","input":{"equations":["x^2 - 4 = 0"]}}' | cargo run --release -- stdin

# List available operations
cargo run --release -- list-ops

# Display version info
cargo run --release -- info
```

JSON request format:
```json
{
  "tool": "solve",
  "input": {
    "equations": ["x^2 - 4 = 0"]
  }
}
```

### 3. WebAssembly (Browser/Node.js)

Use in JavaScript/TypeScript for on-device computation:

```javascript
import init, { ComputationalEngine } from '@brainwires/brainwires-compute-engine';

await init();
const engine = new ComputationalEngine();

const result = engine.solve({
    equations: ["x^2 - 4 = 0"]
});
console.log(result);
```

See [WASM.md](WASM.md) for complete WebAssembly documentation.

## 🏗️ Architecture

Clean layered architecture with 10 unified tools:

```
JSON/WASM API → ToolDispatcher → Tool Traits → Unified Implementations → Domain Modules
```

- **Tool Traits**: Define the 10-tool interface (`src/engine/traits.rs`)
- **Unified Implementations**: Wire tools to domain modules (`src/implementations/`)
- **Domain Modules**: Specialized algorithms (`src/mathematics/`, `src/physics/`, etc.)

See [ARCHITECTURE.md](ARCHITECTURE.md) for detailed documentation.

## 🔬 Building & Testing

### Native Build

```bash
# Build library and CLI
cargo build --release

# Run the MCP server
./target/release/brainwires-compute-engine stdin

# Run all tests (523 unit + integration tests)
cargo test

# Run specific test suites
cargo test --test linear_algebra_comprehensive_tests
cargo test --test tensor_calculus_comprehensive_tests
cargo test mathematics::special_functions

# Generate documentation
cargo doc --open

# Code quality
cargo check
cargo fmt
cargo clippy
```

### WebAssembly Build

```bash
# Build for all WASM targets
./build-wasm.sh

# Or build specific targets
npm run build:web        # Browser ES modules
npm run build:nodejs     # Node.js
npm run build:bundler    # Webpack/Vite/Rollup
npm run build:no-modules # Classic <script> tags

# Run browser example
cd examples/wasm && npm run dev

# Run Node.js example
cd examples/wasm && npm run node

# Test WASM
wasm-pack test --headless --firefox
```

See [WASM.md](WASM.md) for detailed WASM build instructions.

## 📊 Test Coverage

Comprehensive test suite with **100% pass rate**:

- **523 unit tests** in library modules
- **40 integration test files** with 1000+ additional tests
- **Coverage**: Science (181 tests), Physics (178 tests), Math (113+ tests), Tools (142 tests), Specialized (83 tests)

### Test Commands

```bash
# Run all tests
cargo test

# Run with coverage report
cargo test --test '*_comprehensive_tests'

# Test specific modules
cargo test chemistry
cargo test physics::electromagnetism
cargo test mathematics::linear_algebra
```

## 🎯 Use Cases

- **Scientific Computing**: Research calculations, simulations, data analysis
- **AI/ML Applications**: Mathematical operations for AI agents and models
- **Educational Tools**: Interactive mathematics and physics education
- **MCP Servers**: Computational backend for Claude and other AI assistants
- **Web Services**: Server-side or client-side mathematical computation
- **Research**: Rapid prototyping of mathematical models
- **Engineering**: CAD, simulation, optimization problems

## 📖 Documentation

- [ARCHITECTURE.md](ARCHITECTURE.md) - Detailed 10-tool architecture
- [API.md](API.md) - Complete JSON API reference
- [WASM.md](WASM.md) - WebAssembly build and usage guide
- [OPERATIONS_COMPLETE_LIST.md](OPERATIONS_COMPLETE_LIST.md) - All 194+ operations
- [MATHEMATICA_COMPETITION_PROGRESS.md](MATHEMATICA_COMPETITION_PROGRESS.md) - Feature parity tracking
- [STDOUT_STDERR_POLICY.md](STDOUT_STDERR_POLICY.md) - MCP server compatibility
- API documentation: `cargo doc --open`

## 🛠️ Development

### Prerequisites
- Rust 1.75+ (2024 edition)
- wasm-pack (for WebAssembly builds)
- Node.js 18+ (for WASM examples)

### Development Commands
```bash
# Check code
cargo check

# Run tests
cargo test

# Format code
cargo fmt

# Lint
cargo clippy

# Build release
cargo build --release

# Run benchmarks
cargo bench
```

## 🌟 Design Philosophy

- **Simplicity**: 10 intuitive tools instead of 916 individual functions
- **Comprehensiveness**: 194+ operations covering mathematics, physics, and science
- **Self-Contained**: Custom CAS with no Python dependencies
- **Type Safety**: Full Rust type checking with zero-cost abstractions
- **Performance**: Release builds optimized with LTO and single codegen unit
- **Testability**: 100% test pass rate with comprehensive coverage
- **Flexibility**: Rich input schemas with sensible defaults
- **Interoperability**: Native Rust, JSON, and WebAssembly interfaces

## 📝 License

MIT OR Apache-2.0

## 🤝 Contributing

Contributions welcome! The 10-tool architecture makes it easy to:

1. Add new solving methods to existing tools
2. Implement new analysis operations
3. Create new simulation models
4. Add scientific formula modules
5. Integrate specialized libraries

See the implementation files in `src/implementations/` for examples.

## 🙏 Acknowledgments

This project aims to provide comprehensive computational capabilities comparable to Mathematica/Wolfram Alpha, but as free and open-source software.

## 📞 Support

- Issues: [GitHub Issues](https://github.com/nightness/brainwires-compute-engine/issues)
- Documentation: See `docs/` directory and inline documentation

---

**Built with Rust 🦀 | Powered by Mathematics 📐 | Tested with Rigor ✅**
