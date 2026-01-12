# Brainwires Compute Engine

A comprehensive computational engine for mathematical and scientific computing in Rust, providing **400+ mathematical operations** through a clean **8-tool API**.

[![Tests](https://img.shields.io/badge/tests-2067%20passing-brightgreen)](https://github.com/nightness/brainwires-compute-engine)
[![Coverage](https://img.shields.io/badge/coverage-79.99%25-brightgreen)](https://github.com/nightness/brainwires-compute-engine)
[![Rust](https://img.shields.io/badge/rust-2024%20edition-orange)](https://www.rust-lang.org/)

## 🚀 Quick Start

```rust
use computational_engine::{ToolRequest, SolveInput, create_default_dispatcher};

// Create dispatcher
let dispatcher = create_default_dispatcher();

// Solve a quadratic equation
let request = ToolRequest::Solve(SolveInput {
    equations: vec!["x^2 + 2*x - 8 = 0".to_string()],
    variables: None,
    initial_guess: None,
    domain: None,
    method: None,
    ..Default::default()
});

let response = dispatcher.dispatch(request).unwrap();
println!("{:#?}", response);
```

## 🧮 The 8 Core Tools

| Tool | Purpose | Example Operations |
|------|---------|-------------------|
| **Solve** | Equations & optimization | Polynomial, linear systems, differential equations, Einstein field equations, curve fitting, minimization |
| **Compute** | Math & physics computation | Matrix ops, calculus (differentiate/integrate), transforms (FFT, Laplace), field theory (EM, gravity), sampling, tensor calculus, special functions |
| **Analyze** | Expression analysis | Series expansion, limits, roots, extrema, stability analysis, simplification |
| **Simulate** | Dynamic systems | ODEs, PDEs, stochastic processes, Monte Carlo, fluid dynamics |
| **ML** | Machine learning | Clustering (K-means), regression, neural networks, PCA, dimensionality reduction |
| **Chaos** | Chaos theory | Fractals (Mandelbrot, Julia), attractors (Lorenz), Lyapunov exponents, bifurcation |
| **Units** | Dimensional analysis | Unit conversion, dimensional consistency checking |
| **Validate** | Equation validation | Physics compliance, conservation laws, symmetry checks |

## ✨ Key Features

### Core Capabilities
- ✅ **2,067 tests** with 100% pass rate (1,089 integration + 960 unit + 18 doc tests)
- ✅ **~80% production code coverage** (measured with cargo-tarpaulin)
- ✅ **400+ mathematical operations** across 8 unified tools
- ✅ **Adaptive integration** with automatic subdivision and error estimation
- ✅ **24+ traditional mathematical subjects** including calculus, linear algebra, statistics, differential equations, graph theory, optimization, and more
- ✅ **Custom Computer Algebra System (CAS)** - no Python dependencies
- ✅ **Multiple interfaces**: Native Rust, JSON/MCP, WebAssembly
- ✅ **Rust 2024 edition** with full type safety

### Mathematics
- **Calculus**: Symbolic differentiation, integration, series analysis
- **Linear Algebra**: Matrix operations, decompositions (SVD, LU, QR, Cholesky), eigenvalues, PCA
- **Tensor Calculus**: Einstein summation, Christoffel symbols, Riemann curvature, metric tensors
- **Symbolic CAS**: Expression parsing, simplification, symbolic differentiation/integration
- **Special Functions**: Bessel, gamma, beta, error functions, elliptic integrals, orthogonal polynomials
- **Numerical Methods**: Root finding, interpolation, ODE/PDE solvers, adaptive numerical integration (Simpson's rule with automatic subdivision)

### Physics
- **Quantum Mechanics**: Wave functions, operators, perturbation theory
- **Relativity**: Special and general relativity calculations
- **Electromagnetism**: Maxwell's equations, EM waves, antennas, waveguides
- **Nuclear Physics**: Radioactive decay, binding energy, fission/fusion
- **Fluid Dynamics**: Navier-Stokes solvers, boundary conditions, flow analysis
- **Control Systems**: Transfer functions, stability analysis, PID tuning

### Science Formulas
- **Chemistry**: Gas laws, pH calculations, molar mass, thermochemistry
- **Biology**: Michaelis-Menten kinetics, Hardy-Weinberg equilibrium, pharmacokinetics
- **Thermodynamics**: Heat transfer, entropy, thermodynamic cycles
- **Optics**: Thin lens, diffraction, interference, polarization
- **Engineering**: Acoustics, materials science, fluid mechanics, control theory
- **Geophysics**: Seismology, atmospheric physics, radiometric dating, planetary science
- **DateTime**: Date arithmetic, business days, leap years, time zones

### Specialized Modules
- **Signal Processing**: FFT, filters, spectrograms, wavelets, window functions
- **Statistics**: Distributions, hypothesis testing, regression, MCMC
- **Optimization**: Gradient descent, Nelder-Mead, curve fitting, symbolic regression
- **Graph Theory**: Shortest paths, MST, connected components, topological sort
- **Information Theory**: Entropy, mutual information, channel capacity, Huffman coding
- **Cryptography**: RSA, prime generation, modular arithmetic
- **Computational Geometry**: Convex hull, Delaunay triangulation, Voronoi diagrams

## 📚 Traditional Mathematical Subjects Covered

The engine spans **24+ traditional mathematical subjects** familiar to educators and students:

### Core Mathematics
- **Differential Calculus** - Symbolic differentiation, partial derivatives, chain rule, gradient, divergence, curl
- **Integral Calculus** - Definite/indefinite integrals, multiple integrals, contour integration
- **Advanced Calculus** - Taylor series, limits, L'Hôpital's rule, convergence tests
- **Linear Algebra** - Matrix operations, decompositions (SVD, QR, LU, Cholesky), eigenvalues, PCA
- **Symbolic Algebra** - Expression parsing, simplification, factorization, polynomial operations
- **Differential Equations** - ODEs, PDEs (heat, wave, Laplace), stiff equations

### Advanced Pure Mathematics
- **Tensor Calculus** - Christoffel symbols, Riemann tensors, Einstein equations, metric tensors
- **Number Theory** - Primes, modular arithmetic, GCD/LCM, Chinese remainder theorem
- **Special Functions** - Bessel, Gamma, Error functions, Elliptic integrals, orthogonal polynomials
- **Computational Geometry** - Convex hull, Delaunay triangulation, Voronoi diagrams

### Applied Mathematics
- **Statistics & Probability** - Distributions, hypothesis testing, regression, ANOVA, correlation
- **Optimization** - Gradient descent, Nelder-Mead, curve fitting, symbolic regression
- **Stochastic Processes** - Brownian motion, Poisson processes, Lévy processes, MCMC
- **Graph Theory** - Shortest paths, MST, topological sorting, connected components
- **Information Theory** - Entropy, mutual information, KL divergence, channel capacity

### Signal & Transform Analysis
- **Fourier Analysis** - FFT, Fourier/Laplace transforms, Fourier series
- **Wavelet Analysis** - Haar, Daubechies, Morlet wavelets, wavelet transforms
- **Signal Processing** - Digital filters, spectral analysis, autocorrelation, power spectrum

### Physics Mathematics
- **Quantum Mechanics** - Schrödinger equation, perturbation theory, wave functions
- **Relativity** - Lorentz transformations, Schwarzschild metric, gravitational effects
- **Statistical Physics** - Partition functions, Boltzmann/Fermi-Dirac/Bose-Einstein distributions
- **Control Theory** - Transfer functions, Bode plots, stability analysis, state-space
- **Nuclear Physics** - Radioactive decay, binding energy, fission/fusion calculations

### Cryptography & Security
- **Cryptographic Mathematics** - RSA, hashing (SHA256/SHA3), discrete logarithm, elliptic curves

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
use computational_engine::{ToolRequest, ComputeInput, ComputeOp, DifferentiateType, create_default_dispatcher};

let dispatcher = create_default_dispatcher();

// Differentiate x^2*sin(x) symbolically
let request = ToolRequest::Compute(ComputeInput {
    operation: ComputeOp::Differentiate(DifferentiateType::Symbolic),
    data: serde_json::json!({
        "expression": "x^2*sin(x)",
        "variable": "x"
    }),
    parameters: Default::default(),
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
import init, { ComputationalEngine } from '@brainwires/compute-engine';

await init();
const engine = new ComputationalEngine();

const result = engine.solve({
    equations: ["x^2 - 4 = 0"]
});
console.log(result);
```

See [WASM.md](WASM.md) for complete WebAssembly documentation.

## 🏗️ Architecture

Clean layered architecture with 8 unified tools:

```
JSON/WASM API → ToolDispatcher → Tool Traits → Unified Implementations → Domain Modules
```

### Source Structure
```
src/
├── engine/      # Core types, traits, dispatcher
├── solve/       # SOLVE tool + optimization/, specialized/
├── compute/     # COMPUTE tool + 27 domain submodules
├── analyze/     # ANALYZE tool + symbolic/
├── simulate/    # SIMULATE tool + fluids/, stochastic/, ode/, finance/
├── ml/          # ML tool + clustering/, regression/, neural_network/, etc.
├── chaos/       # CHAOS tool + fractals/, attractors/, lyapunov/, etc.
├── units/       # UNITS tool + dimensional_analysis
└── validate/    # VALIDATE tool + equation validation
```

See [ARCHITECTURE.md](ARCHITECTURE.md) for detailed documentation.

## 🔬 Building & Testing

### Native Build

```bash
# Build library and CLI
cargo build --release

# Run the MCP server
./target/release/brainwires-compute-engine stdin

# Run all tests (2,067 tests total)
cargo test

# Run integration tests only
cargo test --test all_integration_tests

# Run specific test suites
cargo test --test all_integration_tests integration::tools::linear_algebra
cargo test --test all_integration_tests integration::tools::tensor_calculus

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

Comprehensive test suite with **100% pass rate** and **~80% production code coverage**:

### Test Statistics
- **2,067 total tests** (100% pass rate)
  - **1,089 integration tests** in `tests/integration/`
  - **960 unit tests** (inline `#[cfg(test)]` modules in `src/`)
  - **18 doc tests**
- **~80% production code coverage** (79.99% as of January 2025)
  - Coverage measured with `cargo-tarpaulin` which excludes test code
  - Production code only, inline test modules automatically excluded

### Coverage Commands

```bash
# Run all tests
cargo test

# Get production code coverage (recommended - excludes test code)
cargo tarpaulin --lib --release --out Stdout

# Alternative: cargo-llvm-cov (may include test code in metrics)
cargo llvm-cov --lib --release

# Generate HTML coverage report
cargo tarpaulin --lib --release --out Html
```

**Note:** Use `cargo-tarpaulin` for accurate production code coverage. The tool automatically excludes inline `#[cfg(test)]` modules from coverage calculations, giving you the true production code coverage percentage.

### Test Organization

The project uses two test patterns:

1. **Integration Tests** (`tests/integration/` directory)
   - 1,089 integration tests across multiple files
   - Organized by tool: `engine/`, `compute/`, `solve/`, `simulate/`, `tools/`, `coverage/`
   - Test the public API and end-to-end functionality

2. **Unit Tests** (inline `#[cfg(test)]` modules)
   - 960 unit tests embedded in source files
   - Located in `#[cfg(test)]` modules within `src/` files
   - Test individual functions and internal implementation details

See [TESTING.md](TESTING.md) for detailed testing guidelines and best practices.

## 🎯 Use Cases

- **Scientific Computing**: Research calculations, simulations, data analysis
- **AI/ML Applications**: Mathematical operations for AI agents and models
- **Educational Tools**: Interactive mathematics and physics education
- **MCP Servers**: Computational backend for Claude and other AI assistants
- **Web Services**: Server-side or client-side mathematical computation
- **Research**: Rapid prototyping of mathematical models
- **Engineering**: CAD, simulation, optimization problems

## 📖 Documentation

- **[TESTING.md](TESTING.md)** - 🧪 **Complete testing guide and best practices**
- **[CLAUDE.md](CLAUDE.md)** - 🤖 **Development guide for Claude Code and contributors**
- [WASM.md](WASM.md) - WebAssembly build and usage guide
- [ARCHITECTURE.md](ARCHITECTURE.md) - Detailed architecture documentation
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

- **Simplicity**: 8 intuitive tools providing access to 400+ mathematical operations
- **Comprehensiveness**: Operations covering 24+ traditional mathematical subjects from calculus and linear algebra to quantum mechanics and control theory
- **Self-Contained**: Custom CAS with no Python dependencies
- **Type Safety**: Full Rust type checking with zero-cost abstractions
- **Performance**: Release builds optimized with LTO and single codegen unit
- **Testability**: 100% test pass rate with 2,067 tests
- **Flexibility**: Rich input schemas with sensible defaults
- **Interoperability**: Native Rust, JSON, and WebAssembly interfaces

## 🆕 Recent Updates

### January 2025
- ✅ **8-Tool Architecture Consolidation**: Unified architecture with 8 primary tools
  - Cleaner, more intuitive API
- ✅ **Test Reorganization**: Tests now mirror source structure
  - `tests/unit/` organized by tool: `solve/`, `compute/`, `analyze/`, `simulate/`, `ml/`, `chaos/`, `units/`, `validate/`
  - `tests/integration/` organized by tool: `engine/`, `compute/`, `solve/`, `simulate/`, `tools/`, `coverage/`, `other/`
  - Total: 2,067 tests (up from 1,685)

### October 2024
- ✅ **Major Test Coverage Improvement**: Significant increase in test count
  - Achieved **~80% production code coverage**
  - Fixed critical black hole physics bug in spin parameter interpretation
  - Created comprehensive [TESTING.md](TESTING.md) guide
- ✅ **Black Hole Physics Fix**: Corrected Kerr black hole implementation
  - Changed spin parameter from mass units to dimensionless (0-1)
  - Fixed event horizon, ergosphere, and ISCO calculations

## 📝 License

MIT OR Apache-2.0

## 🤝 Contributing

Contributions welcome! The 8-tool architecture makes it easy to:

1. Add new solving methods to existing tools
2. Implement new analysis operations
3. Create new simulation models
4. Add scientific formula modules
5. Integrate specialized libraries

See the tool modules in `src/solve/`, `src/compute/`, etc. for examples.

## 🙏 Acknowledgments

This project aims to provide comprehensive computational capabilities comparable to Mathematica/Wolfram Alpha, but as free and open-source software.

## 📞 Support

- Issues: [GitHub Issues](https://github.com/nightness/brainwires-compute-engine/issues)
- Documentation: See `docs/` directory and inline documentation

---

**Built with Rust 🦀 | Powered by Mathematics 📐 | Tested with Rigor ✅**
