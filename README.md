# Brainwires Compute Engine

A comprehensive computational engine for mathematical and scientific computing in Rust, providing **406 mathematical operations** through a clean **10-tool API**.

[![Tests](https://img.shields.io/badge/tests-1131%20passing-brightgreen)](https://github.com/nightness/brainwires-compute-engine)
[![Coverage](https://img.shields.io/badge/coverage-44.78%25-yellow)](https://github.com/nightness/brainwires-compute-engine)
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
| **Integrate** | Integration | Definite, indefinite, adaptive numeric, improper, multiple, contour |
| **Analyze** | Analysis & transforms | Series, limits, roots, extrema, Fourier, wavelets, stability |
| **Simulate** | Dynamic systems | ODEs, PDEs, stochastic processes, Monte Carlo, cellular automata |
| **Compute** | Matrix & tensor ops | Linear algebra, tensor calculus, special functions, scientific formulas |
| **Transform** | Signal transforms | FFT, Fourier, Laplace, wavelets, filters, window functions |
| **FieldTheory** | Physics fields | Electromagnetic, gravitational, quantum fields, Green's functions |
| **Sample** | Statistical sampling | Monte Carlo, MCMC, distributions, bootstrap, signal analysis |
| **Optimize** | Optimization | Curve fitting, interpolation, minimization, symbolic regression |

## ✨ Key Features

### Core Capabilities
- ✅ **1,139 tests** with 100% pass rate (1,131 passing, 8 ignored for performance)
- ✅ **163 test files** professionally organized (123 unit + 40 integration)
- ✅ **406 mathematical operations** across 10 unified tools
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
- **Fluid Dynamics**: Navier-Stokes solvers, boundary conditions, flow analysis (tests disabled for CI performance)
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

# Run all tests (1,139 tests in 163 files)
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

Comprehensive test suite with **100% pass rate**:

### Test Statistics
- **1,139 total tests** (1,131 passing, 8 ignored for performance)
- **163 test files** organized professionally:
  - **123 unit test files** in `tests/unit/` (mirror `src/` structure)
  - **40 integration test files** in `tests/integration/` (organized by category)

### Test Distribution by Category
- **Tools**: 569 tests (10 unified tool comprehensive tests)
- **API**: 182 tests (API & engine functionality)
- **Physics**: 178 tests (quantum, relativity, EM, control, nuclear, statistical)
- **Modules**: 134 tests (crypto, datetime, graph theory, info theory, etc.)
- **Coverage**: ~76 tests (complete operation coverage verification)

### Test Organization
```
tests/
├── all_integration_tests.rs  # Entry point for all integration tests
├── integration/               # 40 integration test files
│   ├── tools/                # 13 files - Tool comprehensive tests
│   ├── physics/              # 6 files - Physics module tests
│   ├── modules/              # 9 files - Module-specific tests
│   ├── api/                  # 7 files - API & engine tests
│   ├── coverage/             # 4 files - Coverage verification
│   └── other/                # 1 file - Uncategorized
└── unit/                     # 123 unit test files (mirrors src/)
    ├── mathematics/
    ├── physics/
    ├── specialized/
    └── ...
```

### Test Commands

```bash
# Run all tests (excludes 8 ignored slow tests)
cargo test

# Run integration tests
cargo test --test all_integration_tests

# Run specific test categories
cargo test --test all_integration_tests integration::tools::
cargo test --test all_integration_tests integration::physics::
cargo test --test all_integration_tests integration::modules::

# Run ignored tests (slow physics simulations)
cargo test -- --ignored

# Run ALL tests including ignored
cargo test -- --include-ignored
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

- **[COMPARISON_ANALYSIS.md](COMPARISON_ANALYSIS.md)** - 📊 **How we compare to Mathematica & Wolfram Alpha**
- [ARCHITECTURE.md](ARCHITECTURE.md) - Detailed 10-tool architecture
- [API.md](API.md) - Complete JSON API reference
- [WASM.md](WASM.md) - WebAssembly build and usage guide
- [OPERATIONS_COMPLETE_LIST.md](OPERATIONS_COMPLETE_LIST.md) - All 406 operations inventory
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

- **Simplicity**: 10 intuitive tools providing access to 406 mathematical operations
- **Comprehensiveness**: 406 operations covering 24+ traditional mathematical subjects from calculus and linear algebra to quantum mechanics and control theory
- **Self-Contained**: Custom CAS with no Python dependencies
- **Type Safety**: Full Rust type checking with zero-cost abstractions
- **Performance**: Release builds optimized with LTO and single codegen unit
- **Testability**: 100% test pass rate with 1,139 tests across 163 professionally organized files
- **Flexibility**: Rich input schemas with sensible defaults
- **Interoperability**: Native Rust, JSON, and WebAssembly interfaces

## 🆕 Recent Updates

### October 2024
- ✅ **Test Reorganization**: Complete refactoring of test structure
  - Extracted all 123 inline tests to separate files mirroring `src/` structure
  - Organized 40 integration tests by category (tools, physics, modules, api, coverage)
  - Professional test organization following Rust best practices
- ✅ **Documentation Updates**: Accurate operation counts from source code (406 operations)
- ✅ **Bug Fixes**: Fixed 8 critical bugs (Spearman correlation, optimization flags, variance calculation)
- ✅ **Test Coverage**: 1,139 tests with 100% pass rate (8 ignored for performance)

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
