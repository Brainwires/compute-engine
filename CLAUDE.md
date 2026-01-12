# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

A unified computational engine for mathematical and scientific computing in Rust, providing **400+ mathematical operations** accessible through a **consolidated 8-tool API**:

**Primary Tools (8):**
- **Solve**: Equations, systems, optimization, root finding
- **Compute**: Calculus, transforms, field theory, sampling, matrix ops
- **Analyze**: Series, limits, stability analysis, simplification
- **Simulate**: Time evolution, stochastic processes, fluid dynamics
- **ML**: Machine learning (clustering, regression, neural networks)
- **Chaos**: Chaos theory (fractals, attractors, Lyapunov exponents)
- **Units**: Dimensional analysis and unit conversion
- **Validate**: Equation and physics validation

**Legacy tool names** (differentiate, integrate, transform, fieldtheory, sample, optimize) are maintained for backward compatibility and route internally to the primary tools.

Supports native Rust, WebAssembly (browser/Node.js), and MCP server interfaces.

**Comprehensive Test Coverage**: 2,067 tests passing (1,089 integration tests + 960 unit tests + 18 doc tests) with 100% pass rate.

## Build Commands

### Native Rust Build
```bash
# Build library and CLI binary
cargo build --release

# Run the MCP server (JSON over stdin/stdout)
cargo run --release -- stdin

# Build with quantum physics support (requires GPU)
cargo build --release --features quantum

# Run tests
cargo test

# Generate documentation
cargo doc --open

# Code quality checks
cargo check
cargo fmt
cargo clippy
```

### WebAssembly Build
```bash
# Build all WASM targets (requires wasm-pack)
./build-wasm.sh

# Or build specific targets
npm run build:bundler    # For webpack/rollup/vite
npm run build:nodejs     # For Node.js
npm run build:web        # For ES modules in browser
npm run build:no-modules # For classic <script> tags

# Test WASM builds
wasm-pack test --headless --firefox
```

### CLI Usage Examples
```bash
# List all available operations
cargo run --release -- list-ops

# Display version and module info
cargo run --release -- info

# Execute JSON request from file
cargo run --release -- json --request '@request.json'

# Execute inline JSON (use single quotes)
cargo run --release -- json --request '{"tool":"solve","input":{"equations":["x^2 - 4 = 0"]}}'

# MCP server mode (reads JSON from stdin, outputs JSON to stdout)
echo '{"tool":"solve","input":{"equations":["x^2 - 4 = 0"]}}' | cargo run --release -- stdin
```

## Architecture Overview

### Consolidated 8-Tool Architecture

**Public API Layer**: 8 primary tools + 6 legacy tools for backward compatibility
- Defined in `src/engine/traits.rs` and `src/engine/dispatcher.rs`
- **Primary Tools**: Solve, Compute, Analyze, Simulate, ML, Chaos, Units, Validate
- **Legacy Tools** (route to primary): Differentiate, Integrate, Transform, FieldTheory, Sample, Optimize
- Clean trait-based interface with JSON-serializable types
- Request/response handled by `ToolDispatcher` in `src/engine/dispatcher.rs`

**Tool Routing**:
- `differentiate` → `compute` with `operation: {differentiate: ...}`
- `integrate` → `compute` with `operation: {integrate: ...}`
- `transform` → `compute` with `operation: {transform: ...}`
- `fieldtheory` → `compute` with `operation: {field: ...}`
- `sample` → `compute` with `operation: {sample: ...}`
- `optimize` → `solve` with `equation_type: {optimize: ...}`

**8 Tool Modules**: Each tool in its own folder with unified implementation + domain submodules
- `src/solve/` - UnifiedSolver + equations, optimization, physics solvers, specialized (game_theory, linear_programming)
- `src/compute/` - UnifiedComputer + 27 domain submodules (matrix, tensor, physics, chemistry, biology, etc.)
- `src/analyze/` - UnifiedAnalyzer + symbolic CAS, series, stability analysis
- `src/simulate/` - UnifiedSimulator + ODE solvers, stochastic, fluids, finance
- `src/ml/` - UnifiedML + clustering, regression, neural networks, dim reduction
- `src/chaos/` - UnifiedChaos + fractals, attractors, Lyapunov, bifurcation
- `src/units/` - UnifiedUnits + dimensional analysis
- `src/validate/` - UnifiedValidator + equation validation

### Module Organization

```
src/
├── engine/              # 8-tool public API
│   ├── traits.rs       # Tool traits (Solve, Compute, etc.)
│   ├── types.rs        # Input/Output types for each tool
│   ├── dispatcher.rs   # Routes requests to implementations
│   └── equations.rs    # Equation type definitions
│
├── solve/               # SOLVE TOOL
│   ├── mod.rs          # UnifiedSolver
│   ├── equations.rs    # Root finding, linear systems
│   ├── differential.rs # ODE/PDE solving
│   ├── optimization/   # Curve fitting, minimization
│   ├── physics/        # Einstein, EM, fluid equation solvers
│   └── specialized/    # Game theory, linear programming
│
├── compute/             # COMPUTE TOOL (largest - 27 submodules)
│   ├── mod.rs          # UnifiedComputer
│   ├── matrix/         # Linear algebra, decompositions
│   ├── tensor/         # Tensor calculus, Christoffel symbols
│   ├── calculus/       # Differentiation, integration
│   ├── transforms/     # FFT, wavelets, filters
│   ├── special_functions/ # Bessel, gamma, elliptic
│   ├── physics/        # 14 physics submodules (quantum, relativity, etc.)
│   ├── sampling/       # Monte Carlo, MCMC, statistics
│   ├── geometry/       # Computational geometry
│   ├── graph/          # Graph algorithms
│   ├── information/    # Information theory
│   ├── biology/        # Biology formulas
│   ├── chemistry/      # Chemistry formulas
│   ├── thermodynamics/ # Thermodynamics
│   └── ...             # More domain modules
│
├── analyze/             # ANALYZE TOOL
│   ├── mod.rs          # UnifiedAnalyzer
│   ├── symbolic/       # Custom CAS (simplify, expand, factor)
│   ├── series/         # Taylor, Laurent series
│   ├── stability/      # Stability analysis
│   └── validation/     # Expression validation
│
├── simulate/            # SIMULATE TOOL
│   ├── mod.rs          # UnifiedSimulator
│   ├── ode/            # Euler, Runge-Kutta solvers
│   ├── stochastic/     # Brownian motion, Monte Carlo
│   ├── fluids/         # Navier-Stokes, flow analysis
│   └── finance/        # Black-Scholes, Heston
│
├── ml/                  # ML TOOL
│   ├── mod.rs          # UnifiedML
│   ├── clustering/     # K-means, DBSCAN
│   ├── regression/     # Linear, logistic, ridge
│   ├── neural_network/ # Dense layers, training
│   ├── classification/ # Decision trees, SVM
│   └── dim_reduction/  # PCA, t-SNE
│
├── chaos/               # CHAOS TOOL
│   ├── mod.rs          # UnifiedChaos
│   ├── fractals/       # Mandelbrot, Julia
│   ├── attractors/     # Lorenz, Rossler
│   ├── lyapunov/       # Lyapunov exponents
│   ├── bifurcation/    # Bifurcation diagrams
│   └── dimension/      # Box-counting, Kaplan-Yorke
│
├── units/               # UNITS TOOL
│   ├── mod.rs          # UnifiedUnits
│   └── dimensional_analysis.rs
│
├── validate/            # VALIDATE TOOL
│   ├── mod.rs          # UnifiedValidator
│   └── equation.rs     # Equation validation
│
├── wasm.rs             # WASM bindings (--features wasm)
├── main.rs             # CLI binary entry point
└── lib.rs              # Library root + create_default_dispatcher()
```

### Data Flow

```
JSON Request → ToolDispatcher → Unified Implementation → Domain Modules → Result
```

**Example:**
1. User sends JSON: `{"tool": "solve", "input": {"equations": ["x^2 - 4 = 0"]}}`
2. Parsed into `ToolRequest::Solve(SolveInput { ... })`
3. Dispatcher calls `UnifiedSolver::solve()` implementation
4. `UnifiedSolver` uses domain modules (symbolic_cas, linear_algebra, etc.)
5. Returns `ToolResponse::Solve(SolveOutput { ... })`
6. Serialized back to JSON

## Critical Rules

### ALWAYS Fix Implementation Bugs, NEVER Weaken Tests (CRITICAL)

**THIS PROJECT REQUIRES PERFECTLY ACCURATE AND DETERMINISTIC BEHAVIOR.**

When a test fails, you MUST:
1. ✅ **Investigate the implementation** - Assume the test is correct and the implementation is wrong
2. ✅ **Fix the root cause** - Correct the algorithm, formula, or logic error in the implementation
3. ✅ **Verify the physics/math** - Ensure the fix matches the correct scientific/mathematical behavior

You MUST NEVER:
1. ❌ **Weaken test assertions** - Don't change expected values to match buggy output
2. ❌ **Add tolerance when not needed** - Don't add `approx_eq` if exact equality should work
3. ❌ **Comment out failing tests** - Don't hide failures by disabling tests
4. ❌ **Make assumptions** - Don't assume the test is wrong without investigating

**Example of the WRONG approach (NEVER DO THIS):**
```rust
// ❌ WRONG - Weakening test to match buggy implementation
assert!(r_ergo >= r_h);  // Changed from > to >= to make test pass
```

**Example of the CORRECT approach:**
```rust
// ✅ CORRECT - Fix the implementation bug
// In implementation file:
let a_geom = self.spin * m_geom;  // Fixed: was using wrong formula
```

**Real incident from 2024-10-26:**
- Black hole ergosphere tests were failing
- Initial response: Weakened test assertions from `>` to `>=`
- **This was WRONG** - It hid a fundamental bug in spin parameter interpretation
- Correct response: Fixed implementation to use dimensionless spin (0-1) instead of mass units
- Result: Tests pass with correct physics, no compromises made

**Remember:** This is a computational engine competing with Mathematica. Wrong results are unacceptable. A failing test is a gift - it shows you where the implementation is incorrect. Fix the implementation, never the test.

### stdout/stderr Policy (CRITICAL for MCP Compatibility)

**NEVER use `println!` in library code.** The engine runs as an MCP server where:
- **stdout** = JSON responses ONLY
- **stderr** = All logs, debug output, errors

```rust
// ❌ WRONG - Corrupts JSON output
println!("Processing...");

// ✅ CORRECT - Use stderr for logs
eprintln!("Processing...");

// ✅ BETTER - Use logging crate (if available)
log::info!("Processing...");
```

**Only `src/main.rs` may write to stdout:**
- CLI info commands (`info`, `list-ops`)
- JSON responses (`json`, `stdin` commands)

**Library modules must use stderr for all debug output.**

See `STDOUT_STDERR_POLICY.md` for complete details.

### Rust Edition

This project uses **Rust 2024 edition** (`edition = "2024"` in Cargo.toml). Be aware of 2024-specific features and syntax.

### Feature Flags

```toml
default = []                  # Minimal build
quantum = ["tokio", "gpu"]   # Quantum physics with GPU acceleration
gpu-acceleration = ["wgpu"]   # GPU compute support
wasm = ["wasm-bindgen"]      # WebAssembly support
```

## Key Dependencies

### Core Math Libraries
- `nalgebra` (0.33) - Linear algebra and matrix operations
- `ndarray` (0.15) - N-dimensional arrays
- `num-complex` (0.4) - Complex number support
- `special` (0.10) - Special functions (Bessel, gamma, etc.)
- `quadrature` (0.1) - Numerical integration

### Signal Processing
- `rustfft` (6.0) - Fast Fourier Transform
- `apodize` (1.0) - Window functions

### Random/Stochastic
- `rand` (0.8) - Random number generation
- `rand_distr` (0.4) - Probability distributions

### Cryptography
- `num-bigint` (0.4) - Arbitrary precision integers
- `sha2` (0.10), `sha3` (0.10) - Hashing algorithms

### Parallelism
- `rayon` (1.7) - Data parallelism
- `tokio` (1.0) - Async runtime (optional, for quantum features)

### WebAssembly
- `wasm-bindgen` (0.2) - JavaScript interop
- `js-sys`, `web-sys` - Browser APIs

## Development Workflow

### Adding a New Operation

**Step 1: Identify which tool it belongs to**
- Is it solving equations/optimization? → `Solve`
- Computing derivatives/integrals/transforms/matrices? → `Compute`
- Simplifying/analyzing expressions? → `Analyze`
- Time evolution/simulation? → `Simulate`
- Machine learning? → `ML`
- Fractals/chaos theory? → `Chaos`
- Unit conversion? → `Units`
- Validation? → `Validate`

**Step 2: Add to the tool's input type**
1. Update the appropriate input enum in `src/engine/types.rs`
2. Add parameters needed for the operation

**Step 3: Implement in the tool module**
1. Edit the corresponding file in `src/{tool}/mod.rs` (e.g., `src/solve/mod.rs` for Solve)
2. Add match arm to handle new operation
3. Call existing domain submodule or create new one

**Step 4: Add domain submodule implementation (if needed)**
1. Create/update submodule in the tool folder (e.g., `src/compute/physics/`)
2. Implement the core algorithm
3. Add tests

### Testing Strategy

```bash
# Run all tests (2,067 total)
cargo test

# Run integration tests only
cargo test --test all_integration_tests

# Run comprehensive integration test suites
cargo test --test '*_comprehensive_tests'

# Run specific comprehensive test suites
cargo test --test differentiate_comprehensive_tests
cargo test --test integrate_comprehensive_tests
cargo test --test analyze_comprehensive_tests
cargo test --test simulate_comprehensive_tests
cargo test --test compute_scientific_comprehensive_tests
cargo test --test transform_comprehensive_tests
cargo test --test fieldtheory_comprehensive_tests
cargo test --test sample_comprehensive_tests
cargo test --test optimize_comprehensive_tests

# Run unit tests only
cargo test --lib

# Test specific module
cargo test tensor_calculus

# Test with all features
cargo test --all-features

# WASM tests (headless browser)
wasm-pack test --headless --firefox
```

### Test Organization

Tests mirror the 8-tool source structure:

```
tests/
├── unit/                    # Unit tests (960 tests)
│   ├── analyze/            # Analyze tool tests
│   │   └── symbolic/
│   ├── chaos/              # Chaos tool tests
│   │   ├── attractors/
│   │   ├── fractals/
│   │   └── lyapunov/
│   ├── compute/            # Compute tool tests
│   │   ├── physics/        # 14 physics subdirectories
│   │   ├── special_functions/
│   │   ├── tensor/
│   │   └── ...
│   ├── ml/                 # ML tool tests
│   ├── simulate/           # Simulate tool tests
│   ├── solve/              # Solve tool tests
│   ├── units/              # Units tool tests
│   └── validate/           # Validate tool tests
│
├── integration/             # Integration tests (1,089 tests)
│   ├── engine/             # 8-tool dispatcher tests
│   ├── compute/            # Compute integration tests
│   ├── solve/              # Solve integration tests
│   ├── simulate/           # Simulate integration tests
│   ├── tools/              # Legacy tool name tests
│   └── coverage/           # Coverage verification tests
```

### WebAssembly Development

WASM builds generate multiple targets for different environments:
- `pkg/bundler/` - Webpack/Rollup/Vite (client-side)
- `pkg/nodejs/` - Node.js (server-side)
- `pkg/web/` - ES modules for browser
- `pkg/no-modules/` - Classic `<script>` tags

**When modifying WASM bindings:**
1. Edit `src/wasm/mod.rs`
2. Rebuild with `./build-wasm.sh`
3. Test in browser example: `cd examples/wasm && npm run dev`
4. Test in Node.js: `cd examples/wasm && npm run node`

## Important Documentation Files

- **ARCHITECTURE.md** - Detailed 8-tool architecture explanation
- **API.md** - Complete JSON API reference
- **WASM.md** - WebAssembly build and usage guide
- **STDOUT_STDERR_POLICY.md** - Critical MCP server compatibility rules
- **OPERATIONS_COMPLETE_LIST.md** - Complete inventory of all 400+ operations with test coverage details

## Common Patterns

### Using the 8-Tool API (Recommended)

```rust
use computational_engine::{ToolRequest, SolveInput, create_default_dispatcher};

let dispatcher = create_default_dispatcher();

let request = ToolRequest::Solve(SolveInput {
    equations: vec!["x^2 + 2*x - 8 = 0".to_string()],
    variables: None,
    initial_guess: None,
    domain: None,
    method: None,
    ..Default::default()
});

let response = dispatcher.dispatch(request).unwrap();
```

### Custom Computer Algebra System

This engine includes a custom CAS (not using external libraries like SymPy):
- Located in `src/analyze/symbolic/`
- Components: `expr.rs`, `parser.rs`, `simplify.rs`, `differentiate.rs`
- Supports symbolic differentiation, simplification, integration
- Matrix symbolic operations in `symbolic_matrix.rs`

## Performance Considerations

### Release Builds
Always use `--release` for production:
```bash
cargo build --release  # ~10-100x faster than debug builds
```

Release profile configured for maximum optimization:
```toml
[profile.release]
lto = true              # Link-time optimization
codegen-units = 1       # Single codegen unit for better optimization
opt-level = 3           # Maximum optimization
```

### Parallelism
Many operations use `rayon` for automatic parallelization. Large computations benefit from multi-core systems.

### GPU Acceleration
Enable with `--features quantum,gpu-acceleration` for quantum physics simulations using `wgpu`.

## Common Pitfalls

1. **Using `println!` in library code** - Breaks MCP server JSON output (use `eprintln!` instead)
2. **Modifying old `main()` functions** - These are legacy and should not be called
3. **Forgetting `--release` flag** - Debug builds are extremely slow for math operations
4. **Not enabling features** - Quantum operations require `--features quantum`
5. **WASM target mismatch** - Use correct target (bundler vs nodejs vs web) for your environment

## Project Goals

This engine aims to compete with Mathematica/Wolfram Alpha by providing:
- Comprehensive mathematical operations (**400+ operations** across 8 tools)
- Extensive test coverage with **2,067 tests passing** (1,089 integration + 960 unit + 18 doc tests)
- High-performance Rust implementation
- Multiple interfaces (Rust API, JSON/MCP, WebAssembly)
- Self-contained (custom CAS, no Python dependencies)
- Free and open-source

## Recent Updates

### January 2025
- ✅ **8-Tool Architecture Consolidation**: Unified from 10 legacy tools to 8 primary tools
  - Primary: Solve, Compute, Analyze, Simulate, ML, Chaos, Units, Validate
  - Legacy names (differentiate, integrate, transform, fieldtheory, sample, optimize) route to primary tools
- ✅ **Test Reorganization**: Tests now mirror source structure
  - `tests/unit/` organized by tool
  - `tests/integration/` organized by tool
  - Total: 2,067 tests (up from 1,685)

### October 2024
- ✅ **Major Test Coverage Improvement**: Achieved ~80% production code coverage
- ✅ **Black Hole Physics Fix**: Corrected Kerr black hole implementation
  - Changed spin parameter from mass units to dimensionless (0-1)
  - Fixed event horizon, ergosphere, and ISCO calculations
- ✅ New signal processing operations (Haar, Daubechies, Mexican hat wavelets)
- ✅ New curve fitting models (exponential, logarithmic, power law, rational)

## Version Information

Current version: 0.1.0 (see Cargo.toml)

Project structure is stable. The 8-tool unified API is the recommended interface for new code.
