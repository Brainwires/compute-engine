# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

A unified computational engine for mathematical and scientific computing in Rust, providing **194+ mathematical operations** accessible through a **consolidated 4-tool API**:

- **Solve**: Equations, systems, optimization, root finding
- **Compute**: Calculus, transforms, field theory, sampling, matrix ops
- **Analyze**: Series, limits, stability analysis, simplification
- **Simulate**: Time evolution, stochastic processes, fluid dynamics

Legacy tool names (differentiate, integrate, transform, fieldtheory, sample, optimize) are maintained for backward compatibility and route internally to the 4 primary tools.

Supports native Rust, WebAssembly (browser/Node.js), and MCP server interfaces.

**Comprehensive Test Coverage**: 2939 tests passing (1222 comprehensive integration tests + 1717 unit tests) with 100% pass rate.

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

### Consolidated 4-Tool Architecture

**Public API Layer**: 4 primary tools + 6 legacy tools for backward compatibility
- Defined in `src/engine/traits.rs` and `src/engine/dispatcher.rs`
- **Primary Tools**: Solve, Compute, Analyze, Simulate
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

**Implementation Layer**: Unified tool implementations
- Located in `src/implementations/`
- Each tool has a `Unified*` implementation (e.g., `UnifiedSolver`, `UnifiedComputer`)
- These wire the 4-tool API to the 245 domain-specific functions
- Legacy tool implementations maintained for backward compatibility

**Domain Modules**: Implementation details (not a separate API)
- Located in `src/mathematics/`, `src/physics/`, `src/tools/`, `src/specialized/`
- Contains the actual math/physics algorithms
- Called by the Unified implementations
- Not meant to be used directly by end users

### Module Organization

```
src/
├── engine/              # 10-tool public API
│   ├── traits.rs       # Tool traits (Solve, Differentiate, etc.)
│   ├── types.rs        # Input/Output types for each tool
│   ├── dispatcher.rs   # Routes requests to implementations
│   └── equations.rs    # Equation type definitions
│
├── implementations/     # Unified tool implementations
│   ├── solver.rs       # UnifiedSolver
│   ├── differentiator.rs
│   ├── integrator.rs
│   ├── analyzer.rs
│   ├── simulator.rs
│   ├── computer.rs
│   ├── transformer.rs
│   ├── field_solver.rs
│   ├── sampler.rs
│   ├── optimizer.rs
│   └── mod.rs          # Dispatcher factory
│
├── mathematics/         # Math implementation modules
│   ├── tensor_calculus/
│   ├── advanced_calculus/
│   ├── linear_algebra/
│   ├── symbolic_regression/
│   ├── special_functions/
│   └── symbolic_cas/   # Custom computer algebra system
│
├── physics/            # Physics implementation modules
│   ├── fluid_dynamics/
│   ├── quantum_physics/
│   └── electromagnetism/
│
├── tools/              # Engineering implementation modules
│   ├── signal_processing/
│   ├── dimensional_analysis/
│   ├── equation_validation/
│   ├── computational_geometry/
│   └── numerical_methods/
│
├── specialized/        # Specialized implementation modules
│   ├── stochastic_processes/
│   ├── cryptographic_mathematics/
│   ├── statistics/
│   ├── optimization/
│   ├── graph_theory/
│   └── information_theory/
│
├── chemistry/          # Chemistry formulas (2025 expansion)
├── biology/            # Biology formulas (2025 expansion)
├── thermodynamics/     # Thermodynamics (2025 expansion)
├── datetime/           # Date/time calculations (2025 expansion)
├── optics/             # Optics (2025 expansion)
├── geophysics/         # Geophysics (2025 expansion)
├── engineering/        # Engineering (2025 expansion)
│
├── api/                # Backwards-compatible JSON API (old format)
├── wasm/               # WASM bindings (--features wasm)
├── main.rs             # CLI binary entry point
└── lib.rs              # Library root
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

**Real incident from 2025-10-26:**
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
- Is it solving equations? → `Solve`
- Computing derivatives? → `Differentiate`
- Computing integrals? → `Integrate`
- Simplifying/analyzing? → `Analyze`
- Time evolution/simulation? → `Simulate`
- Matrix/tensor operations? → `Compute`
- Fourier/Laplace/filters? → `Transform`
- EM/gravity/quantum fields? → `FieldTheory`
- Statistical sampling? → `Sample`
- Optimization/fitting? → `Optimize`

**Step 2: Add to the tool's input type**
1. Update the appropriate input enum in `src/engine/types.rs`
2. Add parameters needed for the operation

**Step 3: Implement in the Unified implementation**
1. Edit the corresponding file in `src/implementations/` (e.g., `solver.rs` for Solve)
2. Add match arm to handle new operation
3. Call existing domain module or create new one

**Step 4: Add domain module implementation (if needed)**
1. Create/update module in `src/mathematics/`, `src/physics/`, etc.
2. Implement the core algorithm
3. Add tests

**Step 5: Update old API (for backwards compatibility)**
1. Add handler in `src/api/handlers/` if needed
2. Register in `src/api/mod.rs`

### Testing Strategy

```bash
# Run all tests (1685 total)
cargo test

# Run comprehensive integration tests (202 tests)
cargo test --test '*_comprehensive_tests'

# Run specific comprehensive test suites
cargo test --test differentiate_comprehensive_tests    # 21 tests
cargo test --test integrate_comprehensive_tests        # 19 tests
cargo test --test analyze_comprehensive_tests          # 25 tests
cargo test --test simulate_comprehensive_tests         # 21 tests
cargo test --test compute_scientific_comprehensive_tests # 29 tests
cargo test --test transform_comprehensive_tests        # 22 tests
cargo test --test fieldtheory_comprehensive_tests      # 15 tests
cargo test --test sample_comprehensive_tests           # 23 tests
cargo test --test optimize_comprehensive_tests         # 27 tests

# Run unit tests (1483 tests)
cargo test --lib

# Test specific module
cargo test tensor_calculus

# Test with all features
cargo test --all-features

# WASM tests (headless browser)
wasm-pack test --headless --firefox
```

**Comprehensive Test Coverage:**
- ✅ **differentiate**: Symbolic, Numeric, Partial, Chain rule, Implicit, Parametric (21 tests)
- ✅ **integrate**: Definite, Indefinite, Numeric, Improper, Multiple, Contour (19 tests)
- ✅ **analyze**: Series, Limits, Roots, Extrema, Stability, Fourier, Wavelets (25 tests)
- ✅ **simulate**: ODEs, PDEs, Stochastic, Monte Carlo, Cellular automata (21 tests)
- ✅ **compute_scientific**: Chemistry, Biology, Thermodynamics, Optics, Geophysics, Engineering (29 tests)
- ✅ **transform**: FFT, Fourier, Laplace, Wavelets, Filters, Windows, Conformal (22 tests)
- ✅ **fieldtheory**: EM Fields, Green's Functions, Quantum Fields (15 tests)
- ✅ **sample**: Monte Carlo, MCMC, Statistical Methods, Signal Analysis (23 tests)
- ✅ **optimize**: Curve Fitting, Interpolation, Minimization, Symbolic Regression (27 tests)

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

- **ARCHITECTURE.md** - Detailed 10-tool architecture explanation
- **API.md** - Complete JSON API reference
- **WASM.md** - WebAssembly build and usage guide
- **STDOUT_STDERR_POLICY.md** - Critical MCP server compatibility rules
- **OPERATIONS_COMPLETE_LIST.md** - Complete inventory of all 194+ operations with test coverage details
- **MATHEMATICA_COMPETITION_PROGRESS.md** - Feature parity tracking

## Common Patterns

### Using the 10-Tool API (Recommended)

```rust
use computational_engine::engine::{ToolRequest, SolveInput};
use computational_engine::implementations::create_default_dispatcher;

let dispatcher = create_default_dispatcher();

let request = ToolRequest::Solve(SolveInput {
    equations: vec!["x^2 + 2*x - 8 = 0".to_string()],
    variables: None,
    initial_guess: None,
    domain: None,
    method: None,
});

let response = dispatcher.dispatch(request).unwrap();
```

### Using the Old JSON API Format (For Backwards Compatibility)

The old `ComputationRequest`/`ComputationResponse` format from `src/api/` is still supported but deprecated:

```rust
use computational_engine::api::{ComputationRequest, process_request};
use serde_json::json;

let request = ComputationRequest {
    module: "advanced_calculus".to_string(),
    operation: "riemann_zeta".to_string(),
    parameters: json!({"s": {"real": 2.0}}),
};

let response = process_request(&request);
```

**Note**: Migrate to the 10-tool API for new code. This old format may be removed in future versions.

### Custom Computer Algebra System

This engine includes a custom CAS (not using external libraries like SymPy):
- Located in `src/mathematics/symbolic_cas/`
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
- Comprehensive mathematical operations (**194+ operations** across 10 tools, see OPERATIONS_COMPLETE_LIST.md)
- Extensive test coverage with **1685 tests passing** (202 comprehensive integration + 1483 unit tests)
- High-performance Rust implementation
- Multiple interfaces (Rust API, JSON/MCP, WebAssembly)
- Self-contained (custom CAS, no Python dependencies)
- Free and open-source

## Recent Implementations (2025-10-14)

### New Signal Processing Operations
- ✅ Haar wavelet transform
- ✅ Daubechies 4 wavelet transform
- ✅ Mexican hat (Ricker) wavelet transform
- ✅ Kaiser window function with Modified Bessel I₀

### New Curve Fitting Models
- ✅ Exponential fitting: y = a·exp(bx)
- ✅ Logarithmic fitting: y = a + b·ln(x)
- ✅ Power law fitting: y = a·x^b
- ✅ Rational function fitting: y = (a+bx)/(1+cx)

### Bug Fixes
- ✅ Fixed 13 scientific computing parameter name mismatches
- ✅ Fixed type ambiguity in mexican_hat_wavelet implementation

## Version Information

Current version: 0.1.0 (see Cargo.toml)

Project structure is stable but API is still evolving. The 10-tool unified API is the recommended interface for new code.
