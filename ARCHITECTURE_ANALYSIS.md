# Architecture Analysis for Refactoring

**Date:** 2026-01-11
**Purpose:** Document the current architectural mess before major refactoring

---

## Executive Summary

The codebase has **THREE SEPARATE API LAYERS** that evolved over time, resulting in:
- Duplicated math implementations
- Inconsistent wiring patterns
- Orphaned/unreachable code
- Confusing entry points for users

**Total Domain Code:** ~40,000 lines across 100+ files

---

## Current Architecture (The Mess)

### Layer 1: Unified Tool API (`src/engine/` + `src/implementations/`)

**Status:** RECOMMENDED, but incomplete

**Location:**
- `src/engine/` - Traits, types, dispatcher
- `src/implementations/` - Unified implementations

**10 Tools:**
| Tool | Implementation | Purpose |
|------|----------------|---------|
| Solve | `UnifiedSolver` | Equations, systems, root finding |
| Differentiate | `UnifiedDifferentiator` | Derivatives, gradients, Jacobians |
| Integrate | `UnifiedIntegrator` | Definite/indefinite integrals |
| Analyze | `UnifiedAnalyzer` | Simplify, expand, factor, transform |
| Simulate | `UnifiedSimulator` | ODEs, PDEs, physics simulations |
| Compute | `UnifiedComputer` | Matrix ops, special functions |
| Transform | `UnifiedTransformer` | FFT, wavelets, Laplace |
| FieldTheory | `UnifiedFieldSolver` | EM fields, quantum fields |
| Sample | `UnifiedSampler` | Monte Carlo, MCMC |
| Optimize | `UnifiedOptimizer` | Gradient descent, curve fitting |

**How it works:**
```
JSON Request → ToolDispatcher → Unified* Implementation → Domain Module → Result
```

**Problems:**
- Not all domain functions are wired
- Some Unified* implementations are incomplete
- Legacy tool names (differentiate, integrate, etc.) route to primary tools

---

### Layer 2: Module-Based JSON API (`src/api/`)

**Status:** MIXED - some real, some stubs

**Location:** `src/api/handlers/` (24+ handler files)

**File Analysis:**

| File | Lines | Status | Notes |
|------|-------|--------|-------|
| `advanced_calculus.rs` | 220 | **REAL** | Inline polynomial derivatives/integrals/limits |
| `linear_algebra.rs` | 218 | **REAL** | Inline matrix operations |
| `cryptographic_mathematics.rs` | 173 | **REAL** | Calls domain module |
| `electromagnetism.rs` | 169 | **REAL** | Calls domain module |
| `information_theory.rs` | 151 | **REAL** | Calls domain module |
| `numerical_methods.rs` | 151 | **REAL** | Calls domain module |
| `special_functions.rs` | 133 | **REAL** | Calls domain module |
| `statistics.rs` | 130 | **REAL** | Calls domain module |
| `chemistry.rs` | 128 | **REAL** | Calls domain module |
| `function_approximator.rs` | 126 | **REAL** | Inline implementations |
| `computational_geometry.rs` | 122 | **REAL** | Calls domain module |
| `graph_theory.rs` | 115 | **REAL** | Calls domain module |
| `stochastic_processes.rs` | 112 | **REAL** | Calls domain module |
| `symbolic_regression.rs` | 106 | **PARTIAL** | Some implemented |
| `signal_processing.rs` | 93 | **PARTIAL** | Some calls domain |
| `equation_validation.rs` | 82 | **PARTIAL** | Limited |
| `dimensional_analysis.rs` | 79 | **PARTIAL** | Limited |
| `optimization.rs` | 64 | **PARTIAL** | Basic |
| `engineering.rs` | 46 | **STUB** | Minimal |
| `geophysics.rs` | 46 | **STUB** | Minimal |
| `thermodynamics.rs` | 45 | **STUB** | Minimal |
| `datetime.rs` | 45 | **STUB** | Minimal |
| `biology.rs` | 44 | **STUB** | Minimal |
| `optics.rs` | 44 | **STUB** | Minimal |
| `tensor_calculus.rs` | 29 | **MOCK** | Returns fake data |
| `quantum_physics.rs` | 29 | **MOCK** | Returns fake data |
| `fluid_dynamics.rs` | 29 | **MOCK** | Returns fake data |

**Entry point:**
```rust
pub fn process_request(request: &ComputationRequest) -> ComputationResponse
```

**Problems:**
- **DUPLICATED MATH**: `advanced_calculus.rs` has inline implementations instead of calling domain modules
- **INCONSISTENT**: Some call domain modules, some have inline code, some are stubs
- **MOCKS**: Several handlers return fake data instead of real computations

---

### Layer 3: Legacy Module-Based API (`src/api_legacy/`)

**Status:** DEPRECATED but has real wiring

**Location:** `src/api_legacy/handlers/` (17 handler files)

**File Analysis:**

| File | Lines | Status | Notes |
|------|-------|--------|-------|
| `signal.rs` | 179 | **REAL** | Calls `crate::signal_processing::*` |
| `information_theory.rs` | 153 | **REAL** | Calls domain module |
| `linear_algebra.rs` | 138 | **REAL** | Calls domain module |
| `statistics.rs` | 135 | **REAL** | Calls domain module |
| `graph_theory.rs` | 117 | **REAL** | Calls domain module |
| `geometry.rs` | 117 | **REAL** | Calls domain module |
| `stochastic.rs` | 116 | **REAL** | Calls domain module |
| `calculus.rs` | 48 | **REAL** | Calls `crate::advanced_calculus::*` |
| `optimization.rs` | 45 | **REAL** | Calls domain module |
| `tensor.rs` | 44 | **REAL** | Calls domain module |
| `quantum.rs` | 12 | **STUB** | Minimal |
| `regression.rs` | 12 | **STUB** | Minimal |
| `validator.rs` | 12 | **STUB** | Minimal |
| `dimensional.rs` | 12 | **STUB** | Minimal |
| `crypto.rs` | 12 | **STUB** | Minimal |
| `fluid.rs` | 11 | **STUB** | Minimal |

**Key Discovery:**
The `calculus.rs` handler shows the **correct pattern**:
```rust
use crate::advanced_calculus::*;

match request.operation.as_str() {
    "fractional_derivative" => handle_op!(handle_fractional_derivative(&request.parameters)),
    "riemann_zeta" => handle_op!(handle_riemann_zeta(&request.parameters)),
    "elliptic_integral" => handle_op!(handle_elliptic_integral(&request.parameters)),
    // ... 18 more operations
}
```

This calls into domain module functions that ARE implemented but NOT wired to the unified API!

---

## Domain Modules (The Real Math)

### src/mathematics/ (~15,000 lines)

```
src/mathematics/
├── calculus/
│   ├── mod.rs
│   ├── fractional.rs        # Fractional calculus (handle_fractional_*)
│   ├── special_functions.rs # Bessel, Gamma, etc.
│   ├── stochastic.rs        # Ito, Stratonovich integrals
│   ├── symbolic_integration.rs
│   └── variational.rs       # Euler-Lagrange, variational calc
│
├── chaos/
│   ├── mod.rs
│   ├── attractors.rs        # Lorenz, Rossler attractors
│   ├── fractals.rs          # Mandelbrot, Julia, Koch
│   └── lyapunov.rs          # Lyapunov exponents
│
├── linear_algebra/
│   └── mod.rs               # SVD, PCA, eigenvalues
│
├── special_functions/
│   ├── mod.rs
│   ├── airy.rs              # Airy functions
│   ├── bessel.rs            # Bessel functions
│   ├── elliptic.rs          # Elliptic integrals
│   ├── error.rs             # Error functions
│   ├── gamma.rs             # Gamma, digamma, beta
│   └── polynomials.rs       # Legendre, Hermite, Laguerre
│
├── symbolic_cas/            # Custom CAS (NOT SymPy!)
│   ├── mod.rs
│   ├── expr.rs              # Expression types
│   ├── parser.rs            # Math parser
│   ├── simplify.rs          # Simplification rules
│   ├── differentiate.rs     # Symbolic differentiation
│   ├── integrate.rs         # Symbolic integration
│   ├── symbolic_matrix.rs   # Symbolic matrices
│   ├── symbolic_tensor.rs   # Symbolic tensors
│   ├── quantum.rs           # Quantum mechanics
│   ├── quantum_advanced.rs  # Advanced QM
│   ├── christoffel.rs       # Christoffel symbols
│   ├── metric_tensors.rs    # Common metrics
│   └── fluid_dynamics.rs    # Fluid CAS
│
├── symbolic_regression/
│   ├── mod.rs
│   └── lib.rs               # Equation discovery
│
├── tensor_calculus/
│   ├── mod.rs
│   ├── tensor.rs            # Christoffel, Riemann, Ricci
│   ├── einstein.rs          # Einstein field equations
│   ├── symbolic.rs          # Symbolic tensor ops
│   └── quantum_tensors.rs   # Quantum tensor fields
│
└── numerical.rs             # Numerical methods
```

### src/physics/ (~10,000 lines)

```
src/physics/
├── mod.rs
│
├── black_holes/
│   ├── mod.rs
│   ├── orbits.rs            # Schwarzschild, ISCO, ergosphere
│   ├── penrose.rs           # Penrose process
│   └── hawking.rs           # Hawking radiation
│
├── cosmology/
│   ├── mod.rs
│   ├── friedmann.rs         # Hubble, distances, evolution
│   ├── cmb.rs               # CMB physics
│   └── dark_energy.rs       # Dark energy models
│
├── electromagnetism/
│   └── mod.rs               # Maxwell, EM waves, antennas
│
├── fluid_dynamics/
│   ├── mod.rs
│   ├── navier_stokes.rs     # 2D Navier-Stokes
│   └── quantum_fluids.rs    # Quantum Navier-Stokes
│
├── quantum_mechanics/
│   ├── mod.rs               # Schrodinger, hydrogen atom
│   ├── bohm_potential.rs    # Bohmian mechanics
│   └── decoherence.rs       # Decoherence scale
│
├── quantum_physics/
│   └── mod.rs               # General quantum physics
│
└── wormholes/
    ├── mod.rs
    ├── metric.rs            # Morris-Thorne metric
    ├── traversal.rs         # Traversability analysis
    └── energy.rs            # Exotic matter requirements
```

### src/tools/ (~8,000 lines)

```
src/tools/
├── mod.rs
│
├── computational_geometry/
│   ├── mod.rs
│   ├── advanced.rs          # Convex hull, Delaunay, Voronoi
│   └── spatial_3d.rs        # 3D geometry
│
├── dimensional_analysis/
│   └── lib.rs               # Unit analysis, conversions
│
├── equation_validation/
│   └── lib.rs               # Equation verification
│
├── numerical_methods/
│   └── mod.rs               # ODE, PDE, root finding, interpolation
│
└── signal_processing/
    ├── mod.rs
    ├── lib.rs               # FFT, filters, spectrograms
    ├── additional.rs        # Peaks, Fourier series
    └── wavelets.rs          # Haar, Daubechies, Mexican hat
```

### src/specialized/ (~12,000 lines)

```
src/specialized/
├── mod.rs
│
├── chemistry/
│   └── mod.rs               # Chemical equations (old, replaced)
│
├── control_theory/
│   ├── mod.rs
│   ├── pid_controller.rs    # PID tuning
│   ├── analysis.rs          # Controllability, observability
│   ├── lqr.rs               # LQR optimal control
│   └── state_space.rs       # State-space models
│
├── cryptographic_mathematics/
│   └── lib.rs               # RSA, ECC, number theory
│
├── game_theory/
│   ├── mod.rs
│   ├── normal_form.rs       # Nash equilibrium
│   ├── extensive_form.rs    # Game trees
│   └── evolutionary.rs      # Evolutionary game theory
│
├── graph_theory/
│   └── mod.rs               # Shortest path, MST, components
│
├── information_theory/
│   └── mod.rs               # Entropy, mutual information
│
├── linear_programming/
│   ├── mod.rs
│   ├── simplex.rs           # Simplex algorithm
│   └── dual.rs              # Dual problems, sensitivity
│
├── machine_learning/
│   ├── mod.rs
│   ├── clustering.rs        # K-means, DBSCAN
│   ├── dimensionality_reduction.rs  # PCA, t-SNE
│   ├── neural_network.rs    # Basic neural networks
│   ├── optimization.rs      # SGD, Adam
│   └── regression.rs        # Linear, polynomial
│
├── optimization/
│   └── mod.rs               # Gradient descent, Nelder-Mead
│
├── statistics/
│   └── mod.rs               # Stats, Monte Carlo, MCMC
│
└── stochastic_processes/
    ├── mod.rs
    └── lib.rs               # Brownian motion, Markov chains
```

### New Scientific Modules (2025 Expansion)

```
src/
├── biology/mod.rs           # Hardy-Weinberg, population genetics
├── chemistry/mod.rs         # Stoichiometry, thermodynamics
├── datetime/mod.rs          # Date/time calculations
├── electrical/mod.rs        # Circuit analysis, NEC
├── engineering/mod.rs       # Structural, fluid, thermal
├── geophysics/mod.rs        # Seismic, gravity, magnetic
├── materials_science/mod.rs # Crystal structures, band theory
├── optics/mod.rs            # Lens equations, diffraction
└── thermodynamics/mod.rs    # Heat transfer, cycles
```

---

## The Problems (Detailed)

### Problem 1: Duplicated Math Implementations

**Example:** Polynomial derivatives are implemented in THREE places:

1. **`src/api/handlers/advanced_calculus.rs:21-63`** - Inline implementation
2. **`src/mathematics/symbolic_cas/differentiate.rs`** - CAS implementation
3. **`src/implementations/differentiator.rs`** - Unified tool implementation

### Problem 2: Orphaned Functions in Domain Modules

**Example:** `src/api_legacy/handlers/calculus.rs` references these operations:
- `handle_fractional_derivative`
- `handle_riemann_zeta`
- `handle_elliptic_integral`
- `handle_hypergeometric`
- `handle_jacobi_theta`
- `handle_bessel_function`
- `handle_legendre_polynomial`
- `handle_euler_lagrange`
- `handle_variational_calculus`
- `handle_ito_integral`
- `handle_stratonovich_integral`
- `handle_sde_solution`
- `handle_symbolic_integral`
- `handle_definite_integral`
- `handle_improper_integral`

These exist in `src/mathematics/calculus/` but are NOT wired to the unified 10-tool API!

### Problem 3: Stubs/Mocks Returning Fake Data

**Example:** `src/api/handlers/quantum_physics.rs`:
```rust
"wavefunction" | "operator" | "operators" | "entanglement" => {
    // Return mock success for API compatibility
    ComputationResponse::success(
        request.module.clone(),
        request.operation.clone(),
        json!({
            "result": "computed",
            "value": 0.0  // FAKE!
        }),
    )
}
```

Meanwhile, REAL quantum implementations exist in:
- `src/physics/quantum_mechanics/mod.rs` (7 functions)
- `src/mathematics/symbolic_cas/quantum.rs`
- `src/mathematics/symbolic_cas/quantum_advanced.rs`

### Problem 4: Inconsistent Entry Points

Users can call the same functionality through:

1. **Unified API:** `{"tool": "compute", "input": {...}}`
2. **Module API:** `{"module": "advanced_calculus", "operation": "derivative", ...}`
3. **Legacy API:** `process_calculus_request({"operation": "fractional_derivative", ...})`
4. **Direct Rust:** `crate::advanced_calculus::handle_fractional_derivative(&params)`

### Problem 5: Not All Chaos/Fractals Exposed

The chaos module (`src/mathematics/chaos/`) has:
- Mandelbrot, Julia, Burning Ship fractals
- Lorenz, Rossler attractors
- Lyapunov exponents
- Bifurcation diagrams

**NONE of these are wired to ANY API!**

### Problem 6: Machine Learning Not Exposed

The ML module (`src/specialized/machine_learning/`) has:
- K-means clustering
- DBSCAN clustering
- Neural networks
- t-SNE dimensionality reduction
- SGD, Adam optimizers

**NONE of these are wired to ANY API!**

### Problem 7: Wormholes Module Not Exposed

The wormholes module (`src/physics/wormholes/`) has:
- Morris-Thorne metric
- Traversability analysis
- Exotic matter calculations

**NONE of these are wired to ANY API!**

---

## Wiring Gaps Summary

### Completely Unwired Modules

| Module | Functions | Priority |
|--------|-----------|----------|
| `mathematics/chaos/` | Fractals, attractors, Lyapunov | HIGH |
| `specialized/machine_learning/` | Clustering, NN, optimization | HIGH |
| `physics/wormholes/` | Metrics, traversal, energy | MEDIUM |
| `specialized/game_theory/` | Nash equilibrium, evolutionary | MEDIUM |
| `tools/dimensional_analysis/` | Unit conversions, checking | HIGH |
| `tools/equation_validation/` | Equation verification | MEDIUM |

### Partially Wired Modules

| Module | Wired | Unwired | Notes |
|--------|-------|---------|-------|
| `mathematics/calculus/` | ~20% | ~80% | Legacy API has most wiring |
| `mathematics/symbolic_cas/` | ~40% | ~60% | Many advanced features unwired |
| `physics/quantum_mechanics/` | ~50% | ~50% | Mock handlers hide real code |
| `specialized/linear_programming/` | Simplex only | Dual, sensitivity | |

---

## Refactoring Recommendations

### Phase 1: Consolidate APIs

1. **Remove `src/api_legacy/`** - Port any unique wiring to unified API
2. **Fix `src/api/` stubs** - Either wire to domain modules or remove
3. **Remove inline implementations** - All math should be in domain modules

### Phase 2: Complete Unified API Wiring

1. Wire all `mathematics/chaos/` functions to `Compute` or `Simulate`
2. Wire all `specialized/machine_learning/` to new `ML` tool or `Compute`
3. Wire all `physics/wormholes/` to `Compute` or `Simulate`
4. Wire `tools/dimensional_analysis/` to `Analyze` or new `Units` tool
5. Wire remaining `mathematics/calculus/` functions

### Phase 3: Simplify Architecture

**Target Architecture:**
```
User Request
    ↓
ToolDispatcher (4 primary + legacy routing)
    ↓
Unified Implementations (10 tools)
    ↓
Domain Modules (mathematics/, physics/, tools/, specialized/)
    ↓
Result
```

**Remove:**
- `src/api/` handlers with inline implementations
- `src/api_legacy/` entirely
- Mock/stub handlers

**Keep:**
- `src/engine/` (tool traits and types)
- `src/implementations/` (enhanced with full wiring)
- Domain modules (as-is, they're the real implementations)

---

## Code Statistics

| Layer | Files | Lines | Status |
|-------|-------|-------|--------|
| Domain Modules | 100+ | ~40,000 | REAL MATH |
| Unified API | 12 | ~3,000 | INCOMPLETE |
| Module API (`src/api/`) | 28 | ~2,800 | MIXED |
| Legacy API (`src/api_legacy/`) | 17 | ~1,200 | DEPRECATED |
| **Total** | **157+** | **~47,000** | |

---

## Next Steps

1. **Audit `src/api_legacy/`** - Extract any unique wiring patterns
2. **Audit `src/api/`** - Identify which handlers have real vs fake implementations
3. **Create wiring plan** - Map each domain function to a tool
4. **Update `FUNCTION_INVENTORY.md`** - Include discovered unwired functions
5. **Begin incremental refactoring** - One module at a time
