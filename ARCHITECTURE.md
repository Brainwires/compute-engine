# Computational Engine - 5-Tool Architecture

## Overview

The Computational Engine has been refactored around **5 core tools** that provide a clean, unified interface for mathematical and scientific computing. This design balances **usability** (small number of intuitive tools) with **power** (deep symbolic and numerical capabilities).

## Philosophy

- **Few tools, deep capability**: Each tool accepts flexible schemas and dispatches to specialized sub-engines internally
- **Declarative interface**: Describe *what* to compute, not *how*
- **Composable results**: Rich JSON output with metadata (symbolic forms, numeric values, LaTeX, etc.)
- **Extensible**: New domains can be added as modules without changing the core API

---

## The 5 Core Tools

### 1. **Solve**
Handles equation solving, systems of equations, root finding, and optimization.

**Capabilities:**
- Symbolic and numeric equation solving
- Systems of linear/nonlinear equations
- Root finding algorithms
- Optimization problems

**Input Schema:**
```rust
SolveInput {
    equations: Vec<String>,           // e.g., ["x^2 + 2x - 8 = 0"]
    variables: Option<Vec<String>>,   // Variables to solve for
    initial_guess: Option<HashMap<String, f64>>,
    domain: Option<Domain>,           // Real or Complex
    method: Option<Method>,           // Symbolic, Numeric, or Auto
}
```

**Example:**
```rust
let request = ToolRequest::Solve(SolveInput {
    equations: vec!["x^2 - 4 = 0".to_string()],
    variables: None,
    initial_guess: None,
    domain: None,
    method: None,
});
```

---

### 2. **Differentiate**
Computes derivatives, gradients, Jacobians, and Hessians.

**Capabilities:**
- Partial derivatives
- Total derivatives
- Gradients and Jacobians
- Higher-order derivatives
- Evaluation at specific points

**Input Schema:**
```rust
DifferentiateInput {
    expression: String,              // e.g., "sin(x*y)"
    variables: Vec<String>,          // e.g., ["x", "y"]
    order: Option<Vec<usize>>,       // Order per variable
    evaluate_at: Option<HashMap<String, f64>>,
}
```

**Example:**
```rust
let request = ToolRequest::Differentiate(DifferentiateInput {
    expression: "x^2*y".to_string(),
    variables: vec!["x".to_string(), "y".to_string()],
    order: None,
    evaluate_at: None,
});
```

---

### 3. **Integrate**
Performs symbolic and numeric integration (single and multivariable).

**Capabilities:**
- Definite and indefinite integrals
- Single and multivariable integration
- Numeric quadrature
- Symbolic antiderivatives

**Input Schema:**
```rust
IntegrateInput {
    expression: String,
    variables: Vec<String>,          // Order matters
    limits: Option<Vec<[f64; 2]>>,  // Integration bounds
    method: Option<Method>,
}
```

**Example:**
```rust
let request = ToolRequest::Integrate(IntegrateInput {
    expression: "sin(x)".to_string(),
    variables: vec!["x".to_string()],
    limits: Some(vec![[0.0, std::f64::consts::PI]]),
    method: None,
});
```

---

### 4. **Analyze**
Performs various mathematical operations: simplification, expansion, factoring, transforms, etc.

**Operations:**
- Simplify
- Expand
- Factor
- Substitute
- Series expansion
- Matrix operations
- Fourier/Laplace transforms
- Limit computation
- Expression evaluation

**Input Schema:**
```rust
AnalyzeInput {
    operation: AnalyzeOperation,     // Simplify, Expand, Factor, etc.
    expression: String,
    options: HashMap<String, Value>, // Context-specific options
}
```

**Example:**
```rust
let request = ToolRequest::Analyze(AnalyzeInput {
    operation: AnalyzeOperation::Simplify,
    expression: "(x^2 - 1)/(x - 1)".to_string(),
    options: HashMap::new(),
});
```

---

### 5. **Simulate**
Runs numerical simulations for dynamic systems, differential equations, and physics problems.

**Models:**
- ODE (Ordinary Differential Equations)
- PDE (Partial Differential Equations)
- Tensor calculus
- Classical mechanics
- Optimization
- Quantum systems
- Fluid dynamics
- Electromagnetic fields

**Input Schema:**
```rust
SimulateInput {
    model: SimulationModel,          // ODE, PDE, Tensor, etc.
    equations: Vec<String>,
    variables: Vec<String>,
    parameters: HashMap<String, f64>,
    initial_conditions: Option<HashMap<String, f64>>,
    range: Option<[f64; 2]>,        // Time/parameter range
    steps: Option<usize>,
    method: Option<String>,          // Solver method
}
```

**Example:**
```rust
let request = ToolRequest::Simulate(SimulateInput {
    model: SimulationModel::ODE,
    equations: vec!["dy/dt = -k*y".to_string()],
    variables: vec!["y".to_string()],
    parameters: {
        let mut p = HashMap::new();
        p.insert("k".to_string(), 0.1);
        p
    },
    initial_conditions: Some({
        let mut ic = HashMap::new();
        ic.insert("y".to_string(), 10.0);
        ic
    }),
    range: Some([0.0, 10.0]),
    steps: Some(100),
    method: None,
});
```

---

## Architecture Layers

```
┌─────────────────────────────────────────────┐
│          JSON/User Interface                │
│     (ToolRequest/ToolResponse enums)        │
└─────────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────────┐
│           ToolDispatcher                    │
│  (Routes requests to trait implementations) │
└─────────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────────┐
│        5 Core Traits                        │
│  Solve | Differentiate | Integrate |        │
│        Analyze | Simulate                   │
└─────────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────────┐
│      Concrete Implementations               │
│  DefaultSolver, DefaultDifferentiator, etc. │
└─────────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────────┐
│       Domain-Specific Modules               │
│  mathematics | physics | tools | specialized│
└─────────────────────────────────────────────┘
```

---

## Usage Patterns

### Pattern 1: Direct Rust API

```rust
use computational_engine::implementations::create_default_dispatcher;
use computational_engine::core::{ToolRequest, SolveInput};

let dispatcher = create_default_dispatcher();
let request = ToolRequest::Solve(SolveInput { /* ... */ });
let response = dispatcher.dispatch(request)?;
```

### Pattern 2: JSON API

```rust
let json = r#"{
    "tool": "solve",
    "input": {
        "equations": ["x^2 - 4 = 0"]
    }
}"#;

let response_json = dispatcher.dispatch_json(json);
```

### Pattern 3: Custom Implementations

Implement the core traits with your own logic:

```rust
use computational_engine::core::{Solve, SolveInput, SolveOutput, ToolResult};

struct MySolver;

impl Solve for MySolver {
    fn solve(&self, input: &SolveInput) -> ToolResult<SolveOutput> {
        // Your custom implementation
    }
}
```

---

## Migration from Legacy API

The legacy domain-specific API is still supported for backwards compatibility:

```rust
// OLD (still works)
use computational_engine::api::{ComputationRequest, process_request};

let request = ComputationRequest {
    module: "linear_algebra".to_string(),
    operation: "compute_qr".to_string(),
    parameters: /* ... */,
};

// NEW (recommended)
use computational_engine::core::{ToolRequest, AnalyzeInput};

let request = ToolRequest::Analyze(AnalyzeInput {
    operation: AnalyzeOperation::Matrix,
    expression: "qr_decomposition".to_string(),
    options: /* ... */,
});
```

---

## Extending the Engine

### Adding a New Solver Method

1. Implement the `Solve` trait
2. Register with the dispatcher
3. Optionally add to domain-specific modules

### Adding a New Analysis Operation

1. Add variant to `AnalyzeOperation` enum
2. Implement in `DefaultAnalyzer` or custom analyzer
3. Update documentation

### Adding a New Simulation Model

1. Add variant to `SimulationModel` enum
2. Implement in `DefaultSimulator` or custom simulator
3. Wire to existing domain modules (e.g., fluid dynamics, quantum physics)

---

## Benefits of This Architecture

1. **Simplicity**: 5 tools instead of dozens of module-specific functions
2. **Consistency**: All tools follow the same request/response pattern
3. **Discoverability**: Easy to understand what the engine can do
4. **Flexibility**: Inputs are flexible schemas with sensible defaults
5. **Extensibility**: New capabilities can be added without changing the interface
6. **Composability**: Results from one tool can feed into another
7. **Type Safety**: Full Rust type checking with serde serialization
8. **JSON API**: Easy integration with web services and AI agents

---

## Examples

See `examples/basic_usage.rs` for complete working examples.

Run with:
```bash
cargo run --example basic_usage
```

---

## Future Enhancements

- [ ] Symbolic math library integration (SymPy-like for Rust)
- [ ] GPU acceleration for simulations
- [ ] More sophisticated ODE/PDE solvers
- [ ] LaTeX rendering
- [ ] Plotting and visualization
- [ ] Multi-precision arithmetic
- [ ] Parallel execution
- [ ] Caching and memoization
- [ ] WebAssembly support
