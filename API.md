# Computational Engine API Documentation

## Overview

The Computational Engine provides a unified JSON-based API to access ~189 mathematical and physical computation functions across 6+ modules. This MCP-style interface allows for programmatic access to all computational capabilities.

## Quick Start

### List All Available Operations

```bash
cargo run --release -- list-ops
```

Output shows all ~42+ operation types across all modules (each operation may have multiple variants).

### Execute via JSON

Create a JSON request file:

```json
{
  "module": "advanced_calculus",
  "operation": "riemann_zeta",
  "parameters": {
    "s": {
      "real": 2.0
    }
  }
}
```

Execute:

```bash
cargo run --release -- json --request '@request.json'
```

Response:

```json
{
  "success": true,
  "module": "advanced_calculus",
  "operation": "riemann_zeta",
  "result": {
    "s": 2.0,
    "zeta_value": 1.643935066848228
  },
  "error": null
}
```

## Request Format

All requests follow this structure:

```json
{
  "module": "module_name",
  "operation": "operation_name",
  "parameters": {
    // Operation-specific parameters
  }
}
```

### Module Names

- `tensor_calculus` (or `tensor`)
- `advanced_calculus` (or `calculus`)
- `fluid_dynamics` (or `fluid`)
- `signal_processing` (or `signal`)
- `stochastic_processes` (or `stochastic`)
- `cryptographic_mathematics` (or `crypto`)
- `symbolic_regression` (or `regression`)
- `quantum_physics` (or `quantum`) - *requires `--features quantum`*

## Available Operations by Module

### Tensor Calculus (1 primary operation)

#### solve_vacuum_einstein

Solve vacuum Einstein field equations in various coordinate systems.

```json
{
  "module": "tensor_calculus",
  "operation": "solve_vacuum_einstein",
  "parameters": {
    "coordinates": ["t", "r", "theta", "phi"],
    "system": "spherical"
  }
}
```

**Returns:** Schwarzschild and Reissner-Nordström solutions

### Advanced Calculus (19 operations)

#### Fractional Calculus

1. **fractional_derivative** - Compute fractional derivatives
2. **fractional_integral** - Compute fractional integrals
3. **fractional_calculus** - General fractional calculus operations

```json
{
  "module": "advanced_calculus",
  "operation": "fractional_derivative",
  "parameters": {
    "function": "x^2",
    "order": 0.5,
    "x": 1.0
  }
}
```

#### Special Functions

4. **riemann_zeta** - Riemann zeta function ζ(s)
5. **elliptic_integral** - Complete elliptic integrals K(k), E(k)
6. **hypergeometric** - Hypergeometric functions
7. **jacobi_theta** - Jacobi theta functions
8. **bessel_function** - Bessel functions Jₙ(x)
9. **legendre_polynomial** - Legendre polynomials Pₙ(x)
10. **special_functions** - General special function operations

```json
{
  "module": "advanced_calculus",
  "operation": "riemann_zeta",
  "parameters": {
    "s": { "real": 2.0 }
  }
}
```

```json
{
  "module": "advanced_calculus",
  "operation": "bessel_function",
  "parameters": {
    "n": 0,
    "x": 1.0
  }
}
```

#### Variational Calculus

11. **euler_lagrange** - Euler-Lagrange equations
12. **variational_calculus** - Variational calculus operations

```json
{
  "module": "advanced_calculus",
  "operation": "euler_lagrange",
  "parameters": {
    "lagrangian": "x'^2 - x^2",
    "variable": "x",
    "domain": [0.0, 1.0]
  }
}
```

#### Stochastic Calculus

13. **ito_integral** - Itô stochastic integrals
14. **stratonovich_integral** - Stratonovich stochastic integrals
15. **sde_solution** - Stochastic differential equation solutions
16. **stochastic_calculus** - General stochastic calculus operations

```json
{
  "module": "advanced_calculus",
  "operation": "ito_integral",
  "parameters": {
    "function": "x",
    "initial_value": 0.0,
    "time_steps": 1000,
    "time_horizon": 1.0
  }
}
```

#### Symbolic Integration

17. **symbolic_integral** - Symbolic indefinite integrals
18. **definite_integral** - Definite integrals with bounds
19. **improper_integral** - Improper integrals

```json
{
  "module": "advanced_calculus",
  "operation": "symbolic_integral",
  "parameters": {
    "integrand": "x^2 + 2*x + 1",
    "variable": "x"
  }
}
```

### Fluid Dynamics (2+ operations)

#### cavity_flow

Run lid-driven cavity flow simulations.

```json
{
  "module": "fluid_dynamics",
  "operation": "cavity_flow",
  "parameters": {
    "grid_size": 64,
    "reynolds_number": 100.0,
    "iterations": 1000
  }
}
```

#### navier_stokes

Solve Navier-Stokes equations.

### Signal Processing (5 operations)

1. **fft** - Fast Fourier Transform
2. **lowpass_filter** - Low-pass digital filter
3. **highpass_filter** - High-pass digital filter
4. **spectrogram** - Time-frequency spectrogram
5. **psd** - Power spectral density estimation

```json
{
  "module": "signal_processing",
  "operation": "fft",
  "parameters": {
    "signal": [1.0, 2.0, 3.0, 4.0],
    "sample_rate": 1000.0,
    "window_type": "hann"
  }
}
```

### Stochastic Processes (5 operations)

1. **brownian_motion** - Brownian motion simulation
2. **markov_chain** - Markov chain simulation
3. **stochastic_integral** - Monte Carlo stochastic integrals
4. **ornstein_uhlenbeck** - Ornstein-Uhlenbeck process
5. **poisson_process** - Poisson process simulation

```json
{
  "module": "stochastic_processes",
  "operation": "brownian_motion",
  "parameters": {
    "initial_value": 0.0,
    "drift": 0.0,
    "volatility": 1.0,
    "time_steps": 1000,
    "time_horizon": 1.0
  }
}
```

### Cryptographic Mathematics (10 operations)

1. **mod_exp** - Modular exponentiation
2. **extended_gcd** - Extended Euclidean algorithm
3. **mod_inverse** - Modular multiplicative inverse
4. **chinese_remainder** - Chinese remainder theorem
5. **miller_rabin** - Miller-Rabin primality test
6. **generate_prime** - Generate random prime numbers
7. **rsa_keygen** - RSA key generation
8. **discrete_log** - Discrete logarithm (BSGS)
9. **sha256** - SHA-256 hash
10. **sha3_256** - SHA3-256 hash

```json
{
  "module": "cryptographic_mathematics",
  "operation": "generate_prime",
  "parameters": {
    "bits": 256
  }
}
```

### Symbolic Regression (planned)

- **discover_equations** - Automated equation discovery
- **evaluate_fitness** - Fitness evaluation for symbolic expressions

### Quantum Physics (optional, requires `--features quantum`)

- **quantum_simulation** - Quantum state evolution
- **particle_physics** - Particle physics simulations
- **tensor_computation** - Tensor-based quantum computations
- **symbolic_mathematics** - Symbolic quantum mathematics

## Response Format

All responses follow this structure:

```json
{
  "success": boolean,
  "module": "module_name",
  "operation": "operation_name",
  "result": {
    // Operation-specific result data
  } | null,
  "error": "error message" | null
}
```

### Success Response

```json
{
  "success": true,
  "module": "advanced_calculus",
  "operation": "riemann_zeta",
  "result": {
    "s": 2.0,
    "zeta_value": 1.643935066848228
  },
  "error": null
}
```

### Error Response

```json
{
  "success": false,
  "module": "advanced_calculus",
  "operation": "unknown_operation",
  "result": null,
  "error": "Unknown calculus operation: unknown_operation"
}
```

## Usage Examples

### Via Command Line

```bash
# List all operations
cargo run --release -- list-ops

# Execute JSON from file
cargo run --release -- json --request '@request.json'

# Execute inline JSON (use single quotes)
cargo run --release -- json --request '{"module":"tensor","operation":"solve_vacuum_einstein","parameters":{"system":"spherical"}}'
```

### Via Rust Library

```rust
use computational_engine::{ComputationRequest, process_request};
use std::collections::HashMap;

let mut params = HashMap::new();
params.insert("system".to_string(), json!("spherical"));

let request = ComputationRequest {
    module: "tensor_calculus".to_string(),
    operation: "solve_vacuum_einstein".to_string(),
    parameters: params,
};

let response = process_request(&request);
println!("{}", serde_json::to_string_pretty(&response).unwrap());
```

### Via JSON String

```rust
use computational_engine::process_json_request;

let json_request = r#"{
  "module": "advanced_calculus",
  "operation": "riemann_zeta",
  "parameters": {"s": {"real": 2.0}}
}"#;

let response = process_json_request(json_request);
println!("{}", response);
```

## Integration Patterns

### MCP Server Integration

The JSON API is designed to be compatible with Model Context Protocol (MCP) servers:

1. Each operation can be exposed as an MCP tool
2. Parameters map directly to tool arguments
3. Results are JSON-serializable for easy transport

### REST API Wrapper

Wrap the engine in a REST API:

```rust
// Example with actix-web
#[post("/compute")]
async fn compute(req: web::Json<ComputationRequest>) -> impl Responder {
    let response = process_request(&req);
    web::Json(response)
}
```

### Python Binding (via PyO3)

```python
import computational_engine

request = {
    "module": "advanced_calculus",
    "operation": "riemann_zeta",
    "parameters": {"s": {"real": 2.0}}
}

result = computational_engine.process_request(request)
print(result)
```

## Performance Considerations

- **Compilation**: Use `--release` for ~10-100x performance improvement
- **Caching**: Results can be cached externally if deterministic
- **Parallelism**: Many operations use Rayon for parallel computation
- **GPU**: Enable `--features quantum` for GPU-accelerated quantum operations

## Error Handling

Common error types:

1. **Parse Error**: Invalid JSON format
2. **Unknown Module**: Module name not recognized
3. **Unknown Operation**: Operation not available in module
4. **Invalid Parameters**: Missing or incorrect parameters
5. **Computation Error**: Error during calculation

All errors include descriptive messages in the `error` field.

## Feature Flags

- `quantum` - Enable quantum physics module with GPU support
- `gpu-acceleration` - GPU support for intensive computations
- `wasm` - WebAssembly compilation support

```bash
# Build with quantum support
cargo build --release --features quantum

# Build all features
cargo build --release --features quantum,wasm
```

## Total Operation Count

- **Tensor Calculus**: 1 primary operation
- **Advanced Calculus**: 19 operations
- **Fluid Dynamics**: 2+ operations
- **Signal Processing**: 5 operations
- **Stochastic Processes**: 5 operations
- **Cryptographic Mathematics**: 10 operations
- **Symbolic Regression**: ~5 operations (planned)
- **Quantum Physics**: ~10 operations (optional)

**Total**: ~42+ exposed operations, with ~189 underlying functions

Many operations have multiple modes and variants, making the actual capability count much higher.

## Support

For issues or questions:
- See [README.md](README.md) for general usage
- See [MIGRATION.md](MIGRATION.md) for migration from standalone projects
- Check function signatures in source code for detailed parameter documentation

## Future Enhancements

- [ ] Complete signal processing API
- [ ] Complete stochastic processes API
- [ ] Complete cryptographic mathematics API
- [ ] Add fluid dynamics API
- [ ] Add symbolic regression API
- [ ] Add dimensional analysis operations
- [ ] Add equation validation operations
- [ ] Add computational geometry operations
- [ ] WebSocket streaming for long-running computations
- [ ] Batch request processing
- [ ] Request validation schema
