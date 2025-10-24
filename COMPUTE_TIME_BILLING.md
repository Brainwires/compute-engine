# Compute Time Tracking for Billing

**Added**: 2025-10-24
**Version**: 0.1.0

## Overview

The Computational Engine now includes **precise compute time tracking** in all MCP responses, enabling accurate customer billing based on actual computation time used.

## Features

### TimedResponse Structure

Every JSON response now includes timing metadata:

```json
{
  "success": true,
  "response": { /* tool-specific output */ },
  "compute_time_ms": 2.456,
  "compute_time_us": 2456,
  "compute_time_ns": 2456789
}
```

### Timing Precision Levels

1. **compute_time_ms** (f64): Floating-point milliseconds
   - Use for: Display to users, logging
   - Precision: ~microsecond level
   - Example: 2.456 ms

2. **compute_time_us** (u128): Integer microseconds
   - Use for: Moderate-precision billing
   - Precision: 1 microsecond
   - Example: 2456 μs

3. **compute_time_ns** (u128): Integer nanoseconds
   - Use for: High-precision billing, performance analytics
   - Precision: 1 nanosecond (system-dependent)
   - Example: 2456789 ns

## API Changes

### New Method (Recommended)
```rust
pub fn dispatch_json(&self, json_str: &str) -> String
```
Returns `TimedResponse` with timing information.

### Legacy Method (Deprecated)
```rust
#[deprecated]
pub fn dispatch_json_legacy(&self, json_str: &str) -> String
```
Returns old format without timing (backward compatibility only).

## Response Format

### Successful Response
```json
{
  "success": true,
  "response": {
    "tool": "compute",
    "output": {
      "result": [2.5, -2.5],
      "additional": null,
      "metadata": null
    }
  },
  "compute_time_ms": 1.234,
  "compute_time_us": 1234,
  "compute_time_ns": 1234567
}
```

### Error Response (with timing)
```json
{
  "success": false,
  "error": "Invalid equation type: unsupported",
  "compute_time_ms": 0.123,
  "compute_time_us": 123,
  "compute_time_ns": 123456
}
```

### Parse Error Response (no compute)
```json
{
  "success": false,
  "error": "Invalid JSON request: expected value at line 1 column 1",
  "compute_time_ms": 0.0,
  "compute_time_us": 0,
  "compute_time_ns": 0
}
```

## Billing Use Cases

### Example 1: Per-Millisecond Billing
```python
# Pricing: $0.001 per millisecond
compute_time_ms = response["compute_time_ms"]
cost = compute_time_ms * 0.001

# Example: 2.5 ms = $0.0025
```

### Example 2: Per-Second Billing with Microsecond Precision
```python
# Pricing: $0.10 per CPU-second
compute_time_us = response["compute_time_us"]
cost = (compute_time_us / 1_000_000) * 0.10

# Example: 2,500 μs (2.5 ms) = $0.00025
```

### Example 3: Tiered Pricing
```python
def calculate_cost(compute_time_ms):
    if compute_time_ms < 10:
        rate = 0.001  # $0.001/ms for fast queries
    elif compute_time_ms < 100:
        rate = 0.0008  # $0.0008/ms for medium queries
    else:
        rate = 0.0005  # $0.0005/ms for long queries (bulk discount)

    return compute_time_ms * rate
```

### Example 4: Compute Unit Billing
```python
# Define 1 Compute Unit = 1 millisecond
# Pricing: $0.50 per 1000 compute units

compute_units = response["compute_time_ms"]
cost = (compute_units / 1000) * 0.50

# Example: 2.5 ms = 2.5 units = $0.00125
```

## Performance Benchmarks

### Typical Compute Times by Operation

| Operation | Typical Time | Range |
|-----------|--------------|-------|
| Solve Linear System (3x3) | 50 μs | 20-100 μs |
| FFT (256 samples) | 100 μs | 50-200 μs |
| Matrix Eigenvalues (10x10) | 200 μs | 100-500 μs |
| ODE Simulation (100 steps) | 500 μs | 200-1000 μs |
| Symbolic Differentiation | 10 μs | 5-50 μs |
| Navier-Stokes (100x100 grid, 10 steps) | 50 ms | 10-200 ms |
| Monte Carlo (10,000 samples) | 5 ms | 1-20 ms |
| RSA Encryption (2048-bit) | 1 ms | 500 μs - 2 ms |

*Note: Times are for release builds on modern hardware. Debug builds are 10-100x slower.*

## Implementation Details

### Timing Measurement
- Uses `std::time::Instant` for high-resolution timing
- Timing starts immediately before `dispatch()` call
- Timing stops immediately after `dispatch()` returns
- Includes all computation time (parsing, execution, serialization of internal results)
- Excludes JSON parsing of request (measured separately for parse errors)
- Excludes JSON serialization of final response (negligible, ~1-10 μs)

### Accuracy
- **System-dependent**: Typical resolution is 1-100 nanoseconds on modern systems
- **Overhead**: Timing overhead is < 100 ns (negligible for billing)
- **Jitter**: Expect ±10% variation due to system load, CPU throttling, etc.

### Thread Safety
- `Instant` is thread-safe
- Multiple concurrent requests each get accurate individual timing
- No timing interference between concurrent requests

## Migration Guide

### For Existing Integrations

**Before** (old format):
```json
{
  "success": true,
  "response": { /* ... */ }
}
```

**After** (new format):
```json
{
  "success": true,
  "response": { /* ... */ },
  "compute_time_ms": 1.234,
  "compute_time_us": 1234,
  "compute_time_ns": 1234567
}
```

### Backward Compatibility

The old API format is still available via `dispatch_json_legacy()` but is deprecated:

```rust
// Old code (still works, but deprecated)
let response = dispatcher.dispatch_json_legacy(request);

// New code (recommended)
let response = dispatcher.dispatch_json(request);
```

All new code should use the timing-enabled `dispatch_json()` method.

## Example: MCP Server Integration

```rust
use computational_engine::implementations::create_default_dispatcher;

fn handle_mcp_request(json_request: &str) -> String {
    let dispatcher = create_default_dispatcher();

    // Returns JSON with timing information
    let response = dispatcher.dispatch_json(json_request);

    // Log timing for analytics
    if let Ok(parsed) = serde_json::from_str::<serde_json::Value>(&response) {
        if let Some(time_ms) = parsed["compute_time_ms"].as_f64() {
            eprintln!("Request completed in {:.3} ms", time_ms);
        }
    }

    response
}
```

## Example: Command-Line Tool

```bash
# Using the MCP stdin mode
echo '{"tool":"solve","input":{"equation_type":"linear_system","equations":["x+y=5","x-y=1"]}}' \
  | cargo run --release -- stdin \
  | jq '.compute_time_ms'

# Output: 0.156
```

## Best Practices

### For Billing Systems
1. **Use `compute_time_us` or `compute_time_ns`** for billing calculations (integer, no rounding errors)
2. **Store compute time** in database for audit trails
3. **Apply rate limits** based on cumulative compute time
4. **Round up** to next billing unit (e.g., round up to nearest millisecond)

### For Performance Monitoring
1. **Use `compute_time_ms`** for dashboards and logs (human-readable)
2. **Track percentiles** (p50, p95, p99) across requests
3. **Alert on anomalies** (e.g., requests > 1000ms)

### For Optimization
1. **Profile slow operations** (> 100ms)
2. **Batch similar requests** to amortize overhead
3. **Use release builds** in production (10-100x faster)
4. **Cache results** for repeated identical requests

## Security Considerations

### Timing Attacks
- Compute time may reveal information about input complexity
- For cryptographic operations, consider constant-time implementations
- Current timing is **not constant-time** (varies with input)

### Rate Limiting
- Use compute time for fair rate limiting
- Example: "1000 ms of compute per minute" instead of "100 requests per minute"
- Prevents abuse via many cheap requests

## Future Enhancements

### Planned Features
1. **Detailed timing breakdown** (parse, compute, serialize)
2. **Memory usage tracking** (peak memory used)
3. **Operation-specific metrics** (e.g., "iterations", "grid size")
4. **Batch request timing** (total and per-operation)
5. **Async timing** (for I/O-bound operations)

### Billing Integration Hooks
1. **Pre/post request hooks** for custom billing logic
2. **Quota management** (reject requests exceeding quota)
3. **Cost estimation** (predict cost before execution)

## Changelog

### v0.1.0 (2025-10-24)
- ✅ Added `TimedResponse` struct with three timing fields
- ✅ Updated `dispatch_json()` to return timed responses
- ✅ Added `dispatch_json_legacy()` for backward compatibility
- ✅ Implemented timing for all 10 tools
- ✅ Added timing to error responses
- ✅ Added example: `examples/test_timing.rs`
- ✅ Documentation: COMPUTE_TIME_BILLING.md

## Support

For questions about billing integration:
- See: `examples/test_timing.rs` for code examples
- See: `FEATURE_AUDIT_2025.md` for operation performance characteristics
- See: `API.md` for complete API reference

---

**Note**: All timing measurements are for the computation itself. Network latency, request queueing, and JSON serialization/deserialization on the client side are not included in compute_time_*.
