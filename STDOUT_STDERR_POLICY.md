# stdout/stderr Policy for MCP Server Compatibility

## Critical Issue ⚠️

This computational engine is designed to be used as an **MCP (Model Context Protocol) server**, where:
- **stdout** = JSON responses ONLY
- **stderr** = All logs, debug output, errors

## Current Status

### ✅ Correct Usage

#### Main Binary (`src/main.rs`)
```rust
// ✓ CORRECT - User-facing output for CLI commands
println!("Computational Engine v{}", VERSION);  // Info command
println!("📦 {} ({} operations):", module);      // ListOps command

// ✓ CORRECT - JSON response to stdout
println!("{}", response);  // JSON command

// ✓ CORRECT - Errors to stderr
eprintln!("Error reading file '{}'", file_path);
```

#### API Module (`src/api.rs`)
```rust
// ✓ CORRECT - Returns JSON via data structure, no direct printing
ComputationResponse { success: true, result: Some(value), ... }
```

### ⚠️ Issues Found

#### Library Modules with Old `main()` Functions

These modules were originally standalone binaries and still have `main()` functions that print to stdout:

**🔴 NEVER CALL THESE MAIN FUNCTIONS:**

1. `src/signal_processing/lib.rs` - Has `main()` with `println!(result)`
2. `src/stochastic_processes/lib.rs` - Has `main()` with `println!(output)`
3. `src/cryptographic_mathematics/lib.rs` - Has `main()` with `println!(output)`
4. `src/symbolic_regression/lib.rs` - Has `main()` with `println!(result)`
5. `src/dimensional_analyzer/lib.rs` - Has `main()` with `println!(result)`
6. `src/equation_validator/lib.rs` - Has `main()` with `println!(result)`

**🔴 DEBUG OUTPUT IN PRODUCTION:**

1. `src/quantum_physics/quantum_simulation.rs` - Uses `println!` for progress
2. `src/quantum_physics/symbolic_mathematics.rs` - Uses `println!` for status

## Policy

### Rule 1: Main Binary stdout
**ONLY** the main binary (`src/main.rs`) may write to stdout, and ONLY for:
- CLI commands that intentionally show output (`info`, `list-ops`)
- JSON responses (`json` command)

### Rule 2: Library Code - NO stdout
**NEVER** use `println!` in library code:
```rust
// ❌ WRONG - Corrupts JSON output
println!("Processing...");
println!("Result: {}", value);

// ✅ CORRECT - Use stderr for logs
eprintln!("Processing...");
eprintln!("Result: {}", value);

// ✅ BETTER - Use logging crate
log::info!("Processing...");
log::debug!("Result: {}", value);
```

### Rule 3: Test Code
Tests may use `println!` for debugging, but should prefer `eprintln!`:
```rust
#[test]
fn test_something() {
    // ✅ ACCEPTABLE in tests
    println!("Debug: {}", value);

    // ✅ BETTER - won't corrupt test output
    eprintln!("Debug: {}", value);
}
```

### Rule 4: Old main() Functions
Old `main()` functions in library modules should be:
1. Removed, OR
2. Marked with `#[cfg(test)]` or `#[cfg(feature = "standalone")]`, OR
3. Converted to use stderr

## How to Fix

### Option 1: Remove Old main() Functions
```rust
// Delete these entirely - they're not used anymore
fn main() {
    // ...
}
```

### Option 2: Mark as Feature-Gated
```rust
#[cfg(feature = "standalone")]
fn main() {
    // Only compiled when building standalone binary
}
```

### Option 3: Convert to stderr
```rust
fn main() {
    let output = process_request(&input);
    eprintln!("{}", output);  // stderr instead of stdout
}
```

## MCP Server Usage

When running as an MCP server:

```rust
// Read JSON request from stdin
let request: ComputationRequest = serde_json::from_reader(std::io::stdin())?;

// Process (all internal logging goes to stderr)
let response = process_request(&request);

// Write JSON response to stdout - THIS IS THE ONLY STDOUT!
println!("{}", serde_json::to_string(&response)?);
```

Any `println!` in library code will corrupt this JSON response.

## Quick Audit Commands

### Find all println in source (excluding tests)
```bash
grep -r "println!" src/ --include="*.rs" | grep -v "test\|//"
```

### Find main functions in library modules
```bash
grep -r "fn main()" src/ --include="*.rs" | grep -v "src/main.rs"
```

### Check for debug output in production code
```bash
grep -r "println!\|print!" src/ --include="*.rs" | \
  grep -v "test\|eprintln\|src/main.rs\|//"
```

## Action Items

### Immediate (Critical for MCP)
- [ ] Remove or feature-gate all `main()` functions in library modules
- [ ] Convert all `println!` in quantum_physics to `eprintln!`
- [ ] Audit all library code for stdout usage

### Short Term
- [ ] Add lint rule to prevent `println!` in library code
- [ ] Set up logging framework (env_logger or tracing)
- [ ] Add CI check for stdout violations

### Long Term
- [ ] Remove all old standalone entry points
- [ ] Unified logging strategy
- [ ] Structured logging for better debugging

## Testing

To verify stdout/stderr separation:

```bash
# Should output ONLY JSON
cargo run -- json --request '@test.json' 2>/dev/null

# Should show ONLY logs/errors
cargo run -- json --request '@test.json' 1>/dev/null
```

## Summary

**Critical Rule**: Library code must NEVER write to stdout. Only the main binary's JSON response handler should write to stdout. Everything else (logs, errors, debug output, progress messages) must go to stderr.

This is non-negotiable for MCP server compatibility.
