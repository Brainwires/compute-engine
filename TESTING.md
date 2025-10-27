# Testing Guide

Complete guide to testing the Brainwires Compute Engine.

## Test Statistics

- **1,685 total tests** (100% pass rate)
  - 202 integration tests
  - 1,483 unit tests
- **83.81% production code coverage** (4,204/5,016 lines covered)

## Quick Reference

```bash
# Run all tests
cargo test

# Run with release optimizations (faster for math-heavy tests)
cargo test --release

# Run specific test
cargo test test_schwarzschild_radius

# Run tests matching pattern
cargo test black_hole

# Get accurate production code coverage
cargo tarpaulin --lib --release --out Stdout

# Generate HTML coverage report
cargo tarpaulin --lib --release --out Html
```

## Test Organization

The project uses two complementary testing approaches:

### 1. Integration Tests (`tests/` directory)

**Location:** `tests/integration/` and `tests/unit/`

**Purpose:** Test public APIs and end-to-end functionality

**Count:** 202 integration tests

**Structure:**
```
tests/
├── integration/
│   ├── tools/          # Tool API tests (solve, differentiate, etc.)
│   ├── physics/        # Physics module tests
│   ├── modules/        # Domain module tests
│   ├── api/            # JSON API tests
│   └── coverage/       # Coverage verification tests
└── unit/               # External unit tests (organized by module)
```

**Example:**
```rust
// tests/integration/tools/solve_tests.rs
#[test]
fn test_solve_quadratic() {
    let dispatcher = create_default_dispatcher();
    let request = ToolRequest::Solve(SolveInput {
        equations: vec!["x^2 - 4 = 0".to_string()],
        variables: None,
        initial_guess: None,
        domain: None,
        method: None,
    });
    let response = dispatcher.dispatch(request).unwrap();
    // Assertions...
}
```

### 2. Unit Tests (Inline `#[cfg(test)]` modules)

**Location:** Inline within `src/` files

**Purpose:** Test individual functions and internal implementation

**Count:** 1,483 unit tests

**Structure:**
```rust
// src/physics/black_holes/mod.rs
pub fn schwarzschild_radius(mass: f64) -> f64 {
    2.0 * G * mass / C2
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_schwarzschild_radius() {
        let mass = SOLAR_MASS;
        let r_s = schwarzschild_radius(mass);
        assert!((r_s - 2953.0).abs() < 1.0);
    }
}
```

**Advantages:**
- Test code lives next to implementation
- Easy to test private functions
- No need for `pub` visibility hacks
- Follows Rust best practices

## Code Coverage

### Getting Accurate Coverage

**Problem:** Standard `cargo-llvm-cov --lib` includes inline test code in coverage calculations, which artificially lowers the percentage.

**Solution:** Use `cargo-tarpaulin` which automatically excludes test code:

```bash
# Install tarpaulin (one-time)
cargo install cargo-tarpaulin

# Get production code coverage
cargo tarpaulin --lib --release --out Stdout

# Generate HTML report
cargo tarpaulin --lib --release --out Html

# Open HTML report
open tarpaulin-report.html  # macOS
xdg-open tarpaulin-report.html  # Linux
```

### Why Tarpaulin?

**cargo-llvm-cov behavior:**
- Includes `#[cfg(test)]` modules in line counts
- Test code doesn't "cover" other test code
- Result: Artificially low coverage percentage

**cargo-tarpaulin behavior:**
- Automatically excludes `#[cfg(test)]` modules
- Only measures production code coverage
- Result: Accurate production code coverage (83.81%)

### Coverage Comparison

```bash
# ❌ WRONG - Includes test code (shows ~35%)
cargo llvm-cov --lib --release

# ✅ CORRECT - Excludes test code (shows 83.81%)
cargo tarpaulin --lib --release --out Stdout
```

## Writing Tests

### Unit Test Guidelines

**DO:**
- ✅ Use inline `#[cfg(test)]` modules for unit tests
- ✅ Test edge cases and error conditions
- ✅ Use descriptive test names: `test_schwarzschild_radius_for_solar_mass`
- ✅ Keep tests focused on single functionality
- ✅ Use `approx_eq` for floating-point comparisons when needed

**DON'T:**
- ❌ Test implementation details that may change
- ❌ Write tests that depend on test execution order
- ❌ Use hardcoded paths or external dependencies
- ❌ Weaken test assertions to make tests pass (fix implementation instead!)

### Integration Test Guidelines

**DO:**
- ✅ Test public API contracts
- ✅ Test end-to-end workflows
- ✅ Organize by feature/module
- ✅ Use meaningful test data
- ✅ Document test purpose in comments

**DON'T:**
- ❌ Duplicate unit test coverage
- ❌ Test internal implementation details
- ❌ Create brittle tests dependent on exact output format

### Testing Best Practices

1. **Test Naming Convention**
   ```rust
   #[test]
   fn test_<function_name>_<scenario>() {
       // Examples:
       // test_schwarzschild_radius_solar_mass()
       // test_event_horizon_kerr_rotating()
       // test_isco_radius_extremal_black_hole()
   }
   ```

2. **Floating-Point Comparisons**
   ```rust
   // ❌ WRONG - Direct equality
   assert_eq!(result, 2.9530);

   // ✅ CORRECT - Tolerance-based comparison
   assert!((result - 2.9530).abs() < 1e-4);

   // ✅ ALSO GOOD - Relative tolerance
   let relative_error = ((result - expected) / expected).abs();
   assert!(relative_error < 1e-6);
   ```

3. **Test Organization**
   ```rust
   #[cfg(test)]
   mod tests {
       use super::*;

       // Group related tests
       mod schwarzschild_tests {
           use super::*;

           #[test]
           fn test_radius() { /* ... */ }

           #[test]
           fn test_horizon() { /* ... */ }
       }

       mod kerr_tests {
           use super::*;

           #[test]
           fn test_ergosphere() { /* ... */ }
       }
   }
   ```

4. **Test Constants**
   ```rust
   #[cfg(test)]
   mod tests {
       use super::*;

       const SOLAR_MASS: f64 = 1.989e30;
       const TEST_TOLERANCE: f64 = 1e-10;

       #[test]
       fn test_something() {
           // Use constants for readability
       }
   }
   ```

## When Tests Fail

### Critical Rule: ALWAYS FIX IMPLEMENTATION, NEVER WEAKEN TESTS

**This project requires perfectly accurate and deterministic behavior.**

When a test fails:

1. ✅ **Investigate the implementation** - Assume the test is correct
2. ✅ **Fix the root cause** - Correct the algorithm or formula
3. ✅ **Verify the physics/math** - Ensure correctness

**NEVER:**
1. ❌ Weaken test assertions
2. ❌ Add tolerance when exact equality should work
3. ❌ Comment out failing tests
4. ❌ Make assumptions about test correctness

**Example - WRONG approach:**
```rust
// ❌ WRONG - Weakening test to make it pass
assert!(r_ergo >= r_h);  // Changed from > to >= to make test pass
```

**Example - CORRECT approach:**
```rust
// ✅ CORRECT - Fix the implementation bug
// In implementation file:
let a_geom = self.spin * m_geom;  // Fixed: was using wrong formula
```

See [CLAUDE.md](CLAUDE.md) section "ALWAYS Fix Implementation Bugs, NEVER Weaken Tests" for the full story.

## Running Specific Tests

### By Module

```bash
# All physics tests
cargo test physics::

# All black hole tests
cargo test physics::black_holes::

# All tensor calculus tests
cargo test mathematics::tensor_calculus::
```

### By Pattern

```bash
# All tests with "schwarzschild" in name
cargo test schwarzschild

# All tests with "integrate" in name
cargo test integrate

# All tests with "matrix" in name
cargo test matrix
```

### By File

```bash
# Run integration tests only
cargo test --test all_integration_tests

# Run specific integration test file
cargo test --test all_integration_tests integration::physics::black_holes_tests
```

### Test Output Options

```bash
# Show test output (println!, etc.)
cargo test -- --nocapture

# Show only failing tests
cargo test -- --quiet

# Run tests in single thread (for debugging)
cargo test -- --test-threads=1

# Stop on first failure
cargo test -- --fail-fast
```

## Performance Testing

### Release Mode

**Always use `--release` for math-heavy tests:**

```bash
# Much faster for numerical computations
cargo test --release

# Debug builds can be 10-100x slower
cargo test  # Avoid for performance-critical tests
```

### Benchmarking

```bash
# Run benchmarks
cargo bench

# Benchmark specific function
cargo bench schwarzschild
```

## Continuous Integration

### GitHub Actions

The project uses GitHub Actions for CI. Tests run on:
- Linux (Ubuntu latest)
- macOS (latest)
- Windows (latest)

### CI Commands

```bash
# What CI runs:
cargo check
cargo test --release
cargo clippy -- -D warnings
cargo fmt -- --check
```

## Troubleshooting

### Test Compilation Errors

```bash
# Check test compilation without running
cargo test --no-run

# Check specific test
cargo test test_name --no-run
```

### Test Hangs

```bash
# Run with timeout
timeout 300 cargo test

# Identify slow tests
cargo test -- --nocapture --test-threads=1
```

### Coverage Issues

```bash
# Clean coverage data
cargo clean

# Regenerate coverage
cargo tarpaulin --lib --release --out Html
```

## Best Practices Summary

1. **Use inline `#[cfg(test)]` for unit tests** - Keep tests close to code
2. **Use `tests/integration/` for integration tests** - Test public APIs
3. **Use `cargo-tarpaulin` for coverage** - Excludes test code automatically
4. **Always test in release mode for math** - 10-100x faster
5. **Fix implementation, never weaken tests** - Maintain code quality
6. **Use descriptive test names** - Make failures easy to understand
7. **Test edge cases** - Zero, negative, NaN, infinity, etc.
8. **Use appropriate tolerances** - Floating-point arithmetic isn't exact
9. **Keep tests independent** - No shared mutable state
10. **Document tricky tests** - Explain why the test is important

## Additional Resources

- [Rust Book - Testing](https://doc.rust-lang.org/book/ch11-00-testing.html)
- [cargo-tarpaulin Documentation](https://github.com/xd009642/tarpaulin)
- [Project CLAUDE.md](CLAUDE.md) - Critical testing guidelines
- [Rust API Guidelines - Testing](https://rust-lang.github.io/api-guidelines/documentation.html)
