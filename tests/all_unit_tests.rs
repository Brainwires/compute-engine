//! Unit test suite
//!
//! Unit tests are embedded in their respective modules via #[cfg(test)] and #[path]
//! declarations in src/. This file exists only to enable `cargo test --test all_unit_tests`.
//!
//! Test files are located in tests/unit/ but are compiled as part of their
//! source modules, giving them access to private functions and types.
