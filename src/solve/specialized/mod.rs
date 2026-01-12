//! Specialized Solvers
//!
//! Domain-specific solvers for:
//! - Number theory (primality, factorization, discrete log)
//! - Game theory (Nash equilibrium, cooperative games, evolutionary dynamics)
//! - Linear programming (Simplex algorithm, duality, sensitivity analysis)

pub mod game_theory;
pub mod linear_programming;
pub mod number_theory;

pub use number_theory::solve_number_theory;
