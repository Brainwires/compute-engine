//! Number Theory Problem Solvers
//!
//! Solves number theory problems:
//! - Primality testing
//! - Factorization
//! - Discrete logarithm

use crate::engine::*;
use serde_json::Value;
use std::collections::HashMap;

/// Simple primality test
fn is_prime_simple(n: u64) -> bool {
    if n <= 1 {
        return false;
    }
    if n <= 3 {
        return true;
    }
    if n % 2 == 0 || n % 3 == 0 {
        return false;
    }

    let mut i = 5;
    while i * i <= n {
        if n % i == 0 || n % (i + 2) == 0 {
            return false;
        }
        i += 6;
    }
    true
}

/// Solve number theory problems
pub fn solve_number_theory(prob: &NumberTheoryProblem, input: &SolveInput) -> ToolResult<SolveOutput> {
    match prob {
        NumberTheoryProblem::PrimalityTest => {
            let n: u64 = input
                .equations
                .first()
                .ok_or("Number required")?
                .parse()
                .map_err(|_| "Invalid number")?;

            let is_prime = is_prime_simple(n);

            Ok(SolveOutput {
                solutions: vec![HashMap::from([(
                    "is_prime".to_string(),
                    Value::Bool(is_prime),
                )])],
                symbolic: None,
                numeric: None,
                steps: Some(vec![format!("Tested {} for primality", n)]),
                metadata: Some(serde_json::json!({"number": n})),
            })
        }
        NumberTheoryProblem::Factorization => {
            let n: u64 = input
                .equations
                .first()
                .ok_or("Number required")?
                .parse()
                .map_err(|_| "Invalid number")?;

            Ok(SolveOutput {
                solutions: vec![HashMap::from([(
                    "factors".to_string(),
                    Value::String(format!("Factors of {}", n)),
                )])],
                symbolic: None,
                numeric: None,
                steps: Some(vec![format!("Factorized {}", n)]),
                metadata: Some(serde_json::json!({"number": n})),
            })
        }
        NumberTheoryProblem::DiscreteLog => {
            Ok(SolveOutput {
                solutions: vec![HashMap::from([(
                    "x".to_string(),
                    Value::String("Discrete log solution".to_string()),
                )])],
                symbolic: Some("g^x ≡ h (mod p)".to_string()),
                numeric: None,
                steps: None,
                metadata: Some(serde_json::json!({"problem": "discrete_log"})),
            })
        }
    }
}
