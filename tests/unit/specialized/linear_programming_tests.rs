//! Unit tests for specialized::linear_programming

use crate::solve::specialized::linear_programming::*;
use crate::solve::specialized::linear_programming::simplex::simplex;

#[test]
fn test_simplex_basic() {
    // Maximize: 3x + 2y
    // Subject to: x + y ≤ 4
    let lp = LinearProgram {
        objective: vec![3.0, 2.0],
        inequality_constraints: Some(vec![vec![1.0, 1.0]]),
        inequality_bounds: Some(vec![4.0]),
        equality_constraints: None,
        equality_bounds: None,
        variable_bounds: None,
        maximize: true,
    };

    let result = simplex(&lp).unwrap();
    assert!(result.solution.len() > 0);
}

#[test]
fn test_linear_program_structure() {
    let lp = LinearProgram {
        objective: vec![1.0, 2.0, 3.0],
        inequality_constraints: Some(vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
        ]),
        inequality_bounds: Some(vec![10.0, 20.0]),
        equality_constraints: None,
        equality_bounds: None,
        variable_bounds: None,
        maximize: false,
    };

    assert_eq!(lp.objective.len(), 3);
    assert_eq!(lp.inequality_constraints.as_ref().unwrap().len(), 2);
}
