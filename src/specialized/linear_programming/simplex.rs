//! Simplex Algorithm Implementation
//!
//! The Simplex algorithm solves linear programming problems in standard form:
//! Maximize c^T · x subject to A · x ≤ b, x ≥ 0

use super::{LinearProgram, LinearProgramSolution, SolutionStatus};

const EPSILON: f64 = 1e-10;
const MAX_ITERATIONS: usize = 10000;

/// Solve a linear program using the Simplex algorithm
///
/// # Example
/// ```
/// use computational_engine::specialized::linear_programming::*;
///
/// // Maximize: 3x + 2y
/// // Subject to: x + y ≤ 4, x + 3y ≤ 6, x ≥ 0, y ≥ 0
/// let lp = LinearProgram::new(vec![3.0, 2.0], true)
///     .with_inequality_constraints(
///         vec![vec![1.0, 1.0], vec![1.0, 3.0]],
///         vec![4.0, 6.0]
///     );
///
/// let solution = simplex(&lp).unwrap();
/// assert_eq!(solution.status, SolutionStatus::Optimal);
/// ```
pub fn simplex(lp: &LinearProgram) -> Result<LinearProgramSolution, String> {
    // Convert to standard form (maximize with ≤ constraints)
    let mut tableau = build_initial_tableau(lp)?;

    let mut iterations = 0;

    loop {
        if iterations >= MAX_ITERATIONS {
            return Ok(LinearProgramSolution {
                solution: vec![],
                objective_value: 0.0,
                status: SolutionStatus::MaxIterations,
                iterations,
                dual_variables: None,
                slack_variables: None,
            });
        }

        // Find entering variable (most positive coefficient in objective row for maximization)
        let entering = find_entering_variable(&tableau);

        if entering.is_none() {
            // Optimal solution found
            return extract_solution(&tableau, lp, iterations);
        }

        let entering_col = entering.unwrap();

        // Find leaving variable (minimum ratio test)
        let leaving = find_leaving_variable(&tableau, entering_col);

        if leaving.is_none() {
            // Unbounded solution
            return Ok(LinearProgramSolution {
                solution: vec![],
                objective_value: f64::INFINITY,
                status: SolutionStatus::Unbounded,
                iterations,
                dual_variables: None,
                slack_variables: None,
            });
        }

        let leaving_row = leaving.unwrap();

        // Perform pivot operation
        pivot(&mut tableau, leaving_row, entering_col);

        iterations += 1;
    }
}

/// Build the initial simplex tableau
///
/// Tableau format:
/// [A | I | b]  <- constraint rows
/// [c | 0 | 0]  <- objective row
fn build_initial_tableau(lp: &LinearProgram) -> Result<Vec<Vec<f64>>, String> {
    let n_vars = lp.num_variables();
    let n_constraints = lp.num_inequality_constraints();

    if n_constraints == 0 {
        return Err("Must have at least one constraint".to_string());
    }

    // Check constraint dimensions
    if let Some(ref constraints) = lp.inequality_constraints {
        for (i, row) in constraints.iter().enumerate() {
            if row.len() != n_vars {
                return Err(format!(
                    "Constraint {} has {} variables, expected {}",
                    i,
                    row.len(),
                    n_vars
                ));
            }
        }
    }

    // Tableau dimensions: (n_constraints + 1) rows, (n_vars + n_constraints + 1) columns
    // Columns: [original vars | slack vars | RHS]
    let n_cols = n_vars + n_constraints + 1;
    let n_rows = n_constraints + 1;

    let mut tableau = vec![vec![0.0; n_cols]; n_rows];

    // Fill constraint rows
    if let (Some(a), Some(b)) = (&lp.inequality_constraints, &lp.inequality_bounds) {
        for i in 0..n_constraints {
            // Check if RHS is negative and multiply entire constraint by -1 if needed
            let sign = if b[i] < 0.0 { -1.0 } else { 1.0 };

            // Original variables
            for j in 0..n_vars {
                tableau[i][j] = sign * a[i][j];
            }
            // Slack variable (identity matrix)
            tableau[i][n_vars + i] = sign * 1.0;
            // RHS (made non-negative)
            tableau[i][n_cols - 1] = sign * b[i];
        }
    }

    // Fill objective row
    // We always convert to maximization:
    // - For maximize f(x): use -f(x) in tableau (standard Simplex for max)
    // - For minimize f(x): maximize -f(x), so use -(-f(x)) = f(x) in tableau
    // Then we look for most negative coefficient to enter
    let obj_sign = if lp.maximize { -1.0 } else { 1.0 };
    for j in 0..n_vars {
        tableau[n_constraints][j] = obj_sign * lp.objective[j];
    }

    Ok(tableau)
}

/// Find the entering variable (most negative coefficient in objective row)
///
/// We use negative coefficients in the tableau for maximization problems,
/// so we look for the most negative value. When it's >= 0, we've reached optimality.
fn find_entering_variable(tableau: &[Vec<f64>]) -> Option<usize> {
    let obj_row = tableau.len() - 1;
    let n_cols = tableau[0].len();

    let mut min_coeff = 0.0;
    let mut entering_col = None;

    // Check all columns except RHS
    // Look for the most negative coefficient
    for j in 0..(n_cols - 1) {
        if tableau[obj_row][j] < min_coeff - EPSILON {
            min_coeff = tableau[obj_row][j];
            entering_col = Some(j);
        }
    }

    entering_col
}

/// Find the leaving variable using minimum ratio test
fn find_leaving_variable(tableau: &[Vec<f64>], entering_col: usize) -> Option<usize> {
    let n_rows = tableau.len() - 1; // Exclude objective row
    let rhs_col = tableau[0].len() - 1;

    let mut min_ratio = f64::INFINITY;
    let mut leaving_row = None;

    for i in 0..n_rows {
        let pivot_element = tableau[i][entering_col];

        if pivot_element > EPSILON {
            let ratio = tableau[i][rhs_col] / pivot_element;

            if ratio >= 0.0 && ratio < min_ratio {
                min_ratio = ratio;
                leaving_row = Some(i);
            }
        }
    }

    leaving_row
}

/// Perform pivot operation
fn pivot(tableau: &mut [Vec<f64>], pivot_row: usize, pivot_col: usize) {
    let n_rows = tableau.len();
    let n_cols = tableau[0].len();

    let pivot_element = tableau[pivot_row][pivot_col];

    // Normalize pivot row
    for j in 0..n_cols {
        tableau[pivot_row][j] /= pivot_element;
    }

    // Eliminate pivot column in other rows
    for i in 0..n_rows {
        if i != pivot_row {
            let factor = tableau[i][pivot_col];
            for j in 0..n_cols {
                tableau[i][j] -= factor * tableau[pivot_row][j];
            }
        }
    }
}

/// Extract solution from final tableau
fn extract_solution(
    tableau: &[Vec<f64>],
    lp: &LinearProgram,
    iterations: usize,
) -> Result<LinearProgramSolution, String> {
    let n_vars = lp.num_variables();
    let n_constraints = lp.num_inequality_constraints();
    let n_cols = tableau[0].len();
    let n_rows = tableau.len() - 1; // Exclude objective row
    let rhs_col = n_cols - 1;
    let obj_row = tableau.len() - 1;

    let mut solution = vec![0.0; n_vars];
    let mut slack_variables = vec![0.0; n_constraints];

    // Extract basic variable values
    for j in 0..(n_vars + n_constraints) {
        // Check if this column is a basic variable (has exactly one 1 and rest 0s)
        let mut one_count = 0;
        let mut one_row = None;

        for i in 0..n_rows {
            if (tableau[i][j] - 1.0).abs() < EPSILON {
                one_count += 1;
                one_row = Some(i);
            } else if tableau[i][j].abs() > EPSILON {
                one_count = 2; // Not a basic variable
                break;
            }
        }

        if one_count == 1 {
            let value = tableau[one_row.unwrap()][rhs_col];
            if j < n_vars {
                solution[j] = value;
            } else {
                slack_variables[j - n_vars] = value;
            }
        }
    }

    // Objective value
    // The RHS of the objective row represents:
    // - For maximization: the maximum value (use directly)
    // - For minimization: we maximized -f(x), so RHS is max(-f(x)) = -min(f(x)),
    //   therefore min(f(x)) = -RHS
    let obj_value = if lp.maximize {
        tableau[obj_row][rhs_col]
    } else {
        -tableau[obj_row][rhs_col]
    };

    Ok(LinearProgramSolution {
        solution,
        objective_value: obj_value,
        status: SolutionStatus::Optimal,
        iterations,
        dual_variables: None,
        slack_variables: Some(slack_variables),
    })
}

/// Check if a solution is feasible
pub fn is_feasible(lp: &LinearProgram, solution: &[f64]) -> bool {
    if solution.len() != lp.num_variables() {
        return false;
    }

    // Check non-negativity
    for &x in solution {
        if x < -EPSILON {
            return false;
        }
    }

    // Check inequality constraints
    if let (Some(a), Some(b)) = (&lp.inequality_constraints, &lp.inequality_bounds) {
        for i in 0..a.len() {
            let lhs: f64 = a[i].iter().zip(solution).map(|(a_ij, x_j)| a_ij * x_j).sum();
            if lhs > b[i] + EPSILON {
                return false;
            }
        }
    }

    // Check equality constraints
    if let (Some(a_eq), Some(b_eq)) = (&lp.equality_constraints, &lp.equality_bounds) {
        for i in 0..a_eq.len() {
            let lhs: f64 = a_eq[i]
                .iter()
                .zip(solution)
                .map(|(a_ij, x_j)| a_ij * x_j)
                .sum();
            if (lhs - b_eq[i]).abs() > EPSILON {
                return false;
            }
        }
    }

    true
}

/// Calculate objective value for a given solution
pub fn objective_value(lp: &LinearProgram, solution: &[f64]) -> f64 {
    let value: f64 = lp
        .objective
        .iter()
        .zip(solution)
        .map(|(c, x)| c * x)
        .sum();

    if lp.maximize {
        value
    } else {
        value
    }
}

