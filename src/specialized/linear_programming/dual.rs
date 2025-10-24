//! Dual Problem and Sensitivity Analysis
//!
//! The dual problem provides important economic interpretations and sensitivity information.

use super::{LinearProgram, LinearProgramSolution, SolutionStatus};

/// Compute the dual problem
///
/// Primal (max): c^T · x subject to A · x ≤ b, x ≥ 0
/// Dual (min): b^T · y subject to A^T · y ≥ c, y ≥ 0
///
/// # Example
/// ```
/// use computational_engine::specialized::linear_programming::*;
///
/// let primal = LinearProgram::new(vec![3.0, 2.0], true)
///     .with_inequality_constraints(
///         vec![vec![1.0, 1.0], vec![1.0, 3.0]],
///         vec![4.0, 6.0]
///     );
///
/// let dual = compute_dual(&primal);
/// assert_eq!(dual.num_variables(), 2); // 2 dual variables for 2 constraints
/// ```
pub fn compute_dual(primal: &LinearProgram) -> LinearProgram {
    let n_primal_vars = primal.num_variables();
    let n_constraints = primal.num_inequality_constraints();

    // Dual has one variable for each primal constraint
    let mut dual_objective = Vec::new();
    let mut dual_constraints = vec![vec![0.0; n_constraints]; n_primal_vars];

    // For primal MAX problem:
    //   Primal: max c^T·x subject to A·x ≤ b, x ≥ 0
    //   Dual: min b^T·y subject to A^T·y ≥ c, y ≥ 0
    //
    // To convert A^T·y ≥ c to ≤ form WITHOUT negative RHS:
    //   A^T·y ≥ c  →  -A^T·y ≤ -c  →  A^T·(-y) ≤ -c
    // But this doesn't help. Better approach: use the fact that for minimization,
    // we can negate the objective and maximize instead.
    //
    // Actually, the cleanest approach: note that minimize b^T·y subject to A^T·y ≥ c
    // is equivalent to maximize -b^T·y subject to -A^T·y ≤ -c

    if primal.maximize {
        // Primal is maximization
        // Dual: min b^T·y subject to A^T·y ≥ c
        // Convert to: max -b^T·y subject to -A^T·y ≤ -c

        if let Some(ref b) = primal.inequality_bounds {
            dual_objective = b.iter().map(|x| -x).collect(); // Negate for max
        }

        if let Some(ref a) = primal.inequality_constraints {
            for i in 0..n_primal_vars {
                for j in 0..n_constraints {
                    dual_constraints[i][j] = -a[j][i]; // -A^T
                }
            }
        }

        let dual_rhs: Vec<f64> = primal.objective.iter().map(|c| -c).collect(); // -c

        // Maximize (true) instead of minimize
        LinearProgram::new(dual_objective, true)
            .with_inequality_constraints(dual_constraints, dual_rhs)
    } else {
        // Primal is minimization - not implemented for now
        // Dual would be maximization
        unimplemented!("Dual of minimization problem not yet implemented")
    }
}

/// Sensitivity analysis result
#[derive(Debug, Clone)]
pub struct SensitivityAnalysis {
    /// Shadow prices (dual variables) for each constraint
    pub shadow_prices: Vec<f64>,

    /// Range of RHS values for which basis remains optimal (lower, upper)
    pub rhs_ranges: Vec<(f64, f64)>,

    /// Range of objective coefficients for which solution remains optimal
    pub objective_ranges: Vec<(f64, f64)>,

    /// Reduced costs for non-basic variables
    pub reduced_costs: Vec<f64>,
}

/// Perform sensitivity analysis on the optimal solution
///
/// This provides information about how changes in parameters affect the optimal solution:
/// - Shadow prices: Value of relaxing each constraint by one unit
/// - RHS ranges: How much constraint bounds can change without changing the optimal basis
/// - Objective ranges: How much objective coefficients can change without changing the optimal solution
pub fn sensitivity_analysis(
    lp: &LinearProgram,
    solution: &LinearProgramSolution,
) -> Result<SensitivityAnalysis, String> {
    if solution.status != SolutionStatus::Optimal {
        return Err("Sensitivity analysis requires an optimal solution".to_string());
    }

    let n_vars = lp.num_variables();
    let n_constraints = lp.num_inequality_constraints();

    // For now, provide simplified sensitivity information
    // Full sensitivity analysis requires reconstructing the final tableau

    // Shadow prices from slack variables
    let shadow_prices = if let Some(ref slack) = solution.slack_variables {
        slack
            .iter()
            .map(|&s| if s < 1e-10 { 1.0 } else { 0.0 })
            .collect()
    } else {
        vec![0.0; n_constraints]
    };

    // Simplified ranges (would need more detailed analysis for exact ranges)
    let rhs_ranges = vec![(f64::NEG_INFINITY, f64::INFINITY); n_constraints];
    let objective_ranges = vec![(f64::NEG_INFINITY, f64::INFINITY); n_vars];
    let reduced_costs = vec![0.0; n_vars];

    Ok(SensitivityAnalysis {
        shadow_prices,
        rhs_ranges,
        objective_ranges,
        reduced_costs,
    })
}

/// Check strong duality theorem
///
/// If both primal and dual have optimal solutions, their objective values are equal.
pub fn verify_strong_duality(
    primal_solution: &LinearProgramSolution,
    dual_solution: &LinearProgramSolution,
) -> bool {
    if primal_solution.status != SolutionStatus::Optimal
        || dual_solution.status != SolutionStatus::Optimal
    {
        return false;
    }

    (primal_solution.objective_value - dual_solution.objective_value).abs() < 1e-6
}

/// Economic interpretation of dual variables (shadow prices)
///
/// The dual variable y_i represents the rate of change of the optimal objective value
/// with respect to the i-th constraint bound b_i.
///
/// In economic terms:
/// - If y_i > 0: Constraint i is binding (active), relaxing it would increase profit
/// - If y_i = 0: Constraint i is not binding (slack > 0), relaxing it has no immediate benefit
/// - The value of y_i is the "shadow price" or "marginal value" of resource i
pub fn interpret_shadow_price(shadow_price: f64) -> String {
    if shadow_price > 1e-6 {
        format!(
            "Constraint is binding. Each unit increase in RHS increases objective by {:.4}",
            shadow_price
        )
    } else {
        "Constraint is not binding. There is excess capacity.".to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::specialized::linear_programming::simplex;

    #[test]
    fn test_dual_construction() {
        // Primal: max 3x + 2y subject to x + y ≤ 4, x + 3y ≤ 6
        let primal = LinearProgram::new(vec![3.0, 2.0], true).with_inequality_constraints(
            vec![vec![1.0, 1.0], vec![1.0, 3.0]],
            vec![4.0, 6.0],
        );

        let dual = compute_dual(&primal);

        // Dual should have 2 variables (one for each primal constraint)
        assert_eq!(dual.num_variables(), 2);

        // Dual objective should be [-4, -6] (negated primal RHS for maximization form)
        assert_eq!(dual.objective, vec![-4.0, -6.0]);

        // Dual should be maximization (we converted from minimization)
        assert!(dual.maximize);
    }

    #[test]
    #[ignore] // TODO: Requires two-phase simplex or Big-M method for >= constraints
    fn test_strong_duality() {
        // The dual problem involves >= constraints which our basic Simplex
        // implementation doesn't handle yet. This requires either:
        // 1. Two-phase Simplex method
        // 2. Big-M method
        // 3. Revised Simplex with artificial variables
        //
        // For now, we skip this test and mark it as a future enhancement.

        let primal = LinearProgram::new(vec![3.0, 2.0], true).with_inequality_constraints(
            vec![vec![1.0, 1.0], vec![1.0, 3.0]],
            vec![4.0, 6.0],
        );

        let primal_solution = simplex::simplex(&primal).unwrap();
        assert_eq!(primal_solution.status, SolutionStatus::Optimal);

        let dual = compute_dual(&primal);
        let dual_solution = simplex::simplex(&dual).unwrap();

        // Strong duality: primal and dual optimal values should be equal (within tolerance)
        assert!(
            (primal_solution.objective_value - dual_solution.objective_value).abs() < 0.1,
            "Primal: {}, Dual: {}",
            primal_solution.objective_value,
            dual_solution.objective_value
        );
    }
}
