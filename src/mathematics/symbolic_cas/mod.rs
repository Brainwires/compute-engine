//! Custom Symbolic Computer Algebra System
//!
//! A lightweight symbolic mathematics library built in-house
//! providing expression manipulation, simplification, differentiation,
//! and basic polynomial operations.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use thiserror::Error;

pub mod christoffel;
pub mod differentiate;
pub mod expr;
pub mod fluid_dynamics;
pub mod integrate;
pub mod metric_tensors;
pub mod parser;
pub mod quantum;
pub mod quantum_advanced;
pub mod simplify;
pub mod statistical_mechanics;
pub mod symbolic_eigenvalues;
pub mod symbolic_matrix;
pub mod symbolic_tensor;

pub use christoffel::{
    christoffel_first_kind, christoffel_symbols, einstein_tensor, geodesic_coefficients,
    ricci_scalar, ricci_tensor, riemann_tensor,
};
pub use expr::{Expr, Rational};
pub use fluid_dynamics::{
    bernoulli_equation, continuity_equation_incompressible, drag_coefficient,
    euler_equation_symbolic, navier_stokes_momentum_symbolic, poiseuille_flow_velocity,
    reynolds_number, stokes_drag_force, stokes_flow_equation, stream_function_velocity,
    vorticity_2d,
};
pub use metric_tensors::*;
pub use parser::parse;
pub use quantum::{
    angular_momentum_x, angular_momentum_y, angular_momentum_z, annihilation_operator_symbolic,
    anticommutator, commutator, creation_operator_symbolic, dirac_gamma_0, dirac_gamma_1,
    dirac_gamma_2, dirac_gamma_3, dirac_gamma_matrices, expectation_value, pauli_matrices, pauli_x,
    pauli_y, pauli_z, time_evolution_operator, uncertainty_commutator, verify_pauli_properties,
};
pub use quantum_advanced::{
    bell_state_phi_minus, bell_state_phi_plus, bell_state_psi_minus, bell_state_psi_plus,
    cnot_gate, density_matrix_mixed, density_matrix_pure_state, density_matrix_trace,
    hadamard_gate, is_potentially_entangled, maximally_mixed_state, partial_trace_qubit,
    pauli_x_gate, pauli_y_gate, pauli_z_gate, phase_gate, rotation_x_gate, rotation_y_gate,
    rotation_z_gate, state_fidelity, swap_gate, t_gate, toffoli_gate, von_neumann_entropy_symbolic,
};
pub use statistical_mechanics::{
    boltzmann_distribution, bose_einstein_distribution, carnot_efficiency,
    entropy_from_partition_function, fermi_dirac_distribution, gibbs_free_energy,
    heat_capacity_constant_pressure, heat_capacity_constant_volume, helmholtz_free_energy,
    ideal_gas_law, maxwell_boltzmann_speed_distribution, partition_function_discrete,
    partition_function_ideal_gas, planck_distribution, stefan_boltzmann_law,
    van_der_waals_equation,
};
pub use symbolic_eigenvalues::{characteristic_polynomial, eigenvalues_2x2, matrix_inverse};
pub use symbolic_matrix::SymbolicMatrix;
pub use symbolic_tensor::{IndexType, SymbolicTensor};

// Re-export main functions
pub use differentiate::differentiate as diff;
pub use integrate::integrate as integrate_expr;
pub use simplify::{expand as expand_expr, simplify as simplify_expr};

/// Error types for symbolic operations
#[derive(Error, Debug)]
pub enum SymbolicError {
    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Simplification error: {0}")]
    SimplificationError(String),

    #[error("Differentiation error: {0}")]
    DifferentiationError(String),

    #[error("Integration error: {0}")]
    IntegrationError(String),

    #[error("Invalid operation: {0}")]
    InvalidOperation(String),

    #[error("Evaluation error: {0}")]
    EvaluationError(String),
}

pub type SymbolicResult<T> = Result<T, SymbolicError>;

/// Output structure for symbolic operations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SymbolicOutput {
    /// The resulting expression as a string
    pub expression: String,

    /// LaTeX representation (optional)
    pub latex: Option<String>,

    /// MathML representation (optional)
    pub mathml: Option<String>,

    /// Additional metadata about the operation
    pub metadata: Option<serde_json::Value>,
}

impl SymbolicOutput {
    pub fn new(expression: String) -> Self {
        SymbolicOutput {
            expression,
            latex: None,
            mathml: None,
            metadata: None,
        }
    }

    pub fn with_latex(mut self, latex: String) -> Self {
        self.latex = Some(latex);
        self
    }

    pub fn with_mathml(mut self, mathml: String) -> Self {
        self.mathml = Some(mathml);
        self
    }

    pub fn with_metadata(mut self, metadata: serde_json::Value) -> Self {
        self.metadata = Some(metadata);
        self
    }
}

/// Simplify an expression
pub fn simplify(expr_str: &str) -> SymbolicResult<SymbolicOutput> {
    let expr = parse(expr_str)?;
    let simplified = simplify_expr(&expr);

    let result_str = format!("{}", simplified);
    let latex = format!("${}", result_str.replace("*", " ")) + "$";

    Ok(SymbolicOutput::new(result_str)
        .with_latex(latex)
        .with_metadata(serde_json::json!({
            "operation": "simplify",
            "original": expr_str,
        })))
}

/// Expand an expression
pub fn expand(expr_str: &str) -> SymbolicResult<SymbolicOutput> {
    let expr = parse(expr_str)?;
    let expanded = expand_expr(&expr);

    let result_str = format!("{}", expanded);
    let latex = format!("${}", result_str.replace("*", " ")) + "$";

    Ok(SymbolicOutput::new(result_str)
        .with_latex(latex)
        .with_metadata(serde_json::json!({
            "operation": "expand",
            "original": expr_str,
        })))
}

/// Compute symbolic derivative
pub fn differentiate(
    expr_str: &str,
    var: &str,
    order: Option<usize>,
) -> SymbolicResult<SymbolicOutput> {
    let expr = parse(expr_str)?;
    let order = order.unwrap_or(1);

    let mut derivative = expr;
    for _ in 0..order {
        derivative = diff(&derivative, var);
    }

    let result_str = format!("{}", derivative);
    let latex = format!("${}", result_str.replace("*", " ")) + "$";

    Ok(SymbolicOutput::new(result_str)
        .with_latex(latex)
        .with_metadata(serde_json::json!({
            "operation": "differentiate",
            "variable": var,
            "order": order,
            "original": expr_str,
        })))
}

/// Substitute variables in an expression
pub fn substitute(
    expr_str: &str,
    rules: &HashMap<String, String>,
) -> SymbolicResult<SymbolicOutput> {
    let expr = parse(expr_str)?;
    let mut result = expr;

    // Apply each substitution rule
    for (from, to) in rules {
        let replacement = parse(to)?;
        result = substitute_internal(&result, from, &replacement);
    }

    let result_str = format!("{}", result);
    let latex = format!("${}", result_str.replace("*", " ")) + "$";

    Ok(SymbolicOutput::new(result_str)
        .with_latex(latex)
        .with_metadata(serde_json::json!({
            "operation": "substitute",
            "original": expr_str,
            "rules": rules,
        })))
}

fn substitute_internal(expr: &Expr, var: &str, replacement: &Expr) -> Expr {
    match expr {
        Expr::Symbol(s) if s == var => replacement.clone(),
        Expr::Add(a, b) => Expr::add(
            substitute_internal(a, var, replacement),
            substitute_internal(b, var, replacement),
        ),
        Expr::Mul(a, b) => Expr::mul(
            substitute_internal(a, var, replacement),
            substitute_internal(b, var, replacement),
        ),
        Expr::Pow(base, exp) => Expr::pow(
            substitute_internal(base, var, replacement),
            substitute_internal(exp, var, replacement),
        ),
        Expr::Function(name, args) => {
            let new_args: Vec<Expr> = args
                .iter()
                .map(|arg| substitute_internal(arg, var, replacement))
                .collect();
            Expr::func(name.clone(), new_args)
        }
        _ => expr.clone(),
    }
}

/// Evaluate an expression at specific values
pub fn evaluate_at(expr_str: &str, values: &HashMap<String, f64>) -> SymbolicResult<f64> {
    let expr = parse(expr_str)?;
    expr.evaluate(values)
}

/// Factor an expression (basic implementation)
pub fn factor(expr_str: &str) -> SymbolicResult<SymbolicOutput> {
    // For now, just parse and simplify
    // Full factorization requires polynomial GCD and other advanced algorithms
    let expr = parse(expr_str)?;
    let simplified = simplify_expr(&expr);

    let result_str = format!("{}", simplified);
    let latex = format!("${}", result_str.replace("*", " \\cdot ")) + "$";

    Ok(SymbolicOutput::new(result_str)
        .with_latex(latex)
        .with_metadata(serde_json::json!({
            "operation": "factor",
            "original": expr_str,
            "note": "Basic factorization - advanced polynomial factoring in progress"
        })))
}

/// Collect terms by a variable (basic implementation)
pub fn collect(expr_str: &str, var: &str) -> SymbolicResult<SymbolicOutput> {
    let expr = parse(expr_str)?;
    let simplified = simplify_expr(&expr);

    let result_str = format!("{}", simplified);
    let latex = format!("${}", result_str.replace("*", " ")) + "$";

    Ok(SymbolicOutput::new(result_str)
        .with_latex(latex)
        .with_metadata(serde_json::json!({
            "operation": "collect",
            "original": expr_str,
            "variable": var,
        })))
}

/// Compute polynomial GCD (basic implementation)
pub fn gcd_poly(poly1_str: &str, poly2_str: &str) -> SymbolicResult<SymbolicOutput> {
    // Placeholder - full implementation requires polynomial division algorithm
    let result_str = "1".to_string();
    let latex = format!("$\\gcd({}, {}) = 1$", poly1_str, poly2_str);

    Ok(SymbolicOutput::new(result_str)
        .with_latex(latex)
        .with_metadata(serde_json::json!({
            "operation": "gcd",
            "polynomial1": poly1_str,
            "polynomial2": poly2_str,
            "note": "Advanced polynomial GCD in progress"
        })))
}

/// Compute polynomial LCM (basic implementation)
pub fn lcm_poly(poly1_str: &str, poly2_str: &str) -> SymbolicResult<SymbolicOutput> {
    // Placeholder - LCM = (p1 * p2) / GCD(p1, p2)
    let expr1 = parse(poly1_str)?;
    let expr2 = parse(poly2_str)?;
    let product = Expr::mul(expr1, expr2);

    let result_str = format!("{}", product);
    let latex = format!(
        "$\\text{{lcm}}({}, {}) = {}$",
        poly1_str, poly2_str, result_str
    );

    Ok(SymbolicOutput::new(result_str)
        .with_latex(latex)
        .with_metadata(serde_json::json!({
            "operation": "lcm",
            "polynomial1": poly1_str,
            "polynomial2": poly2_str,
            "note": "Advanced polynomial LCM in progress"
        })))
}

/// Compute partial derivatives (gradient)
pub fn gradient(expr_str: &str, vars: &[String]) -> SymbolicResult<Vec<SymbolicOutput>> {
    let mut results = Vec::new();

    for var in vars {
        let partial = differentiate(expr_str, var, None)?;
        results.push(partial);
    }

    Ok(results)
}

/// Compute symbolic integration
pub fn integrate(expr_str: &str, var: &str) -> SymbolicResult<SymbolicOutput> {
    let expr = parse(expr_str)?;
    let integral = integrate_expr(&expr, var);

    let result_str = format!("{}", integral);
    let latex = format!("${}", result_str.replace("*", " \\cdot ")) + "$";

    Ok(SymbolicOutput::new(result_str)
        .with_latex(latex)
        .with_metadata(serde_json::json!({
            "operation": "integrate",
            "variable": var,
            "original": expr_str,
            "method": "symbolic_integration",
        })))
}

/// Compute Taylor/Maclaurin series expansion (basic implementation)
pub fn series_expansion(
    expr_str: &str,
    var: &str,
    point: f64,
    order: usize,
) -> SymbolicResult<SymbolicOutput> {
    let expr = parse(expr_str)?;

    // Compute derivatives
    let mut terms = Vec::new();
    let mut current = expr.clone();

    for n in 0..=order {
        // Evaluate derivative at point
        let mut values = HashMap::new();
        values.insert(var.to_string(), point);

        if let Ok(coeff) = current.evaluate(&values) {
            if coeff.abs() > 1e-10 {
                // factorial(n)
                let factorial: i64 = if n == 0 { 1 } else { (1..=n as i64).product() };

                // For exact Taylor series, use rational arithmetic
                // Convert f^(n)(point) to rational and divide by n!
                // For now, approximate with higher precision
                let numerator = (coeff * 1_000_000.0).round() as i64;
                let term_coeff_rational =
                    Expr::rational_unchecked(numerator, 1_000_000 * factorial);

                // Create term: coeff * (x - point)^n
                let term = if n == 0 {
                    term_coeff_rational
                } else if point == 0.0 {
                    Expr::mul(
                        term_coeff_rational,
                        Expr::pow(Expr::sym(var), Expr::num(n as i64)),
                    )
                } else {
                    Expr::mul(
                        term_coeff_rational,
                        Expr::pow(
                            Expr::add(Expr::sym(var), Expr::num(-(point as i64))),
                            Expr::num(n as i64),
                        ),
                    )
                };
                terms.push(term);
            }
        }

        // Next derivative
        current = diff(&current, var);
    }

    let result = terms
        .into_iter()
        .reduce(|acc, term| Expr::add(acc, term))
        .unwrap_or(Expr::num(0));

    let result_str = format!("{}", result);
    let latex = format!("${}", result_str.replace("*", " ")) + "$";

    Ok(SymbolicOutput::new(result_str)
        .with_latex(latex)
        .with_metadata(serde_json::json!({
            "operation": "series",
            "variable": var,
            "point": point,
            "order": order,
            "original": expr_str,
        })))
}

/// Compute limit
pub fn limit(
    expr_str: &str,
    var: &str,
    point: &str,
    direction: Option<&str>,
) -> SymbolicResult<SymbolicOutput> {
    let expr = parse(expr_str)?;

    // Handle special cases for point
    let limit_result = if point == "inf" || point == "∞" || point == "infinity" {
        compute_limit_at_infinity(&expr, var, true)
    } else if point == "-inf" || point == "-∞" || point == "-infinity" {
        compute_limit_at_infinity(&expr, var, false)
    } else {
        // Numeric point
        let point_value = point
            .parse::<f64>()
            .map_err(|_| SymbolicError::ParseError(format!("Invalid point: {}", point)))?;
        compute_limit_at_point(&expr, var, point_value, direction)
    };

    let result_str = match limit_result {
        LimitValue::Finite(val) => val,
        LimitValue::Infinity => "∞".to_string(),
        LimitValue::NegativeInfinity => "-∞".to_string(),
        LimitValue::Undefined => "undefined".to_string(),
        LimitValue::DNE => "does not exist".to_string(),
    };

    let latex = format!(
        "$\\lim_{{{}\\to{}}} {} = {}$",
        var,
        point,
        expr_str.replace("*", " \\cdot "),
        result_str
    );

    Ok(SymbolicOutput::new(result_str.clone())
        .with_latex(latex)
        .with_metadata(serde_json::json!({
            "operation": "limit",
            "variable": var,
            "point": point,
            "direction": direction,
            "original": expr_str,
            "result": result_str,
        })))
}

#[derive(Debug, Clone)]
enum LimitValue {
    Finite(String),
    Infinity,
    NegativeInfinity,
    Undefined,
    DNE, // Does Not Exist
}

/// Compute limit at a specific point
fn compute_limit_at_point(
    expr: &Expr,
    var: &str,
    point: f64,
    direction: Option<&str>,
) -> LimitValue {
    // First try direct substitution
    let mut values = HashMap::new();
    values.insert(var.to_string(), point);

    match expr.evaluate(&values) {
        Ok(val) if val.is_finite() => {
            return LimitValue::Finite(format_float(val));
        }
        Ok(val) if val.is_infinite() && val.is_sign_positive() => {
            return LimitValue::Infinity;
        }
        Ok(val) if val.is_infinite() && val.is_sign_negative() => {
            return LimitValue::NegativeInfinity;
        }
        _ => {
            // Direct substitution failed, try L'Hôpital's rule or algebraic simplification
        }
    }

    // Try algebraic simplification first
    let simplified = simplify_expr(expr);
    match simplified.evaluate(&values) {
        Ok(val) if val.is_finite() => {
            return LimitValue::Finite(format_float(val));
        }
        _ => {}
    }

    // Try L'Hôpital's rule for 0/0 or ∞/∞ forms
    if let Some(limit_val) = try_lhopital_rule(expr, var, point) {
        return limit_val;
    }

    // Try approaching from the specified direction
    let epsilon = 1e-6;
    match direction {
        Some("left") | Some("-") => {
            values.insert(var.to_string(), point - epsilon);
        }
        Some("right") | Some("+") => {
            values.insert(var.to_string(), point + epsilon);
        }
        _ => {
            // Check both sides
            values.insert(var.to_string(), point - epsilon);
            let left_val = expr.evaluate(&values).ok();
            values.insert(var.to_string(), point + epsilon);
            let right_val = expr.evaluate(&values).ok();

            match (left_val, right_val) {
                (Some(l), Some(r)) if (l - r).abs() < 1e-4 && l.is_finite() => {
                    return LimitValue::Finite(format_float(l));
                }
                (Some(l), Some(r))
                    if l.is_infinite()
                        && r.is_infinite()
                        && l.is_sign_positive() == r.is_sign_positive() =>
                {
                    return if l.is_sign_positive() {
                        LimitValue::Infinity
                    } else {
                        LimitValue::NegativeInfinity
                    };
                }
                _ => return LimitValue::DNE,
            }
        }
    }

    match expr.evaluate(&values) {
        Ok(val) if val.is_finite() => LimitValue::Finite(format_float(val)),
        Ok(val) if val.is_infinite() && val.is_sign_positive() => LimitValue::Infinity,
        Ok(val) if val.is_infinite() && val.is_sign_negative() => LimitValue::NegativeInfinity,
        _ => LimitValue::Undefined,
    }
}

/// Compute limit at infinity
fn compute_limit_at_infinity(expr: &Expr, var: &str, positive: bool) -> LimitValue {
    // Simplify the expression first
    let simplified = simplify_expr(expr);

    // Try evaluating at very large values
    let test_value = if positive { 1e6 } else { -1e6 };
    let mut values = HashMap::new();
    values.insert(var.to_string(), test_value);

    match simplified.evaluate(&values) {
        Ok(val) if val.abs() < 1e-3 => {
            // Approaching zero
            LimitValue::Finite("0".to_string())
        }
        Ok(val) if val.is_finite() => {
            // Appears to converge to a finite value
            // Test with larger value to confirm
            values.insert(var.to_string(), if positive { 1e8 } else { -1e8 });
            match simplified.evaluate(&values) {
                Ok(val2) if (val - val2).abs() < val.abs() * 0.1 => {
                    LimitValue::Finite(format_float(val))
                }
                _ => LimitValue::Infinity,
            }
        }
        Ok(val) if val.is_infinite() && val.is_sign_positive() => LimitValue::Infinity,
        Ok(val) if val.is_infinite() && val.is_sign_negative() => LimitValue::NegativeInfinity,
        _ => LimitValue::Undefined,
    }
}

/// Try L'Hôpital's rule for indeterminate forms
fn try_lhopital_rule(expr: &Expr, var: &str, point: f64) -> Option<LimitValue> {
    // Check if expression is a quotient (a * b^(-1) form)
    if let Expr::Mul(a, b) = expr {
        if let Expr::Pow(base, exp) = b.as_ref() {
            // Check if exponent is -1 (division)
            let mut test_values = HashMap::new();
            test_values.insert("_test_".to_string(), 1.0);
            if let Ok(exp_val) = exp.evaluate(&test_values) {
                if (exp_val + 1.0).abs() < 1e-10 {
                    // This is a / b form
                    let numerator = a.as_ref();
                    let denominator = base.as_ref();

                    // Evaluate both at the point
                    let mut values = HashMap::new();
                    values.insert(var.to_string(), point);

                    let num_val = numerator.evaluate(&values).ok()?;
                    let den_val = denominator.evaluate(&values).ok()?;

                    // Check for 0/0 form
                    if num_val.abs() < 1e-6 && den_val.abs() < 1e-6 {
                        // Apply L'Hôpital's rule: take derivatives
                        let num_deriv = diff(numerator, var);
                        let den_deriv = diff(denominator, var);

                        let num_deriv_val = num_deriv.evaluate(&values).ok()?;
                        let den_deriv_val = den_deriv.evaluate(&values).ok()?;

                        if den_deriv_val.abs() > 1e-10 {
                            let result = num_deriv_val / den_deriv_val;
                            if result.is_finite() {
                                return Some(LimitValue::Finite(format_float(result)));
                            }
                        }
                    }
                }
            }
        }
    }

    None
}

/// Format float with reasonable precision
fn format_float(val: f64) -> String {
    if val.abs() < 1e-10 {
        "0".to_string()
    } else if val.fract().abs() < 1e-10 {
        format!("{}", val.round() as i64)
    } else {
        format!("{:.6}", val)
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string()
    }
}

/// Rationalize an expression (placeholder)
pub fn rationalize(expr_str: &str) -> SymbolicResult<SymbolicOutput> {
    simplify(expr_str)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expand() {
        let result = expand("(x + 1)^2").unwrap();
        assert!(result.expression.len() > 0);
    }

    #[test]
    fn test_differentiate() {
        let result = differentiate("x^2 + 3*x", "x", None).unwrap();
        assert!(result.expression.len() > 0);
    }

    #[test]
    fn test_evaluate() {
        let mut values = HashMap::new();
        values.insert("x".to_string(), 2.0);
        let result = evaluate_at("x^2 + 3*x + 2", &values).unwrap();
        assert!((result - 12.0).abs() < 1e-5);
    }

    #[test]
    fn test_substitute() {
        let mut rules = HashMap::new();
        rules.insert("x".to_string(), "y".to_string());
        let result = substitute("x + 1", &rules).unwrap();
        assert!(result.expression.contains("y"));
    }
}
