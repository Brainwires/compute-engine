use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;

#[derive(Debug, Serialize, Deserialize)]
pub struct SymbolicIntegrationResult {
    pub success: bool,
    pub operation: String,
    pub result: Option<Value>,
    pub error: Option<String>,
}

// Basic symbolic integration patterns (simplified symbolic math)
pub fn symbolic_integrate(function_expr: &str, variable: &str) -> String {
    match function_expr.trim() {
        // Power functions
        "1" => format!("{}", variable),
        expr if expr == variable => format!("{}^2/2", variable),
        expr if expr == &format!("{}^2", variable) => format!("{}^3/3", variable),
        expr if expr == &format!("{}^3", variable) => format!("{}^4/4", variable),
        expr if expr.starts_with(&format!("{}^", variable)) => {
            if let Some(power_str) = expr.strip_prefix(&format!("{}^", variable)) {
                if let Ok(power) = power_str.parse::<f64>() {
                    if power != -1.0 {
                        return format!("{}^{}/{}", variable, power + 1.0, power + 1.0);
                    } else {
                        return format!("ln({})", variable);
                    }
                }
            }
            format!("∫{} d{}", expr, variable)
        }

        // Exponential functions
        expr if expr == &format!("exp({})", variable) => format!("exp({})", variable),
        expr if expr == &format!("e^{}", variable) => format!("e^{}", variable),
        expr if expr == &format!("exp(2*{})", variable) => format!("exp(2*{})/2", variable),
        expr if expr == &format!("exp(-{})", variable) => format!("-exp(-{})", variable),
        expr if expr.starts_with(&format!("exp(")) && expr.contains(&format!("*{}", variable)) => {
            // exp(a*x) -> exp(a*x)/a
            if let Some(coef) = extract_coefficient_exp(expr, variable) {
                format!("exp({}*{})/({:.6})", coef, variable, coef)
            } else {
                format!("∫{} d{}", expr, variable)
            }
        }

        // Trigonometric functions
        expr if expr == &format!("sin({})", variable) => format!("-cos({})", variable),
        expr if expr == &format!("cos({})", variable) => format!("sin({})", variable),
        expr if expr == &format!("tan({})", variable) => format!("-ln|cos({})|", variable),
        expr if expr == &format!("sec({})", variable) => {
            format!("ln|sec({}) + tan({})|", variable, variable)
        }
        expr if expr == &format!("csc({})", variable) => {
            format!("-ln|csc({}) + cot({})|", variable, variable)
        }
        expr if expr == &format!("cot({})", variable) => format!("ln|sin({})|", variable),
        expr if expr == &format!("sec({})^2", variable) => format!("tan({})", variable),
        expr if expr == &format!("csc({})^2", variable) => format!("-cot({})", variable),
        expr if expr == &format!("sec({})*tan({})", variable, variable) => {
            format!("sec({})", variable)
        }
        expr if expr == &format!("csc({})*cot({})", variable, variable) => {
            format!("-csc({})", variable)
        }

        // Trigonometric products
        expr if expr == &format!("sin({})*cos({})", variable, variable) => {
            format!("sin({}^2)/2", variable)
        }
        expr if expr == &format!("sin({})^2", variable) => {
            format!("{}/2 - sin(2*{})/4", variable, variable)
        }
        expr if expr == &format!("cos({})^2", variable) => {
            format!("{}/2 + sin(2*{})/4", variable, variable)
        }
        expr if expr == &format!("sin({})^3", variable) => {
            format!("-cos({}) + cos({}^3)/3", variable, variable)
        }
        expr if expr == &format!("cos({})^3", variable) => {
            format!("sin({}) - sin({}^3)/3", variable, variable)
        }

        // Inverse trigonometric functions
        expr if expr == &format!("1/(1+{}^2)", variable) => format!("arctan({})", variable),
        expr if expr == &format!("1/sqrt(1-{}^2)", variable) => format!("arcsin({})", variable),
        expr if expr == &format!("-1/sqrt(1-{}^2)", variable) => format!("arccos({})", variable),
        expr if expr == &format!("arcsin({})", variable) => {
            format!("{}*arcsin({}) + sqrt(1-{}^2)", variable, variable, variable)
        }
        expr if expr == &format!("arccos({})", variable) => {
            format!("{}*arccos({}) - sqrt(1-{}^2)", variable, variable, variable)
        }
        expr if expr == &format!("arctan({})", variable) => {
            format!("{}*arctan({}) - ln(1+{}^2)/2", variable, variable, variable)
        }

        // Logarithmic functions
        expr if expr == &format!("1/{}", variable) => format!("ln|{}|", variable),
        expr if expr == &format!("ln({})", variable) => {
            format!("{}*ln({}) - {}", variable, variable, variable)
        }
        expr if expr == &format!("ln({})^2", variable) => format!(
            "{}*ln({}^2) - 2*{}*ln({}) + 2*{}",
            variable, variable, variable, variable, variable
        ),
        expr if expr == &format!("{}*ln({})", variable, variable) => {
            format!("{}^2*ln({})/2 - {}^2/4", variable, variable, variable)
        }

        // Hyperbolic functions
        expr if expr == &format!("sinh({})", variable) => format!("cosh({})", variable),
        expr if expr == &format!("cosh({})", variable) => format!("sinh({})", variable),
        expr if expr == &format!("tanh({})", variable) => format!("ln(cosh({}))", variable),
        expr if expr == &format!("sech({})", variable) => format!("arctan(sinh({}))", variable),
        expr if expr == &format!("csch({})", variable) => format!("ln|tanh({}/2)|", variable),
        expr if expr == &format!("coth({})", variable) => format!("ln|sinh({})|", variable),
        expr if expr == &format!("sinh({})^2", variable) => {
            format!("sinh(2*{})/4 - {}/2", variable, variable)
        }
        expr if expr == &format!("cosh({})^2", variable) => {
            format!("sinh(2*{})/4 + {}/2", variable, variable)
        }

        // Square root and radicals
        expr if expr == &format!("sqrt({})", variable) => format!("2/3 * {}^(3/2)", variable),
        expr if expr == &format!("1/sqrt({})", variable) => format!("2*sqrt({})", variable),
        expr if expr == &format!("sqrt({}^2 + a^2)", variable) => {
            format!(
                "{}/2*sqrt({}^2 + a^2) + a^2/2*ln|{} + sqrt({}^2 + a^2)|",
                variable, variable, variable, variable
            )
        }
        expr if expr == &format!("sqrt({}^2 - a^2)", variable) => {
            format!(
                "{}/2*sqrt({}^2 - a^2) - a^2/2*ln|{} + sqrt({}^2 - a^2)|",
                variable, variable, variable, variable
            )
        }
        expr if expr == &format!("1/sqrt(a^2 - {}^2)", variable) => {
            format!("arcsin({}/a)", variable)
        }
        expr if expr == &format!("1/sqrt({}^2 + a^2)", variable) => {
            format!("ln|{} + sqrt({}^2 + a^2)|", variable, variable)
        }

        // Rational functions (simple cases)
        expr if expr == &format!("1/{}^2", variable) => format!("-1/{}", variable),
        expr if expr == &format!("1/{}^3", variable) => format!("-1/(2*{}^2)", variable),
        expr if expr == &format!("1/{}^4", variable) => format!("-1/(3*{}^3)", variable),
        expr if expr == &format!("1/({}^2 + a^2)", variable) => format!("arctan({}/a)/a", variable),
        expr if expr == &format!("1/({}^2 - a^2)", variable) => {
            format!("ln|({} - a)/({} + a)|/(2*a)", variable, variable)
        }
        expr if expr == &format!("{}/({}^2 + a^2)", variable, variable) => {
            format!("ln({}^2 + a^2)/2", variable)
        }

        // Exponential times polynomial
        expr if expr == &format!("{}*exp({})", variable, variable) => {
            format!("({} - 1)*exp({})", variable, variable)
        }
        expr if expr == &format!("{}^2*exp({})", variable, variable) => {
            format!("({}^2 - 2*{} + 2)*exp({})", variable, variable, variable)
        }

        // Exponential times trigonometric
        expr if expr == &format!("exp({})*sin({})", variable, variable) => format!(
            "exp({})/2 * (sin({}) - cos({}))",
            variable, variable, variable
        ),
        expr if expr == &format!("exp({})*cos({})", variable, variable) => format!(
            "exp({})/2 * (sin({}) + cos({}))",
            variable, variable, variable
        ),

        _ => {
            // Try to handle linear combinations and constants
            if let Some(result) = handle_linear_combination(function_expr, variable) {
                result
            } else {
                format!("∫{} d{} (not implemented)", function_expr, variable)
            }
        }
    }
}

// Extract coefficient from exp(a*x)
fn extract_coefficient_exp(expr: &str, variable: &str) -> Option<f64> {
    // Try to extract coefficient from "exp(a*x)" or "exp(x*a)"
    if !expr.starts_with("exp(") || !expr.ends_with(')') {
        return None;
    }

    let inner = &expr[4..expr.len() - 1]; // Remove "exp(" and ")"

    // Handle "a*x" or "x*a"
    if let Some(_star_pos) = inner.find('*') {
        let parts: Vec<&str> = inner.split('*').collect();
        if parts.len() == 2 {
            if parts[0] == variable {
                return parts[1].trim().parse::<f64>().ok();
            } else if parts[1] == variable {
                return parts[0].trim().parse::<f64>().ok();
            }
        }
    }

    None
}

// Handle linear combinations like "2*x + 3" or "a*sin(x) + b*cos(x)"
fn handle_linear_combination(expr: &str, variable: &str) -> Option<String> {
    // Simple coefficient extraction for basic cases
    if expr.contains('+') {
        let parts: Vec<&str> = expr.split('+').collect();
        let integrated_parts: Vec<String> = parts
            .iter()
            .map(|part| symbolic_integrate(part.trim(), variable))
            .collect();
        return Some(integrated_parts.join(" + "));
    }

    if expr.contains('-') && !expr.starts_with('-') {
        let parts: Vec<&str> = expr.split('-').collect();
        if parts.len() == 2 {
            let first = symbolic_integrate(parts[0].trim(), variable);
            let second = symbolic_integrate(parts[1].trim(), variable);
            return Some(format!("{} - {}", first, second));
        }
    }

    // Handle simple constants multiplied by functions
    if expr.contains('*') {
        let parts: Vec<&str> = expr.split('*').collect();
        if parts.len() == 2 {
            // Try to identify constant and function
            if let Ok(constant) = parts[0].trim().parse::<f64>() {
                let function_integral = symbolic_integrate(parts[1].trim(), variable);
                if !function_integral.starts_with('∫') {
                    return Some(format!("{} * ({})", constant, function_integral));
                }
            } else if let Ok(constant) = parts[1].trim().parse::<f64>() {
                let function_integral = symbolic_integrate(parts[0].trim(), variable);
                if !function_integral.starts_with('∫') {
                    return Some(format!("{} * ({})", constant, function_integral));
                }
            }
        }
    }

    None
}

// Numerical definite integration using Simpson's rule
pub fn definite_integral_numerical(
    function_expr: &str,
    variable: &str,
    lower_bound: f64,
    upper_bound: f64,
    num_intervals: usize,
) -> f64 {
    let h = (upper_bound - lower_bound) / num_intervals as f64;
    let mut sum = evaluate_function(function_expr, variable, lower_bound)
        + evaluate_function(function_expr, variable, upper_bound);

    // Simpson's rule coefficients
    for i in 1..num_intervals {
        let x = lower_bound + i as f64 * h;
        let coefficient = if i % 2 == 0 { 2.0 } else { 4.0 };
        sum += coefficient * evaluate_function(function_expr, variable, x);
    }

    sum * h / 3.0
}

// Improper integral evaluation (basic cases)
pub fn improper_integral_evaluate(
    function_expr: &str,
    variable: &str,
    lower_bound: Option<f64>,
    upper_bound: Option<f64>,
) -> f64 {
    // Handle some specific improper integrals
    match (lower_bound, upper_bound) {
        (Some(a), None) => {
            // Integral from a to infinity
            if function_expr == &format!("1/{}^2", variable) && a > 0.0 {
                return 1.0 / a;
            }
            if function_expr == &format!("exp(-{})", variable) && a >= 0.0 {
                return (-a).exp();
            }
            // For other cases, approximate with large upper bound
            definite_integral_numerical(function_expr, variable, a, 1000.0, 1000)
        }
        (None, Some(b)) => {
            // Integral from negative infinity to b
            if function_expr == &format!("exp({})", variable) && b <= 0.0 {
                return b.exp();
            }
            // Approximate with large negative lower bound
            definite_integral_numerical(function_expr, variable, -1000.0, b, 1000)
        }
        (None, None) => {
            // Integral from negative infinity to positive infinity
            if function_expr == &format!("exp(-{}^2)", variable) {
                return std::f64::consts::PI.sqrt();
            }
            // Gaussian integral and similar cases
            0.0 // Default for unsupported cases
        }
        (Some(a), Some(b)) => {
            // Regular definite integral
            definite_integral_numerical(function_expr, variable, a, b, 1000)
        }
    }
}

// Simple function evaluation (limited symbolic parser)
fn evaluate_function(expr: &str, variable: &str, value: f64) -> f64 {
    let expression = expr.replace(variable, &value.to_string());

    match expr {
        e if e == variable => value,
        e if e == &format!("{}^2", variable) => value * value,
        e if e == &format!("{}^3", variable) => value * value * value,
        e if e == &format!("sin({})", variable) => value.sin(),
        e if e == &format!("cos({})", variable) => value.cos(),
        e if e == &format!("tan({})", variable) => value.tan(),
        e if e == &format!("exp({})", variable) => value.exp(),
        e if e == &format!("ln({})", variable) => value.ln(),
        e if e == &format!("sqrt({})", variable) => value.sqrt(),
        e if e == &format!("1/{}", variable) => 1.0 / value,
        e if e == &format!("1/{}^2", variable) => 1.0 / (value * value),
        e if e == "1" => 1.0,
        _ => {
            // Try to parse as a number
            expression.parse::<f64>().unwrap_or(0.0)
        }
    }
}

// Integration by parts (for simple cases)
pub fn integration_by_parts(u: &str, dv: &str, variable: &str) -> String {
    let v = symbolic_integrate(dv, variable);
    let du = symbolic_derivative(u, variable);
    let remaining_integral = symbolic_integrate(&format!("({}) * ({})", v, du), variable);

    format!("({}) * ({}) - ({})", u, v, remaining_integral)
}

// Simple symbolic differentiation (helper for integration by parts)
fn symbolic_derivative(expr: &str, variable: &str) -> String {
    match expr {
        e if e == variable => "1".to_string(),
        e if e == &format!("{}^2", variable) => format!("2*{}", variable),
        e if e == &format!("{}^3", variable) => format!("3*{}^2", variable),
        e if e == &format!("sin({})", variable) => format!("cos({})", variable),
        e if e == &format!("cos({})", variable) => format!("-sin({})", variable),
        e if e == &format!("exp({})", variable) => format!("exp({})", variable),
        e if e == &format!("ln({})", variable) => format!("1/{}", variable),
        _ => {
            if expr.parse::<f64>().is_ok() {
                "0".to_string() // Derivative of constant is zero
            } else {
                format!("d/d{} [{}]", variable, expr) // Unsupported
            }
        }
    }
}

// Handler functions
pub fn handle_symbolic_integral(params: &HashMap<String, Value>) -> SymbolicIntegrationResult {
    let function = params
        .get("function")
        .and_then(|v| v.as_str())
        .unwrap_or("x");
    let variable = params
        .get("variable")
        .and_then(|v| v.as_str())
        .unwrap_or("x");

    let integral = symbolic_integrate(function, variable);

    SymbolicIntegrationResult {
        success: true,
        operation: "symbolic_integral".to_string(),
        result: Some(serde_json::json!({
            "function": function,
            "variable": variable,
            "integral": integral
        })),
        error: None,
    }
}

pub fn handle_definite_integral(params: &HashMap<String, Value>) -> SymbolicIntegrationResult {
    let function = params
        .get("function")
        .and_then(|v| v.as_str())
        .unwrap_or("x");
    let variable = params
        .get("variable")
        .and_then(|v| v.as_str())
        .unwrap_or("x");
    let lower_bound = params
        .get("lower_bound")
        .and_then(|v| v.as_f64())
        .unwrap_or(0.0);
    let upper_bound = params
        .get("upper_bound")
        .and_then(|v| v.as_f64())
        .unwrap_or(1.0);
    let num_intervals = params
        .get("num_intervals")
        .and_then(|v| v.as_u64())
        .unwrap_or(1000) as usize;

    let numerical_value =
        definite_integral_numerical(function, variable, lower_bound, upper_bound, num_intervals);
    let symbolic_integral = symbolic_integrate(function, variable);

    SymbolicIntegrationResult {
        success: true,
        operation: "definite_integral".to_string(),
        result: Some(serde_json::json!({
            "function": function,
            "variable": variable,
            "lower_bound": lower_bound,
            "upper_bound": upper_bound,
            "symbolic_integral": symbolic_integral,
            "numerical_value": numerical_value,
            "num_intervals": num_intervals
        })),
        error: None,
    }
}

pub fn handle_improper_integral(params: &HashMap<String, Value>) -> SymbolicIntegrationResult {
    let function = params
        .get("function")
        .and_then(|v| v.as_str())
        .unwrap_or("1/x^2");
    let variable = params
        .get("variable")
        .and_then(|v| v.as_str())
        .unwrap_or("x");
    let lower_bound = params.get("lower_bound").and_then(|v| v.as_f64());
    let upper_bound = params.get("upper_bound").and_then(|v| v.as_f64());

    let convergence_value =
        improper_integral_evaluate(function, variable, lower_bound, upper_bound);

    let bounds_description = match (lower_bound, upper_bound) {
        (Some(a), Some(b)) => format!("from {} to {}", a, b),
        (Some(a), None) => format!("from {} to ∞", a),
        (None, Some(b)) => format!("from -∞ to {}", b),
        (None, None) => "from -∞ to ∞".to_string(),
    };

    SymbolicIntegrationResult {
        success: true,
        operation: "improper_integral".to_string(),
        result: Some(serde_json::json!({
            "function": function,
            "variable": variable,
            "lower_bound": lower_bound,
            "upper_bound": upper_bound,
            "bounds_description": bounds_description,
            "convergence_value": convergence_value
        })),
        error: None,
    }
}
