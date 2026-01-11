use anyhow::Result;
use serde::Serialize;
use std::collections::HashMap;

#[derive(Debug, Serialize)]
pub struct SymbolicResult {
    pub operation: String,
    pub input_expression: String,
    pub output_expression: String,
    pub steps: Vec<String>,
    pub variables_used: Vec<String>,
    pub complexity_analysis: ComplexityAnalysis,
    pub computation_time_ms: u128,
}

#[derive(Debug, Serialize)]
pub struct ComplexityAnalysis {
    pub input_terms: usize,
    pub output_terms: usize,
    pub polynomial_degree: Option<i32>,
    pub function_types: Vec<String>,
    pub computational_complexity: String,
}

/// Symbolic differentiation
pub fn symbolic_differentiate(
    expression: &str,
    variables: &HashMap<String, f64>,
) -> Result<SymbolicResult> {
    let start = std::time::Instant::now();

    eprintln!("📊 Performing symbolic differentiation...");
    eprintln!("   Expression: {}", expression);

    let mut steps = vec![format!("Input: {}", expression)];

    // Parse expression and identify differentiation variable
    let diff_var = variables.keys().next().cloned().unwrap_or("x".to_string());
    steps.push(format!("Differentiating with respect to: {}", diff_var));

    // Basic differentiation rules
    let mut output = match expression.trim() {
        // Constant rule: d/dx(c) = 0
        expr if is_constant(expr, &diff_var) => {
            steps.push("Applied constant rule: d/dx(c) = 0".to_string());
            "0".to_string()
        }

        // Power rule: d/dx(x^n) = n*x^(n-1)
        expr if is_power_of_x(expr, &diff_var) => {
            let (base, exp) = parse_power(expr);
            if exp == "1" {
                steps.push("Applied power rule: d/dx(x) = 1".to_string());
                "1".to_string()
            } else if exp.parse::<i32>().is_ok() {
                let n: i32 = exp.parse().unwrap();
                if n == 2 {
                    steps.push(format!(
                        "Applied power rule: d/dx({}^{}) = {}*{}",
                        base, n, n, base
                    ));
                    format!("{}*{}", n, base)
                } else {
                    steps.push(format!(
                        "Applied power rule: d/dx({}^{}) = {}*{}^{}",
                        base,
                        n,
                        n,
                        base,
                        n - 1
                    ));
                    if n - 1 == 1 {
                        format!("{}*{}", n, base)
                    } else {
                        format!("{}*{}^{}", n, base, n - 1)
                    }
                }
            } else {
                steps.push("General power rule applied".to_string());
                format!("{}*{}^({}-1)*d/d{}[{}]", exp, base, exp, diff_var, exp)
            }
        }

        // Trigonometric functions
        expr if expr.starts_with("sin(") => {
            let inner = extract_function_argument(expr);
            steps.push(format!(
                "Applied trigonometric rule: d/dx(sin(u)) = cos(u)*du/dx"
            ));
            if inner == diff_var {
                format!("cos({})", inner)
            } else {
                format!("cos({})*d/d{}[{}]", inner, diff_var, inner)
            }
        }

        expr if expr.starts_with("cos(") => {
            let inner = extract_function_argument(expr);
            steps.push(format!(
                "Applied trigonometric rule: d/dx(cos(u)) = -sin(u)*du/dx"
            ));
            if inner == diff_var {
                format!("-sin({})", inner)
            } else {
                format!("-sin({})*d/d{}[{}]", inner, diff_var, inner)
            }
        }

        // Exponential functions
        expr if expr.starts_with("exp(") => {
            let inner = extract_function_argument(expr);
            steps.push(format!(
                "Applied exponential rule: d/dx(exp(u)) = exp(u)*du/dx"
            ));
            if inner == diff_var {
                format!("exp({})", inner)
            } else {
                format!("exp({})*d/d{}[{}]", inner, diff_var, inner)
            }
        }

        // Logarithmic functions
        expr if expr.starts_with("log(") || expr.starts_with("ln(") => {
            let inner = extract_function_argument(expr);
            steps.push(format!(
                "Applied logarithmic rule: d/dx(ln(u)) = (1/u)*du/dx"
            ));
            if inner == diff_var {
                format!("1/{}", inner)
            } else {
                format!("(1/{})*d/d{}[{}]", inner, diff_var, inner)
            }
        }

        // Sum rule: d/dx(f + g) = df/dx + dg/dx
        expr if expr.contains(" + ") => {
            let parts: Vec<&str> = expr.split(" + ").collect();
            steps.push(format!("Applied sum rule: d/dx(f + g) = df/dx + dg/dx"));

            let mut diff_parts = vec![];
            for part in parts {
                let part_diff = symbolic_differentiate(part, variables)?;
                diff_parts.push(part_diff.output_expression);
            }
            diff_parts.join(" + ")
        }

        // Product rule: d/dx(f*g) = f'*g + f*g'
        expr if expr.contains("*") && !expr.starts_with("sin") && !expr.starts_with("cos") => {
            let parts: Vec<&str> = expr.split("*").collect();
            if parts.len() == 2 {
                steps.push(format!("Applied product rule: d/dx(f*g) = f'*g + f*g'"));
                let f = parts[0].trim();
                let g = parts[1].trim();

                let f_prime = symbolic_differentiate(f, variables)?.output_expression;
                let g_prime = symbolic_differentiate(g, variables)?.output_expression;

                format!("({})*({}) + ({})*({}))", f_prime, g, f, g_prime)
            } else {
                // Multiple factors - generalized product rule
                steps.push("Applied generalized product rule".to_string());
                "d/dx[product_of_multiple_terms]".to_string()
            }
        }

        // Variable itself
        expr if expr.trim() == diff_var => {
            steps.push(format!("d/dx({}) = 1", diff_var));
            "1".to_string()
        }

        // Default case
        _ => {
            steps.push("Expression not recognized, treating as general function".to_string());
            format!("d/d{}[{}]", diff_var, expression)
        }
    };

    // Simplify the result
    output = simplify_symbolic(&output);
    steps.push(format!("Simplified result: {}", output));

    let variables_used = extract_variables(&format!("{} {}", expression, output));
    let complexity = analyze_expression_complexity(expression, &output);

    let computation_time = start.elapsed().as_millis();
    eprintln!("✅ Differentiation completed in {} ms", computation_time);

    Ok(SymbolicResult {
        operation: "differentiation".to_string(),
        input_expression: expression.to_string(),
        output_expression: output,
        steps,
        variables_used,
        complexity_analysis: complexity,
        computation_time_ms: computation_time,
    })
}

/// Symbolic integration
pub fn symbolic_integrate(
    expression: &str,
    variables: &HashMap<String, f64>,
) -> Result<SymbolicResult> {
    let start = std::time::Instant::now();

    eprintln!("∫ Performing symbolic integration...");
    eprintln!("   Expression: {}", expression);

    let mut steps = vec![format!("Input: ∫ {} dx", expression)];

    let int_var = variables.keys().next().cloned().unwrap_or("x".to_string());
    steps.push(format!("Integrating with respect to: {}", int_var));

    // Basic integration rules
    let mut output = match expression.trim() {
        // Constant rule: ∫c dx = c*x + C
        expr if is_constant(expr, &int_var) => {
            steps.push("Applied constant rule: ∫c dx = c*x + C".to_string());
            format!("{}*{} + C", expr, int_var)
        }

        // Power rule: ∫x^n dx = x^(n+1)/(n+1) + C
        expr if expr == int_var => {
            steps.push("Applied power rule: ∫x dx = x²/2 + C".to_string());
            format!("{}^2/2 + C", int_var)
        }

        expr if is_power_of_x(expr, &int_var) => {
            let (base, exp) = parse_power(expr);
            if let Ok(n) = exp.parse::<i32>() {
                if n == -1 {
                    steps.push("Applied logarithmic rule: ∫(1/x) dx = ln|x| + C".to_string());
                    format!("ln|{}| + C", base)
                } else {
                    let new_exp = n + 1;
                    steps.push(format!(
                        "Applied power rule: ∫x^{} dx = x^{}/{} + C",
                        n, new_exp, new_exp
                    ));
                    format!("{}^{}/{} + C", base, new_exp, new_exp)
                }
            } else {
                steps.push("General power rule applied".to_string());
                format!("{}^({}+1)/({}+1) + C", base, exp, exp)
            }
        }

        // Exponential functions
        expr if expr.starts_with("exp(") => {
            let inner = extract_function_argument(expr);
            steps.push("Applied exponential rule: ∫exp(u) du = exp(u) + C".to_string());
            if inner == int_var {
                format!("exp({}) + C", inner)
            } else {
                format!("∫exp({}) d{}", inner, int_var)
            }
        }

        // Trigonometric functions
        expr if expr.starts_with("sin(") => {
            let inner = extract_function_argument(expr);
            steps.push("Applied trigonometric rule: ∫sin(u) du = -cos(u) + C".to_string());
            if inner == int_var {
                format!("-cos({}) + C", inner)
            } else {
                format!("∫sin({}) d{}", inner, int_var)
            }
        }

        expr if expr.starts_with("cos(") => {
            let inner = extract_function_argument(expr);
            steps.push("Applied trigonometric rule: ∫cos(u) du = sin(u) + C".to_string());
            if inner == int_var {
                format!("sin({}) + C", inner)
            } else {
                format!("∫cos({}) d{}", inner, int_var)
            }
        }

        // Sum rule: ∫(f + g) dx = ∫f dx + ∫g dx
        expr if expr.contains(" + ") => {
            let parts: Vec<&str> = expr.split(" + ").collect();
            steps.push("Applied sum rule: ∫(f + g) dx = ∫f dx + ∫g dx".to_string());

            let mut int_parts = vec![];
            for part in parts {
                let part_int = symbolic_integrate(part, variables)?;
                int_parts.push(part_int.output_expression.replace(" + C", ""));
            }
            format!("{} + C", int_parts.join(" + "))
        }

        // Default case
        _ => {
            steps.push(
                "Complex expression - integration may require advanced techniques".to_string(),
            );
            format!("∫({}) d{}", expression, int_var)
        }
    };

    output = simplify_symbolic(&output);
    steps.push(format!("Final result: {}", output));

    let variables_used = extract_variables(&format!("{} {}", expression, output));
    let complexity = analyze_expression_complexity(expression, &output);

    let computation_time = start.elapsed().as_millis();
    eprintln!("✅ Integration completed in {} ms", computation_time);

    Ok(SymbolicResult {
        operation: "integration".to_string(),
        input_expression: expression.to_string(),
        output_expression: output,
        steps,
        variables_used,
        complexity_analysis: complexity,
        computation_time_ms: computation_time,
    })
}

/// Symbolic simplification
pub fn symbolic_simplify(
    expression: &str,
    _variables: &HashMap<String, f64>,
) -> Result<SymbolicResult> {
    let start = std::time::Instant::now();

    eprintln!("🔧 Simplifying expression...");
    eprintln!("   Input: {}", expression);

    let mut steps = vec![format!("Input: {}", expression)];
    let mut output = expression.to_string();

    // Basic simplification rules

    // Arithmetic simplifications
    output = output.replace("+ 0", "");
    output = output.replace("0 +", "");
    output = output.replace("- 0", "");
    output = output.replace("* 1", "");
    output = output.replace("1 *", "");
    output = output.replace("* 0", "* 0");
    output = output.replace("0 *", "0");

    if output.contains("* 0") || output.contains("0 *") {
        steps.push("Applied zero property: anything × 0 = 0".to_string());
        output = "0".to_string();
    }

    // Trigonometric identities
    if output.contains("sin^2") && output.contains("cos^2") {
        steps.push("Applied Pythagorean identity: sin²x + cos²x = 1".to_string());
        output = output.replace("sin^2(x) + cos^2(x)", "1");
        output = output.replace("cos^2(x) + sin^2(x)", "1");
    }

    // Exponential and logarithmic simplifications
    output = output.replace("exp(0)", "1");
    output = output.replace("ln(1)", "0");
    output = output.replace("log(1)", "0");
    output = output.replace("exp(ln(", "(");
    output = output.replace("ln(exp(", "(");

    if output != expression {
        steps.push(format!("Applied algebraic simplifications: {}", output));
    }

    // Factor common terms
    if output.contains("+") || output.contains("-") {
        let factored = attempt_factoring(&output);
        if factored != output {
            steps.push(format!("Factored common terms: {}", factored));
            output = factored;
        }
    }

    // Combine like terms
    let combined = combine_like_terms(&output);
    if combined != output {
        steps.push(format!("Combined like terms: {}", combined));
        output = combined;
    }

    let variables_used = extract_variables(&format!("{} {}", expression, output));
    let complexity = analyze_expression_complexity(expression, &output);

    let computation_time = start.elapsed().as_millis();
    eprintln!("✅ Simplification completed in {} ms", computation_time);

    Ok(SymbolicResult {
        operation: "simplification".to_string(),
        input_expression: expression.to_string(),
        output_expression: output,
        steps,
        variables_used,
        complexity_analysis: complexity,
        computation_time_ms: computation_time,
    })
}

/// Symbolic expansion
pub fn symbolic_expand(
    expression: &str,
    _variables: &HashMap<String, f64>,
) -> Result<SymbolicResult> {
    let start = std::time::Instant::now();

    eprintln!("📈 Expanding expression...");
    eprintln!("   Input: {}", expression);

    let mut steps = vec![format!("Input: {}", expression)];
    let mut output = expression.to_string();

    // Expand common patterns

    // (a + b)^2 = a^2 + 2ab + b^2
    if let Some(squared_binomial) = extract_squared_binomial(&output) {
        let expanded = expand_squared_binomial(&squared_binomial);
        steps.push(format!("Expanded (a + b)² = a² + 2ab + b²"));
        output = output.replace(&squared_binomial, &expanded);
    }

    // (a + b)(c + d) = ac + ad + bc + bd
    if let Some(product) = extract_binomial_product(&output) {
        let expanded = expand_binomial_product(&product);
        steps.push(format!("Expanded (a + b)(c + d) = ac + ad + bc + bd"));
        output = output.replace(&product.0, &expanded);
    }

    // Distribute multiplication over addition: a(b + c) = ab + ac
    output = distribute_multiplication(&output);
    if output != expression {
        steps.push("Distributed multiplication over addition".to_string());
    }

    let variables_used = extract_variables(&format!("{} {}", expression, output));
    let complexity = analyze_expression_complexity(expression, &output);

    let computation_time = start.elapsed().as_millis();
    eprintln!("✅ Expansion completed in {} ms", computation_time);

    Ok(SymbolicResult {
        operation: "expansion".to_string(),
        input_expression: expression.to_string(),
        output_expression: output,
        steps,
        variables_used,
        complexity_analysis: complexity,
        computation_time_ms: computation_time,
    })
}

// Helper functions

fn is_constant(expr: &str, var: &str) -> bool {
    !expr.contains(var)
        && expr
            .chars()
            .all(|c| c.is_ascii_digit() || c == '.' || c == '-')
}

fn is_power_of_x(expr: &str, var: &str) -> bool {
    expr.contains(var) && (expr.contains("^") || expr == var)
}

fn parse_power(expr: &str) -> (String, String) {
    if let Some(pos) = expr.find('^') {
        let base = expr[..pos].to_string();
        let exp = expr[pos + 1..].to_string();
        (base, exp)
    } else {
        (expr.to_string(), "1".to_string())
    }
}

fn extract_function_argument(expr: &str) -> String {
    if let Some(start) = expr.find('(') {
        if let Some(end) = expr.rfind(')') {
            return expr[start + 1..end].to_string();
        }
    }
    "x".to_string()
}

fn simplify_symbolic(expr: &str) -> String {
    let mut result = expr.to_string();

    // Remove redundant operations
    result = result.replace("+ -", "- ");
    result = result.replace("- -", "+ ");
    result = result.replace("  ", " ");
    result = result.trim().to_string();

    // Remove unnecessary parentheses around single terms
    if result.starts_with('(') && result.ends_with(')') && result.matches('(').count() == 1 {
        result = result[1..result.len() - 1].to_string();
    }

    result
}

fn extract_variables(expr: &str) -> Vec<String> {
    let mut vars = std::collections::HashSet::new();

    for word in expr.split_whitespace() {
        for ch in word.chars() {
            if ch.is_ascii_alphabetic() && ch != 'C' {
                vars.insert(ch.to_string());
            }
        }
    }

    let mut result: Vec<String> = vars.into_iter().collect();
    result.sort();
    result
}

fn analyze_expression_complexity(input: &str, output: &str) -> ComplexityAnalysis {
    let input_terms = input.matches('+').count() + input.matches('-').count() + 1;
    let output_terms = output.matches('+').count() + output.matches('-').count() + 1;

    let function_types = vec![
        if input.contains("sin") || input.contains("cos") {
            Some("trigonometric")
        } else {
            None
        },
        if input.contains("exp") || input.contains("ln") {
            Some("exponential")
        } else {
            None
        },
        if input.contains("^") {
            Some("polynomial")
        } else {
            None
        },
    ]
    .into_iter()
    .flatten()
    .map(|s| s.to_string())
    .collect();

    let polynomial_degree = if input.contains("^") {
        // Try to extract highest power
        Some(2) // Simplified
    } else if input.contains("*") {
        Some(1)
    } else {
        Some(0)
    };

    ComplexityAnalysis {
        input_terms,
        output_terms,
        polynomial_degree,
        function_types,
        computational_complexity: "O(n)".to_string(),
    }
}

fn attempt_factoring(expr: &str) -> String {
    // Very basic factoring - in practice would need sophisticated algorithms
    expr.to_string()
}

fn combine_like_terms(expr: &str) -> String {
    // Very basic combination - in practice would need term analysis
    expr.to_string()
}

fn extract_squared_binomial(expr: &str) -> Option<String> {
    // Look for patterns like (a + b)^2
    if expr.contains(")^2") {
        // Very simplified detection
        None
    } else {
        None
    }
}

fn expand_squared_binomial(expr: &str) -> String {
    // Expand (a + b)^2 = a^2 + 2ab + b^2
    expr.to_string()
}

fn extract_binomial_product(_expr: &str) -> Option<(String, String, String)> {
    // Look for patterns like (a + b)(c + d)
    None
}

fn expand_binomial_product(product: &(String, String, String)) -> String {
    // Expand (a + b)(c + d) = ac + ad + bc + bd
    product.0.clone()
}

fn distribute_multiplication(expr: &str) -> String {
    // Distribute a(b + c) = ab + ac
    expr.to_string()
}
