//! Unified Validator implementation
//!
//! Routes validation requests to tools/equation_validation/

use crate::engine::*;
use std::collections::HashMap;

pub struct UnifiedValidator;

impl UnifiedValidator {
    pub fn new() -> Self {
        Self
    }
}

impl Default for UnifiedValidator {
    fn default() -> Self {
        Self::new()
    }
}

impl Validate for UnifiedValidator {
    fn validate(&self, input: &ValidateInput) -> ToolResult<ValidateOutput> {
        use crate::tools::equation_validation;

        match &input.operation {
            ValidateOp::Equation => {
                // Full equation validation
                let domain = input.parameters
                    .get("domain")
                    .and_then(|v| v.as_str())
                    .unwrap_or("general")
                    .to_string();

                let conservation_laws: Option<Vec<String>> = input.parameters
                    .get("conservation_laws")
                    .and_then(|v| v.as_array())
                    .map(|arr| arr.iter().filter_map(|v| v.as_str().map(String::from)).collect());

                let symmetries: Option<Vec<String>> = input.parameters
                    .get("symmetries")
                    .and_then(|v| v.as_array())
                    .map(|arr| arr.iter().filter_map(|v| v.as_str().map(String::from)).collect());

                let units = if input.variable_units.is_empty() {
                    None
                } else {
                    Some(input.variable_units.clone())
                };

                let result = equation_validation::validate_equation(
                    input.expression.clone(),
                    domain,
                    units,
                    conservation_laws,
                    symmetries,
                ).map_err(|e| format!("Equation validation failed: {}", e))?;

                Ok(ValidateOutput {
                    is_valid: result.is_valid,
                    result: serde_json::json!({
                        "dimensional_consistency": result.dimensional_consistency,
                        "physics_compliance": result.physics_compliance,
                        "mathematical_correctness": result.mathematical_correctness,
                        "unit_analysis": result.unit_analysis,
                        "confidence": result.confidence
                    }),
                    errors: result.violations.clone(),
                    warnings: vec![],
                    suggestions: if result.violations.is_empty() {
                        None
                    } else {
                        Some(vec!["Review the violations and correct the equation".to_string()])
                    },
                    metadata: Some(serde_json::json!({
                        "operation": "equation",
                        "expression": input.expression
                    })),
                })
            }

            ValidateOp::Dimensions => {
                // Check dimensional consistency only
                let units_map: HashMap<String, String> = input.variable_units.clone();

                let (consistent, unit_analysis) = equation_validation::check_dimensional_consistency(
                    &input.expression,
                    &units_map,
                ).map_err(|e| format!("Dimensional check failed: {}", e))?;

                Ok(ValidateOutput {
                    is_valid: consistent,
                    result: serde_json::json!({
                        "consistent": consistent,
                        "unit_analysis": unit_analysis
                    }),
                    errors: if consistent {
                        vec![]
                    } else {
                        vec!["Dimensional inconsistency detected".to_string()]
                    },
                    warnings: vec![],
                    suggestions: None,
                    metadata: Some(serde_json::json!({
                        "operation": "dimensions",
                        "expression": input.expression
                    })),
                })
            }

            ValidateOp::Conservation => {
                // Check conservation laws
                let domain = input.parameters
                    .get("domain")
                    .and_then(|v| v.as_str())
                    .unwrap_or("general");

                let conservation_laws: Vec<String> = input.parameters
                    .get("conservation_laws")
                    .and_then(|v| v.as_array())
                    .map(|arr| arr.iter().filter_map(|v| v.as_str().map(String::from)).collect())
                    .unwrap_or_else(|| vec!["energy".to_string(), "momentum".to_string(), "charge".to_string()]);

                let violations = equation_validation::check_conservation_laws(
                    &input.expression,
                    domain,
                    &conservation_laws,
                ).map_err(|e| format!("Conservation check failed: {}", e))?;

                let is_valid = violations.is_empty();

                Ok(ValidateOutput {
                    is_valid,
                    result: serde_json::json!({
                        "laws_checked": conservation_laws,
                        "violations": violations
                    }),
                    errors: violations,
                    warnings: vec![],
                    suggestions: None,
                    metadata: Some(serde_json::json!({
                        "operation": "conservation",
                        "domain": domain
                    })),
                })
            }

            ValidateOp::Symmetry => {
                // Check symmetry properties
                let domain = input.parameters
                    .get("domain")
                    .and_then(|v| v.as_str())
                    .unwrap_or("general");

                let symmetries: Vec<String> = input.parameters
                    .get("symmetries")
                    .and_then(|v| v.as_array())
                    .map(|arr| arr.iter().filter_map(|v| v.as_str().map(String::from)).collect())
                    .unwrap_or_else(|| vec![
                        "time_translation".to_string(),
                        "space_translation".to_string(),
                        "rotation".to_string(),
                    ]);

                let violations = equation_validation::check_symmetries(
                    &input.expression,
                    domain,
                    &symmetries,
                ).map_err(|e| format!("Symmetry check failed: {}", e))?;

                let is_valid = violations.is_empty();

                Ok(ValidateOutput {
                    is_valid,
                    result: serde_json::json!({
                        "symmetries_checked": symmetries,
                        "violations": violations
                    }),
                    errors: violations,
                    warnings: vec![],
                    suggestions: None,
                    metadata: Some(serde_json::json!({
                        "operation": "symmetry",
                        "domain": domain
                    })),
                })
            }

            ValidateOp::Physics => {
                // Check physics compliance
                let domain = input.parameters
                    .get("domain")
                    .and_then(|v| v.as_str())
                    .unwrap_or("general");

                let (compliant, violations) = equation_validation::check_physics_compliance(
                    &input.expression,
                    domain,
                ).map_err(|e| format!("Physics compliance check failed: {}", e))?;

                Ok(ValidateOutput {
                    is_valid: compliant,
                    result: serde_json::json!({
                        "compliant": compliant,
                        "domain": domain,
                        "violations": violations
                    }),
                    errors: violations,
                    warnings: vec![],
                    suggestions: None,
                    metadata: Some(serde_json::json!({
                        "operation": "physics",
                        "domain": domain
                    })),
                })
            }

            ValidateOp::Bounds => {
                // Check mathematical bounds and correctness
                let is_correct = equation_validation::check_mathematical_correctness(&input.expression)
                    .map_err(|e| format!("Mathematical correctness check failed: {}", e))?;

                // Additional bounds checking
                let mut warnings = Vec::new();
                let mut errors = Vec::new();

                // Check for division by potential zeros
                if input.expression.contains("/") {
                    warnings.push("Expression contains division - ensure denominator is non-zero".to_string());
                }

                // Check for logarithms
                if input.expression.contains("log") || input.expression.contains("ln") {
                    warnings.push("Expression contains logarithm - ensure argument is positive".to_string());
                }

                // Check for square roots
                if input.expression.contains("sqrt") {
                    warnings.push("Expression contains square root - ensure argument is non-negative".to_string());
                }

                if !is_correct {
                    errors.push("Mathematical bounds or syntax error detected".to_string());
                }

                Ok(ValidateOutput {
                    is_valid: is_correct,
                    result: serde_json::json!({
                        "mathematically_correct": is_correct,
                        "has_division": input.expression.contains("/"),
                        "has_logarithm": input.expression.contains("log") || input.expression.contains("ln"),
                        "has_sqrt": input.expression.contains("sqrt")
                    }),
                    errors,
                    warnings,
                    suggestions: if is_correct { None } else {
                        Some(vec![
                            "Check parentheses matching".to_string(),
                            "Ensure no division by zero".to_string(),
                            "Ensure no sqrt of negative numbers".to_string(),
                        ])
                    },
                    metadata: Some(serde_json::json!({
                        "operation": "bounds"
                    })),
                })
            }

            ValidateOp::Singularities => {
                // Check for singularities in the expression
                let mut singularities = Vec::new();
                let mut warnings = Vec::new();

                // Parse the equation
                let (left, right) = equation_validation::parse_equation(&input.expression)
                    .map_err(|e| format!("Failed to parse equation: {}", e))?;

                // Check for division patterns that could lead to singularities
                let mut check_singularities = |expr: &str, singularities: &mut Vec<String>| {
                    // Division by variable
                    if expr.contains("/x") || expr.contains("/y") || expr.contains("/z") {
                        singularities.push("Potential singularity: division by variable".to_string());
                    }

                    // 1/0 patterns
                    if expr.contains("/0") {
                        singularities.push("Singularity: division by zero".to_string());
                    }

                    // 1/(x-a) patterns (poles)
                    if expr.contains("/(") {
                        singularities.push("Potential pole: division by expression".to_string());
                    }

                    // tan, cot, sec, csc have singularities
                    if expr.contains("tan(") {
                        singularities.push("Potential singularity: tan function has poles at odd multiples of pi/2".to_string());
                    }
                    if expr.contains("cot(") {
                        singularities.push("Potential singularity: cot function has poles at multiples of pi".to_string());
                    }
                    if expr.contains("sec(") {
                        singularities.push("Potential singularity: sec function has poles at odd multiples of pi/2".to_string());
                    }
                    if expr.contains("csc(") {
                        singularities.push("Potential singularity: csc function has poles at multiples of pi".to_string());
                    }

                    // Log at zero
                    if expr.contains("log(") || expr.contains("ln(") {
                        warnings.push("Logarithm undefined at zero".to_string());
                    }

                    // Gamma function singularities
                    if expr.contains("gamma(") {
                        singularities.push("Gamma function has poles at non-positive integers".to_string());
                    }
                };

                check_singularities(&left, &mut singularities);
                check_singularities(&right, &mut singularities);

                let is_valid = singularities.iter().all(|s| s.contains("Potential"));

                Ok(ValidateOutput {
                    is_valid,
                    result: serde_json::json!({
                        "singularities": singularities,
                        "left_side": left,
                        "right_side": right
                    }),
                    errors: singularities.iter()
                        .filter(|s| !s.contains("Potential"))
                        .cloned()
                        .collect(),
                    warnings: singularities.iter()
                        .filter(|s| s.contains("Potential"))
                        .cloned()
                        .collect(),
                    suggestions: if singularities.is_empty() {
                        None
                    } else {
                        Some(vec![
                            "Consider the domain of the function carefully".to_string(),
                            "Add constraints to avoid singular points".to_string(),
                        ])
                    },
                    metadata: Some(serde_json::json!({
                        "operation": "singularities"
                    })),
                })
            }
        }
    }
}
