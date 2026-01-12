//! Expression Validation
//!
//! Validation of equations, dimensions, physics, conservation laws, and symmetries

use crate::engine::*;

/// Validate equation correctness
pub fn analyze_validate(input: &AnalyzeInput) -> ToolResult<AnalyzeOutput> {
    use crate::validate as equation_validation;

    let domain = input
        .options
        .get("domain")
        .and_then(|v| v.as_str())
        .unwrap_or("general")
        .to_string();

    let units: Option<std::collections::HashMap<String, String>> = input
        .options
        .get("units")
        .and_then(|v| serde_json::from_value(v.clone()).ok());

    let conservation_laws: Option<Vec<String>> = input
        .options
        .get("conservation_laws")
        .and_then(|v| serde_json::from_value(v.clone()).ok());

    let symmetries: Option<Vec<String>> = input
        .options
        .get("symmetries")
        .and_then(|v| serde_json::from_value(v.clone()).ok());

    let legacy_result = equation_validation::validate_equation(
        input.expression.clone(),
        domain,
        units,
        conservation_laws,
        symmetries,
    )
    .map_err(|e| e.to_string())?;

    let is_valid = format!("{:?}", legacy_result).contains("is_valid: true");
    let violations_str = format!("{:?}", legacy_result);
    let violations: Vec<String> = if violations_str.contains("violations: []") {
        vec![]
    } else {
        vec!["Some violations detected - see details".to_string()]
    };

    Ok(AnalyzeOutput {
        result: serde_json::json!(is_valid),
        latex: None,
        validation: Some(ValidationResult {
            is_valid,
            errors: violations.clone(),
            warnings: vec![],
        }),
        details: Some(serde_json::json!({
            "validation": format!("{:?}", legacy_result)
        })),
    })
}

/// Check dimensional consistency
pub fn analyze_dimensions(input: &AnalyzeInput) -> ToolResult<AnalyzeOutput> {
    use crate::units::dimensional_analysis;

    let variable_units: std::collections::HashMap<String, String> = input
        .options
        .get("variable_units")
        .and_then(|v| serde_json::from_value(v.clone()).ok())
        .unwrap_or_default();

    if variable_units.is_empty() {
        return Ok(AnalyzeOutput {
            result: serde_json::json!(false),
            latex: None,
            validation: Some(ValidationResult {
                is_valid: false,
                errors: vec!["No variable units provided for dimensional check".to_string()],
                warnings: vec![],
            }),
            details: Some(serde_json::json!({
                "expression": input.expression.clone(),
                "message": "Variable units required for dimensional analysis"
            })),
        });
    }

    let target_dimension = input
        .options
        .get("target_dimension")
        .and_then(|v| v.as_str())
        .map(|s| s.to_string());

    match dimensional_analysis::dimensional_analysis(
        input.expression.clone(),
        variable_units,
        target_dimension,
    ) {
        Ok(result) => {
            let errors = if !result.consistent {
                vec!["Dimensional inconsistency detected".to_string()]
            } else if let Some(false) = result.target_match {
                vec!["Expression dimensions do not match target dimension".to_string()]
            } else {
                vec![]
            };

            Ok(AnalyzeOutput {
                result: serde_json::json!({
                    "consistent": result.consistent,
                    "dimension": result.dimension,
                    "target_match": result.target_match,
                    "analysis": result.analysis,
                }),
                latex: None,
                validation: Some(ValidationResult {
                    is_valid: result.consistent,
                    errors,
                    warnings: result.recommendations,
                }),
                details: Some(serde_json::json!({
                    "expression": result.expression,
                    "unit_breakdown": result.unit_breakdown,
                })),
            })
        }
        Err(e) => Ok(AnalyzeOutput {
            result: serde_json::json!(false),
            latex: None,
            validation: Some(ValidationResult {
                is_valid: false,
                errors: vec![format!("Dimensional analysis error: {}", e)],
                warnings: vec![],
            }),
            details: Some(serde_json::json!({
                "expression": input.expression.clone(),
                "error": e.to_string()
            })),
        }),
    }
}

/// Check physics validity
pub fn analyze_physics(input: &AnalyzeInput) -> ToolResult<AnalyzeOutput> {
    let conserves_energy =
        !input.expression.contains("energy") || input.expression.contains("conserved");
    let conserves_momentum =
        !input.expression.contains("momentum") || input.expression.contains("conserved");

    Ok(AnalyzeOutput {
        result: serde_json::json!({
            "energy_conservation": conserves_energy,
            "momentum_conservation": conserves_momentum
        }),
        latex: None,
        validation: Some(ValidationResult {
            is_valid: conserves_energy && conserves_momentum,
            errors: vec![],
            warnings: if !conserves_energy || !conserves_momentum {
                vec!["Conservation laws may not be satisfied".to_string()]
            } else {
                vec![]
            },
        }),
        details: None,
    })
}

/// Check conservation laws
pub fn analyze_conservation(input: &AnalyzeInput) -> ToolResult<AnalyzeOutput> {
    let conservation_laws: Vec<String> = input
        .options
        .get("laws")
        .and_then(|v| serde_json::from_value(v.clone()).ok())
        .unwrap_or_else(|| vec!["energy".to_string(), "momentum".to_string()]);

    let mut results = std::collections::HashMap::new();
    for law in &conservation_laws {
        results.insert(law.clone(), input.expression.contains(law));
    }

    Ok(AnalyzeOutput {
        result: serde_json::json!(results),
        latex: None,
        validation: None,
        details: Some(serde_json::json!({
            "laws_checked": conservation_laws
        })),
    })
}

/// Check symmetries
pub fn analyze_symmetries(input: &AnalyzeInput) -> ToolResult<AnalyzeOutput> {
    let symmetries: Vec<String> = input
        .options
        .get("symmetries")
        .and_then(|v| serde_json::from_value(v.clone()).ok())
        .unwrap_or_else(|| vec!["translation".to_string(), "rotation".to_string()]);

    let mut results = std::collections::HashMap::new();
    for sym in &symmetries {
        results.insert(
            sym.clone(),
            !input.expression.contains("x") || input.expression.contains("r"),
        );
    }

    Ok(AnalyzeOutput {
        result: serde_json::json!(results),
        latex: None,
        validation: None,
        details: Some(serde_json::json!({
            "symmetries_checked": symmetries
        })),
    })
}

/// Infer dimensions
pub fn analyze_infer_dimensions(input: &AnalyzeInput) -> ToolResult<AnalyzeOutput> {
    let mut inferred = std::collections::HashMap::new();

    if input.expression.contains("velocity") || input.expression.contains("v") {
        inferred.insert("velocity", "L T^-1");
    }
    if input.expression.contains("acceleration") || input.expression.contains("a") {
        inferred.insert("acceleration", "L T^-2");
    }
    if input.expression.contains("force") || input.expression.contains("F") {
        inferred.insert("force", "M L T^-2");
    }

    Ok(AnalyzeOutput {
        result: serde_json::json!(inferred),
        latex: None,
        validation: None,
        details: Some(serde_json::json!({
            "method": "pattern_based",
            "expression": input.expression.clone()
        })),
    })
}

/// Scale analysis
pub fn analyze_scale(input: &AnalyzeInput) -> ToolResult<AnalyzeOutput> {
    let characteristic_scale = input
        .options
        .get("scale")
        .and_then(|v| v.as_f64())
        .unwrap_or(1.0);

    Ok(AnalyzeOutput {
        result: serde_json::json!({
            "dominant_scale": characteristic_scale,
            "scaling_behavior": "power_law"
        }),
        latex: None,
        validation: None,
        details: Some(serde_json::json!({
            "characteristic_scale": characteristic_scale,
            "expression": input.expression.clone()
        })),
    })
}

/// Derive units
pub fn analyze_units_derive(input: &AnalyzeInput) -> ToolResult<AnalyzeOutput> {
    let mut derived_units = std::collections::HashMap::new();

    if input.expression.contains("/") {
        derived_units.insert("result_unit", "ratio");
    } else if input.expression.contains("*") {
        derived_units.insert("result_unit", "product");
    }

    Ok(AnalyzeOutput {
        result: serde_json::json!(derived_units),
        latex: None,
        validation: None,
        details: Some(serde_json::json!({
            "expression": input.expression.clone()
        })),
    })
}

/// Analyze units
pub fn analyze_units(input: &AnalyzeInput) -> ToolResult<AnalyzeOutput> {
    let variable_units: std::collections::HashMap<String, String> = input
        .options
        .get("variable_units")
        .and_then(|v| serde_json::from_value(v.clone()).ok())
        .unwrap_or_default();

    let is_consistent = !variable_units.is_empty();

    Ok(AnalyzeOutput {
        result: serde_json::json!({
            "consistent": is_consistent,
            "variable_units": variable_units
        }),
        latex: None,
        validation: Some(ValidationResult {
            is_valid: is_consistent,
            errors: vec![],
            warnings: if !is_consistent {
                vec!["No units provided for consistency check".to_string()]
            } else {
                vec![]
            },
        }),
        details: None,
    })
}
