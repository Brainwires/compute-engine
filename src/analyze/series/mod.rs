//! Series Expansions
//!
//! Taylor series, Laurent series, and limit computations

use crate::engine::*;

/// Series expansion (Taylor/Maclaurin)
pub fn analyze_series_expansion(input: &AnalyzeInput) -> ToolResult<AnalyzeOutput> {
    use crate::analyze::symbolic as symbolic_cas;

    let order = input
        .options
        .get("order")
        .and_then(|v| v.as_u64())
        .unwrap_or(3) as usize;
    let variable = input
        .options
        .get("variable")
        .and_then(|v| v.as_str())
        .unwrap_or("x");
    let point = input
        .options
        .get("point")
        .and_then(|v| v.as_f64())
        .unwrap_or(0.0);

    let result = symbolic_cas::series_expansion(&input.expression, variable, point, order)
        .map_err(|e| format!("Series expansion failed: {}", e))?;

    Ok(AnalyzeOutput {
        result: serde_json::json!({
            "expansion": result.expression,
            "order": order,
            "variable": variable,
            "point": point
        }),
        latex: result.latex.clone(),
        validation: None,
        details: Some(serde_json::json!({
            "expansion_type": "taylor",
            "around": point,
            "metadata": result.metadata
        })),
    })
}

/// Laurent series expansion
pub fn analyze_laurent_series(input: &AnalyzeInput) -> ToolResult<AnalyzeOutput> {
    use crate::analyze::symbolic as symbolic_cas;

    let order = input
        .options
        .get("order")
        .and_then(|v| v.as_i64())
        .unwrap_or(3) as usize;

    let variable = input
        .options
        .get("variable")
        .and_then(|v| v.as_str())
        .unwrap_or("z");

    let point = input
        .options
        .get("point")
        .and_then(|v| v.as_f64())
        .unwrap_or(0.0);

    let taylor_part = symbolic_cas::series_expansion(&input.expression, variable, point, order)
        .map_err(|e| format!("Laurent series (analytic part) failed: {}", e))?;

    let mut eval_map = std::collections::HashMap::new();
    eval_map.insert(variable.to_string(), point);
    let has_singularity = symbolic_cas::evaluate_at(&input.expression, &eval_map).is_err();

    Ok(AnalyzeOutput {
        result: serde_json::json!({
            "analytic_part": taylor_part.expression,
            "expansion_point": point,
            "order": order,
            "variable": variable,
            "has_singularity_at_point": has_singularity,
            "note": "Full Laurent series with principal part (negative powers) computation in progress"
        }),
        latex: Some(format!(
            "$\\sum_{{n=-{}}}^{{{}}} a_n ({}-{})^n$",
            order, order, variable, point
        )),
        validation: None,
        details: Some(serde_json::json!({
            "expansion_type": "laurent",
            "singularity": point,
            "analytic_part_metadata": taylor_part.metadata,
            "method": "taylor_based_with_singularity_detection"
        })),
    })
}

/// Compute limit using symbolic CAS
pub fn analyze_limit(input: &AnalyzeInput) -> ToolResult<AnalyzeOutput> {
    use crate::analyze::symbolic as symbolic_cas;

    let point = input
        .options
        .get("point")
        .and_then(|v| v.as_str())
        .unwrap_or("0");
    let variable = input
        .options
        .get("variable")
        .and_then(|v| v.as_str())
        .unwrap_or("x");
    let direction = input.options.get("direction").and_then(|v| v.as_str());

    let result = symbolic_cas::limit(&input.expression, variable, point, direction)
        .map_err(|e| format!("Limit computation failed: {}", e))?;

    Ok(AnalyzeOutput {
        result: serde_json::json!({
            "limit_value": result.expression,
            "variable": variable,
            "point": point,
            "direction": direction.unwrap_or("both")
        }),
        latex: result.latex.clone(),
        validation: None,
        details: Some(serde_json::json!({
            "point": point,
            "direction": direction,
            "metadata": result.metadata
        })),
    })
}
