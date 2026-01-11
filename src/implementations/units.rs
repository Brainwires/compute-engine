//! Unified Units implementation
//!
//! Routes unit and dimensional analysis requests to tools/dimensional_analysis/

use crate::engine::*;

pub struct UnifiedUnits;

impl UnifiedUnits {
    pub fn new() -> Self {
        Self
    }
}

impl Default for UnifiedUnits {
    fn default() -> Self {
        Self::new()
    }
}

impl Units for UnifiedUnits {
    fn units(&self, input: &UnitsInput) -> ToolResult<UnitsOutput> {
        use crate::tools::dimensional_analysis;

        match &input.operation {
            UnitsOp::Convert => {
                // Extract parameters for unit conversion
                let value = input.value.unwrap_or(1.0);
                let from_unit = input.from_unit.as_ref()
                    .ok_or("from_unit required for Convert operation")?;
                let to_unit = input.to_unit.as_ref()
                    .ok_or("to_unit required for Convert operation")?;

                let result = dimensional_analysis::convert_units(from_unit, to_unit, value)
                    .map_err(|e| format!("Unit conversion failed: {}", e))?;

                let dims = std::collections::HashMap::new();
                // Parse dimension string to extract powers
                // (simplified - the actual dimension is a string like "M⋅L⋅T^-2")

                Ok(UnitsOutput {
                    result: serde_json::json!({
                        "converted_value": result.factor,
                        "from": from_unit,
                        "to": to_unit,
                        "original_value": value,
                        "dimension": result.dimension
                    }),
                    converted_value: Some(result.factor),
                    si_units: Some(result.dimension.clone()),
                    dimensions: Some(dims),
                    compatible: Some(true),
                    metadata: Some(serde_json::json!({
                        "operation": "convert",
                        "from_unit": from_unit,
                        "to_unit": to_unit
                    })),
                })
            }

            UnitsOp::Analyze => {
                // Full dimensional analysis of an expression
                let expression = input.expression.as_ref()
                    .ok_or("expression required for Analyze operation")?;

                let variable_units: std::collections::HashMap<String, String> = input.variable_units
                    .iter()
                    .map(|(k, v)| (k.clone(), v.clone()))
                    .collect();

                let target = input.parameters
                    .get("target")
                    .and_then(|v| v.as_str())
                    .map(String::from);

                let result = dimensional_analysis::dimensional_analysis(
                    expression.clone(),
                    variable_units,
                    target,
                ).map_err(|e| format!("Dimensional analysis failed: {}", e))?;

                Ok(UnitsOutput {
                    result: serde_json::json!({
                        "expression": result.expression,
                        "dimension": result.dimension,
                        "consistent": result.consistent,
                        "unit_breakdown": result.unit_breakdown,
                        "analysis": result.analysis,
                        "recommendations": result.recommendations,
                        "target_match": result.target_match
                    }),
                    converted_value: None,
                    si_units: Some(result.dimension),
                    dimensions: None,
                    compatible: Some(result.consistent),
                    metadata: Some(serde_json::json!({
                        "operation": "analyze",
                        "expression": expression
                    })),
                })
            }

            UnitsOp::CheckCompatibility => {
                // Check if two units are compatible
                let from_unit = input.from_unit.as_ref()
                    .ok_or("from_unit required for CheckCompatibility operation")?;
                let to_unit = input.to_unit.as_ref()
                    .ok_or("to_unit required for CheckCompatibility operation")?;

                let compatible = dimensional_analysis::are_units_compatible(from_unit, to_unit)
                    .map_err(|e| format!("Compatibility check failed: {}", e))?;

                Ok(UnitsOutput {
                    result: serde_json::json!({
                        "compatible": compatible,
                        "unit1": from_unit,
                        "unit2": to_unit
                    }),
                    converted_value: None,
                    si_units: None,
                    dimensions: None,
                    compatible: Some(compatible),
                    metadata: Some(serde_json::json!({
                        "operation": "check_compatibility",
                        "units": [from_unit, to_unit]
                    })),
                })
            }

            UnitsOp::GetBase => {
                // Get SI base units for a given unit
                let unit = input.from_unit.as_ref()
                    .ok_or("from_unit required for GetBase operation")?;

                let si_base = dimensional_analysis::get_si_base_units(unit)
                    .map_err(|e| format!("Get SI base units failed: {}", e))?;

                Ok(UnitsOutput {
                    result: serde_json::json!({
                        "unit": unit,
                        "si_base": si_base
                    }),
                    converted_value: None,
                    si_units: Some(si_base.clone()),
                    dimensions: None,
                    compatible: None,
                    metadata: Some(serde_json::json!({
                        "operation": "get_base",
                        "input_unit": unit
                    })),
                })
            }

            UnitsOp::Derive => {
                // Derive units for a physical quantity
                let quantity = input.parameters
                    .get("quantity")
                    .and_then(|v| v.as_str())
                    .or_else(|| input.expression.as_deref())
                    .ok_or("quantity parameter or expression required for Derive operation")?;

                let units = dimensional_analysis::derive_units_for_quantity(quantity)
                    .map_err(|e| format!("Derive units failed: {}", e))?;

                Ok(UnitsOutput {
                    result: serde_json::json!({
                        "quantity": quantity,
                        "units": units
                    }),
                    converted_value: None,
                    si_units: units.first().cloned(),
                    dimensions: None,
                    compatible: None,
                    metadata: Some(serde_json::json!({
                        "operation": "derive",
                        "quantity": quantity
                    })),
                })
            }

            UnitsOp::Parse => {
                // Parse a unit string and extract its components
                let unit = input.from_unit.as_ref()
                    .or_else(|| input.expression.as_ref())
                    .ok_or("from_unit or expression required for Parse operation")?;

                let tokens = dimensional_analysis::tokenize_unit(unit);
                let si_base = dimensional_analysis::get_si_base_units(unit)
                    .unwrap_or_else(|_| "unknown".to_string());

                Ok(UnitsOutput {
                    result: serde_json::json!({
                        "unit": unit,
                        "tokens": tokens,
                        "si_base": si_base
                    }),
                    converted_value: None,
                    si_units: Some(si_base),
                    dimensions: None,
                    compatible: None,
                    metadata: Some(serde_json::json!({
                        "operation": "parse",
                        "input": unit
                    })),
                })
            }

            UnitsOp::Simplify => {
                // Simplify a unit expression
                let unit = input.from_unit.as_ref()
                    .or_else(|| input.expression.as_ref())
                    .ok_or("from_unit or expression required for Simplify operation")?;

                let si_base = dimensional_analysis::get_si_base_units(unit)
                    .map_err(|e| format!("Simplify failed: {}", e))?;

                Ok(UnitsOutput {
                    result: serde_json::json!({
                        "original": unit,
                        "simplified": si_base,
                        "is_simplified": unit == &si_base
                    }),
                    converted_value: None,
                    si_units: Some(si_base),
                    dimensions: None,
                    compatible: None,
                    metadata: Some(serde_json::json!({
                        "operation": "simplify"
                    })),
                })
            }
        }
    }
}
