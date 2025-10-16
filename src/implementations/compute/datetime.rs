//! DateTime computation operations
//!
//! Handles date and time calculations including interval arithmetic, date differences,
//! age calculations, business day operations, and calendar functions.

use crate::engine::*;

/// Compute datetime operations
pub fn compute_datetime(op: &DateTimeOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
    use crate::datetime::*;

    let dt_op = match op {
        DateTimeOp::AddInterval => DateTimeOperation::AddInterval,
        DateTimeOp::SubtractInterval => DateTimeOperation::SubtractInterval,
        DateTimeOp::DateDifference => DateTimeOperation::DateDifference,
        DateTimeOp::AgeCurrent => DateTimeOperation::AgeCurrent,
        DateTimeOp::AgeAtDate => DateTimeOperation::AgeAtDate,
        DateTimeOp::BusinessDays => DateTimeOperation::BusinessDays,
        DateTimeOp::AddBusinessDays => DateTimeOperation::AddBusinessDays,
        DateTimeOp::IsLeapYear => DateTimeOperation::IsLeapYear,
        DateTimeOp::DaysInMonth => DateTimeOperation::DaysInMonth,
        DateTimeOp::WeekNumber => DateTimeOperation::WeekNumber,
        DateTimeOp::DayOfWeek => DateTimeOperation::DayOfWeek,
    };

    let params: DateTimeParams = serde_json::from_value(
        serde_json::to_value(&input.parameters)
            .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
    )
    .map_err(|e| format!("Failed to parse datetime parameters: {}", e))?;

    let dt_input = DateTimeInput {
        operation: dt_op,
        parameters: params,
    };

    let result =
        calculate_datetime(dt_input).map_err(|e| format!("DateTime calculation error: {}", e))?;

    Ok(ComputeOutput {
        result: serde_json::json!({
            "value": result.value,
            "numeric_value": result.numeric_value,
            "unit": result.unit,
            "operation_used": result.operation_used,
            "interpretation": result.interpretation,
            "additional_info": result.additional_info
        }),
        additional: None,
        metadata: None,
    })
}
