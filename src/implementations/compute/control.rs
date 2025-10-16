//! Control systems computation operations
//!
//! This module handles control theory and systems engineering calculations including:
//! - Transfer functions and frequency response
//! - Stability analysis (poles, zeros, margins)
//! - State-space representations
//! - Time-domain analysis (step response, impulse response)

use crate::engine::*;

/// Compute control systems operations
///
/// Handles control systems analysis and design:
///
/// **Transfer Functions:**
/// - Evaluation at specific frequencies
/// - Series and parallel combination
/// - Frequency response (magnitude and phase)
///
/// **Frequency Domain Analysis:**
/// - Bode plots: Magnitude and phase vs frequency
/// - Nyquist plots: Complex frequency response
/// - Root locus: Pole movement with varying gain
///
/// **Stability Analysis:**
/// - Pole-zero analysis: System poles and zeros
/// - Routh-Hurwitz criterion: Stability test
/// - Gain margin: How much gain increase until instability
/// - Phase margin: Phase shift tolerance before instability
///
/// **State-Space Methods:**
/// - State-space representation: A, B, C, D matrices
/// - Controllability: Can all states be controlled?
/// - Observability: Can all states be observed?
/// - Conversion to/from transfer function
///
/// **Time-Domain Analysis:**
/// - Step response: System response to unit step input
/// - Impulse response: System response to delta function
pub fn compute_control_systems(
    op: &ControlSystemsOp,
    input: &ComputeInput,
) -> ToolResult<ComputeOutput> {
    use crate::physics::control_systems::*;

    let result_json = match op {
        ControlSystemsOp::TransferFunction => {
            let numerator: Vec<f64> = input
                .parameters
                .get("numerator")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("numerator coefficients required")?;
            let denominator: Vec<f64> = input
                .parameters
                .get("denominator")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("denominator coefficients required")?;
            let operation = input
                .parameters
                .get("operation")
                .and_then(|v| v.as_str())
                .unwrap_or("evaluate")
                .to_string();
            let frequency = input.parameters.get("frequency").and_then(|v| v.as_f64());
            let second_tf: Option<(Vec<f64>, Vec<f64>)> = input
                .parameters
                .get("second_tf")
                .and_then(|v| serde_json::from_value(v.clone()).ok());

            let result = transfer_function(TransferFunctionRequest {
                numerator,
                denominator,
                operation,
                frequency,
                second_tf,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        ControlSystemsOp::PoleZeroAnalysis => {
            let numerator: Vec<f64> = input
                .parameters
                .get("numerator")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("numerator required")?;
            let denominator: Vec<f64> = input
                .parameters
                .get("denominator")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("denominator required")?;

            let result = pole_zero_analysis(PoleZeroRequest {
                numerator,
                denominator,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        ControlSystemsOp::BodePlot => {
            let numerator: Vec<f64> = input
                .parameters
                .get("numerator")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("numerator required")?;
            let denominator: Vec<f64> = input
                .parameters
                .get("denominator")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("denominator required")?;
            let freq_min = input
                .parameters
                .get("freq_min")
                .and_then(|v| v.as_f64())
                .unwrap_or(0.1);
            let freq_max = input
                .parameters
                .get("freq_max")
                .and_then(|v| v.as_f64())
                .unwrap_or(100.0);
            let num_points = input
                .parameters
                .get("num_points")
                .and_then(|v| v.as_u64())
                .unwrap_or(50) as usize;

            let result = bode_plot(BodePlotRequest {
                numerator,
                denominator,
                freq_min,
                freq_max,
                num_points,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        ControlSystemsOp::NyquistPlot => {
            let numerator: Vec<f64> = input
                .parameters
                .get("numerator")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("numerator required")?;
            let denominator: Vec<f64> = input
                .parameters
                .get("denominator")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("denominator required")?;
            let freq_min = input
                .parameters
                .get("freq_min")
                .and_then(|v| v.as_f64())
                .unwrap_or(0.1);
            let freq_max = input
                .parameters
                .get("freq_max")
                .and_then(|v| v.as_f64())
                .unwrap_or(100.0);
            let num_points = input
                .parameters
                .get("num_points")
                .and_then(|v| v.as_u64())
                .unwrap_or(50) as usize;

            let result = nyquist_plot(NyquistPlotRequest {
                numerator,
                denominator,
                freq_min,
                freq_max,
                num_points,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        ControlSystemsOp::RootLocus => {
            let numerator: Vec<f64> = input
                .parameters
                .get("numerator")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("numerator required")?;
            let denominator: Vec<f64> = input
                .parameters
                .get("denominator")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("denominator required")?;
            let gain_min = input
                .parameters
                .get("gain_min")
                .and_then(|v| v.as_f64())
                .unwrap_or(0.0);
            let gain_max = input
                .parameters
                .get("gain_max")
                .and_then(|v| v.as_f64())
                .unwrap_or(10.0);
            let num_points = input
                .parameters
                .get("num_points")
                .and_then(|v| v.as_u64())
                .unwrap_or(50) as usize;

            let result = root_locus(RootLocusRequest {
                numerator,
                denominator,
                gain_min,
                gain_max,
                num_points,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        ControlSystemsOp::StateSpace => {
            let a_matrix: Vec<Vec<f64>> = input
                .parameters
                .get("a_matrix")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("a_matrix required")?;
            let b_matrix: Vec<Vec<f64>> = input
                .parameters
                .get("b_matrix")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("b_matrix required")?;
            let c_matrix: Vec<Vec<f64>> = input
                .parameters
                .get("c_matrix")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("c_matrix required")?;
            let d_matrix: Vec<Vec<f64>> = input
                .parameters
                .get("d_matrix")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("d_matrix required")?;
            let operation = input
                .parameters
                .get("operation")
                .and_then(|v| v.as_str())
                .unwrap_or("to_transfer_function")
                .to_string();
            let time = input.parameters.get("time").and_then(|v| v.as_f64());

            let result = state_space(StateSpaceRequest {
                a_matrix,
                b_matrix,
                c_matrix,
                d_matrix,
                operation,
                time,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        ControlSystemsOp::Controllability => {
            let a_matrix: Vec<Vec<f64>> = input
                .parameters
                .get("a_matrix")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("a_matrix required")?;
            let b_matrix: Vec<Vec<f64>> = input
                .parameters
                .get("b_matrix")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("b_matrix required")?;

            let result = controllability(ControllabilityRequest { a_matrix, b_matrix })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        ControlSystemsOp::Observability => {
            let a_matrix: Vec<Vec<f64>> = input
                .parameters
                .get("a_matrix")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("a_matrix required")?;
            let c_matrix: Vec<Vec<f64>> = input
                .parameters
                .get("c_matrix")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("c_matrix required")?;

            let result = observability(ObservabilityRequest { a_matrix, c_matrix })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        ControlSystemsOp::RouthHurwitz => {
            let characteristic_polynomial: Vec<f64> = input
                .parameters
                .get("characteristic_polynomial")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("characteristic_polynomial required")?;

            let result = routh_hurwitz(RouthHurwitzRequest {
                characteristic_polynomial,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        ControlSystemsOp::GainMargin => {
            let numerator: Vec<f64> = input
                .parameters
                .get("numerator")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("numerator required")?;
            let denominator: Vec<f64> = input
                .parameters
                .get("denominator")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("denominator required")?;

            let result = gain_margin(GainMarginRequest {
                numerator,
                denominator,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        ControlSystemsOp::PhaseMargin => {
            let numerator: Vec<f64> = input
                .parameters
                .get("numerator")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("numerator required")?;
            let denominator: Vec<f64> = input
                .parameters
                .get("denominator")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("denominator required")?;

            let result = phase_margin(PhaseMarginRequest {
                numerator,
                denominator,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        ControlSystemsOp::StepResponse => {
            let numerator: Vec<f64> = input
                .parameters
                .get("numerator")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("numerator required")?;
            let denominator: Vec<f64> = input
                .parameters
                .get("denominator")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("denominator required")?;
            let time_span = input
                .parameters
                .get("time_span")
                .and_then(|v| v.as_f64())
                .unwrap_or(10.0);
            let num_points = input
                .parameters
                .get("num_points")
                .and_then(|v| v.as_u64())
                .unwrap_or(100) as usize;

            let result = step_response(StepResponseRequest {
                numerator,
                denominator,
                time_span,
                num_points,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
    };

    Ok(ComputeOutput {
        result: result_json,
        additional: None,
        metadata: None,
    })
}
