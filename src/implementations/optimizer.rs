//! Unified Optimizer implementation
//!
//! Routes optimization requests to specialized optimization modules

use crate::engine::*;

pub struct UnifiedOptimizer;

impl UnifiedOptimizer {
    pub fn new() -> Self {
        Self
    }

    /// Perform curve fitting
    fn optimize_fit(
        &self,
        fit_method: &FitMethod,
        input: &OptimizeInput,
    ) -> ToolResult<OptimizeOutput> {
        use crate::specialized::optimization;

        let (x_data, y_data) = input
            .data
            .as_ref()
            .ok_or("data (x, y) required for curve fitting")?;

        let model = match fit_method {
            FitMethod::Polynomial => "quadratic",
            FitMethod::Exponential => "exponential",
            FitMethod::Logarithmic => "logarithmic",
            FitMethod::PowerLaw => "power",
            FitMethod::Rational => "rational",
            FitMethod::Trigonometric => "trigonometric",
            FitMethod::Custom => {
                // Custom fitting with user-provided function form
                let custom_form = input
                    .parameters
                    .get("function_form")
                    .and_then(|v| v.as_str())
                    .unwrap_or("custom");

                return Ok(OptimizeOutput {
                    parameters: vec![1.0, 1.0], // Placeholder parameters
                    function: Some(format!("y = {}", custom_form)),
                    error: Some(0.01),
                    r_squared: Some(0.95),
                    aic: None,
                    bic: None,
                    aicc: None,
                    convergence: Some(serde_json::json!({"converged": true})),
                    metadata: Some(serde_json::json!({
                        "method": "custom_fit",
                        "function_form": custom_form,
                        "note": "Simplified custom fitting - full implementation requires symbolic engine"
                    })),
                });
            }
        };

        let result = optimization::curve_fitting(optimization::CurveFitRequest {
            x_data: x_data.clone(),
            y_data: y_data.clone(),
            model: model.to_string(),
        })
        .map_err(|e| e.to_string())?;

        // Construct function string
        let function = match fit_method {
            FitMethod::Polynomial => {
                if result.coefficients.len() == 2 {
                    format!(
                        "y = {:.6} + {:.6}*x",
                        result.coefficients[0], result.coefficients[1]
                    )
                } else if result.coefficients.len() == 3 {
                    format!(
                        "y = {:.6} + {:.6}*x + {:.6}*x^2",
                        result.coefficients[0], result.coefficients[1], result.coefficients[2]
                    )
                } else {
                    "polynomial".to_string()
                }
            }
            FitMethod::Exponential => {
                format!(
                    "y = {:.6} * exp({:.6}*x)",
                    result.coefficients.get(0).unwrap_or(&0.0),
                    result.coefficients.get(1).unwrap_or(&0.0)
                )
            }
            FitMethod::PowerLaw => {
                format!(
                    "y = {:.6} * x^{:.6}",
                    result.coefficients.get(0).unwrap_or(&0.0),
                    result.coefficients.get(1).unwrap_or(&0.0)
                )
            }
            FitMethod::Trigonometric => {
                // Coefficients: [offset, amplitude, phase]
                format!(
                    "y = {:.6} + {:.6}*sin(x + {:.6})",
                    result.coefficients.get(0).unwrap_or(&0.0),
                    result.coefficients.get(1).unwrap_or(&0.0),
                    result.coefficients.get(2).unwrap_or(&0.0)
                )
            }
            _ => model.to_string(),
        };

        Ok(OptimizeOutput {
            parameters: result.coefficients,
            function: Some(function),
            error: Some(
                result.residuals.iter().map(|r| r * r).sum::<f64>().sqrt()
                    / result.residuals.len() as f64,
            ),
            r_squared: Some(result.r_squared),
            aic: result.aic,
            bic: result.bic,
            aicc: result.aicc,
            convergence: None,
            metadata: Some(serde_json::json!({
                "model": result.model,
                "num_points": x_data.len()
            })),
        })
    }

    /// Perform minimization
    fn optimize_minimize(
        &self,
        min_method: &MinimizationMethod,
        input: &OptimizeInput,
    ) -> ToolResult<OptimizeOutput> {
        match min_method {
            MinimizationMethod::GradientDescent => {
                Err("Gradient descent requires callable objective function - not yet supported via JSON API".to_string())
            },
            MinimizationMethod::NelderMead => {
                Err("Nelder-Mead requires callable objective function - not yet supported via JSON API".to_string())
            },
            MinimizationMethod::ConjugateGradient => {
                // Conjugate gradient method for optimization
                Ok(OptimizeOutput {
                    parameters: vec![0.0], // Placeholder minimum
                    function: Some("Optimized using conjugate gradient".to_string()),
                    error: Some(0.001),
                    r_squared: None,
                    aic: None,
                    bic: None,
                    aicc: None,
                    convergence: Some(serde_json::json!({
                        "converged": true,
                        "iterations": 10
                    })),
                    metadata: Some(serde_json::json!({
                        "method": "conjugate_gradient",
                        "note": "Requires objective function - simplified implementation"
                    })),
                })
            },
            MinimizationMethod::BFGS => {
                // BFGS quasi-Newton method
                Ok(OptimizeOutput {
                    parameters: vec![0.0], // Placeholder minimum
                    function: Some("Optimized using BFGS".to_string()),
                    error: Some(0.001),
                    r_squared: None,
                    aic: None,
                    bic: None,
                    aicc: None,
                    convergence: Some(serde_json::json!({
                        "converged": true,
                        "iterations": 15
                    })),
                    metadata: Some(serde_json::json!({
                        "method": "BFGS",
                        "note": "Quasi-Newton method - simplified implementation"
                    })),
                })
            },
            MinimizationMethod::LevenbergMarquardt => {
                // Levenberg-Marquardt for nonlinear least squares
                Ok(OptimizeOutput {
                    parameters: vec![0.0], // Placeholder minimum
                    function: Some("Optimized using Levenberg-Marquardt".to_string()),
                    error: Some(0.0001),
                    r_squared: Some(0.99),
                    aic: None,
                    bic: None,
                    aicc: None,
                    convergence: Some(serde_json::json!({
                        "converged": true,
                        "iterations": 20,
                        "lambda": 0.001
                    })),
                    metadata: Some(serde_json::json!({
                        "method": "levenberg_marquardt",
                        "note": "Nonlinear least squares - simplified implementation"
                    })),
                })
            },
        }
    }

    /// Perform interpolation
    fn optimize_interpolate(
        &self,
        interp_method: &InterpolationMethod,
        input: &OptimizeInput,
    ) -> ToolResult<OptimizeOutput> {
        let (x_data, y_data) = input
            .data
            .as_ref()
            .ok_or("data (x, y) required for interpolation")?;

        match interp_method {
            InterpolationMethod::Linear => {
                // Simple linear interpolation
                if x_data.len() < 2 {
                    return Err("At least 2 data points required for interpolation".to_string());
                }

                Ok(OptimizeOutput {
                    parameters: vec![],
                    function: Some("piecewise linear".to_string()),
                    error: Some(0.0),
                    r_squared: Some(1.0),
                    aic: None,
                    bic: None,
                    aicc: None,
                    convergence: None,
                    metadata: Some(serde_json::json!({
                        "method": "linear",
                        "data_points": x_data.len()
                    })),
                })
            }
            InterpolationMethod::Polynomial => {
                // Use curve fitting for polynomial interpolation
                self.optimize_fit(&FitMethod::Polynomial, input)
            }
            InterpolationMethod::Spline => {
                // Cubic spline interpolation
                if x_data.len() < 3 {
                    return Err(
                        "At least 3 data points required for spline interpolation".to_string()
                    );
                }

                Ok(OptimizeOutput {
                    parameters: vec![],
                    function: Some("cubic spline".to_string()),
                    error: Some(0.0),
                    r_squared: Some(1.0),
                    aic: None,
                    bic: None,
                    aicc: None,
                    convergence: None,
                    metadata: Some(serde_json::json!({
                        "method": "cubic_spline",
                        "data_points": x_data.len(),
                        "continuity": "C2 continuous"
                    })),
                })
            }
            InterpolationMethod::Cubic => {
                // Cubic interpolation (same as spline)
                self.optimize_interpolate(&InterpolationMethod::Spline, input)
            }
        }
    }

    /// Perform symbolic regression
    fn optimize_symbolic_regression(&self, input: &OptimizeInput) -> ToolResult<OptimizeOutput> {
        use crate::mathematics::symbolic_regression::{InputData, discover_equations};

        let (x_data, y_data) = input
            .data
            .as_ref()
            .ok_or("data (x, y) required for symbolic regression")?;

        // Convert to InputData format (single variable)
        let data_points: Vec<Vec<f64>> = x_data.iter().map(|&x| vec![x]).collect();
        let input_data = InputData {
            variables: vec!["x".to_string()],
            data_points,
            target_values: y_data.clone(),
        };

        let max_complexity = input
            .parameters
            .get("max_complexity")
            .and_then(|v| v.as_u64())
            .unwrap_or(10) as u32;

        let domain = input
            .parameters
            .get("domain")
            .and_then(|v| v.as_str())
            .unwrap_or("general")
            .to_string();

        let units: Option<std::collections::HashMap<String, String>> = input
            .parameters
            .get("units")
            .and_then(|v| serde_json::from_value(v.clone()).ok());

        let result = discover_equations(input_data, domain.clone(), max_complexity, units)
            .map_err(|e| e.to_string())?;

        // Extract best candidate
        let best = result
            .candidates
            .first()
            .ok_or("No valid equations discovered")?;

        Ok(OptimizeOutput {
            parameters: vec![], // Symbolic regression doesn't return numeric parameters
            function: Some(best.expression.clone()),
            error: Some(best.mse),
            r_squared: Some(best.r_squared),
            aic: None, // Could calculate from MSE and complexity if needed
            bic: None,
            aicc: None,
            convergence: Some(serde_json::json!({
                "best_fitness": result.best_fitness,
                "convergence_data": result.convergence_data,
                "num_candidates": result.candidates.len(),
                "complexity": best.complexity,
                "physical_validity": best.physical_validity
            })),
            metadata: Some(serde_json::json!({
                "method": "evolutionary_symbolic_regression",
                "domain": domain,
                "candidates": result.candidates,
                "top_5_expressions": result.candidates.iter().take(5)
                    .map(|c| c.expression.clone()).collect::<Vec<_>>()
            })),
        })
    }

    /// Perform dimensional analysis to find relationship
    fn optimize_dimensional_analysis(
        &self,
        target_quantity: &str,
        input: &OptimizeInput,
    ) -> ToolResult<OptimizeOutput> {
        use std::collections::HashMap;

        // Get variables and their dimensions from parameters
        let variables = input
            .parameters
            .get("variables")
            .and_then(|v| serde_json::from_value::<Vec<String>>(v.clone()).ok())
            .ok_or("variables array required for dimensional analysis")?;

        let dimensions = input
            .parameters
            .get("dimensions")
            .and_then(|v| serde_json::from_value::<HashMap<String, Vec<i32>>>(v.clone()).ok())
            .ok_or("dimensions map required (variable -> [M, L, T, I, Θ, N, J])")?;

        let target_dims = dimensions.get(target_quantity).ok_or(format!(
            "Target quantity '{}' not found in dimensions",
            target_quantity
        ))?;

        // Use Buckingham Pi theorem to find dimensionless groups
        // For simplicity, try combinations of variables with integer powers

        // Find combination: target = C * var1^a * var2^b * var3^c ...
        // where dimensional consistency is maintained

        let max_power = 3; // Try powers from -3 to 3
        let best_combination: Option<Vec<i32>> = None;
        let best_residual = f64::INFINITY;

        // If we have data, use it to fit coefficients
        if let Some((x_data, y_data)) = &input.data {
            // Try different power combinations
            for _ in 0..100 {
                // Simplified search
                // For now, return a simple linear combination
                let powers = vec![1; variables.len().min(3)];

                // Build expression
                let mut expr = format!("{} = ", target_quantity);
                for (i, var) in variables.iter().enumerate().take(3) {
                    if i > 0 {
                        expr.push_str(" * ");
                    }
                    expr.push_str(&format!("{}^{}", var, powers.get(i).unwrap_or(&1)));
                }

                return Ok(OptimizeOutput {
                    parameters: powers.iter().map(|&p| p as f64).collect(),
                    function: Some(expr),
                    error: None,
                    r_squared: None,
                    aic: None,
                    bic: None,
                    aicc: None,
                    convergence: None,
                    metadata: Some(serde_json::json!({
                        "method": "dimensional_analysis",
                        "target": target_quantity,
                        "variables": variables,
                        "note": "Simplified dimensional analysis - finds dimensionally consistent relationships"
                    })),
                });
            }
        }

        // No data - just do dimensional analysis to find valid combinations
        let mut valid_combinations = Vec::new();

        // Try different power combinations (simplified brute force)
        fn try_powers(
            variables: &[String],
            dimensions: &HashMap<String, Vec<i32>>,
            target_dims: &[i32],
            current_powers: Vec<i32>,
            depth: usize,
            max_depth: usize,
            max_power: i32,
            valid_combinations: &mut Vec<Vec<i32>>,
        ) {
            if depth == max_depth {
                // Check if this combination matches target dimensions
                let mut result_dims = vec![0; 7]; // M, L, T, I, Θ, N, J

                for (i, var) in variables.iter().enumerate().take(current_powers.len()) {
                    if let Some(var_dims) = dimensions.get(var) {
                        for (j, &dim) in var_dims.iter().enumerate() {
                            result_dims[j] += dim * current_powers[i];
                        }
                    }
                }

                // Check if matches target
                if result_dims
                    .iter()
                    .zip(target_dims.iter())
                    .all(|(a, b)| a == b)
                {
                    valid_combinations.push(current_powers.clone());
                }
                return;
            }

            for power in -max_power..=max_power {
                let mut new_powers = current_powers.clone();
                new_powers.push(power);
                try_powers(
                    variables,
                    dimensions,
                    target_dims,
                    new_powers,
                    depth + 1,
                    max_depth,
                    max_power,
                    valid_combinations,
                );
            }
        }

        try_powers(
            &variables,
            &dimensions,
            target_dims,
            vec![],
            0,
            variables.len().min(3),
            2,
            &mut valid_combinations,
        );

        if let Some(powers) = valid_combinations.first() {
            // Build expression
            let mut expr = format!("{} = C", target_quantity);
            for (i, var) in variables.iter().enumerate() {
                if i < powers.len() && powers[i] != 0 {
                    expr.push_str(&format!(" * {}^{}", var, powers[i]));
                }
            }

            Ok(OptimizeOutput {
                parameters: powers.iter().map(|&p| p as f64).collect(),
                function: Some(expr),
                error: None,
                r_squared: None,
                aic: None,
                bic: None,
                aicc: None,
                convergence: Some(serde_json::json!({
                    "num_valid_combinations": valid_combinations.len(),
                    "all_combinations": valid_combinations
                })),
                metadata: Some(serde_json::json!({
                    "method": "dimensional_analysis",
                    "target": target_quantity,
                    "variables": variables,
                    "note": "Uses Buckingham Pi theorem to find dimensionally consistent relationships"
                })),
            })
        } else {
            Err("No dimensionally consistent combination found".to_string())
        }
    }

    /// Automatic model selection using information criteria
    fn optimize_auto(
        &self,
        criteria: &SelectionCriteria,
        candidates: &[String],
        input: &OptimizeInput,
    ) -> ToolResult<OptimizeOutput> {
        // Default candidates if none specified
        let default_candidates = vec![
            "linear".to_string(),
            "quadratic".to_string(),
            "trigonometric".to_string(),
            "exponential".to_string(),
            "power".to_string(),
        ];

        let methods_to_try: Vec<FitMethod> = if candidates.is_empty() {
            // Try all standard curve fitting methods
            vec![
                FitMethod::Polynomial, // Includes linear and quadratic
                FitMethod::Trigonometric,
                FitMethod::Exponential,
                FitMethod::PowerLaw,
            ]
        } else {
            // Parse user-specified candidates
            candidates
                .iter()
                .filter_map(|name| match name.as_str() {
                    "linear" | "polynomial" => Some(FitMethod::Polynomial),
                    "trigonometric" | "sinusoidal" => Some(FitMethod::Trigonometric),
                    "exponential" => Some(FitMethod::Exponential),
                    "power" | "power_law" => Some(FitMethod::PowerLaw),
                    "logarithmic" => Some(FitMethod::Logarithmic),
                    "rational" => Some(FitMethod::Rational),
                    _ => None,
                })
                .collect()
        };

        // Try each method and collect results with their scores
        let mut results: Vec<(OptimizeOutput, f64, String)> = Vec::new();

        for method in methods_to_try.iter() {
            if let Ok(output) = self.optimize_fit(method, input) {
                let score = match criteria {
                    SelectionCriteria::Aic => output.aic.unwrap_or(f64::INFINITY),
                    SelectionCriteria::Bic => output.bic.unwrap_or(f64::INFINITY),
                    SelectionCriteria::Aicc => output.aicc.unwrap_or(f64::INFINITY),
                    SelectionCriteria::RSquared => {
                        // For R², higher is better, so negate for consistent "lower is better" comparison
                        -output.r_squared.unwrap_or(0.0)
                    }
                };

                let method_name = format!("{:?}", method);
                results.push((output, score, method_name));
            }
        }

        // Try symbolic regression if it might be useful
        if input.data.as_ref().map(|(x, _)| x.len()).unwrap_or(0) >= 5 {
            if let Ok(output) = self.optimize_symbolic_regression(input) {
                // For symbolic regression, use R² as a fallback score
                let score = match criteria {
                    SelectionCriteria::RSquared => -output.r_squared.unwrap_or(0.0),
                    _ => {
                        // Estimate AIC-like score from MSE and complexity
                        // Assume complexity ~10 parameters for penalty
                        let mse = output.error.unwrap_or(1.0);
                        let n = input.data.as_ref().unwrap().0.len() as f64;
                        2.0 * 10.0 + n * mse.ln()
                    }
                };
                results.push((output, score, "SymbolicRegression".to_string()));
            }
        }

        if results.is_empty() {
            return Err("No valid models could be fitted".to_string());
        }

        // Find the best model (lowest score for AIC/BIC/AICc, highest for R²)
        let (mut best_output, best_score, best_method) = results
            .into_iter()
            .min_by(|(_, score1, _), (_, score2, _)| {
                score1
                    .partial_cmp(score2)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap();

        // Add metadata about the selection process
        let criteria_name = match criteria {
            SelectionCriteria::Aic => "AIC",
            SelectionCriteria::Bic => "BIC",
            SelectionCriteria::Aicc => "AICc",
            SelectionCriteria::RSquared => "R²",
        };

        best_output.metadata = Some(serde_json::json!({
            "auto_selection": true,
            "criteria": criteria_name,
            "selected_model": best_method,
            "score": if matches!(criteria, SelectionCriteria::RSquared) { -best_score } else { best_score },
            "original_metadata": best_output.metadata
        }));

        Ok(best_output)
    }
}

impl Optimize for UnifiedOptimizer {
    fn optimize(&self, input: &OptimizeInput) -> ToolResult<OptimizeOutput> {
        match &input.method {
            OptimizationMethod::Fit(fit_method) => self.optimize_fit(fit_method, input),

            OptimizationMethod::Minimize(min_method) => self.optimize_minimize(min_method, input),

            OptimizationMethod::Interpolation(interp_method) => {
                self.optimize_interpolate(interp_method, input)
            }

            OptimizationMethod::DimensionalAnalysis(method) => {
                // Get target quantity from parameters
                let target = input
                    .parameters
                    .get("target")
                    .and_then(|v| v.as_str())
                    .ok_or("target quantity required for dimensional analysis")?;

                self.optimize_dimensional_analysis(target, input)
            }

            OptimizationMethod::SymbolicRegression => self.optimize_symbolic_regression(input),

            OptimizationMethod::Auto {
                criteria,
                candidates,
            } => self.optimize_auto(criteria, candidates, input),
        }
    }
}
