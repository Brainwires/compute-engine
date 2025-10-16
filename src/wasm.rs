//! WebAssembly bindings for the computational engine
//!
//! This module provides JavaScript-friendly bindings for browser/Node.js environments.
//!
//! # Usage
//!
//! ```javascript
//! import init, { ComputationalEngine } from './computational_engine.js';
//!
//! await init();
//! const engine = new ComputationalEngine();
//!
//! // Solve an equation
//! const result = engine.solve({
//!   equations: ["x^2 - 4 = 0"]
//! });
//! console.log(result);
//!
//! // Process raw JSON
//! const jsonResult = engine.processJson(JSON.stringify({
//!   tool: "solve",
//!   input: { equations: ["x^2 - 4 = 0"] }
//! }));
//! ```

use crate::{
    api::process_json_request,
    engine::{ToolRequest, ToolResponse},
    implementations::create_default_dispatcher,
};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

/// Set panic hook for better error messages in WASM
#[wasm_bindgen(start)]
pub fn start() {
    #[cfg(feature = "console_error_panic_hook")]
    console_error_panic_hook::set_once();
}

/// Main computational engine interface for WebAssembly
#[wasm_bindgen]
pub struct ComputationalEngine {
    // We don't store the dispatcher because it's stateless
}

#[wasm_bindgen]
impl ComputationalEngine {
    /// Create a new computational engine instance
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {}
    }

    /// Process a JSON request and return a JSON response
    ///
    /// # Arguments
    /// * `json_request` - JSON string containing the tool request
    ///
    /// # Returns
    /// JSON string with the computation result
    ///
    /// # Example
    /// ```javascript
    /// const result = engine.processJson(JSON.stringify({
    ///   tool: "solve",
    ///   input: { equations: ["x^2 - 4 = 0"] }
    /// }));
    /// ```
    #[wasm_bindgen(js_name = processJson)]
    pub fn process_json(&self, json_request: &str) -> String {
        // Try new unified API first
        match serde_json::from_str::<ToolRequest>(json_request) {
            Ok(tool_request) => {
                let dispatcher = create_default_dispatcher();
                match dispatcher.dispatch(tool_request) {
                    Ok(response) => serde_json::to_string(&response).unwrap_or_else(|e| {
                        format!(
                            r#"{{"success":false,"error":"Serialization error: {}"}}"#,
                            e
                        )
                    }),
                    Err(e) => {
                        format!(r#"{{"success":false,"error":"{}"}}"#, e)
                    }
                }
            }
            Err(_) => {
                // Fall back to legacy API
                process_json_request(json_request)
            }
        }
    }

    /// Solve equations
    ///
    /// # Arguments
    /// * `input` - JsValue containing SolveInput structure
    ///
    /// # Returns
    /// JsValue with the solution result
    #[wasm_bindgen(js_name = solve)]
    pub fn solve(&self, input: JsValue) -> Result<JsValue, JsValue> {
        let solve_input: crate::engine::SolveInput = serde_wasm_bindgen::from_value(input)
            .map_err(|e| JsValue::from_str(&format!("Invalid input: {}", e)))?;

        let request = ToolRequest::Solve(solve_input);
        let dispatcher = create_default_dispatcher();

        match dispatcher.dispatch(request) {
            Ok(response) => serde_wasm_bindgen::to_value(&response)
                .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e))),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    }

    /// Differentiate expressions
    #[wasm_bindgen(js_name = differentiate)]
    pub fn differentiate(&self, input: JsValue) -> Result<JsValue, JsValue> {
        let diff_input: crate::engine::DifferentiateInput =
            serde_wasm_bindgen::from_value(input)
                .map_err(|e| JsValue::from_str(&format!("Invalid input: {}", e)))?;

        let request = ToolRequest::Differentiate(diff_input);
        let dispatcher = create_default_dispatcher();

        match dispatcher.dispatch(request) {
            Ok(response) => serde_wasm_bindgen::to_value(&response)
                .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e))),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    }

    /// Integrate expressions
    #[wasm_bindgen(js_name = integrate)]
    pub fn integrate(&self, input: JsValue) -> Result<JsValue, JsValue> {
        let int_input: crate::engine::IntegrateInput = serde_wasm_bindgen::from_value(input)
            .map_err(|e| JsValue::from_str(&format!("Invalid input: {}", e)))?;

        let request = ToolRequest::Integrate(int_input);
        let dispatcher = create_default_dispatcher();

        match dispatcher.dispatch(request) {
            Ok(response) => serde_wasm_bindgen::to_value(&response)
                .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e))),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    }

    /// Analyze expressions (simplify, expand, factor, etc.)
    #[wasm_bindgen(js_name = analyze)]
    pub fn analyze(&self, input: JsValue) -> Result<JsValue, JsValue> {
        let analyze_input: crate::engine::AnalyzeInput = serde_wasm_bindgen::from_value(input)
            .map_err(|e| JsValue::from_str(&format!("Invalid input: {}", e)))?;

        let request = ToolRequest::Analyze(analyze_input);
        let dispatcher = create_default_dispatcher();

        match dispatcher.dispatch(request) {
            Ok(response) => serde_wasm_bindgen::to_value(&response)
                .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e))),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    }

    /// Simulate systems (ODEs, PDEs, physics, etc.)
    #[wasm_bindgen(js_name = simulate)]
    pub fn simulate(&self, input: JsValue) -> Result<JsValue, JsValue> {
        let sim_input: crate::engine::SimulateInput = serde_wasm_bindgen::from_value(input)
            .map_err(|e| JsValue::from_str(&format!("Invalid input: {}", e)))?;

        let request = ToolRequest::Simulate(sim_input);
        let dispatcher = create_default_dispatcher();

        match dispatcher.dispatch(request) {
            Ok(response) => serde_wasm_bindgen::to_value(&response)
                .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e))),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    }

    /// Compute (tensors, matrices, special functions)
    #[wasm_bindgen(js_name = compute)]
    pub fn compute(&self, input: JsValue) -> Result<JsValue, JsValue> {
        let compute_input: crate::engine::ComputeInput = serde_wasm_bindgen::from_value(input)
            .map_err(|e| JsValue::from_str(&format!("Invalid input: {}", e)))?;

        let request = ToolRequest::Compute(compute_input);
        let dispatcher = create_default_dispatcher();

        match dispatcher.dispatch(request) {
            Ok(response) => serde_wasm_bindgen::to_value(&response)
                .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e))),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    }

    /// Transform (FFT, Fourier, Laplace, wavelets)
    #[wasm_bindgen(js_name = transform)]
    pub fn transform(&self, input: JsValue) -> Result<JsValue, JsValue> {
        let transform_input: crate::engine::TransformInput = serde_wasm_bindgen::from_value(input)
            .map_err(|e| JsValue::from_str(&format!("Invalid input: {}", e)))?;

        let request = ToolRequest::Transform(transform_input);
        let dispatcher = create_default_dispatcher();

        match dispatcher.dispatch(request) {
            Ok(response) => serde_wasm_bindgen::to_value(&response)
                .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e))),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    }

    /// Field theory (EM fields, quantum fields)
    #[wasm_bindgen(js_name = fieldtheory)]
    pub fn fieldtheory(&self, input: JsValue) -> Result<JsValue, JsValue> {
        let field_input: crate::engine::FieldTheoryInput = serde_wasm_bindgen::from_value(input)
            .map_err(|e| JsValue::from_str(&format!("Invalid input: {}", e)))?;

        let request = ToolRequest::FieldTheory(field_input);
        let dispatcher = create_default_dispatcher();

        match dispatcher.dispatch(request) {
            Ok(response) => serde_wasm_bindgen::to_value(&response)
                .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e))),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    }

    /// Sample (Monte Carlo, statistics)
    #[wasm_bindgen(js_name = sample)]
    pub fn sample(&self, input: JsValue) -> Result<JsValue, JsValue> {
        let sample_input: crate::engine::SampleInput = serde_wasm_bindgen::from_value(input)
            .map_err(|e| JsValue::from_str(&format!("Invalid input: {}", e)))?;

        let request = ToolRequest::Sample(sample_input);
        let dispatcher = create_default_dispatcher();

        match dispatcher.dispatch(request) {
            Ok(response) => serde_wasm_bindgen::to_value(&response)
                .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e))),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    }

    /// Optimize (curve fitting, minimization)
    #[wasm_bindgen(js_name = optimize)]
    pub fn optimize(&self, input: JsValue) -> Result<JsValue, JsValue> {
        let optimize_input: crate::engine::OptimizeInput = serde_wasm_bindgen::from_value(input)
            .map_err(|e| JsValue::from_str(&format!("Invalid input: {}", e)))?;

        let request = ToolRequest::Optimize(optimize_input);
        let dispatcher = create_default_dispatcher();

        match dispatcher.dispatch(request) {
            Ok(response) => serde_wasm_bindgen::to_value(&response)
                .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e))),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    }

    /// Get version information
    #[wasm_bindgen(js_name = version)]
    pub fn version(&self) -> String {
        crate::VERSION.to_string()
    }

    /// List all available operations
    #[wasm_bindgen(js_name = listOperations)]
    pub fn list_operations(&self) -> JsValue {
        let ops = crate::api::list_all_operations();
        serde_wasm_bindgen::to_value(&ops).unwrap_or(JsValue::NULL)
    }
}

impl Default for ComputationalEngine {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_engine_creation() {
        let engine = ComputationalEngine::new();
        assert_eq!(engine.version(), crate::VERSION);
    }

    #[test]
    fn test_json_processing() {
        let engine = ComputationalEngine::new();
        let request = r#"{"tool":"solve","input":{"equations":["x^2 - 4 = 0"]}}"#;
        let response = engine.process_json(request);
        assert!(!response.is_empty());
    }
}
