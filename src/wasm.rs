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
    create_default_dispatcher,
    engine::ToolRequest,
};
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
            Err(e) => {
                format!(r#"{{"success":false,"error":"Parse error: {}"}}"#, e)
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

    /// Compute operations (calculus, matrices, transforms, field theory, sampling)
    ///
    /// # Arguments
    /// * `input` - JsValue containing ComputeInput structure
    ///
    /// # Returns
    /// JsValue with the computation result
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

    /// Analyze expressions (simplify, expand, factor, series, etc.)
    ///
    /// # Arguments
    /// * `input` - JsValue containing AnalyzeInput structure
    ///
    /// # Returns
    /// JsValue with the analysis result
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

    /// Simulate systems (ODEs, PDEs, physics, stochastic, finance)
    ///
    /// # Arguments
    /// * `input` - JsValue containing SimulateInput structure
    ///
    /// # Returns
    /// JsValue with the simulation result
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

    /// Machine learning operations (clustering, regression, neural networks)
    ///
    /// # Arguments
    /// * `input` - JsValue containing MLInput structure
    ///
    /// # Returns
    /// JsValue with the ML result
    #[wasm_bindgen(js_name = ml)]
    pub fn ml(&self, input: JsValue) -> Result<JsValue, JsValue> {
        let ml_input: crate::engine::MLInput = serde_wasm_bindgen::from_value(input)
            .map_err(|e| JsValue::from_str(&format!("Invalid input: {}", e)))?;

        let request = ToolRequest::Ml(ml_input);
        let dispatcher = create_default_dispatcher();

        match dispatcher.dispatch(request) {
            Ok(response) => serde_wasm_bindgen::to_value(&response)
                .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e))),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    }

    /// Chaos theory operations (fractals, attractors, Lyapunov exponents)
    ///
    /// # Arguments
    /// * `input` - JsValue containing ChaosInput structure
    ///
    /// # Returns
    /// JsValue with the chaos analysis result
    #[wasm_bindgen(js_name = chaos)]
    pub fn chaos(&self, input: JsValue) -> Result<JsValue, JsValue> {
        let chaos_input: crate::engine::ChaosInput = serde_wasm_bindgen::from_value(input)
            .map_err(|e| JsValue::from_str(&format!("Invalid input: {}", e)))?;

        let request = ToolRequest::Chaos(chaos_input);
        let dispatcher = create_default_dispatcher();

        match dispatcher.dispatch(request) {
            Ok(response) => serde_wasm_bindgen::to_value(&response)
                .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e))),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    }

    /// Unit conversion and dimensional analysis
    ///
    /// # Arguments
    /// * `input` - JsValue containing UnitsInput structure
    ///
    /// # Returns
    /// JsValue with the conversion result
    #[wasm_bindgen(js_name = units)]
    pub fn units(&self, input: JsValue) -> Result<JsValue, JsValue> {
        let units_input: crate::engine::UnitsInput = serde_wasm_bindgen::from_value(input)
            .map_err(|e| JsValue::from_str(&format!("Invalid input: {}", e)))?;

        let request = ToolRequest::Units(units_input);
        let dispatcher = create_default_dispatcher();

        match dispatcher.dispatch(request) {
            Ok(response) => serde_wasm_bindgen::to_value(&response)
                .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e))),
            Err(e) => Err(JsValue::from_str(&e)),
        }
    }

    /// Validate equations and physics
    ///
    /// # Arguments
    /// * `input` - JsValue containing ValidateInput structure
    ///
    /// # Returns
    /// JsValue with the validation result
    #[wasm_bindgen(js_name = validate)]
    pub fn validate(&self, input: JsValue) -> Result<JsValue, JsValue> {
        let validate_input: crate::engine::ValidateInput = serde_wasm_bindgen::from_value(input)
            .map_err(|e| JsValue::from_str(&format!("Invalid input: {}", e)))?;

        let request = ToolRequest::Validate(validate_input);
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
        let ops = crate::engine::list_all_operations();
        serde_wasm_bindgen::to_value(&ops).unwrap_or(JsValue::NULL)
    }
}

impl Default for ComputationalEngine {
    fn default() -> Self {
        Self::new()
    }
}
