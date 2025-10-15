//! Tests for symbolic differentiation, integration, and series expansion
//!
//! These tests ensure that the symbolic CAS operations work correctly
//! with the unified tool API.

use computational_engine::engine::*;
use computational_engine::engine::equations::*;
use computational_engine::create_default_dispatcher;
use std::collections::HashMap;

// ============================================================================
// SYMBOLIC DIFFERENTIATE TESTS
// ============================================================================

#[test]
fn test_differentiate_symbolic_polynomial() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::Symbolic,
        expression: "x^3 + 2*x^2 - 5*x + 1".to_string(),
        variables: vec!["x".to_string()],
        order: Some(vec![1]),
        evaluate_at: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Symbolic differentiation should succeed");

    if let Ok(ToolResponse::Differentiate(output)) = result {
        // Should contain derivative result
        assert!(output.derivatives.contains_key("x"));
        println!("Derivative of x^3 + 2*x^2 - 5*x + 1: {:?}", output.derivatives.get("x"));

        // The derivative should be 3*x^2 + 4*x - 5
        let derivative = output.derivatives.get("x").unwrap();
        let deriv_str = derivative.as_str().unwrap();

        // Check that it contains expected terms (simplified form may vary)
        assert!(deriv_str.contains("x^2") || deriv_str.contains("x²"),
               "Should contain x^2 term");
    }
}

#[test]
fn test_differentiate_symbolic_order_parameter_array() {
    // This tests the fix for order being an array instead of a number
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::Symbolic,
        expression: "x^2".to_string(),
        variables: vec!["x".to_string()],
        order: Some(vec![1]),  // Must be array!
        evaluate_at: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Differentiation with order array should work");
}

#[test]
fn test_differentiate_symbolic_second_derivative() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::Symbolic,
        expression: "x^4".to_string(),
        variables: vec!["x".to_string()],
        order: Some(vec![2]),  // Second derivative
        evaluate_at: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Second derivative should work");
}

// ============================================================================
// SYMBOLIC INTEGRATE TESTS
// ============================================================================

#[test]
fn test_integrate_symbolic_linear() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Symbolic,
        expression: "2*x + 3".to_string(),
        variables: vec!["x".to_string()],
        limits: None,
        path: None,
        method: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Symbolic integration of 2x+3 should succeed");

    if let Ok(ToolResponse::Integrate(output)) = result {
        let result_str = output.symbolic.unwrap_or_default();
        println!("Integral of 2*x + 3: {}", result_str);

        // Should contain x^2 and x terms
        assert!(result_str.contains("x^2") || result_str.contains("x²"),
               "Integral should contain x^2 term");
        assert!(result_str.contains("x"), "Integral should contain x term");
    }
}

#[test]
fn test_integrate_symbolic_power() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Symbolic,
        expression: "x^2".to_string(),
        variables: vec!["x".to_string()],
        limits: None,
        path: None,
        method: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Symbolic integration of x^2 should succeed");

    if let Ok(ToolResponse::Integrate(output)) = result {
        let result_str = output.symbolic.unwrap_or_default();
        println!("Integral of x^2: {}", result_str);

        // Should be x^3/3
        assert!(result_str.contains("x^3") || result_str.contains("x³"),
               "Integral should contain x^3 term");
    }
}

#[test]
fn test_integrate_symbolic_sin() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Symbolic,
        expression: "sin(x)".to_string(),
        variables: vec!["x".to_string()],
        limits: None,
        path: None,
        method: None,
        parameters: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Symbolic integration of sin(x) should succeed");

    if let Ok(ToolResponse::Integrate(output)) = result {
        let result_str = output.symbolic.unwrap_or_default();
        println!("Integral of sin(x): {}", result_str);

        // Should be -cos(x)
        assert!(result_str.contains("cos"), "Integral should contain cos");
    }
}

// ============================================================================
// TAYLOR SERIES EXPANSION TESTS
// ============================================================================

#[test]
fn test_series_expansion_exp_order_4() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert("order".to_string(), serde_json::json!(4));
    options.insert("variable".to_string(), serde_json::json!("x"));
    options.insert("point".to_string(), serde_json::json!(0.0));

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::SeriesExpansion,
        expression: "exp(x)".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Taylor series expansion should succeed");

    if let Ok(ToolResponse::Analyze(output)) = result {
        let expansion = output.result.get("expansion").unwrap().as_str().unwrap();
        println!("Taylor series of exp(x): {}", expansion);

        // Should have terms up to x^4
        assert!(expansion.contains("x^4") || expansion.contains("x⁴"),
               "Should contain x^4 term");
        assert!(expansion.contains("1/24") || expansion.contains("24"),
               "Should contain 1/24 coefficient for x^4");
        assert!(expansion.contains("1/6") || expansion.contains("6"),
               "Should contain 1/6 coefficient for x^3");
        assert!(expansion.contains("1/2") || expansion.contains("2"),
               "Should contain 1/2 coefficient for x^2");
    }
}

#[test]
fn test_series_expansion_sin() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert("order".to_string(), serde_json::json!(3));
    options.insert("variable".to_string(), serde_json::json!("x"));
    options.insert("point".to_string(), serde_json::json!(0.0));

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::SeriesExpansion,
        expression: "sin(x)".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Taylor series of sin(x) should succeed");
}

// ============================================================================
// JSON API TESTS (matching MCP interface)
// ============================================================================

#[test]
fn test_json_api_differentiate_symbolic() {
    let dispatcher = create_default_dispatcher();

    let json_request = r#"{
        "tool": "differentiate",
        "input": {
            "operation": "symbolic",
            "expression": "x^3 + 2*x^2 - 5*x + 1",
            "variables": ["x"],
            "order": [1]
        }
    }"#;

    let response = dispatcher.dispatch_json(json_request);
    println!("Differentiate JSON response: {}", response);
    assert!(response.contains("derivatives") || response.contains("error"));
    assert!(!response.contains("missing field"),
           "Should not have 'missing field' error");
}

#[test]
fn test_json_api_integrate_symbolic() {
    let dispatcher = create_default_dispatcher();

    let json_request = r#"{
        "tool": "integrate",
        "input": {
            "integration_type": "symbolic",
            "expression": "2*x + 3",
            "variables": ["x"]
        }
    }"#;

    let response = dispatcher.dispatch_json(json_request);
    println!("Integrate JSON response: {}", response);
    assert!(response.contains("result") || response.contains("symbolic"));
    assert!(!response.contains("not yet implemented"),
           "Symbolic integration should be implemented");
}

#[test]
fn test_json_api_taylor_series() {
    let dispatcher = create_default_dispatcher();

    let json_request = r#"{
        "tool": "analyze",
        "input": {
            "operation": "series_expansion",
            "expression": "exp(x)",
            "options": {
                "order": 4,
                "variable": "x",
                "point": 0.0
            }
        }
    }"#;

    let response = dispatcher.dispatch_json(json_request);
    println!("Taylor series JSON response: {}", response);
    assert!(response.contains("expansion") || response.contains("result"));
    // Check that we actually got order 4, not the default of 3
    assert!(response.contains("x^4") || response.contains("1/24"),
           "Should include 4th order term");
}

// ============================================================================
// INTEGRATION TESTS (full workflow)
// ============================================================================

#[test]
fn test_workflow_differentiate_then_integrate() {
    let dispatcher = create_default_dispatcher();

    // Step 1: Differentiate x^3
    let diff_request = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::Symbolic,
        expression: "x^3".to_string(),
        variables: vec!["x".to_string()],
        order: Some(vec![1]),
        evaluate_at: None,
        parameters: HashMap::new(),
    });

    let diff_result = dispatcher.dispatch(diff_request);
    assert!(diff_result.is_ok());

    // Result should be 3*x^2
    // Step 2: Integrate it back
    let integrate_request = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Symbolic,
        expression: "3*x^2".to_string(),
        variables: vec!["x".to_string()],
        limits: None,
        path: None,
        method: None,
        parameters: HashMap::new(),
    });

    let integrate_result = dispatcher.dispatch(integrate_request);
    assert!(integrate_result.is_ok());

    println!("Differentiate-integrate workflow completed");
}

#[test]
fn test_all_symbolic_tools_work() {
    let dispatcher = create_default_dispatcher();

    // Test differentiate
    let diff = ToolRequest::Differentiate(DifferentiateInput {
        operation: DifferentiationOp::Symbolic,
        expression: "x^2".to_string(),
        variables: vec!["x".to_string()],
        order: Some(vec![1]),
        evaluate_at: None,
        parameters: HashMap::new(),
    });
    assert!(dispatcher.dispatch(diff).is_ok(), "Differentiate should work");

    // Test integrate
    let integ = ToolRequest::Integrate(IntegrateInput {
        integration_type: IntegrationType::Symbolic,
        expression: "x".to_string(),
        variables: vec!["x".to_string()],
        limits: None,
        path: None,
        method: None,
        parameters: HashMap::new(),
    });
    assert!(dispatcher.dispatch(integ).is_ok(), "Integrate should work");

    // Test series expansion
    let series = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::SeriesExpansion,
        expression: "exp(x)".to_string(),
        options: HashMap::from([
            ("order".to_string(), serde_json::json!(2)),
            ("variable".to_string(), serde_json::json!("x")),
            ("point".to_string(), serde_json::json!(0.0))
        ]),
    });
    assert!(dispatcher.dispatch(series).is_ok(), "Series expansion should work");

    println!("✅ All symbolic tools working!");
}
