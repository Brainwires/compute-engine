//! Comprehensive analyzer test suite
//!
//! Tests for all ANALYZE tool operations including:
//! - Field analysis (vector, scalar, tensor)
//! - Series operations (Taylor, Laurent, partial fractions)
//! - Validation operations (dimensional, physics, conservation)
//! - Simplification and parsing
//! - Limit computation

use computational_engine::engine::*;
use computational_engine::engine::equations::*;
use computational_engine::create_default_dispatcher;
use std::collections::HashMap;

// ============================================================================
// FIELD ANALYSIS TESTS
// ============================================================================

#[test]
fn test_vector_field_analysis() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert("components".to_string(), serde_json::json!(["x", "y", "z"]));
    options.insert("variables".to_string(), serde_json::json!(["x", "y", "z"]));

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::FieldAnalysis(FieldAnalysisType::Vector),
        expression: "vector_field".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should analyze vector field: {:?}", result);

    if let Ok(ToolResponse::Analyze(output)) = result {
        let res = output.result.as_object().unwrap();
        assert_eq!(res.get("field_type").unwrap().as_str().unwrap(), "vector");
        assert!(res.contains_key("divergence"), "Should compute divergence");
        assert!(res.contains_key("curl"), "Should compute curl");
    }
}

#[test]
fn test_scalar_field_analysis() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert("variables".to_string(), serde_json::json!(["x", "y", "z"]));

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::FieldAnalysis(FieldAnalysisType::Scalar),
        expression: "x^2 + y^2 + z^2".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should analyze scalar field: {:?}", result);

    if let Ok(ToolResponse::Analyze(output)) = result {
        let res = output.result.as_object().unwrap();
        assert_eq!(res.get("field_type").unwrap().as_str().unwrap(), "scalar");
        assert!(res.contains_key("gradient"), "Should compute gradient");
        assert!(res.contains_key("laplacian"), "Should compute Laplacian");
    }
}

#[test]
fn test_tensor_field_analysis() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert("dimension".to_string(), serde_json::json!(4));
    options.insert("variables".to_string(), serde_json::json!(["t", "x", "y", "z"]));

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::FieldAnalysis(FieldAnalysisType::Tensor),
        expression: "metric_tensor".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should analyze tensor field: {:?}", result);

    if let Ok(ToolResponse::Analyze(output)) = result {
        let res = output.result.as_object().unwrap();
        assert_eq!(res.get("field_type").unwrap().as_str().unwrap(), "tensor");
    }
}

#[test]
fn test_tensor_field_with_metric() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert("dimension".to_string(), serde_json::json!(2));
    options.insert("variables".to_string(), serde_json::json!(["x", "y"]));
    options.insert("metric".to_string(), serde_json::json!([[1.0, 0.0], [0.0, 1.0]]));

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::FieldAnalysis(FieldAnalysisType::Tensor),
        expression: "metric".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should analyze tensor with metric: {:?}", result);
}

// ============================================================================
// SERIES OPERATIONS TESTS
// ============================================================================

#[test]
fn test_series_expansion_polynomial() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert("variable".to_string(), serde_json::json!("x"));
    options.insert("point".to_string(), serde_json::json!(0.0));
    options.insert("order".to_string(), serde_json::json!(3));

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::SeriesExpansion,
        expression: "x^2 + 2*x".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should expand series: {:?}", result);

    if let Ok(ToolResponse::Analyze(output)) = result {
        assert!(output.result.as_object().unwrap().contains_key("expansion"));
    }
}

#[test]
fn test_series_expansion_trigonometric() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert("variable".to_string(), serde_json::json!("x"));
    options.insert("point".to_string(), serde_json::json!(0.0));
    options.insert("order".to_string(), serde_json::json!(5));

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::SeriesExpansion,
        expression: "sin(x)".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should expand trigonometric series: {:?}", result);
}

#[test]
fn test_laurent_series() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert("variable".to_string(), serde_json::json!("z"));
    options.insert("point".to_string(), serde_json::json!(1.0));
    options.insert("order".to_string(), serde_json::json!(3));

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::LaurentSeries,
        expression: "z^2 + z".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute Laurent series: {:?}", result);

    if let Ok(ToolResponse::Analyze(output)) = result {
        let res = output.result.as_object().unwrap();
        assert!(res.contains_key("analytic_part"), "Should have analytic part");
        assert!(res.contains_key("has_singularity_at_point"), "Should detect singularity");
    }
}

#[test]
fn test_partial_fraction() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert("variable".to_string(), serde_json::json!("x"));

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::PartialFraction,
        expression: "(x + 1) / (x^2 - 1)".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute partial fractions: {:?}", result);

    if let Ok(ToolResponse::Analyze(output)) = result {
        let res = output.result.as_object().unwrap();
        assert!(res.contains_key("factored"), "Should have factored form");
        assert!(res.contains_key("simplified"), "Should have simplified form");
    }
}

// ============================================================================
// LIMIT TESTS
// ============================================================================

#[test]
fn test_limit_computation() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert("variable".to_string(), serde_json::json!("x"));
    options.insert("point".to_string(), serde_json::json!("0"));

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::Limit,
        expression: "sin(x) / x".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute limit: {:?}", result);

    if let Ok(ToolResponse::Analyze(output)) = result {
        assert!(output.result.as_object().unwrap().contains_key("limit_value"));
    }
}

#[test]
fn test_limit_infinity() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert("variable".to_string(), serde_json::json!("x"));
    options.insert("point".to_string(), serde_json::json!("infinity"));

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::Limit,
        expression: "1 / x".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should compute limit to infinity: {:?}", result);
}

// ============================================================================
// SIMPLIFICATION TESTS
// ============================================================================

#[test]
fn test_simplify_expression() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::Simplify,
        expression: "x + x".to_string(),
        options: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simplify expression: {:?}", result);
}

#[test]
fn test_simplify_complex() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::Simplify,
        expression: "(x + 1)^2 - x^2 - 2*x - 1".to_string(),
        options: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simplify complex expression: {:?}", result);
}

// ============================================================================
// PARSING TESTS
// ============================================================================

#[test]
fn test_parse_equation() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::Parse,
        expression: "x^2 + 3*x + 2 = 0".to_string(),
        options: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should parse equation: {:?}", result);

    if let Ok(ToolResponse::Analyze(output)) = result {
        assert!(output.validation.is_some());
        assert!(output.validation.unwrap().is_valid);
    }
}

#[test]
fn test_extract_variables() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::ExtractVariables,
        expression: "x^2 + y*z + 3*a".to_string(),
        options: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should extract variables: {:?}", result);

    if let Ok(ToolResponse::Analyze(output)) = result {
        let vars = output.result.as_array().unwrap();
        assert!(vars.len() > 0, "Should find variables");
    }
}

// ============================================================================
// VALIDATION TESTS
// ============================================================================

#[test]
fn test_validate_equation() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert("domain".to_string(), serde_json::json!("physics"));

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::Validate,
        expression: "F = m * a".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should validate equation: {:?}", result);
}

#[test]
fn test_check_dimensions() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    let mut units = HashMap::new();
    units.insert("x".to_string(), "m".to_string());
    units.insert("v".to_string(), "m/s".to_string());
    options.insert("variable_units".to_string(), serde_json::to_value(units).unwrap());

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::CheckDimensions,
        expression: "x = v * t".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should check dimensions: {:?}", result);
}

#[test]
fn test_check_physics() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::CheckPhysics,
        expression: "energy is conserved".to_string(),
        options: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should check physics validity: {:?}", result);
}

#[test]
fn test_check_conservation_laws() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert("laws".to_string(), serde_json::json!(["energy", "momentum"]));

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::CheckConservation,
        expression: "energy and momentum are conserved".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should check conservation laws: {:?}", result);
}

#[test]
fn test_check_symmetries() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert("symmetries".to_string(), serde_json::json!(["translation", "rotation"]));

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::CheckSymmetries,
        expression: "f(r)".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should check symmetries: {:?}", result);
}

// ============================================================================
// DIMENSIONAL ANALYSIS TESTS
// ============================================================================

#[test]
fn test_infer_dimensions() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::InferDimensions,
        expression: "velocity * acceleration".to_string(),
        options: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should infer dimensions: {:?}", result);
}

#[test]
fn test_scale_analysis() {
    let dispatcher = create_default_dispatcher();

    let mut options = HashMap::new();
    options.insert("scale".to_string(), serde_json::json!(1.0));

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::ScaleAnalysis,
        expression: "x^2".to_string(),
        options,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should perform scale analysis: {:?}", result);
}

// ============================================================================
// NUMBER THEORY TESTS
// ============================================================================

#[test]
fn test_is_prime_true() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::IsPrime,
        expression: "17".to_string(),
        options: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should check primality: {:?}", result);

    if let Ok(ToolResponse::Analyze(output)) = result {
        assert_eq!(output.result.as_bool().unwrap(), true, "17 is prime");
    }
}

#[test]
fn test_is_prime_false() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::IsPrime,
        expression: "15".to_string(),
        options: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should check primality: {:?}", result);

    if let Ok(ToolResponse::Analyze(output)) = result {
        assert_eq!(output.result.as_bool().unwrap(), false, "15 is not prime");
    }
}

// ============================================================================
// ERROR HANDLING TESTS
// ============================================================================

#[test]
fn test_vector_field_requires_components() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::FieldAnalysis(FieldAnalysisType::Vector),
        expression: "vector_field".to_string(),
        options: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Vector field analysis should require components");
}

#[test]
fn test_invalid_prime_number() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::IsPrime,
        expression: "not_a_number".to_string(),
        options: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Should reject non-numeric input for primality test");
}
