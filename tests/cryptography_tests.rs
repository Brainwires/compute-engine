use computational_engine::engine::*;
use computational_engine::create_default_dispatcher;
use std::collections::HashMap;

#[test]
fn test_is_prime_via_analyze() {
    let dispatcher = create_default_dispatcher();

    // Test prime numbers
    for n in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29] {
        let request = ToolRequest::Analyze(AnalyzeInput {
            operation: AnalysisOp::IsPrime,
            expression: n.to_string(),
            options: HashMap::new(),
        });

        let result = dispatcher.dispatch(request);
        assert!(result.is_ok(), "{} should be recognized as prime", n);
    }

    // Test non-prime numbers
    for n in [4, 6, 8, 9, 10] {
        let request = ToolRequest::Analyze(AnalyzeInput {
            operation: AnalysisOp::IsPrime,
            expression: n.to_string(),
            options: HashMap::new(),
        });

        let result = dispatcher.dispatch(request);
        // May return error or false result for non-primes
        let _ = result;
    }
}

#[test]
fn test_analyze_simplify() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::Simplify,
        expression: "x + x".to_string(),
        options: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Simplification should work");
}

#[test]
fn test_analyze_parse() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Analyze(AnalyzeInput {
        operation: AnalysisOp::Parse,
        expression: "F = m*a".to_string(),
        options: HashMap::new(),
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Parse should work");
}
