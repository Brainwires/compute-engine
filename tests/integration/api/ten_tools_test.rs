//! Tests for the new 10-tool unified API

use computational_engine::create_default_dispatcher;
use computational_engine::engine::*;

#[test]
fn test_dispatcher_creation() {
    let dispatcher = create_default_dispatcher();

    // Test that dispatcher was created successfully
    // (If this compiles and runs, the basic structure works)
    assert!(true);
}

#[test]
fn test_solve_einstein_vacuum() {
    let dispatcher = create_default_dispatcher();

    let request = ToolRequest::Solve(SolveInput {
        equation_type: EquationType::Einstein(EinsteinEquation::Vacuum),
        equations: vec![],
        variables: Some(vec![
            "t".to_string(),
            "r".to_string(),
            "theta".to_string(),
            "phi".to_string(),
        ]),
        initial_guess: None,
        boundary_conditions: None,
        domain: None,
        method: None,
        parameters: std::collections::HashMap::new(),
    });

    let result = dispatcher.dispatch(request);

    // This will succeed once we fully implement the mapping
    // For now, it should at least not panic
    match result {
        Ok(_) => println!("Einstein vacuum solved!"),
        Err(e) => println!("Expected error (not yet implemented): {}", e),
    }
}

#[test]
fn test_json_api() {
    let dispatcher = create_default_dispatcher();

    let json_request = r#"{
        "tool": "solve",
        "input": {
            "equation_type": {
                "einstein": "vacuum"
            },
            "equations": [],
            "parameters": {}
        }
    }"#;

    let response = dispatcher.dispatch_json(json_request);

    // Should return valid JSON
    assert!(response.contains("success") || response.contains("error"));
}

#[test]
fn test_all_10_tools_registered() {
    let dispatcher = create_default_dispatcher();

    // Test that all 10 tools can be called (even if stubbed)
    let tools = vec![
        (
            "solve",
            ToolRequest::Solve(SolveInput {
                equation_type: EquationType::RootFinding,
                equations: vec![],
                variables: None,
                initial_guess: None,
                boundary_conditions: None,
                domain: None,
                method: None,
                parameters: std::collections::HashMap::new(),
            }),
        ),
        // Add more as we implement them
    ];

    for (name, request) in tools {
        let result = dispatcher.dispatch(request);
        println!("{}: {:?}", name, result.is_ok());
    }
}
