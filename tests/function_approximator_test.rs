use computational_engine::api::process_json_request;
use serde_json::json;

#[test]
fn test_function_approximator_basic() {
    // Test approximating x^2 + 1 (from the video: 0->1, 1->2, 3->10, 5->26)
    let request = json!({
        "module": "function_approximator",
        "operation": "approximate",
        "parameters": {
            "inputs": [0.0, 1.0, 3.0, 5.0],
            "outputs": [1.0, 2.0, 10.0, 26.0],
            "population_size": 600,
            "generations": 100,
            "complexity": 2.0,
            "death_probability": 0.5,
            "diversity": 10
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();

    eprintln!(
        "Response: {}",
        serde_json::to_string_pretty(&parsed).unwrap()
    );

    assert!(parsed["success"].as_bool().unwrap_or(false));

    // Check that we got best functions
    let best_functions = &parsed["result"]["best_functions"];
    assert!(best_functions.is_array());
    assert!(best_functions.as_array().unwrap().len() > 0);

    // Check first function has required fields
    let best = &best_functions[0];
    assert!(best["expression"].is_string());
    assert!(best["fitness"].is_number());
    assert!(best["complexity"].is_number());
    assert!(best["predictions"].is_array());

    eprintln!("Best function: {}", best["expression"]);
    eprintln!("Fitness: {}", best["fitness"]);
    eprintln!("Complexity: {}", best["complexity"]);
}

#[test]
fn test_function_approximator_simple() {
    // Test approximating 2x (0->0, 1->2, 2->4, 3->6)
    let request = json!({
        "module": "function_approximator",
        "operation": "approximate",
        "parameters": {
            "inputs": [0.0, 1.0, 2.0, 3.0],
            "outputs": [0.0, 2.0, 4.0, 6.0],
            "population_size": 200,
            "generations": 50,
            "complexity": 2.0
        }
    });

    let response = process_json_request(&request.to_string());
    let parsed: serde_json::Value = serde_json::from_str(&response).unwrap();

    assert!(parsed["success"].as_bool().unwrap_or(false));

    eprintln!(
        "Best function for 2x: {}",
        parsed["result"]["best_functions"][0]["expression"]
    );
}
