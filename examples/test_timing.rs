//! Example demonstrating compute time tracking for billing purposes
//!
//! This example shows how the MCP server now returns timing information
//! with every response, allowing customers to be billed based on compute time.
//!
//! Run with: cargo run --example test_timing

use computational_engine::implementations::create_default_dispatcher;

fn main() {
    println!("=== Computational Engine - Compute Time Tracking Demo ===\n");

    let dispatcher = create_default_dispatcher();

    // Test 1: Simple equation solving
    println!("Test 1: Solve x^2 - 4 = 0");
    let request1 = r#"{
        "tool": "solve",
        "input": {
            "equation_type": "linear_system",
            "equations": ["x^2 - 4 = 0"],
            "variables": ["x"]
        }
    }"#;

    let response1 = dispatcher.dispatch_json(request1);
    println!("Response:");
    println!("{}\n", serde_json::to_string_pretty(&serde_json::from_str::<serde_json::Value>(&response1).unwrap()).unwrap());

    // Test 2: FFT transform (more compute intensive)
    println!("\nTest 2: FFT Transform (256 samples)");
    let data: Vec<f64> = (0..256).map(|i| (i as f64 * 0.1).sin()).collect();
    let request2 = format!(
        r#"{{
        "tool": "transform",
        "input": {{
            "transform_type": {{"fft": "forward"}},
            "data": {:?},
            "sampling_rate": 1000.0
        }}
    }}"#,
        data
    );

    let response2 = dispatcher.dispatch_json(&request2);
    let response_json: serde_json::Value = serde_json::from_str(&response2).unwrap();
    println!("Response (timing info only):");
    println!("  Success: {}", response_json["success"]);
    println!("  Compute Time: {} ms", response_json["compute_time_ms"]);
    println!("  Compute Time: {} μs", response_json["compute_time_us"]);
    println!("  Compute Time: {} ns", response_json["compute_time_ns"]);

    // Test 3: Matrix eigenvalue decomposition
    println!("\nTest 3: Matrix Eigenvalue Decomposition");
    let request3 = r#"{
        "tool": "compute",
        "input": {
            "operation": {"matrix_decomp": "eigen"},
            "data": {
                "matrix": [[4.0, 1.0], [1.0, 3.0]]
            }
        }
    }"#;

    let response3 = dispatcher.dispatch_json(request3);
    let response_json: serde_json::Value = serde_json::from_str(&response3).unwrap();
    println!("Response (timing info only):");
    println!("  Success: {}", response_json["success"]);
    println!("  Compute Time: {} ms", response_json["compute_time_ms"]);
    println!("  Compute Time: {} μs", response_json["compute_time_us"]);
    println!("  Compute Time: {} ns", response_json["compute_time_ns"]);

    // Test 4: Invalid request (should still return timing)
    println!("\nTest 4: Invalid Request (error handling)");
    let request4 = r#"{"invalid": "json"}"#;

    let response4 = dispatcher.dispatch_json(request4);
    let response_json: serde_json::Value = serde_json::from_str(&response4).unwrap();
    println!("Response:");
    println!("  Success: {}", response_json["success"]);
    println!("  Error: {}", response_json["error"].as_str().unwrap_or("None"));
    println!("  Compute Time: {} ms", response_json["compute_time_ms"]);

    println!("\n=== Billing Use Case ===");
    println!("The compute_time_ns field provides nanosecond precision for accurate billing.");
    println!("Example billing calculation:");
    println!("  - Rate: $0.001 per millisecond of compute");
    println!("  - Compute time: 2.5 ms");
    println!("  - Cost: $0.0025");
    println!("\nAll responses now include:");
    println!("  - compute_time_ms: floating-point milliseconds (for display)");
    println!("  - compute_time_us: integer microseconds (for moderate precision)");
    println!("  - compute_time_ns: integer nanoseconds (for highest precision billing)");
}
